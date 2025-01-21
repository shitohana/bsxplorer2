use crate::io::bsx::ipc::IpcFileReader;
use crate::io::report::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::io::report::report_batch_utils::first_position;
use crate::region::{GenomicPosition, RegionCoordinates};
use itertools::Itertools;
use polars::error::PolarsResult;
use polars::export::arrow::array::Array;
use polars::export::arrow::record_batch::RecordBatchT;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::btree_map::Cursor;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::error::Error;
use std::io::{Read, Seek};
use std::ops::Bound::Included;
use std::sync::mpsc::{Receiver, Sender, SyncSender};
use std::sync::{mpsc, Arc, Mutex, RwLock};
use std::thread;
use std::thread::JoinHandle;

/// chromosome -> (start_position -> batch_index)
type BSXIndex = HashMap<String, BTreeMap<u64, usize>>;

/// Structure, that holds chromosome-wise index of batches from BSXplorer
/// IPC file.
#[derive(Clone)]
pub struct BSXBTree(BSXIndex);

impl BSXBTree {
    pub(crate) fn from_file<R: Read + Seek>(
        reader: &mut BsxFileReader<R>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut new = Self::empty();

        for (num, df) in {
            (0..reader.reader.blocks_total())
                .map(|i| reader.get_batch(i))
                .enumerate()
        } {
            let df = df.expect("failed to read df")?;
            let first_pos = df.first_position()?;
            new.insert(first_pos, num);
        }
        Ok(new)
    }

    pub(crate) fn empty() -> Self {
        Self(HashMap::new())
    }

    /// Insert node to a [BTreeMap].
    pub(crate) fn insert(&mut self, pos: GenomicPosition<u64>, idx: usize) {
        match self.0.get_mut(pos.chr()) {
            Some(btree) => {
                btree.insert(pos.position(), idx);
            }
            None => {
                let mut btree = BTreeMap::new();
                btree.insert(pos.position(), idx);
                self.0.insert(pos.chr().to_string(), btree);
            }
        };
    }
    /// Get batch, which starts exactly on specified position.
    pub fn get_exact(&self, pos: GenomicPosition<u64>) -> Option<&usize> {
        self.0.get(pos.chr())?.get(&pos.position())
    }
    /// Get cursor to the first batch, which has lower start position.
    pub fn get_lower_bound(&self, chr: String, start: u64) -> Option<Cursor<'_, u64, usize>> {
        Some(self.0.get(&chr)?.lower_bound(Included(&start)))
    }

    /// Get all batch indexes, which contain data about specified region
    pub fn get_region(&self, coordinates: &RegionCoordinates<u64>) -> Option<Vec<usize>> {
        let mut lower_bound =
            self.get_lower_bound(coordinates.chr.to_string(), coordinates.start())?;
        let mut batches = vec![lower_bound.prev()];
        while let Some((start_val, index)) = lower_bound.next() {
            if start_val > &coordinates.end() {
                break;
            } else {
                batches.push(Some((start_val, index)));
            }
        }
        Some(
            batches
                .iter()
                .flatten()
                .map(|(_key, index)| *(*index))
                .collect(),
        )
    }
}

#[derive(Clone)]
pub struct LinkedReadPlan {
    reg_mapping: HashMap<RegionCoordinates<u64>, HashSet<usize>>,
    batch_mapping: BTreeMap<usize, HashSet<RegionCoordinates<u64>>>,
}

impl LinkedReadPlan {
    fn new() -> Self {
        Self {
            reg_mapping: HashMap::new(),
            batch_mapping: BTreeMap::new(),
        }
    }

    pub(crate) fn from_index_and_regions(
        index: BSXBTree,
        regions: &[RegionCoordinates<u64>],
    ) -> Self {
        let new = Self::new();
        let new_mutex = Mutex::new(new);

        regions.par_iter().for_each(|region| {
            let batches = index.get_region(region);
            if let Some(batches) = batches {
                batches
                    .iter()
                    .for_each(|idx| new_mutex.lock().unwrap().insert(region.clone(), *idx))
            }
        });

        new_mutex.into_inner().unwrap()
    }

    pub fn insert(&mut self, region: RegionCoordinates<u64>, index: usize) {
        match self.reg_mapping.get_mut(&region) {
            Some(indices) => {
                indices.insert(index);
            }
            None => {
                self.reg_mapping
                    .insert(region.clone(), HashSet::from_iter(vec![index]));
            }
        }
        match self.batch_mapping.get_mut(&index) {
            Some(regions) => {
                regions.insert(region.clone());
            }
            None => {
                self.batch_mapping
                    .insert(index, HashSet::from_iter(vec![region.clone()]));
            }
        }
    }

    pub fn remove_region(&mut self, region: &RegionCoordinates<u64>) -> Option<HashSet<usize>> {
        if let Some(batches) = self.reg_mapping.remove(region) {
            for batch in batches.iter() {
                self.batch_mapping
                    .get_mut(batch)
                    .map(|regions| regions.remove(region));
            }
            Some(batches)
        } else {
            None
        }
    }

    pub fn remove_batch(&mut self, index: usize) -> Option<HashSet<RegionCoordinates<u64>>> {
        if let Some(regions) = self.batch_mapping.remove(&index) {
            for region in regions.iter() {
                self.reg_mapping
                    .get_mut(region)
                    .map(|indices| indices.remove(&index));
            }
            Some(regions)
        } else {
            None
        }
    }

    pub fn get_region(&self, region: &RegionCoordinates<u64>) -> Option<&HashSet<usize>> {
        self.reg_mapping.get(region)
    }

    pub fn get_batch(&self, index: usize) -> Option<&HashSet<RegionCoordinates<u64>>> {
        self.batch_mapping.get(&index)
    }
}

type BatchSliceCache = Vec<EncodedBsxBatch>;
type AssembleCache = RwLock<HashMap<RegionCoordinates<u64>, BatchSliceCache>>;

struct RegionAssembler {
    read_plan: Arc<RwLock<LinkedReadPlan>>,
    assemble_cache: Arc<AssembleCache>,
    output_queue: VecDeque<(RegionCoordinates<u64>, EncodedBsxBatch)>,
}

impl RegionAssembler {
    pub fn new(read_plan: LinkedReadPlan) -> Self {
        Self {
            read_plan: Arc::new(RwLock::new(read_plan)),
            assemble_cache: Arc::new(RwLock::new(HashMap::new())),
            output_queue: VecDeque::new(),
        }
    }

    pub fn consume_batch(&mut self, data: EncodedBsxBatch, batch_idx: usize) {
        {
            if let Some(regions) = self.read_plan.read().unwrap().get_batch(batch_idx) {
                let data_ref = Arc::new(data);

                regions.par_iter().for_each(|region| {
                    let trimmed = data_ref.trim_region(region).unwrap();
                    let cache = Arc::clone(&self.assemble_cache);
                    
                    let was_cached= cache.read().unwrap().get(region).is_some();
                    if was_cached {
                        cache.write().unwrap().get_mut(region).unwrap().push(trimmed);
                    } else {
                        cache.write().unwrap().insert(region.clone(), vec![trimmed]);
                    }
                });
            }
        }
        
        let removed_batch_regions = self.read_plan.write().unwrap().remove_batch(batch_idx);
        
        if let Some(regions_to_check) = removed_batch_regions {
            for region in regions_to_check.iter() {
                if self
                    .read_plan
                    .read()
                    .unwrap()
                    .get_region(region)
                    .map_or(false, |indices| indices.is_empty())
                {
                    if let Some(region_batches) =
                        self.assemble_cache.write().unwrap().remove(region)
                    {
                        let assembled = region_batches
                            .into_iter()
                            .sorted_by_cached_key(|b| b.first_position().unwrap().position())
                            .reduce(|mut acc, new| {
                                acc.extend(&new).unwrap();
                                acc
                            })
                            .unwrap();
                        self.output_queue.push_front((region.clone(), assembled));
                    }
                }
            }
        }
    }

    fn try_get(&mut self) -> Option<(RegionCoordinates<u64>, EncodedBsxBatch)> {
        self.output_queue.pop_front()
    }

    fn read_plan(&self) -> Arc<RwLock<LinkedReadPlan>> {
        Arc::clone(&self.read_plan)
    }

    fn is_finished(&self) -> bool {
        self.read_plan.read().unwrap().batch_mapping.is_empty()
    }
}

impl From<LinkedReadPlan> for RegionAssembler {
    fn from(value: LinkedReadPlan) -> Self {
        Self::new(value)
    }
}

pub struct RegionReader {
    assembler: RegionAssembler,
    _read_thread: JoinHandle<()>,
    receiver: Receiver<(usize, EncodedBsxBatch)>,
}

impl RegionReader {
    pub fn try_new<R: Read + Seek + Send + 'static>(handle: R, regions: &[RegionCoordinates<u64>]) -> Result<Self, Box<dyn Error>> {
        let mut reader = BsxFileReader::new(handle);

        let read_plan = LinkedReadPlan::from_index_and_regions(
            BSXBTree::from_file(&mut reader)?,
            regions,
        );

        let assembler = RegionAssembler::new(read_plan.clone());
        let read_plan_copy = Arc::new(read_plan.clone());

        const BATCH_SIZE: usize = 10;
        let (sender, receiver) = mpsc::sync_channel(BATCH_SIZE);
        
        let read_thread =
            thread::spawn(move || region_reader_thread(reader, sender, read_plan_copy));

        Ok(Self {
            assembler,
            _read_thread: read_thread,
            receiver,
        })
    }
}

impl Iterator for RegionReader {
    type Item = (RegionCoordinates<u64>, EncodedBsxBatch);

    fn next(&mut self) -> Option<Self::Item> {
        match self.assembler.try_get() {
            Some(region_data) => Some(region_data),
            None => match self.receiver.recv() {
                Ok((index, batch)) => {
                    self.assembler.consume_batch(batch, index);
                    self.next()
                }
                Err(e) => {
                    if self.assembler.is_finished() {
                        None
                    } else {
                        panic!("{}", e)
                    }
                },
            },
        }
    }
}

pub fn df_to_regions(
    data_frame: &DataFrame, 
    chr_colname: Option<&str>, 
    start_colname: Option<&str>, 
    end_colname: Option<&str>
) -> PolarsResult<Vec<RegionCoordinates<u64>>> {
    let chr_colname = chr_colname.unwrap_or("chr");
    let start_colname = start_colname.unwrap_or("start");
    let end_colname = end_colname.unwrap_or("end");
    
    let target_schema = Schema::from_iter([
        (chr_colname.into(), DataType::String),
        (start_colname.into(), DataType::UInt64),
        (end_colname.into(), DataType::UInt64),
    ]);
    
    let mut data_casted = data_frame
        .select([chr_colname, start_colname, end_colname])?
        .select_with_schema(
            [chr_colname, start_colname, end_colname],
            &SchemaRef::new(target_schema),
        )?;
    data_casted = data_casted.drop_nulls(None::<&[String]>)?;
    
    let chr_col = data_casted.column(chr_colname)?.str()?;
    let start_col = data_casted.column(start_colname)?.u64()?;
    let pos_col = data_casted.column(end_colname)?.u64()?;
    
    Ok(
        itertools::izip!(chr_col.into_iter(), start_col.into_iter(), pos_col.into_iter())
            .map(|(chr, start, end)| RegionCoordinates::new(chr.unwrap().to_string(), start.unwrap(), end.unwrap()))
            .collect_vec()
    )
    
    
}

fn region_reader_thread<R: Read + Seek>(
    mut reader: BsxFileReader<R>,
    sender: SyncSender<(usize, EncodedBsxBatch)>,
    read_plan: Arc<LinkedReadPlan>,
) {
    for batch_idx in read_plan.batch_mapping.keys() {
        match reader.get_batch(*batch_idx) {
            Some(Ok(batch)) => {
                sender.send((*batch_idx, batch)).unwrap();
            }
            None => {
                continue;
            }
            Some(Err(e)) => panic!("{:?}", e),
        }
    }
}

pub struct BsxFileReader<R: Read + Seek> {
    reader: IpcFileReader<R>,
}

impl<R: Read + Seek> BsxFileReader<R> {
    pub fn new(handle: R) -> Self {
        Self {
            reader: IpcFileReader::new(handle, None, None),
        }
    }

    fn process_record_batch(
        &self,
        batch: PolarsResult<RecordBatchT<Box<dyn Array>>>,
    ) -> PolarsResult<EncodedBsxBatch> {
        batch
            .map(|batch| {
                DataFrame::try_from((batch, self.reader.metadata().schema.as_ref()))
                    .expect("Failed to create dataFrame from batch")
            })
            .map(|df| unsafe { EncodedBsxBatch::new_unchecked(df) })
    }

    pub fn get_batch(&mut self, batch_idx: usize) -> Option<PolarsResult<EncodedBsxBatch>> {
        self.reader
            .read_at(batch_idx)
            .map(|res| self.process_record_batch(res))
    }
}

impl<R: Read + Seek> Iterator for BsxFileReader<R> {
    type Item = PolarsResult<EncodedBsxBatch>;
    fn next(&mut self) -> Option<Self::Item> {
        let next = self.reader.next();
        next.map(|res| self.process_record_batch(res))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tree() {
        let mut tree = BSXBTree::empty();
        tree.insert(GenomicPosition::new("chr1".to_string(), 1), 0);
        tree.insert(GenomicPosition::new("chr1".to_string(), 100), 1);
        tree.insert(GenomicPosition::new("chr1".to_string(), 200), 2);

        let test = tree.get_region(&RegionCoordinates::new("chr1".to_string(), 1, 100));
        assert_eq!(test, Some(vec![0, 1]));
        let test = tree.get_region(&RegionCoordinates::new("chr1".to_string(), 300, 400));
        // As we compare only start of batch. This result gives us an information
        // that this batch can not be anywhere else, rather than in 2nd batch
        assert_eq!(test, Some(vec![2, 2]));
    }

    #[test]
    fn test_plan() {
        let mut tree = BSXBTree::empty();
        tree.insert(GenomicPosition::new("chr1".to_string(), 1), 0);
        tree.insert(GenomicPosition::new("chr1".to_string(), 100), 1);
        tree.insert(GenomicPosition::new("chr1".to_string(), 200), 2);

        let reg1 = RegionCoordinates::new("chr1".to_string(), 1, 100);
        let reg2 = RegionCoordinates::new("chr1".to_string(), 1, 200);
        let reg3 = RegionCoordinates::new("chr2".to_string(), 2, 400);

        let plan = LinkedReadPlan::from_index_and_regions(
            tree,
            &[reg1.clone(), reg2.clone(), reg3.clone()],
        );

        assert_eq!(
            plan.batch_mapping,
            BTreeMap::from_iter([
                (0, HashSet::from_iter([reg1.clone(), reg2.clone()])),
                (1, HashSet::from_iter([reg1.clone(), reg2.clone()])),
                (2, HashSet::from_iter([reg2.clone()])),
            ])
        );

        assert_eq!(
            plan.reg_mapping,
            HashMap::from_iter([
                (reg1, HashSet::from_iter([0, 1])),
                (reg2, HashSet::from_iter([0, 1, 2]))
            ])
        )
    }
}
