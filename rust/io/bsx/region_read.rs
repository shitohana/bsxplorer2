use crate::data_structs::bsx_batch::{BsxBatchMethods, EncodedBsxBatch};
use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
use crate::io::bsx::ipc::IpcFileReader;
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
use std::sync::mpsc::{Receiver, SyncSender};
use std::sync::{mpsc, Arc, RwLock};
use std::thread;
use std::thread::JoinHandle;

#[cfg(feature = "python")]
use crate::data_structs::region::PyRegionCoordinates;
use crate::io::bsx::read::{BSXIndex, BsxFileReader};
use crate::utils::types::Strand;
#[cfg(feature = "python")]
use crate::utils::wrap_polars_result;
#[cfg(feature = "python")]
use crate::wrap_box_result;
#[cfg(feature = "python")]
use pyo3::exceptions::PyIOError;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3_polars::error::PyPolarsErr;

/// Structure, that holds chromosome-wise index of batches from BSXplorer
/// IPC file.
#[derive(Clone)]
pub struct BSXBTree(BSXIndex);

impl BSXBTree {
    pub(crate) fn from_file<R: Read + Seek + Send>(
        reader: &mut BsxFileReader<R>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut new = Self::empty();

        for (num, df) in {
            (0..reader.blocks_total())
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
        let pairs: Vec<(RegionCoordinates<u64>, usize)> = regions
            .par_iter()
            .flat_map_iter(|region| {
                index
                    .get_region(region)
                    .into_iter()
                    .flatten()
                    .map(|idx| (region.clone(), idx))
            })
            .collect();

        let mut plan = LinkedReadPlan::new();
        for (region, idx) in pairs {
            plan.insert(region, idx);
        }
        plan
    }

    pub fn insert(&mut self, region: RegionCoordinates<u64>, index: usize) {
        self.reg_mapping
            .entry(region.clone())
            .or_insert_with(|| {
                let mut set = HashSet::with_capacity(1);
                set.insert(index);
                set
            })
            .insert(index);
        self.batch_mapping
            .entry(index)
            .or_insert_with(|| {
                let mut set = HashSet::with_capacity(1);
                set.insert(region.clone());
                set
            })
            .insert(region);
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

struct RegionAssembler {
    read_plan: LinkedReadPlan,
    assemble_cache: Arc<RwLock<HashMap<RegionCoordinates<u64>, Vec<EncodedBsxBatch>>>>,
    output_queue: Arc<RwLock<VecDeque<(RegionCoordinates<u64>, EncodedBsxBatch)>>>,
}

impl RegionAssembler {
    pub fn new(read_plan: LinkedReadPlan) -> Self {
        Self {
            read_plan,
            assemble_cache: Arc::new(RwLock::new(HashMap::new())),
            output_queue: Arc::new(RwLock::new(VecDeque::new())),
        }
    }

    pub fn consume_batch(&mut self, data: EncodedBsxBatch, batch_idx: usize) {
        let data_ref = Arc::new(data);
        let affected_regions = self.read_plan.remove_batch(batch_idx);

        if let Some(affected_regions) = affected_regions {
            // Collect trimmed data in parallel.
            let trimmed_results: Vec<(RegionCoordinates<u64>, EncodedBsxBatch)> = affected_regions
                .par_iter()
                .map(|region_coordinates| {
                    let mut batch = data_ref.trim_region(region_coordinates);
                    batch.map(|mut df| {
                        if matches!(region_coordinates.strand, Strand::None) {
                            df = df.filter(None, Some(region_coordinates.strand))
                        };
                        (region_coordinates.clone(), df)
                    })
                })
                .collect::<PolarsResult<Vec<_>>>()
                .expect("Failed to trim data");

            // Update the cache using a single write lock.
            {
                let mut cache_write = self.assemble_cache.write().unwrap();
                for (region_coordinates, region_data) in trimmed_results.iter() {
                    cache_write
                        .entry(region_coordinates.clone())
                        .or_insert_with(Vec::new)
                        .push(region_data.clone());
                }
            }

            // Check for finished regions (those with no pending batches) and assemble.
            affected_regions.par_iter().for_each(|region| {
                let is_finished = self
                    .read_plan
                    .get_region(region)
                    .map_or(false, |indices| indices.is_empty());

                if is_finished {
                    let region_batches = {
                        let mut cache_write = self.assemble_cache.write().unwrap();
                        cache_write
                            .remove(region)
                            .expect("Region marked as read, but missing in cache")
                    };

                    // Sort and then vertically stack batches.
                    let mut sorted_batches = region_batches
                        .into_iter()
                        .sorted_unstable_by_key(|b| b.first_position().unwrap().position())
                        .collect_vec();

                    let assembled = sorted_batches
                        .into_iter()
                        .reduce(|acc, df| acc.vstack(df).unwrap())
                        .unwrap();

                    self.output_queue
                        .write()
                        .unwrap()
                        .push_front((region.clone(), assembled));
                }
            });
        }
    }

    fn try_get(&mut self) -> Option<(RegionCoordinates<u64>, EncodedBsxBatch)> {
        self.output_queue.write().unwrap().pop_front()
    }

    #[allow(dead_code)]
    fn read_plan(&self) -> &LinkedReadPlan {
        &self.read_plan
    }

    fn is_finished(&self) -> bool {
        self.read_plan.batch_mapping.is_empty()
    }
}

impl From<LinkedReadPlan> for RegionAssembler {
    fn from(value: LinkedReadPlan) -> Self {
        Self::new(value)
    }
}

#[cfg_attr(feature = "python", pyclass)]
pub struct RegionReader {
    assembler: RegionAssembler,
    _read_thread: JoinHandle<()>,
    receiver: Receiver<(usize, EncodedBsxBatch)>,
}

impl RegionReader {
    pub fn try_new<R: Read + Seek + Send + 'static>(
        handle: R,
        regions: &[RegionCoordinates<u64>],
    ) -> Result<Self, Box<dyn Error>> {
        let mut reader = BsxFileReader::new(handle);

        let read_plan =
            LinkedReadPlan::from_index_and_regions(BSXBTree::from_file(&mut reader)?, regions);

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
        loop {
            if let Some(region_data) = self.assembler.try_get() {
                return Some(region_data);
            }
            match self.receiver.recv() {
                Ok((index, batch)) => {
                    self.assembler.consume_batch(batch, index);
                    // Continue looping instead of making a recursive call.
                }
                Err(e) => {
                    if self.assembler.is_finished() {
                        return None;
                    } else {
                        panic!("{}", e)
                    }
                }
            }
        }
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl RegionReader {
    #[new]
    fn new_py(path: String, regions: Vec<PyRegionCoordinates>) -> PyResult<Self> {
        let handle = File::open(path)?;
        let regions = regions
            .into_iter()
            .map(RegionCoordinates::from)
            .collect_vec();
        wrap_box_result!(PyIOError, Self::try_new(handle, &regions))
    }

    fn __next__(&mut self) -> Option<(PyRegionCoordinates, EncodedBsxBatch)> {
        self.next().map(|(index, batch)| (index.into(), batch))
    }
}

pub fn df_to_regions(
    data_frame: &DataFrame,
    chr_colname: Option<&str>,
    start_colname: Option<&str>,
    end_colname: Option<&str>,
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

    Ok(itertools::izip!(
        chr_col.into_iter(),
        start_col.into_iter(),
        pos_col.into_iter()
    )
    .map(|(chr, start, end)| {
        RegionCoordinates::new(
            chr.unwrap().to_string(),
            start.unwrap(),
            end.unwrap(),
            Strand::None,
        )
    })
    .collect_vec())
}

fn region_reader_thread<R: Read + Seek>(
    mut reader: BsxFileReader<R>,
    sender: SyncSender<(usize, EncodedBsxBatch)>,
    read_plan: Arc<LinkedReadPlan>,
) {
    let batch_indices = read_plan.batch_mapping.keys().cloned().collect_vec();
    for batch_idx in batch_indices {
        match reader.get_batch(batch_idx) {
            Some(Ok(batch)) => {
                sender.send((batch_idx, batch)).unwrap();
            }
            None => {
                continue;
            }
            Some(Err(e)) => panic!("{:?}", e),
        }
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

        let test = tree.get_region(&RegionCoordinates::new(
            "chr1".to_string(),
            1,
            100,
            Strand::None,
        ));
        assert_eq!(test, Some(vec![0, 1]));
        let test = tree.get_region(&RegionCoordinates::new(
            "chr1".to_string(),
            300,
            400,
            Strand::None,
        ));
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

        let reg1 = RegionCoordinates::new("chr1".to_string(), 1, 100, Strand::None);
        let reg2 = RegionCoordinates::new("chr1".to_string(), 1, 200, Strand::None);
        let reg3 = RegionCoordinates::new("chr2".to_string(), 2, 400, Strand::None);

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
