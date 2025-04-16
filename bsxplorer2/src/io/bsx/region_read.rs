use anyhow::anyhow;
use itertools::Itertools;
use polars::error::PolarsResult;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::error::Error;
use std::io::{Read, Seek};
use std::ops::Bound::Included;
use std::ops::Deref;
use std::sync::mpsc::{Receiver, SyncSender};
use std::sync::{mpsc, Arc, RwLock};
use std::thread;
use std::thread::JoinHandle;

use crate::data_structs::batch::BsxBatchMethods;
use crate::data_structs::batch::EncodedBsxBatch;
use crate::data_structs::batch::LazyBsxBatch;
use crate::data_structs::region::{GenomicPosition, RegionCoordinates};
use crate::data_structs::region_data::RegionData;
use crate::io::bsx::read::{BSXIndex, BsxFileReader};
use crate::utils::types::{IPCEncodedEnum, Strand};

/// Structure, that holds chromosome-wise index of batches from BSXplorer
/// IPC file.
#[derive(Clone)]
pub struct BSXBTree(BSXIndex);

impl BSXBTree {
    pub(crate) fn from_file<R: Read + Seek + Send>(
        reader: &mut BsxFileReader<R>
    ) -> Result<Self, Box<dyn Error>> {
        let mut new = Self::empty();

        for (num, df) in {
            (0..reader.blocks_total())
                .map(|i| reader.get_batch(i))
                .enumerate()
        } {
            let df = df.expect("failed to read df")?;
            let first_pos = GenomicPosition::new(
                df.chr_val()?.to_string(),
                df.start_pos()
                    .ok_or(anyhow!("Empty batch"))? as u64,
            );
            new.insert(first_pos, num);
        }
        Ok(new)
    }

    pub(crate) fn empty() -> Self {
        Self(HashMap::new())
    }

    /// Insert node to a [BTreeMap].
    pub(crate) fn insert(
        &mut self,
        pos: GenomicPosition<u64>,
        idx: usize,
    ) {
        match self.0.get_mut(pos.chr()) {
            Some(btree) => {
                btree.insert(pos.position(), idx);
            },
            None => {
                let mut btree = BTreeMap::new();
                btree.insert(pos.position(), idx);
                self.0
                    .insert(pos.chr().to_string(), btree);
            },
        };
    }

    /// Get batch, which starts exactly on specified position.
    pub fn get_exact(
        &self,
        pos: GenomicPosition<u64>,
    ) -> Option<&usize> {
        self.0
            .get(pos.chr())?
            .get(&pos.position())
    }

    /// Get cursor to the first batch, which has lower start position.
    pub fn get_lower_bound(
        &self,
        chr: &str,
        start: u64,
    ) -> Option<std::collections::btree_map::Range<'_, u64, usize>> {
        self.0.get(chr).map(|inner_map| {
            inner_map.range((Included(start), Included(u64::MAX))) // Use range instead of lower_bound
        })
    }

    /// Get all batch indexes, which contain data_structs about specified region
    pub fn get_region(
        &self,
        coordinates: &RegionCoordinates<u64>,
    ) -> Option<Vec<usize>> {
        let range =
            self.get_lower_bound(coordinates.chr(), coordinates.start())?;
        let mut batches: Vec<usize> = Vec::new();

        for (start_val, index) in range {
            if *start_val > coordinates.end() {
                break;
            }
            batches.push(*index);
        }

        if batches.is_empty() {
            return if let Some(inner_map) = self.0.get(coordinates.chr()) {
                if !inner_map.is_empty() {
                    let (_start, idx) = inner_map.last_key_value().unwrap();
                    Some(vec![*idx; 2])
                } else {
                    Some(vec![])
                }
            } else {
                None
            };
        }

        Some(batches)
    }
}

#[derive(Clone)]
pub struct LinkedReadPlan {
    reg_mapping: HashMap<RegionData<String, u64, ()>, HashSet<usize>>,
    batch_mapping: BTreeMap<usize, HashSet<RegionData<String, u64, ()>>>,
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
        regions: &[RegionData<String, u64, ()>],
    ) -> Self {
        let pairs: Vec<(RegionData<String, u64, ()>, usize)> = regions
            .par_iter()
            .flat_map_iter(|region| {
                index
                    .get_region(&region.as_region())
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

    pub fn insert(
        &mut self,
        region: RegionData<String, u64, ()>,
        index: usize,
    ) {
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

    pub fn remove_region(
        &mut self,
        region: &RegionData<String, u64, ()>,
    ) -> Option<HashSet<usize>> {
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

    pub fn remove_batch(
        &mut self,
        index: usize,
    ) -> Option<HashSet<RegionData<String, u64, ()>>> {
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

    pub fn get_region(
        &self,
        region: &RegionData<String, u64, ()>,
    ) -> Option<&HashSet<usize>> {
        self.reg_mapping.get(region)
    }

    pub fn get_batch(
        &self,
        index: usize,
    ) -> Option<&HashSet<RegionData<String, u64, ()>>> {
        self.batch_mapping.get(&index)
    }
}

struct RegionAssembler {
    read_plan: LinkedReadPlan,
    assemble_cache:
        Arc<RwLock<HashMap<RegionData<String, u64, ()>, Vec<EncodedBsxBatch>>>>,
    output_queue:
        Arc<RwLock<VecDeque<RegionData<String, u64, EncodedBsxBatch>>>>,
}

impl RegionAssembler {
    pub fn new(read_plan: LinkedReadPlan) -> Self {
        Self {
            read_plan,
            assemble_cache: Arc::new(RwLock::new(HashMap::new())),
            output_queue: Arc::new(RwLock::new(VecDeque::new())),
        }
    }

    pub fn consume_batch(
        &mut self,
        data: EncodedBsxBatch,
        batch_idx: usize,
    ) {
        let data_ref = Arc::new(data);
        let affected_regions = self.read_plan.remove_batch(batch_idx);

        if let Some(affected_regions) = affected_regions {
            // Collect trimmed data_structs in parallel.
            let trimmed_results: Vec<(
                RegionData<String, u64, ()>,
                EncodedBsxBatch,
            )> = affected_regions
                .par_iter()
                .map(|region_data| {
                    let coordinates = region_data.as_region();
                    let mut lazy = LazyBsxBatch::from(data_ref.deref().clone())
                        .filter_pos_gt(coordinates.start() - 1)
                        .filter_pos_lt(coordinates.end);
                    if !matches!(region_data.strand(), Strand::None) {
                        lazy = lazy.filter_strand(region_data.strand());
                    }
                    EncodedBsxBatch::try_from(lazy)
                        .map(|batch| (region_data.clone(), batch))
                })
                .collect::<anyhow::Result<Vec<_>>>()
                .expect("Failed to trim data_structs");

            // Update the cache using a single write lock.
            {
                let mut cache_write = self.assemble_cache.write().unwrap();
                for (region_data, batch_data) in trimmed_results.iter() {
                    cache_write
                        .entry(region_data.clone())
                        .or_insert_with(Vec::new)
                        .push(batch_data.clone());
                }
            }

            // Check for finished regions (those with no pending batches) and
            // assemble.
            affected_regions
                .par_iter()
                .for_each(|region_data| {
                    let is_finished = self
                        .read_plan
                        .get_region(region_data)
                        .map_or(false, |indices| indices.is_empty());

                    if is_finished {
                        let region_batches = {
                            let mut cache_write =
                                self.assemble_cache.write().unwrap();
                            cache_write.remove(region_data).expect(
                                "Region marked as read, but missing in cache",
                            )
                        };

                        // Sort and then vertically stack batches.
                        let sorted_batches = region_batches
                            .into_iter()
                            .sorted_unstable_by_key(|b| b.start_pos())
                            .collect_vec();

                        let assembled = sorted_batches
                            .into_iter()
                            .reduce(|acc, df| acc.vstack(&df).unwrap())
                            .unwrap();

                        self.output_queue
                            .write()
                            .unwrap()
                            .push_front(
                                region_data.clone().with_data(assembled),
                            );
                    }
                });
        }
    }

    fn try_get(&mut self) -> Option<RegionData<String, u64, EncodedBsxBatch>> {
        self.output_queue
            .write()
            .unwrap()
            .pop_front()
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

pub struct RegionReader {
    assembler: RegionAssembler,
    _read_thread: JoinHandle<()>,
    receiver: Receiver<(usize, EncodedBsxBatch)>,
}

impl RegionReader {
    pub fn try_new<R: Read + Seek + Send + 'static>(
        handle: R,
        regions: &[RegionData<String, u64, ()>],
    ) -> Result<Self, Box<dyn Error>> {
        let mut reader = BsxFileReader::new(handle);

        let read_plan = LinkedReadPlan::from_index_and_regions(
            BSXBTree::from_file(&mut reader)?,
            regions,
        );

        let assembler = RegionAssembler::new(read_plan.clone());
        let read_plan_copy = Arc::new(read_plan.clone());

        const BATCH_SIZE: usize = 10;
        let (sender, receiver) = mpsc::sync_channel(BATCH_SIZE);

        let read_thread = thread::spawn(move || {
            region_reader_thread(reader, sender, read_plan_copy)
        });

        Ok(Self {
            assembler,
            _read_thread: read_thread,
            receiver,
        })
    }
}

impl Iterator for RegionReader {
    type Item = RegionData<String, u64, EncodedBsxBatch>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(region_data) = self.assembler.try_get() {
                return Some(region_data);
            }
            match self.receiver.recv() {
                Ok((index, batch)) => {
                    self.assembler
                        .consume_batch(batch, index);
                    // Continue looping instead of making a recursive call.
                },
                Err(e) => {
                    if self.assembler.is_finished() {
                        return None;
                    } else {
                        panic!("{}", e)
                    }
                },
            }
        }
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
    let start_col = data_casted
        .column(start_colname)?
        .u64()?;
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
    let batch_indices = read_plan
        .batch_mapping
        .keys()
        .cloned()
        .collect_vec();
    for batch_idx in batch_indices {
        match reader.get_batch(batch_idx) {
            Some(Ok(batch)) => {
                sender
                    .send((batch_idx, batch))
                    .expect("Failed to send");
            },
            None => {
                continue;
            },
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
        // As we compare only start of batch. This result gives us an
        // information that this batch can not be anywhere else, rather
        // than in 2nd batch
        assert_eq!(test, Some(vec![2, 2]));
    }
}
