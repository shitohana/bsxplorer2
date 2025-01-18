use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::collections::btree_map::Cursor;
use std::error::Error;
use std::io::{Read, Seek, Write};
use std::ops::Bound::Included;
use std::sync::{Arc, Mutex};
use polars::frame::DataFrame;
use rayon::prelude::*;
use crate::io::bsx::ipc::IpcFileReader;
use crate::io::report::bsx_batch::EncodedBsxBatch;
use crate::io::report::report_batch_utils::first_position;
use crate::region::{GenomicPosition, RegionCoordinates};

/// chromosome -> (start_position -> batch_index)
type BSXIndex = HashMap<String, BTreeMap<u64, usize>>;

/// Structure, that holds chromosome-wise index of batches from BSXplorer
/// IPC file.
#[derive(Clone)]
pub struct BSXBTree(BSXIndex);

impl BSXBTree {
    pub(crate) fn from_file<R: Read + Seek>(reader: &mut IpcFileReader<R>) -> Result<Self, Box<dyn Error>> {
        let mut new = Self::empty();

        for (num, df) in {
            (0..reader.blocks_total())
                .map(|i| reader.read_df_at(i))
                .enumerate()
        } {
            let df = df.expect("failed to read df");
            let first_pos = first_position(&df, "chr", "position")?;
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
        Some(
            self.0
                .get(&chr)?
                .lower_bound(Included(&start)),
        )
    }
    
    /// Get all batch indexes, which contain data about specified region
    pub fn get_region(&self, coordinates: &RegionCoordinates<u64>) -> Option<Vec<usize>> {
        let mut lower_bound = self.get_lower_bound(coordinates.chr.to_string(), coordinates.start())?;
        let mut batches = vec![lower_bound.prev()];
        while let Some((start_val, index)) = lower_bound.next() {
            if start_val > &coordinates.end() {
                break;
            }
            else {
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

struct ReadPlan {
    /// Region --> {batch indexes}
    reg_mapping: HashMap<RegionCoordinates<u64>, Vec<usize>>,
    /// Batch index --> {regions}
    batch_mapping: BTreeMap<usize, HashSet<RegionCoordinates<u64>>>
}

impl ReadPlan {
    pub fn new() -> Self {
        Self {
            reg_mapping: HashMap::new(),
            batch_mapping: BTreeMap::new()
        }
    }
    
    pub(crate) fn from_index_and_regions(index: BSXBTree, regions: &[RegionCoordinates<u64>]) -> Self {
        let new = Self::new();
        let new_mutex = Mutex::new(new);
        
        regions.par_iter().for_each(|region| {
            let batches = index.get_region(region);
            if let Some(batches) = batches {
                new_mutex.lock().unwrap().insert(region.clone(), batches);
            }
        });
        
        new_mutex.into_inner().unwrap()
    }
    
    pub fn insert(&mut self, coords: RegionCoordinates<u64>, indices: Vec<usize>) {
        match self.reg_mapping.get_mut(&coords) {
            Some(batches) => { batches.extend(indices.clone()); },
            None => { self.reg_mapping.insert(coords.clone(), Vec::from_iter(indices.clone())); },
        }
        for index in indices {
            match self.batch_mapping.get_mut(&index) {
                Some(regions) => { regions.insert(coords.clone()); },
                None => { self.batch_mapping.insert(index, HashSet::from_iter([coords.clone()])); }
            }
        }
        
    }
    
    pub(crate) fn get_region(&self, coordinates: &RegionCoordinates<u64>) -> Option<&Vec<usize>> {
        self.reg_mapping.get(coordinates)
    }
    
    pub(crate) fn get_batch(&self, index: usize) -> Option<&HashSet<RegionCoordinates<u64>>> {
        self.batch_mapping.get(&index)
    }

    pub fn reg_mapping(&self) -> &HashMap<RegionCoordinates<u64>, Vec<usize>> {
        &self.reg_mapping
    }

    pub fn batch_mapping(&self) -> &BTreeMap<usize, HashSet<RegionCoordinates<u64>>> {
        &self.batch_mapping
    }
}

type BatchSliceCache = HashMap<usize, DataFrame>;
type AssembleCache = HashMap<RegionCoordinates<u64>, BatchSliceCache>;

struct RegionAssembler {
    pending_regions: HashMap<RegionCoordinates<u64>, Vec<usize>>,
    read_plan: ReadPlan,
    assemble_cache: AssembleCache,
    output_queue: VecDeque<(RegionCoordinates<u64>, DataFrame)>,
}

impl RegionAssembler {
    pub fn new(read_plan: ReadPlan) -> Self {
        Self {
            pending_regions: read_plan.reg_mapping().clone(),
            read_plan,
            assemble_cache: HashMap::new(),
            output_queue: VecDeque::new(),
        }
    }
    
    pub fn consume_batch(&mut self, data: DataFrame, batch_idx: usize) {
        let pending_regions = self.read_plan.get_batch(batch_idx);
        let data_ref = Arc::new(data);
        
        if let Some(pending_regions) = pending_regions {
            pending_regions.par_iter().for_each(|region| {
                todo!()
            })
        }
    }
}

impl From<ReadPlan> for RegionAssembler {
    fn from(value: ReadPlan) -> Self {
        Self::new(value)
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
        
        let mut plan = ReadPlan::from_index_and_regions(
            tree,
            &[reg1.clone(), reg2.clone(), reg3.clone()],
        );
        
        assert_eq!(plan.batch_mapping, BTreeMap::from_iter([
            (0, HashSet::from_iter([reg1.clone(), reg2.clone()])),
            (1, HashSet::from_iter([reg1.clone(), reg2.clone()])),
            (2, HashSet::from_iter([reg2.clone()])),
        ]));
        
        assert_eq!(plan.reg_mapping, HashMap::from_iter([
            (reg1, Vec::from_iter([0, 1])),
            (reg2, Vec::from_iter([0, 1, 2]))
        ]))
    }
}