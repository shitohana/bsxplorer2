use std::collections::{BTreeMap, HashMap, VecDeque, HashSet};
use std::collections::btree_map::Cursor;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::ops::Bound::Included;
use std::ops::Not;
use itertools::Itertools;
use log::{debug, info};
use polars::prelude::*;
use polars_arrow::io::ipc::read::FileMetadata;
use rayon::prelude::*;

use crate::io::bsx::ipc::IpcFileReader;
use crate::region::RegionData;
use crate::utils::types::{Context, Strand};

type BSXIndex = HashMap<String, BTreeMap<u32, usize>>;
#[derive(Clone)]
pub struct BSXBTree(BSXIndex);

impl BSXBTree {
    pub(crate) fn new() -> Self {
        Self(HashMap::new())
    }
    pub(crate) fn insert(&mut self, chr: String, start: u32, idx: usize) {
        match self.0.get_mut(&chr) {
            Some(btree) => {btree.insert(start, idx);},
            None => {
                let mut btree = BTreeMap::new();
                btree.insert(start, idx);
                self.0.insert(chr, btree);
            }
        };
    }
    pub fn get_exact(&self, chr: String, start: u32) -> Option<&usize> {
        self.0.get(&chr)?.get(&start)
    }
    pub fn get_lower_bound(&self, chr: String, start: u32) -> Option<Cursor<'_, u32, usize>> {
        Some(self.0.get(&chr)?.lower_bound(Included(&start)))
    }
    pub fn get_region(&self, chr: &str, start: u32, end: u32) -> Option<Vec<usize>> {
        let mut lower_bound = self.get_lower_bound(chr.to_string(), start)?;
        let mut batches = vec![lower_bound.prev()];
        while let Some((start_val, index)) = lower_bound.next() {
            if start_val > &end {break}
            else {
                batches.push(Some((start_val, index)));
            }
        };
        Some(batches.iter().flatten().map(|(_key, index)| (*index).clone()).collect())
    }
}


pub struct BSXReader {
    reader: IpcFileReader<File>,
    metadata: FileMetadata,
    pub index: BSXBTree,
}

impl BSXReader {
    pub fn new(mut handle: File, projection: Option<Vec<usize>>) -> Self {
        let temp_handle = handle.try_clone().unwrap();
        let mut temp_reader = IpcFileReader::new(temp_handle, projection.clone(), None);;
        info!("Opened BSX file");
        let mut index = BSXBTree::new();

        for (num, df) in {
            (0..temp_reader.blocks_total())
                .map(|i| temp_reader.read_df_at(i))
                .enumerate()
        } {
            let df = df.expect("failed to read df");
            
            let start = df.column("position").unwrap()
                .u32().unwrap()
                .first().unwrap();
            
            let chr_col = df.column("chr").unwrap()
                .categorical().unwrap();
            let chr_mapping = chr_col.get_rev_map();
            
            let chr = chr_mapping.get(chr_col.physical().get(0).unwrap()).to_owned();
            
            index.insert(chr, start, num);
        }
        info!("Indexed BSX file with {} batches", temp_reader.blocks_total());
        handle.seek(SeekFrom::Start(0)).unwrap();
        let reader = IpcFileReader::new(handle, projection, None);
        let metadata = reader.metadata().clone();

        Self {reader, metadata, index}
    }
}

impl Iterator for BSXReader {
    type Item = DataFrame;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next() {
            Some(record_batch) => {
                let record_batch = record_batch.expect("Could not create record batch");
                let df = DataFrame::try_from((record_batch, self.metadata.schema.as_ref()))
                    .expect("Could not convert record batch into DataFrame");
                Some(df)
            },
            None => None,
        }
    }

    fn count(self) -> usize
    where
        Self: Sized
    {
        self.reader.blocks_total()
    }

    fn last(mut self) -> Option<Self::Item>
    where
        Self: Sized
    {
        self.nth(self.reader.blocks_total() - 1)
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        Some(self.reader.read_df_at(n).expect("Could not read df"))
    }
}

#[derive(Clone, Hash, PartialEq, Eq, Debug)]
pub struct RegionCoordinates {
    chr: String,
    start: u32,
    end: u32,
}

impl RegionCoordinates {
    pub fn chr(&self) -> &str { self.chr.as_str() }
    pub fn start(&self) -> u32 { self.start }
    pub fn end(&self) -> u32 { self.end }
    pub fn length(&self) -> u32 { self.end - self.start }
}

#[derive(Clone)]
pub struct ReadFilters {
    context: Option<Context>,
    strand: Option<Strand>,
}

impl ReadFilters {
    pub fn new(context: Option<Context>, strand: Option<Strand>) -> Self {
        Self { context, strand }
    }
    
    pub fn empty() -> Self {
        Self::new(None, None)
    }
}

pub struct RegionsDataReader {
    reader: BSXReader,
    regions: DataFrame,
}

type ReadPlan = HashMap<usize, HashSet<usize>>;
type RegionMapping = HashMap<usize, RegionCoordinates>;

impl RegionsDataReader {
    pub fn new(reader: BSXReader, regions: DataFrame) -> Self {
        let regions_cols = regions.get_column_names_str();
        if !["chr", "start", "end"].iter().map(|name| regions_cols.contains(&name)).all(|exists| exists) {
            panic!("Obligatory columns ('chr', 'start', 'end') not present!");
        }
        let regions = regions.lazy().cast(PlHashMap::from_iter([
            ("chr", DataType::String),
            ("start", DataType::UInt32),
            ("end", DataType::UInt32)
        ]), true).collect().unwrap();
        info!("RegionsDataReader created. Reading {} regions.", regions.height());
        Self {reader, regions}
    }
    
    fn make_plan(&mut self) -> Result<(ReadPlan, RegionMapping), PolarsError> {
        debug!("Making plan");
        // Create mapping
        let chr_col = self.regions.column("chr")?.str()?;
        let start_col = self.regions.column("start")?.u32()?;
        let end_col = self.regions.column("end")?.u32()?;
        let gene_idx_mapping: HashMap<usize, RegionCoordinates> = HashMap::from_iter(
            itertools::izip!(chr_col.into_iter(), start_col.into_iter(), end_col.into_iter())
                .enumerate()
                .filter_map(|(idx, values)| match values {
                    (Some(chr), Some(start), Some(end)) => Some((
                        idx,
                        RegionCoordinates {chr: chr.to_string(), start, end}
                    )),
                    _ => None
                })
        );
        
        // Hashmap of batch idx -> gene_idx
        let mut batch_hashmap: ReadPlan = HashMap::new();
        'gene_iter: for (gene_idx, coordinates) in gene_idx_mapping.iter() {
            let indices = self.reader.index.get_region(
                coordinates.chr.as_str(), coordinates.start, coordinates.end
            );
            
            if indices.is_none() { continue 'gene_iter }
            
            for batch_idx in indices.unwrap() {
                match batch_hashmap.get_mut(&batch_idx) {
                    Some(regions) => {
                        regions.insert(gene_idx.clone());
                    },
                    None => {
                        let regions = HashSet::from_iter([gene_idx.clone()]);
                        batch_hashmap.insert(batch_idx, regions);
                    }
                }
            }
        }
        info!("Plan created");
        Ok((batch_hashmap, gene_idx_mapping))
    }
}

impl IntoIterator for RegionsDataReader {
    type Item = RegionData;
    type IntoIter = RegionsDataIterator;

    fn into_iter<'a>(mut self) -> Self::IntoIter {
        let (plan, mapping) = self.make_plan().unwrap();
        RegionsDataIterator::new(self.reader, plan, mapping, None, ReadFilters::empty())
    }
}

pub struct RegionsDataIterator {
    /// Instance of [BSXReader] to read batches with
    reader: BSXReader,
    /// Read plan - mapping of IPC File batch indexes, and regions,
    /// which data they contain
    plan: ReadPlan,
    /// Mapping between region id and [RegionCoordinates]
    region_mapping: HashMap<usize, RegionCoordinates>,
    /// Mapping between region id and set of batches, which contain its
    /// data. Region is considered read when HashSet.len() == 0;
    regions_remaining: HashMap<usize, HashSet<usize>>,
    /// Iterator for plan
    plan_iter: Box<dyn Iterator<Item=usize>>,
    /// Queue of assembled regions to return from
    assembled_queue: VecDeque<(usize, DataFrame)>,
    /// Mapping between ids of regions, which data is not yet completely read
    /// and the data itself
    incomplete: HashMap<usize, Vec<DataFrame>>,
    /// How many RecordBatches to read at once
    batch_per_read: usize,
    /// TODO
    read_filters: ReadFilters 
}

impl RegionsDataIterator {
    pub fn new(
        reader: BSXReader,
        plan: ReadPlan,
        region_mapping: HashMap<usize, RegionCoordinates>,
        batch_per_read: Option<usize>,
        read_filters: ReadFilters,
    ) -> Self {
        info!("RegionsDataIterator created");
        let plan_iter = Box::new(
            plan.clone().into_iter()
                .sorted_by_key(|(batch_idx, _reg_idx)| {batch_idx.clone()})
                .map(|(batch_idx, _reg_idx)| (batch_idx))
                .collect::<Vec<_>>()
                .into_iter()
        );
        let mut region_queue: HashMap<usize, HashSet<usize>> = HashMap::new();
        for (batch_idx, region_idxs) in plan.clone().into_iter() {
            for idx in region_idxs {
                match region_queue.get_mut(&idx) {
                    Some(mapped_batches) => { mapped_batches.insert(batch_idx); },
                    None => {
                        region_queue.insert(idx, HashSet::from_iter([batch_idx.clone()]));
                    }
                }
            }
        };
        let incomplete = HashMap::new();
        let assembled_queue = VecDeque::new();
        Self {
            reader, 
            region_mapping,
            regions_remaining: region_queue,
            plan_iter, 
            plan,
            incomplete,
            assembled_queue, 
            batch_per_read: batch_per_read.unwrap_or(std::thread::available_parallelism().unwrap().into()),
            read_filters,
        }
    }

    pub fn with_batch_per_read(mut self, value: usize) -> Self {
        self.batch_per_read = value;
        self
    }
    
    pub fn with_filters (mut self, filters: ReadFilters) -> Self {
        self.read_filters = filters;
        self
    
    }
    
    /// Trim batches and assign them to regions
    /// 
    /// batches: Vec<(batch_idx, [DataFrame])>
    fn process_batches(&mut self, batches: Vec<(usize, DataFrame)>) {
        // Convert incoming vector into raw batch ids
        let batch_idxs = batches.iter()
            .map(|(idx, _)| idx.clone())
            .collect::<Vec<_>>();
        // Determine, which regions are waiting for this batch data
        let pending_regions = batch_idxs.iter()
            .map(|idx| self.plan.get(idx).unwrap().clone())
            .collect::<Vec<_>>();
        // Remove batch_idx from remaining set for each region
        for (batch_idx, region_indexes) in itertools::izip!(batch_idxs.into_iter(), pending_regions.clone().into_iter()) {
            for region_id in region_indexes {
                self.regions_remaining.get_mut(&region_id).unwrap().remove(&batch_idx);
            }
        }
        
        // Get trimmed region coordinates data
        let regions_set: HashSet<usize> = HashSet::from_iter(
            pending_regions.clone().into_iter()
                .map(|set| set.into_iter().collect::<Vec<_>>())
                .flatten()
        );
        let region_coordinates: HashMap<usize, RegionCoordinates> = HashMap::from_iter(regions_set.iter()
            .map(|reg_id| (*reg_id, self.region_mapping.get(reg_id).unwrap().clone()))
        );
        let filters_copy = self.read_filters.clone();
        
        // Trim dataframes
        let res = batches.into_par_iter().zip(pending_regions.into_par_iter())
            // Filter batches if needed
            .map(|((_batch_idx, mut df), _region_indexes)| {
                'strand_filter: {
                    if let Some(strand_filter) = filters_copy.strand {
                        let condition = match strand_filter {
                            Strand::Forward => df.column("strand").unwrap().bool().unwrap().to_owned(),
                            Strand::Reverse => df.column("strand").unwrap().bool().unwrap().not(),
                            _ => break 'strand_filter
                        };
                        df = df.filter(&condition).unwrap();
                    }
                }

                'context_filter: {
                    if let Some(context_filter) = filters_copy.context {
                        let condition = match context_filter {
                            Context::CG => df.column("context").unwrap().bool().unwrap().to_owned(),
                            Context::CHG => df.column("context").unwrap().bool().unwrap().not(),
                            Context::CHH => df.column("context").unwrap().bool().unwrap().is_null(),
                            Context::ALL => break 'context_filter
                        };
                        df = df.filter(&condition).unwrap();
                    }
                }
                ((_batch_idx, df), _region_indexes)
            })
            // Extract region DataFrames
            .map(
                |((_, df), region_indexes)| {
                    let mut res = Vec::new();
                    for reg_index in region_indexes {
                        let coordinates = region_coordinates.get(&reg_index).unwrap();
    
                        let positions = df.column("position").unwrap().u32().unwrap();
                        let start = positions.iter().position(|x| { coordinates.start < x.unwrap() }).unwrap_or_else(|| 0);
                        let slice_length = positions.iter().skip(start).position(|x| coordinates.end < x.unwrap()).unwrap_or(positions.len() - start - 1);
    
                        let mut df_slice = df.slice(start as i64, slice_length);
    
                        res.push((
                            reg_index,
                            df_slice
                        ))
                    }
                    res
                }).collect::<Vec<_>>();
        
        // Append processed data to self.incomplete
        for (reg_idx, reg_df) in res.into_iter().flatten() {
            match self.incomplete.get_mut(&reg_idx) {
                Some(vec) => vec.push(reg_df),
                None => {
                    self.incomplete.insert(reg_idx, vec![reg_df]);
                }
            }
        }
        // Check if any of affected regions are all completed and can
        // be outputted
        for region_id in regions_set {
            if self.regions_remaining.get(&region_id).unwrap().is_empty() {
                let region_data = match self.incomplete.remove(&region_id) {
                    Some(vec) => vec,
                    None => {unreachable!("Could not remove region")}
                };
                let reg_data_len = region_data.len();
                match self.assemble_region(region_data.clone()) {
                    Some(data) => { 
                        let reg_data = self.region_mapping.get(&region_id).unwrap();
                        debug!("Assembled region {}:{}-{} from {} parts", reg_data.chr, reg_data.start, reg_data.end, reg_data_len);
                        self.assembled_queue.push_back((region_id, data)) 
                    },
                    None => {}
                }
            }
        }
    }

    fn assemble_region(&mut self, region_data: Vec<DataFrame>) -> Option<DataFrame> {
        let mut sorted_iter = region_data.iter()
            .filter(|df| df.height() > 0)
            .sorted_by_key(|df| df.column("position").unwrap().u32().unwrap().first().unwrap());
        let mut res = match sorted_iter.next() {
            Some(df) => df.clone(),
            None => return None
        };
        for df in sorted_iter {
            res.extend(df).unwrap();
        }
        res.align_chunks_par();
        Some(res)
    }
}

impl Iterator for RegionsDataIterator {
    type Item = RegionData;
    
    fn next(&mut self) -> Option<Self::Item> {
        while self.assembled_queue.is_empty() {
            let required_batches = std::iter::repeat(())
                .take(self.batch_per_read)
                .map(|_| {
                    let batch_idx = match self.plan_iter.next() {
                        Some(batch_idx) => batch_idx,
                        None => return None
                    };
                    Some((batch_idx, self.reader.nth(batch_idx).unwrap()))
                })
                .flatten()
                .collect::<Vec<_>>();
            debug!("Read {:?} batches", required_batches.len());
            if required_batches.is_empty() {
                return None;
            }
            self.process_batches(required_batches);
        };
        let (reg_id, data) = self.assembled_queue.pop_front().unwrap();
        Some(RegionData::new(data, self.region_mapping.get(&reg_id).unwrap().clone()))
    }
}