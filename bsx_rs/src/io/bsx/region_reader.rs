use crate::io::bsx::bsx_reader::BSXReader;
use crate::region::{RegionCoordinates, RegionData};
use crate::ubatch::UniversalBatch;
use crate::utils::types::{Context, Strand};
use itertools::Itertools;
use log::{debug, info};
use polars::datatypes::{DataType, PlHashMap};
use polars::error::PolarsError;
use polars::frame::DataFrame;
use polars::prelude::IntoLazy;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use std::collections::{HashMap, HashSet, VecDeque};

/// Configuration for [RegionsDataReader]
#[derive(Clone)]
pub struct ReadFilters {
    /// Which [Context] to filter
    pub(crate) context: Option<Context>,
    /// Which [Strand] to filter
    pub(crate) strand: Option<Strand>,
}

impl ReadFilters {
    pub fn new(context: Option<Context>, strand: Option<Strand>) -> Self {
        Self { context, strand }
    }

    /// Constructor for no filtering
    pub fn empty() -> Self {
        Self::new(None, None)
    }
}

/// High-level interface for reading BSXplorer IPC file, which
/// iterates over regions
pub struct RegionsDataReader {
    reader: BSXReader,
    /// [DataFrame] with obligatory columns: "chr", "start", "end"
    regions: DataFrame,
}

type ReadPlan = HashMap<usize, HashSet<usize>>;
type RegionMapping = HashMap<usize, RegionCoordinates>;

impl RegionsDataReader {
    pub fn new(reader: BSXReader, regions: DataFrame) -> Self {
        let regions_cols = regions.get_column_names_str();
        if !["chr", "start", "end"]
            .iter()
            .map(|name| regions_cols.contains(&name))
            .all(|exists| exists)
        {
            panic!("Obligatory columns ('chr', 'start', 'end') not present!");
        }
        let regions = regions
            .lazy()
            .cast(
                PlHashMap::from_iter([
                    ("chr", DataType::String),
                    ("start", DataType::UInt32),
                    ("end", DataType::UInt32),
                ]),
                true,
            )
            .collect()
            .unwrap();
        info!(
            "RegionsDataReader created. Reading {} regions.",
            regions.height()
        );
        Self { reader, regions }
    }

    fn make_plan(&mut self) -> Result<(ReadPlan, RegionMapping), PolarsError> {
        debug!("Making plan");
        // Create mapping
        let chr_col = self.regions.column("chr")?.str()?;
        let start_col = self.regions.column("start")?.u32()?;
        let end_col = self.regions.column("end")?.u32()?;
        let gene_idx_mapping: HashMap<usize, RegionCoordinates> = HashMap::from_iter(
            itertools::izip!(
                chr_col.into_iter(),
                start_col.into_iter(),
                end_col.into_iter()
            )
            .enumerate()
            .filter_map(|(idx, values)| match values {
                (Some(chr), Some(start), Some(end)) => {
                    Some((idx, RegionCoordinates::new(chr.to_owned(), start, end)))
                }
                _ => None,
            }),
        );

        // Hashmap of batch idx -> gene_idx
        let mut batch_hashmap: ReadPlan = HashMap::new();
        'gene_iter: for (gene_idx, coordinates) in gene_idx_mapping.iter() {
            let indices = self.reader.index.get_region(
                coordinates.chr.as_str(),
                coordinates.start,
                coordinates.end,
            );

            if indices.is_none() {
                continue 'gene_iter;
            }

            for batch_idx in indices.unwrap() {
                match batch_hashmap.get_mut(&batch_idx) {
                    Some(regions) => {
                        regions.insert(gene_idx.clone());
                    }
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
    plan_iter: Box<dyn Iterator<Item = usize>>,
    /// Queue of assembled regions to return from
    assembled_queue: VecDeque<RegionData>,
    /// Mapping between ids of regions, which data is not yet completely read
    /// and the data itself
    incomplete: HashMap<usize, Vec<DataFrame>>,
    /// How many RecordBatches to read at once
    batch_per_read: usize,
    read_filters: ReadFilters,
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
            plan.clone()
                .into_iter()
                .sorted_by_key(|(batch_idx, _reg_idx)| batch_idx.clone())
                .map(|(batch_idx, _reg_idx)| batch_idx)
                .collect::<Vec<_>>()
                .into_iter(),
        );
        let mut region_queue: HashMap<usize, HashSet<usize>> = HashMap::new();
        for (batch_idx, region_idxs) in plan.clone().into_iter() {
            for idx in region_idxs {
                match region_queue.get_mut(&idx) {
                    Some(mapped_batches) => {
                        mapped_batches.insert(batch_idx);
                    }
                    None => {
                        region_queue.insert(idx, HashSet::from_iter([batch_idx.clone()]));
                    }
                }
            }
        }
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
            batch_per_read: batch_per_read
                .unwrap_or(std::thread::available_parallelism().unwrap().into()),
            read_filters,
        }
    }

    pub fn with_batch_per_read(mut self, value: usize) -> Self {
        self.batch_per_read = value;
        self
    }

    pub fn with_filters(mut self, filters: ReadFilters) -> Self {
        self.read_filters = filters;
        self
    }

    /// Trim batches and assign them to regions
    ///
    /// batches: Vec<(batch_idx, [DataFrame])>
    fn process_batches(&mut self, batches: Vec<(usize, UniversalBatch)>) {
        // Convert incoming vector into raw batch ids
        let batch_idxs = batches
            .iter()
            .map(|(idx, _)| idx.clone())
            .collect::<Vec<_>>();
        // Determine, which regions are waiting for this batch data
        let pending_regions = batch_idxs
            .iter()
            .map(|idx| self.plan.get(idx).unwrap().clone())
            .collect::<Vec<_>>();
        // Remove batch_idx from remaining set for each region
        for (batch_idx, region_indexes) in
            itertools::izip!(batch_idxs.into_iter(), pending_regions.clone().into_iter())
        {
            for region_id in region_indexes {
                self.regions_remaining
                    .get_mut(&region_id)
                    .unwrap()
                    .remove(&batch_idx);
            }
        }

        // Get trimmed region coordinates data
        let regions_set: HashSet<usize> = HashSet::from_iter(
            pending_regions
                .clone()
                .into_iter()
                .map(|set| set.into_iter().collect::<Vec<_>>())
                .flatten(),
        );
        let region_coordinates: HashMap<usize, RegionCoordinates> = HashMap::from_iter(
            regions_set
                .iter()
                .map(|reg_id| (*reg_id, self.region_mapping.get(reg_id).unwrap().clone())),
        );
        let filters_copy = self.read_filters.clone();

        // Trim dataframes
        let res = batches
            .into_par_iter()
            .zip(pending_regions.into_par_iter())
            // Filter batches if needed
            .map(|((_batch_idx, batch), _region_indexes)| {
                ((_batch_idx, batch.filter(&filters_copy)), _region_indexes)
            })
            // Extract region DataFrames
            .map(|((_, batch), region_indexes)| {
                let mut res = Vec::new();
                for reg_index in region_indexes {
                    let coordinates = region_coordinates.get(&reg_index).unwrap();
                    let region_data = batch.slice(coordinates).data;
                    res.push((reg_index, region_data))
                }
                res
            })
            .collect::<Vec<_>>();

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
                    None => {
                        unreachable!("Could not remove region")
                    }
                };
                let reg_data_len = region_data.len();

                match RegionData::from_parts(
                    region_data,
                    self.region_mapping.get(&region_id).unwrap().clone(),
                ) {
                    Some(data) => {
                        debug!(
                            "Assembled region {}:{}-{} from {} parts",
                            data.get_coordinates().chr,
                            data.get_coordinates().start,
                            data.get_coordinates().end,
                            reg_data_len
                        );
                        self.assembled_queue.push_back(data)
                    }
                    None => {
                        debug!(
                            "Could not assemble region {:?}. Empty DataFrame",
                            self.region_mapping.get(&region_id).unwrap()
                        );
                    }
                }
            }
        }
    }
}

impl Iterator for RegionsDataIterator {
    type Item = RegionData;

    fn next(&mut self) -> Option<Self::Item> {
        while self.assembled_queue.is_empty() {
            let required_batches = std::iter::repeat(())
                .take(self.batch_per_read)
                .filter_map(|_| {
                    let batch_idx = match self.plan_iter.next() {
                        Some(batch_idx) => batch_idx,
                        None => return None,
                    };
                    Some((batch_idx, self.reader.nth(batch_idx).unwrap()))
                })
                .collect::<Vec<_>>();

            debug!("Read {:?} batches", required_batches.len());
            if required_batches.len() == 0 {
                return None;
            }

            self.process_batches(required_batches);
        }
        let data = self.assembled_queue.pop_front().unwrap();
        Some(data)
    }
}
