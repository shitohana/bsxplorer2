use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::io::bsx::multiple_reader::MultiBsxFileReader;
use crate::io::bsx::read::BsxFileReader;
use crate::tools::dmr::meth_region::{MethylatedRegionOwned, MethylatedRegionView};
use crate::tools::dmr::penalty_segment::PenaltySegmentModel;
use crate::tools::dmr::{utils, DMRegion, DmrModel, FilterConfig, ReaderMetadata};
use crate::utils::types::Context;
use anyhow::anyhow;
use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;
use log::error;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use statrs::statistics::Statistics;
use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::io::{Read, Seek};
use std::sync::mpsc::Receiver;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct DmrConfig {
    pub context: Context,
    pub n_missing: usize,
    pub min_coverage: i16,
    pub diff_threshold: f64,
    pub min_cpgs: usize,
    pub segment_model: PenaltySegmentModel,
}

impl DmrConfig {
    pub fn new(
        context: Context,
        n_missing: usize,
        min_coverage: i16,
        diff_threshold: f64,
        min_cpgs: usize,
        segment_model: PenaltySegmentModel,
    ) -> Self {
        Self {
            context,
            n_missing,
            min_coverage,
            diff_threshold,
            min_cpgs,
            segment_model,
        }
    }

    pub fn filter_config(&self) -> FilterConfig {
        FilterConfig {
            context: self.context,
            n_missing: self.n_missing,
            min_coverage: self.min_coverage,
        }
    }

    pub fn try_finish<F, R>(&self, readers: Vec<(R, F)>) -> anyhow::Result<DmrIterator<R>>
    where
        F: Read + Seek + Send + Sync + 'static,
        R: Display + Eq + Hash + Clone + Default + std::fmt::Debug + Send + 'static,
    {
        let sample_mapping: HashMap<uuid::Uuid, (R, F)> = HashMap::from_iter(
            readers
                .into_iter()
                .map(|(group, handle)| (uuid::Uuid::new_v4(), (group, handle))),
        );

        let group_mapping: HashMap<uuid::Uuid, R> = HashMap::from_iter(
            sample_mapping
                .iter()
                .map(|(id, (group, _))| (id.clone(), group.clone())),
        );
        let group_order = {
            let group_set: HashSet<R> =
                HashSet::from_iter(group_mapping.values().map(|x| x.clone()));
            if group_set.len() != 2 {
                return Err(anyhow!(
                    "There should be only two groups of samples! ({:?})",
                    group_set
                ));
            }
            let groups_vec = group_set.into_iter().collect_vec();
            (groups_vec[0].clone(), groups_vec[1].clone())
        };
        let readers_mapping: HashMap<uuid::Uuid, BsxFileReader<F>> = HashMap::from_iter(
            sample_mapping
                .into_iter()
                .map(|(id, (_, handle))| (id, BsxFileReader::new(handle))),
        );

        let config_copy = self.clone();

        let (sender, receiver) = std::sync::mpsc::sync_channel(10);
        let last_chr = Arc::new(String::new());
        let ref_idset = RefIDSet::new();

        let group_mapping_clone = group_mapping.clone();
        let mut multi_reader =
            MultiBsxFileReader::try_new(readers_mapping).map_err(|e| anyhow!("{:?}", e))?;
        let reader_stat = ReaderMetadata::new(multi_reader.blocks_total());

        let join_handle = std::thread::spawn(move || {
            let group_mapping = group_mapping_clone;

            while let Some(batches) = multi_reader.next() {
                let labels = batches
                    .iter()
                    .map(|(id, _)| group_mapping.get(id).unwrap().clone())
                    .collect();
                let data = batches.into_iter().map(|(_, batch)| batch).collect();
                let group = EncodedBsxBatchGroup::try_new(data, Some(labels)).unwrap();
                sender.send(Some(group)).unwrap();
            }
            sender.send(None).unwrap();
        });

        let out = DmrIterator {
            config: config_copy,
            ref_idset,
            group_pair: group_order,
            leftover: None,
            regions_cache: Vec::new(),
            last_chr,
            receiver,
            reader_stat,
            _join_handle: join_handle,
        };

        Ok(out)
    }
}

impl Default for DmrConfig {
    fn default() -> Self {
        Self {
            context: Context::CG,
            n_missing: 0,
            min_coverage: 5,
            diff_threshold: 0.1,
            min_cpgs: 5,
            segment_model: PenaltySegmentModel::default(),
        }
    }
}

pub struct DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default,
{
    config: DmrConfig,
    group_pair: (R, R),
    ref_idset: RefIDSet<Arc<String>>,
    leftover: Option<(MethylatedRegionOwned, MethylatedRegionOwned)>,
    regions_cache: Vec<DMRegion>,
    last_chr: Arc<String>,
    receiver: Receiver<Option<EncodedBsxBatchGroup<R>>>,
    reader_stat: ReaderMetadata,
    _join_handle: std::thread::JoinHandle<()>,
}

impl<R> DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + std::fmt::Debug,
{
    fn t_test(
        chr: Arc<String>,
        left: MethylatedRegionView,
        right: MethylatedRegionView,
    ) -> DMRegion {
        let p_value = utils::welch_t_test(&left.density, &right.density);
        DMRegion::new(
            chr,
            left.positions.first().cloned().unwrap_or(0),
            left.positions.last().cloned().unwrap_or(0),
            p_value,
            left.density.mean(),
            right.density.mean(),
            left.size(),
        )
    }

    pub fn blocks_total(&self) -> usize {
        self.reader_metadata().blocks_total.clone()
    }
}

impl<R> DmrModel<R> for DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + std::fmt::Debug,
{
    fn reader_metadata(&self) -> &ReaderMetadata {
        &self.reader_stat
    }

    fn reader_metadata_mut(&mut self) -> &mut ReaderMetadata {
        &mut self.reader_stat
    }

    fn filter_config(&self) -> FilterConfig {
        self.config.filter_config()
    }

    fn group_pair(&self) -> &(R, R) {
        &self.group_pair
    }

    fn ref_idset(&mut self) -> &mut RefIDSet<Arc<String>> {
        &mut self.ref_idset
    }

    fn last_chr(&self) -> Arc<String> {
        self.last_chr.clone()
    }

    fn set_last_chr(&mut self, chr: Arc<String>) {
        self.last_chr = chr;
    }

    fn process_leftover(
        &mut self,
        new_chr: Arc<String>,
        region_left: MethylatedRegionOwned,
        region_right: MethylatedRegionOwned,
    ) -> (MethylatedRegionOwned, MethylatedRegionOwned) {
        if let Some((leftover_left, leftover_right)) = self.leftover.take() {
            if new_chr != self.last_chr {
                self.regions_cache.push(Self::t_test(
                    self.last_chr.clone(),
                    leftover_left.as_view(),
                    leftover_right.as_view(),
                ));
                (region_left, region_right)
            } else {
                let region_left =
                    MethylatedRegionView::concat(leftover_left.as_view(), region_left.as_view());
                let region_right =
                    MethylatedRegionView::concat(leftover_right.as_view(), region_right.as_view());
                (region_left, region_right)
            }
        } else {
            (region_left, region_right)
        }
    }

    fn get_intersecting_indices(
        &self,
        region_left: MethylatedRegionView,
        region_right: MethylatedRegionView,
    ) -> anyhow::Result<Vec<(usize, usize)>> {
        let segments_left = region_left.segment_indices_penalty(&self.config.segment_model);
        let segments_right = region_right.segment_indices_penalty(&self.config.segment_model);
        let positions = region_left.positions;
        Ok(utils::filter_intersecting_indices(
            &segments_left,
            &segments_right,
            positions,
            self.config.segment_model.union_threshold,
        ))
    }

    fn update_cache(
        &mut self,
        chr: Arc<String>,
        region_left: MethylatedRegionOwned,
        region_right: MethylatedRegionOwned,
        intersecting_indices: Vec<(usize, usize)>,
    ) -> anyhow::Result<()> {
        let mut left = region_left.slice_multiple(&intersecting_indices);
        let mut right = region_right.slice_multiple(&intersecting_indices);

        self.leftover = Some((
            left.pop().unwrap().to_owned(),
            right.pop().unwrap().to_owned(),
        ));

        let mut new_cache = left
            .clone()
            .into_par_iter()
            .zip(right.clone().into_par_iter())
            .map(|(left, right)| Self::t_test(chr.clone(), left, right))
            .collect::<Vec<DMRegion>>();

        self.regions_cache.append(&mut new_cache);
        Ok(())
    }

    fn region_cache(&mut self) -> &mut Vec<DMRegion> {
        &mut self.regions_cache
    }

    fn receiver(&mut self) -> &Receiver<Option<EncodedBsxBatchGroup<R>>> {
        &self.receiver
    }

    fn process_last_leftover(&mut self) -> Option<DMRegion> {
        if let Some(leftover) = self.leftover.take() {
            let region = Self::t_test(
                self.last_chr.clone(),
                leftover.0.as_view(),
                leftover.1.as_view(),
            );
            Some(region)
        } else {
            None
        }
    }
}

impl<R> Iterator for DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + Debug,
{
    type Item = (usize, DMRegion);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Try to pop from region cache first
            if let Some(region) = self.region_cache().pop() {
                return Some((self.current_block(), region));
            }

            // Try to receive a new group
            match self.receiver().recv() {
                Ok(Some(group)) => {
                    self.reader_metadata_mut().current_block += 1;
                    if let Err(e) = self.process_group(group) {
                        error!("Error processing group: {}", e);
                    }
                }
                Ok(None) => {
                    // No more groups; process the last leftover and break
                    return self
                        .process_last_leftover()
                        .map(|region| (self.current_block(), region));
                }
                Err(_) => {
                    // Handle the case where receiving fails (e.g., channel closed)
                    return self
                        .process_last_leftover()
                        .map(|region| (self.current_block(), region));
                }
            }
        }
    }
}
