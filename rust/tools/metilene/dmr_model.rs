use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::io::bsx::multiple_reader::MultiBsxFileReader;
use crate::io::bsx::read::BsxFileReader;
use crate::tools::dmr::{DMRegion, FilterConfig, ReaderMetadata};
use crate::tools::metilene::ks2d::mann_whitney_u;
use crate::tools::metilene::{recurse_segmentation, PreSegmentOwned};
use crate::utils::types::Context;
use anyhow::anyhow;
use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;
use log::error;
use rayon::prelude::*;
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
    pub valley: f64,
    pub trend_threshold: f64,
    pub max_dist: u64,
}

impl DmrConfig {
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

    pub fn new(
        context: Context,
        n_missing: usize,
        min_coverage: i16,
        diff_threshold: f64,
        min_cpgs: usize,
        valley: f64,
        trend_threshold: f64,
        max_dist: u64,
    ) -> Self {
        Self {
            context,
            n_missing,
            min_coverage,
            diff_threshold,
            min_cpgs,
            valley,
            trend_threshold,
            max_dist,
        }
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
            valley: 0.7,
            trend_threshold: 0.5,
            max_dist: 100,
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
    leftover: Option<PreSegmentOwned>,
    regions_cache: Vec<DMRegion>,
    last_chr: Arc<String>,
    receiver: Receiver<Option<EncodedBsxBatchGroup<R>>>,
    reader_stat: ReaderMetadata,
    _join_handle: std::thread::JoinHandle<()>,
}

impl<R> DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + Debug,
{
    pub fn blocks_total(&self) -> usize {
        self.reader_metadata().blocks_total
    }

    fn current_block(&self) -> usize {
        self.reader_metadata().current_block
    }

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
    fn region_cache(&mut self) -> &mut Vec<DMRegion> {
        &mut self.regions_cache
    }

    fn receiver(&mut self) -> &Receiver<Option<EncodedBsxBatchGroup<R>>> {
        &self.receiver
    }

    fn check_groups(&self, group: &EncodedBsxBatchGroup<R>) -> anyhow::Result<()> {
        if !(group
            .labels()
            .as_ref()
            .map(|groups| groups.iter().unique().count() == 2)
            .unwrap_or(false))
        {
            return Err(anyhow!(
                "There should be EXACTLY two sample groups! Found: {:?}",
                group.labels()
            ));
        }
        Ok(())
    }

    fn divide_groups(
        &self,
        group: EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<(EncodedBsxBatchGroup<R>, EncodedBsxBatchGroup<R>)> {
        let mut individual_groups = group.split_groups().into_iter().collect_vec();
        if individual_groups.len() != 2 {
            return Err(anyhow!("Too many groups"));
        }
        let (group_left, group_right) = if individual_groups[0].0 == self.group_pair().0 {
            (
                individual_groups.pop().unwrap().1,
                individual_groups.pop().unwrap().1,
            )
        } else {
            let first = individual_groups.pop().unwrap().1;
            (individual_groups.pop().unwrap().1, first)
        };

        Ok((group_left, group_right))
    }

    fn apply_filters(
        &self,
        group: EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<EncodedBsxBatchGroup<R>> {
        group
            .filter_context(self.filter_config().context)?
            .mark_low_counts(self.filter_config().min_coverage)?
            .filter_n_missing(self.filter_config().n_missing)
    }

    fn process_group(&mut self, group: EncodedBsxBatchGroup<R>) -> anyhow::Result<()> {
        // Check if the group has exactly two sample groups
        self.check_groups(&group)?;

        // Apply filters to the group
        let group = self.apply_filters(group)?;

        // If the group is empty after filtering, return early
        if group.height() == 0 {
            return Ok(());
        }

        // Get the chromosome of the current group
        let new_chr = self.ref_idset().intern(group.get_chr()?.as_str());

        // Divide the group into left and right groups
        let (group_left, group_right) = self.divide_groups(group)?;
        let (density_left, density_right) = (
            group_left.get_average_density(true)?,
            group_right.get_average_density(true)?,
        );
        let positions = group_left
            .get_positions()?
            .into_iter()
            .map_into::<u64>()
            .collect_vec();

        let mut segment = PreSegmentOwned::new(positions, density_left, density_right);
        let mut presegmented = Vec::new();
        if let Some(leftover) = self.leftover.take() {
            if new_chr != self.last_chr() {
                presegmented.push(leftover);
            } else {
                segment = leftover.concat(segment);
            }
        }
        presegmented.append(&mut segment.split_by_dist(self.config.max_dist));

        self.set_last_chr(new_chr.clone());
        self.update_cache(new_chr.clone(), presegmented)?;
        Ok(())
    }

    fn update_cache(
        &mut self,
        chr: Arc<String>,
        initial_segments: Vec<PreSegmentOwned>,
    ) -> anyhow::Result<()> {
        let mut initial_segments = initial_segments;
        self.leftover = Some(initial_segments.pop().unwrap());

        let mut new_cache = initial_segments
            .into_par_iter()
            .map(|segment| {
                let result = recurse_segmentation(
                    segment.to_view(),
                    self.config.min_cpgs,
                    self.config.trend_threshold,
                    self.config.valley,
                    self.config.diff_threshold,
                );
                result
                    .into_iter()
                    .filter(|s| s.pvalue.is_some())
                    .map(|segment| {
                        let a_mean = segment.group_a().mean();
                        let b_mean = segment.group_b().mean();

                        let pval = if (b_mean - a_mean).abs() > self.config.diff_threshold {
                            let (_, u_test_pval) =
                                mann_whitney_u(segment.group_a(), segment.group_b());
                            u_test_pval
                        } else {
                            segment.pvalue.unwrap()
                        };

                        DMRegion::new(
                            self.last_chr.clone(),
                            segment.start_pos() as u32,
                            segment.end_pos() as u32,
                            pval,
                            a_mean,
                            b_mean,
                            segment.size(),
                        )
                    })
                    .collect_vec()
            })
            .flatten()
            .collect::<Vec<_>>();

        self.regions_cache.append(&mut new_cache);
        Ok(())
    }

    fn process_last_leftover(&mut self) -> anyhow::Result<bool> {
        if let Some(leftover) = self.leftover.take() {
            let result = recurse_segmentation(
                leftover.to_view(),
                self.config.min_cpgs,
                self.config.trend_threshold,
                self.config.valley,
                self.config.diff_threshold,
            );
            let mut new_cache = result
                .into_iter()
                .filter(|s| s.pvalue.is_some())
                .map(|segment| {
                    let a_mean = segment.group_a().mean();
                    let b_mean = segment.group_b().mean();

                    let pval = if (b_mean - a_mean).abs() > self.config.diff_threshold {
                        let (_, u_test_pval) = mann_whitney_u(segment.group_a(), segment.group_b());
                        u_test_pval
                    } else {
                        segment.pvalue.unwrap()
                    };

                    DMRegion::new(
                        self.last_chr.clone(),
                        segment.start_pos() as u32,
                        segment.end_pos() as u32,
                        pval,
                        a_mean,
                        b_mean,
                        segment.size(),
                    )
                })
                .collect_vec();

            self.regions_cache.append(&mut new_cache);
            Ok(true)
        } else {
            Ok(false)
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
                    return if self
                        .process_last_leftover()
                        .expect("Error processing last leftover")
                    {
                        self.next()
                    } else {
                        None
                    };
                }
                Err(_) => {
                    // Handle the case where receiving fails (e.g., channel closed)
                    return if self
                        .process_last_leftover()
                        .expect("Error processing last leftover")
                    {
                        self.next()
                    } else {
                        None
                    };
                }
            }
        }
    }
}
