use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::io::bsx::multiple_reader::MultiBsxFileReader;
use crate::io::bsx::read::BsxFileReader;
use crate::tools::dmr::config::DmrConfig;
use crate::tools::dmr::data_structs::{DMRegion, ReaderMetadata, SegmentOwned};
use crate::tools::dmr::segmentation::{tv_recurse_segment, FilterConfig};
use crate::utils::mann_whitney_u;
use crate::utils::types::Context;
use anyhow::anyhow;
use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;
use log::error;
use rayon::prelude::*;
use statrs::statistics::Statistics;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::io::{Read, Seek};
use std::sync::mpsc::Receiver;
use std::sync::Arc;

pub struct DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default,
{
    pub(crate) config: DmrConfig,
    pub(crate) group_pair: (R, R),
    pub(crate) ref_idset: RefIDSet<Arc<String>>,
    pub(crate) leftover: Option<SegmentOwned>,
    pub(crate) regions_cache: Vec<DMRegion>,
    pub(crate) last_chr: Arc<String>,
    pub(crate) receiver: Receiver<Option<EncodedBsxBatchGroup<R>>>,
    pub(crate) reader_stat: ReaderMetadata,
    pub(crate) _join_handle: std::thread::JoinHandle<()>,
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

        let mut segment = SegmentOwned::new(positions, density_left, density_right);
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
        initial_segments: Vec<SegmentOwned>,
    ) -> anyhow::Result<()> {
        let mut initial_segments = initial_segments;
        self.leftover = Some(initial_segments.pop().unwrap());

        let mut new_cache = initial_segments
            .into_par_iter()
            .map(|segment| {
                let result = tv_recurse_segment(
                    segment.to_view().clone(),
                    self.config.initial_l,
                    self.config.l_min,
                    self.config.l_coef,
                    self.config.min_cpgs,
                    self.config.diff_threshold,
                    self.config.seg_tolerance,
                    self.config.merge_pvalue,
                );
                result
                    .into_iter()
                    // .filter(|s| s.pvalue.is_some())
                    .map(|segment| {
                        let a_mean = segment.group_a().mean();
                        let b_mean = segment.group_b().mean();

                        let pval = if (b_mean - a_mean).abs() > self.config.diff_threshold {
                            let (_, u_test_pval) =
                                mann_whitney_u(segment.group_a(), segment.group_b());
                            u_test_pval
                        } else {
                            segment.pvalue.unwrap_or(1.0)
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
                    .filter(|dmr| dmr.end != dmr.start)
                    .filter(|dmr| {
                        dmr.meth_diff.abs() >= self.config.diff_threshold
                            && dmr.p_value <= self.config.seg_pvalue
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
            let result = tv_recurse_segment(
                leftover.to_view().clone(),
                self.config.initial_l,
                self.config.l_min,
                self.config.l_coef,
                self.config.min_cpgs,
                self.config.diff_threshold,
                self.config.seg_tolerance,
                self.config.merge_pvalue,
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
                .filter(|dmr| dmr.end != dmr.start)
                .filter(|dmr| {
                    dmr.meth_diff.abs() >= self.config.diff_threshold
                        && dmr.p_value <= self.config.seg_pvalue
                })
                .collect_vec();

            self.regions_cache.append(&mut new_cache);
            self.regions_cache.sort_by_key(|d| d.start);
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
            if let Some(region) = if self.regions_cache.len() > 0 {
                Some(self.region_cache().remove(0))
            } else {
                None
            } {
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
