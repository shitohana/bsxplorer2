use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::Arc;

use anyhow::anyhow;
use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;
use log::error;
use rayon::prelude::*;
use statrs::statistics::Statistics;

use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::tools::dmr::config::DmrConfig;
use crate::tools::dmr::data_structs::{DMRegion, ReaderMetadata, SegmentOwned};
use crate::tools::dmr::segmentation::{tv_recurse_segment, FilterConfig};
use crate::utils::mann_whitney_u;

/// An iterator over differentially methylated regions (DMRs).
pub struct DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default, {
    /// Configuration for DMR analysis.
    pub(crate) config:        DmrConfig,
    /// Pair of sample groups being compared.
    pub(crate) group_pair:    (R, R),
    /// Set of reference IDs.
    pub(crate) ref_idset:     RefIDSet<Arc<String>>,
    /// Leftover segment from the previous batch.
    pub(crate) leftover:      Option<SegmentOwned>,
    /// Cache of DMRs.
    pub(crate) regions_cache: Vec<DMRegion>,
    /// Last chromosome processed.
    pub(crate) last_chr:      Arc<String>,
    /// Receiver for encoded BSx batch groups.
    pub(crate) receiver: crossbeam::channel::Receiver<EncodedBsxBatchGroup<R>>,
    /// Metadata about the reader.
    pub(crate) reader_stat:   ReaderMetadata,
    /// Join handle for the reader thread.
    pub(crate) _join_handle:  std::thread::JoinHandle<()>,
}

impl<R> DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + Debug,
{
    /// Returns the total number of blocks processed.
    pub fn blocks_total(&self) -> usize { self.reader_metadata().blocks_total }

    /// Returns the current block being processed.
    fn current_block(&self) -> usize { self.reader_metadata().current_block }

    /// Returns a reference to the reader metadata.
    fn reader_metadata(&self) -> &ReaderMetadata { &self.reader_stat }

    /// Returns a mutable reference to the reader metadata.
    fn reader_metadata_mut(&mut self) -> &mut ReaderMetadata {
        &mut self.reader_stat
    }

    /// Returns the filter configuration.
    fn filter_config(&self) -> FilterConfig { self.config.filter_config() }

    /// Returns a reference to the group pair.
    fn group_pair(&self) -> &(R, R) { &self.group_pair }

    /// Returns a mutable reference to the reference ID set.
    fn ref_idset(&mut self) -> &mut RefIDSet<Arc<String>> {
        &mut self.ref_idset
    }

    /// Returns the last chromosome processed.
    fn last_chr(&self) -> Arc<String> { self.last_chr.clone() }

    /// Sets the last chromosome processed.
    fn set_last_chr(
        &mut self,
        chr: Arc<String>,
    ) {
        self.last_chr = chr;
    }

    /// Returns a mutable reference to the region cache.
    fn region_cache(&mut self) -> &mut Vec<DMRegion> { &mut self.regions_cache }

    /// Returns a reference to the receiver.
    fn receiver(
        &mut self
    ) -> &crossbeam::channel::Receiver<EncodedBsxBatchGroup<R>> {
        &self.receiver
    }

    /// Checks if the group has exactly two sample groups.
    fn check_groups(
        &self,
        group: &EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<()> {
        if !group
            .labels()
            .as_ref()
            .map(|groups| groups.iter().unique().count() == 2)
            .unwrap_or(false)
        {
            return Err(anyhow!(
                "There should be EXACTLY two sample groups! Found: {:?}",
                group.labels()
            ));
        }
        Ok(())
    }

    /// Divides the group into two subgroups based on the group pair.
    fn divide_groups(
        &self,
        group: EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<(EncodedBsxBatchGroup<R>, EncodedBsxBatchGroup<R>)>
    {
        let mut individual_groups = group
            .split_groups()
            .into_iter()
            .collect_vec();
        if individual_groups.len() != 2 {
            return Err(anyhow!("Too many groups"));
        }
        let (group_left, group_right) =
            if individual_groups[0].0 == self.group_pair().0 {
                (
                    individual_groups.pop().unwrap().1,
                    individual_groups.pop().unwrap().1,
                )
            }
            else {
                let first = individual_groups.pop().unwrap().1;
                (individual_groups.pop().unwrap().1, first)
            };

        Ok((group_left, group_right))
    }

    /// Applies filters to the group.
    fn apply_filters(
        &self,
        group: EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<EncodedBsxBatchGroup<R>> {
        group
            .filter_context(self.filter_config().context)?
            .mark_low_counts(self.filter_config().min_coverage)?
            .filter_n_missing(self.filter_config().n_missing)
    }

    /// Processes a single batch group.
    fn process_group(
        &mut self,
        group: EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<()> {
        // Check if the group has exactly two sample groups
        self.check_groups(&group)?;

        // Apply filters to the group
        let group = self.apply_filters(group)?;

        // If the group is empty after filtering, return early
        if group.height() == 0 {
            return Ok(());
        }

        // Get the chromosome of the current group
        let new_chr = self
            .ref_idset()
            .intern(group.get_chr()?.as_str());

        // Divide the group into left and right groups
        let (group_left, group_right) = self.divide_groups(group)?;
        let (density_left, density_right) = (
            group_left.get_average_density(true)?.into_iter().map_into().collect_vec(),
            group_right.get_average_density(true)?.into_iter().map_into().collect_vec(),
        );
        let positions = group_left
            .get_positions()?
            .into_iter()
            .map_into::<u64>()
            .collect_vec();

        let mut segment =
            SegmentOwned::new(positions, density_left, density_right);
        let mut presegmented = Vec::new();
        if let Some(leftover) = self.leftover.take() {
            if new_chr != self.last_chr() {
                presegmented.push(leftover);
            }
            else {
                segment = leftover.concat(segment);
            }
        }
        presegmented.append(&mut segment.split_by_dist(self.config.max_dist));

        self.set_last_chr(new_chr.clone());
        self.update_cache(presegmented)?;
        Ok(())
    }

    /// Updates the cache with new segments.
    fn update_cache(
        &mut self,
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
                    .map(|segment| {
                        DMRegion::from_segment_view(segment, self.last_chr.to_string())
                    })
                    .filter(|dmr| dmr.end != dmr.start)
                    .collect_vec()
            })
            .flatten()
            .collect::<Vec<_>>();

        self.regions_cache
            .append(&mut new_cache);
        Ok(())
    }

    /// Processes the last leftover segment.
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
                .filter(|s| s.pvalue.get().is_some())
                .map(|segment| {
                    DMRegion::from_segment_view(segment, self.last_chr.to_string())
                })
                .filter(|dmr| dmr.end != dmr.start)
                .collect_vec();

            self.regions_cache
                .append(&mut new_cache);
            self.regions_cache
                .sort_by_key(|d| d.start);
            Ok(true)
        }
        else {
            Ok(false)
        }
    }
}

impl<R> Iterator for DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + Debug,
{
    type Item = (usize, DMRegion);

    /// Returns the next DMR in the iterator.
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Try to pop from region cache first
            if let Some(region) = if self.regions_cache.len() > 0 {
                Some(self.region_cache().remove(0))
            }
            else {
                None
            } {
                return Some((self.current_block(), region));
            }

            // Try to receive a new group
            match self.receiver().recv() {
                Ok(group) => {
                    self.reader_metadata_mut().current_block += 1;
                    if let Err(e) = self.process_group(group) {
                        error!("Error processing group: {}", e);
                    }
                },
                Err(_) => {
                    // Handle the case where receiving fails (e.g., channel
                    // closed)
                    return if self
                        .process_last_leftover()
                        .expect("Error processing last leftover")
                    {
                        self.next()
                    }
                    else {
                        None
                    };
                },
            }
        }
    }
}
