pub(crate) mod beta_binom;
pub mod dmr_binom;
pub mod dmr_fast;
pub mod meth_region;
pub mod penalty_segment;
pub mod sure_segment;
pub(crate) mod utils;

use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::utils::types::Context;
use anyhow::anyhow;
use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;
use meth_region::{MethylatedRegionOwned, MethylatedRegionView};
use serde::{Deserialize, Serialize, Serializer};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::sync::Arc;
pub struct ReaderMetadata {
    pub(crate) blocks_total: usize,
    pub(crate) current_block: usize,
}

impl ReaderMetadata {
    pub(crate) fn new(blocks_total: usize) -> Self {
        Self {
            blocks_total,
            current_block: 0,
        }
    }
}

pub struct FilterConfig {
    pub context: Context,
    pub n_missing: usize,
    pub min_coverage: i16,
}

pub trait DmrModel<R>
where
    R: Display + Eq + Hash + Clone + Default + std::fmt::Debug,
{
    fn reader_metadata(&self) -> &ReaderMetadata;
    fn reader_metadata_mut(&mut self) -> &mut ReaderMetadata;
    fn filter_config(&self) -> FilterConfig;
    fn group_pair(&self) -> &(R, R);
    fn ref_idset(&mut self) -> &mut RefIDSet<Arc<String>>;
    fn last_chr(&self) -> Arc<String>;
    fn set_last_chr(&mut self, chr: Arc<String>);

    fn blocks_total(&self) -> usize {
        self.reader_metadata().blocks_total
    }

    fn current_block(&self) -> usize {
        self.reader_metadata().current_block
    }

    /// Processes a group of methylation data, applies filters, and performs statistical tests.
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

        // Create methylated regions for the left and right groups
        let (region_left, region_right) = self.create_regions(&group_left, &group_right)?;

        // Handle leftover regions from the previous group
        let (region_left, region_right) =
            self.process_leftover(new_chr.clone(), region_left, region_right);

        self.set_last_chr(new_chr.clone());

        // Get intersecting indices between the left and right regions
        let intersecting_indices =
            self.get_intersecting_indices(region_left.as_view(), region_right.as_view())?;

        // If there are no intersecting indices, return early
        if intersecting_indices.is_empty() {
            return Ok(());
        }

        // Update the cache with the intersecting regions
        self.update_cache(
            new_chr.clone(),
            region_left,
            region_right,
            intersecting_indices,
        )?;
        Ok(())
    }

    fn process_leftover(
        &mut self,
        new_chr: Arc<String>,
        region_left: MethylatedRegionOwned,
        region_right: MethylatedRegionOwned,
    ) -> (MethylatedRegionOwned, MethylatedRegionOwned);

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

    fn apply_filters(
        &self,
        group: EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<EncodedBsxBatchGroup<R>> {
        group
            .filter_context(self.filter_config().context)?
            .mark_low_counts(self.filter_config().min_coverage)?
            .filter_n_missing(self.filter_config().n_missing)
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

    fn create_regions(
        &self,
        group_left: &EncodedBsxBatchGroup<R>,
        group_right: &EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<(MethylatedRegionOwned, MethylatedRegionOwned)> {
        let positions = group_left.get_positions()?;
        let region_left = MethylatedRegionOwned::new(
            positions.clone(),
            group_left.get_average_density(true)?,
            group_left.get_sum_counts_m()?,
            group_left.get_sum_counts_total()?,
        );
        let region_right = MethylatedRegionOwned::new(
            positions,
            group_right.get_average_density(true)?,
            group_right.get_sum_counts_m()?,
            group_right.get_sum_counts_total()?,
        );
        Ok((region_left, region_right))
    }

    fn get_intersecting_indices(
        &self,
        region_left: MethylatedRegionView,
        region_right: MethylatedRegionView,
    ) -> anyhow::Result<Vec<(usize, usize)>>;

    fn update_cache(
        &mut self,
        chr: Arc<String>,
        region_left: MethylatedRegionOwned,
        region_right: MethylatedRegionOwned,
        intersecting_indices: Vec<(usize, usize)>,
    ) -> anyhow::Result<()>;

    fn region_cache(&mut self) -> &mut Vec<DMRegion>;
    fn receiver(&mut self) -> &std::sync::mpsc::Receiver<Option<EncodedBsxBatchGroup<R>>>;
    fn process_last_leftover(&mut self) -> Option<DMRegion>;
}

fn serialize_scientific<S>(x: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&format!("{:e}", x))
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DMRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    #[serde(serialize_with = "serialize_scientific")]
    pub p_value: f64,
    pub meth_left: f64,
    pub meth_right: f64,
    pub n_cytosines: usize,
    pub meth_diff: f64,
    pub meth_mean: f64,
}

impl DMRegion {
    pub(crate) fn new(
        chr: Arc<String>,
        start: u32,
        end: u32,
        p_value: f64,
        meth_left: f64,
        meth_right: f64,
        n_cytosines: usize,
    ) -> Self {
        DMRegion {
            chr: chr.to_string(),
            start,
            end,
            p_value,
            meth_left,
            meth_right,
            n_cytosines,
            meth_diff: meth_left - meth_right,
            meth_mean: (meth_left + meth_right) / 2.0,
        }
    }

    pub fn meth_diff(&self) -> f64 {
        self.meth_left - self.meth_right
    }

    fn meth_mean(&self) -> f64 {
        (self.meth_left + self.meth_right) / 2.0
    }

    pub fn length(&self) -> u32 {
        self.end - self.start + 1
    }
}
