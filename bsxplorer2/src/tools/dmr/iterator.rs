use std::ops::Div;
use std::rc::Rc;
use std::sync::Arc;

use bio_types::annot::refids::RefIDSet;
use crossbeam::channel::{Receiver, Sender};
use itertools::Itertools;
use log::error;
use polars::prelude::{
    Column, DataType,
};
use rayon::prelude::*;

use crate::data_structs::batch::merge_replicates;
use crate::data_structs::batch::{colnames, BsxBatchMethods};
use crate::tools::dmr::config::DmrConfig;
use crate::tools::dmr::data_structs::{DMRegion, ReaderMetadata, SegmentOwned};
use crate::tools::dmr::segmentation::tv_recurse_segment;

fn merge_counts(columns: Vec<&Column>) -> Column {
    let mut columns = columns;
    let mut agg = Rc::new(columns.pop().unwrap().to_owned());
    for column in columns {
        agg = Rc::new((agg.as_ref() + column).unwrap());
    }
    Rc::unwrap_or_clone(agg)
}

fn merge_density(columns: Vec<&Column>) -> Column {
    let mut columns = columns;
    let length = columns.len();
    let mut agg = Rc::new(columns.pop().unwrap().to_owned());
    for column in columns {
        agg = Rc::new((agg.as_ref() + column).unwrap());
    }
    Rc::unwrap_or_clone(agg).div(length)
}

pub(crate) fn segment_reading<I, B>(
    mut left_readers: Vec<I>,
    mut right_readers: Vec<I>,
    sender: Sender<(String, SegmentOwned)>,
) where
    I: Iterator<Item = B>,
    B: BsxBatchMethods,
{
    loop {
        match (
            left_readers
                .iter_mut()
                .map(|x| x.next())
                .collect::<Option<Vec<B>>>(),
            right_readers
                .iter_mut()
                .map(|x| x.next())
                .collect::<Option<Vec<B>>>(),
        ) {
            (Some(left), Some(right)) => {
                let left_merged =
                    merge_replicates(left, merge_counts, merge_density)
                        .expect("Failed to merge replicates");
                let right_merged =
                    merge_replicates(right, merge_counts, merge_density)
                        .expect("Failed to merge replicates");

                let chr = left_merged
                    .chr_val()
                    .unwrap()
                    .to_string();
                // TODO add filtering
                // TODO change segment position datatype to u32
                let positions = left_merged
                    .data()
                    .column(colnames::POS_NAME)
                    .unwrap()
                    .cast(&DataType::UInt64)
                    .unwrap()
                    .u64()
                    .unwrap()
                    .to_vec_null_aware()
                    .left()
                    .expect("Unexpected nulls");

                let left_density: Vec<f32> = left_merged
                    .data()
                    .column(colnames::DENSITY_NAME)
                    .unwrap()
                    .cast(&DataType::Float32)
                    .unwrap()
                    .f32()
                    .unwrap()
                    .into_iter()
                    .map(|v| v.unwrap_or(f32::NAN))
                    .collect_vec();
                let right_density: Vec<f32> = right_merged
                    .data()
                    .column(colnames::DENSITY_NAME)
                    .unwrap()
                    .cast(&DataType::Float32)
                    .unwrap()
                    .f32()
                    .unwrap()
                    .into_iter()
                    .map(|v| v.unwrap_or(f32::NAN))
                    .collect_vec();

                let segment =
                    SegmentOwned::new(positions, left_density, right_density);

                sender
                    .send((chr, segment))
                    .expect("Failed to send segment");
            },
            (None, Some(_)) => panic!("Unexpected end of input (left readers)"),
            (Some(_), None) => {
                panic!("Unexpected end of input (right readers)")
            },
            (None, None) => break,
        }
    }
}

/// An iterator over differentially methylated regions (DMRs).
pub struct DmrIterator {
    /// Configuration for DMR analysis.
    pub(crate) config: DmrConfig,
    /// Set of reference IDs.
    pub(crate) ref_idset: RefIDSet<Arc<String>>,
    /// Leftover segment from the previous batch.
    pub(crate) leftover: Option<SegmentOwned>,
    /// Cache of DMRs.
    pub(crate) regions_cache: Vec<DMRegion>,
    /// Last chromosome processed.
    pub(crate) last_chr: Arc<String>,
    /// Receiver for encoded BSx batch groups.
    pub(crate) receiver: Receiver<(String, SegmentOwned)>,
    /// Metadata about the reader.
    pub(crate) reader_stat: ReaderMetadata,
    /// Join handle for the reader thread.
    pub(crate) _join_handle: std::thread::JoinHandle<()>,
}

impl DmrIterator {
    /// Returns the total number of blocks processed.
    pub fn blocks_total(&self) -> usize {
        self.reader_stat.blocks_total
    }

    /// Returns the current block being processed.
    fn current_block(&self) -> usize {
        self.reader_stat.current_block
    }

    /// Returns the last chromosome processed.
    fn last_chr(&self) -> Arc<String> {
        self.last_chr.clone()
    }

    fn comp_segment(
        &self,
        segment: SegmentOwned,
    ) -> Vec<DMRegion> {
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
    }

    /// Updates the cache with new segments.
    fn update_cache(
        &mut self,
        mut preseg: SegmentOwned,
        chr: String,
    ) -> anyhow::Result<()> {
        // Get the chromosome of the current group
        let chr = self.ref_idset.intern(&chr);
        let mut initial_segments = Vec::new();

        if let Some(leftover) = self.leftover.take() {
            if chr != self.last_chr() {
                initial_segments.push(leftover);
            } else {
                preseg = leftover.concat(preseg);
            }
        }
        initial_segments
            .append(&mut preseg.split_by_dist(self.config.max_dist));

        self.last_chr = chr.clone();
        self.leftover = Some(initial_segments.pop().unwrap());

        let mut new_cache = initial_segments
            .into_par_iter()
            .map(|seg| self.comp_segment(seg))
            .flatten()
            .collect::<Vec<_>>();

        self.regions_cache
            .append(&mut new_cache);
        Ok(())
    }

    /// Processes the last leftover segment.
    fn process_last_leftover(&mut self) -> anyhow::Result<bool> {
        if let Some(leftover) = self.leftover.take() {
            let mut new_cache = self.comp_segment(leftover);

            self.regions_cache
                .append(&mut new_cache);
            self.regions_cache
                .sort_by_key(|d| d.start);
            Ok(true)
        } else {
            Ok(false)
        }
    }
}

impl Iterator for DmrIterator {
    type Item = (usize, DMRegion);

    /// Returns the next DMR in the iterator.
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Try to pop from region cache first
            if let Some(region) = if !self.regions_cache.is_empty() {
                Some(self.regions_cache.remove(0))
            } else {
                None
            } {
                return Some((self.current_block(), region));
            }

            // Try to receive a new group
            match self.receiver.recv() {
                Ok((chr, segment)) => {
                    self.reader_stat.current_block += 1;
                    if let Err(e) = self.update_cache(segment, chr) {
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
                    } else {
                        None
                    };
                },
            }
        }
    }
}
