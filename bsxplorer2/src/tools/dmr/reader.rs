use std::hash::Hash;
use std::os::fd::AsRawFd;

use anyhow::anyhow;
use arcstr::ArcStr;
use crossbeam::queue::SegQueue;
use itertools::Itertools;
use polars::error::PolarsResult;
use rayon::prelude::*;

use super::{
    tv_recurse_segment,
    DMRegion,
    DmrConfig,
    SegmentOwned,
};
use crate::data_structs::batch::merge_replicates;
use crate::prelude::{
    AggMethod,
    BsxBatch,
    BsxFileReader,
    MultiBsxFileReader,
};
use crate::utils::THREAD_POOL;
use crate::BsxColumns;

fn comp_segment(
    segment: SegmentOwned,
    config: &DmrConfig,
    chr: ArcStr,
) -> Vec<DMRegion> {
    tv_recurse_segment(
        segment.to_view(),
        config.initial_l,
        config.l_min,
        config.l_coef,
        config.min_cpgs,
        config.diff_threshold,
        config.seg_tolerance,
        config.merge_pvalue,
    )
    .into_iter()
    .filter(|s| s.end_pos() != s.start_pos())
    .map(|s| DMRegion::from_segment_view(s, chr.to_string()))
    .collect_vec()
}

fn merge_for_dmr(
    batches: Vec<BsxBatch>,
    config: &DmrConfig,
) -> PolarsResult<BsxBatch> {
    let bitmap = batches
        .iter()
        .map(|b| b.column(BsxColumns::Density))
        .map(|c| c.is_null())
        .fold(vec![0; batches.len()], |agg, v| {
            let arr_iter = v.iter().map(|b| {
                if b.unwrap() {
                    1
                }
                else {
                    0
                }
            });
            agg.iter()
                .zip(arr_iter)
                .map(|(val, b)| b + *val)
                .collect_vec()
        })
        .into_iter()
        .map(|v| v >= config.n_missing)
        .collect_vec();

    let mut result = merge_replicates(
        batches,
        AggMethod::Sum.get_expr(),
        AggMethod::Mean.get_expr(),
    )?
    .lazy()
    .filter_context(config.context)
    .filter_coverage_gt(config.min_coverage)
    .collect()?;

    result.set_validity_bitmap(BsxColumns::Density, bitmap)?;
    Ok(result)
}

pub struct DmrReader {
    config:  DmrConfig,
    readers: (MultiBsxFileReader, MultiBsxFileReader),
}

impl DmrReader {
    pub fn new(
        config: DmrConfig,
        readers: (MultiBsxFileReader, MultiBsxFileReader),
    ) -> Self {
        Self { config, readers }
    }

    pub fn from_readers<R, F>(
        readers: Vec<(R, F)>,
        config: DmrConfig,
    ) -> anyhow::Result<Self>
    where
        F: AsRawFd + 'static + Sync + Send,
        R: Hash + Eq, {
        let mut readers = readers
            .into_iter()
            .into_group_map()
            .into_values()
            .map(|v| {
                v.into_iter()
                    .map(BsxFileReader::try_new)
                    .collect::<anyhow::Result<Vec<_>>>()
            })
            .collect::<anyhow::Result<Vec<_>>>()?
            .into_iter()
            .map(|v| MultiBsxFileReader::from_iter(v))
            .collect_vec();

        readers.iter_mut().try_for_each(|r| r.validate(false))?;

        match readers.len() {
            0 => anyhow::bail!("No readers supplied"),
            2 => {},
            n => anyhow::bail!("Too many readers groups supplied: {}", n),
        }

        let left = readers.pop().unwrap();
        let right = readers.pop().unwrap();

        let out = Self::new(config, (left, right));

        Ok(out)
    }

    pub fn iter(&mut self) -> DmrIterator<'_> {
        let left_iter = Box::new(self.readers.0.iter());
        let right_iter = Box::new(self.readers.1.iter());

        DmrIterator::new(left_iter, right_iter, &self.config)
    }
}

pub struct DmrIterator<'a> {
    right_iter: Box<dyn Iterator<Item = PolarsResult<Vec<BsxBatch>>> + 'a>,
    left_iter:  Box<dyn Iterator<Item = PolarsResult<Vec<BsxBatch>>> + 'a>,
    config:     &'a DmrConfig,
    results:    SegQueue<DMRegion>,
    last_chr:   ArcStr,
    leftover:   Option<SegmentOwned>,
}

impl<'a> DmrIterator<'a> {
    pub fn new(
        left_iter: Box<dyn Iterator<Item = PolarsResult<Vec<BsxBatch>>> + 'a>,
        right_iter: Box<dyn Iterator<Item = PolarsResult<Vec<BsxBatch>>> + 'a>,
        config: &'a DmrConfig,
    ) -> DmrIterator<'a> {
        let leftover = None;
        let last_chr = ArcStr::from("");
        let results = SegQueue::new();
        Self {
            left_iter,
            right_iter,
            config,
            leftover,
            last_chr,
            results,
        }
    }

    fn get_batches(&mut self) -> anyhow::Result<Option<(BsxBatch, BsxBatch)>> {
        let right = self.right_iter.next().transpose()?;
        let left = self.left_iter.next().transpose()?;

        let iterators_ended = right.is_none() && left.is_none();
        let both_batches_empty = right.as_ref().is_some_and(|b| b.is_empty())
            && left.as_ref().is_some_and(|b| b.is_empty());

        if both_batches_empty || iterators_ended {
            return Ok(None);
        }

        let right = right.ok_or(anyhow!("Unexpected end of input (right readers)"))?;
        let left = left.ok_or(anyhow!("Unexpected end of input (left readers)"))?;

        let right = merge_for_dmr(right, &self.config)?;
        let left = merge_for_dmr(left, &self.config)?;

        Ok(Some((left, right)))
    }

    fn read_segment(&mut self) -> anyhow::Result<bool> {
        let (left, right) = match self.get_batches()? {
            Some(v) => v,
            None => return Ok(false),
        };
        let mut segment = SegmentOwned::try_from((&left, &right))?;

        let chr = ArcStr::from(left.seqname().unwrap());

        let mut initial_segments = Vec::new();

        // If some batch left
        if let Some(leftover) = self.leftover.take() {
            // If another chromosome -> treat as independent segment
            if chr != self.last_chr {
                initial_segments.push(leftover);
            }
            // Otherwise merge with previous
            else {
                segment = leftover.concat(segment);
            }
        }

        let mut seg_by_dist = segment.split_by_dist(self.config.max_dist);
        initial_segments.append(&mut seg_by_dist);

        self.last_chr = chr.clone();
        self.leftover = initial_segments.pop();

        THREAD_POOL.install(|| {
            initial_segments
                .into_par_iter()
                .map(|seg| comp_segment(seg, &self.config, chr.clone()))
                .flatten()
                .for_each(|res| self.results.push(res));
        });

        Ok(true)
    }

    fn process_last_leftover(&mut self) -> bool {
        if let Some(leftover) = self.leftover.take() {
            for res in comp_segment(leftover, &self.config, self.last_chr.clone()) {
                self.results.push(res);
            }
            true
        }
        else {
            false
        }
    }
}

impl Iterator for DmrIterator<'_> {
    type Item = anyhow::Result<DMRegion>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(dmr) = self.results.pop() {
            Some(Ok(dmr))
        }
        else {
            match self.read_segment() {
                Ok(true) => self.next(),
                Ok(false) => {
                    if self.process_last_leftover() {
                        self.next()
                    }
                    else {
                        None
                    }
                },
                Err(e) => Some(Err(e)),
            }
        }
    }
}
