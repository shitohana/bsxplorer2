use itertools::{izip, Itertools};
use rayon::prelude::*;
use statrs::statistics::Statistics;
use std::cmp::{Ordering, PartialEq};
use std::collections::BTreeSet;
use std::ops::Not;
use std::sync::Arc;

mod dmr_model;
use crate::tools::dmr::tv1d_clone::condat;

pub use crate::utils::ks2d_2sample;
pub use crate::utils::mann_whitney_u;
pub use dmr_model::*;

/// Numerator    = |S(a, b)| - (b - a) * |S(start, end)| / (t - s)
///
/// Denominator  = (b - a) * \[1 - (b - a) / (end - start)\]
///
/// Z = Numerator^2 / Denominator

fn round(x: f64, decimals: u32) -> f64 {
    let y = 10i32.pow(decimals) as f64;
    (x * y).round() / y
}
fn scoring_function(mds_cumsum: &[f64], start: usize, end: usize, a: usize, b: usize) -> f64 {
    let numerator = (mds_cumsum[b] - mds_cumsum[a]).abs()
        - (b - a) as f64 / (end - start) as f64 * (mds_cumsum[end - 1] - mds_cumsum[start]).abs();

    let denominator = (b - a) as f64 * (1f64 - (b - a) as f64 / (end - start) as f64);

    round(numerator.powi(2) / denominator, 6)
}

/// ### Returns
///
/// Tuple (a, b, score), where a, b - interval [a, b) and score - maximum score value
fn maximize_scoring(mds_cumsum: &[f64], start: usize, end: usize) -> (usize, usize, f64) {
    assert!(end - start > 0);
    (start..end)
        .combinations(2)
        .par_bridge()
        .map(|ab_vec| {
            let score = scoring_function(mds_cumsum, start, end, ab_vec[0], ab_vec[1]);
            assert!(!score.is_nan());
            (ab_vec[0], ab_vec[1], score)
        })
        .max_by(|x, y| x.2.total_cmp(&y.2))
        .unwrap()
}

#[derive(Clone, Debug)]
pub struct PreSegmentView<'a> {
    pvalue: Option<f64>,
    rel_start: usize,
    rel_end: usize,
    parent: Arc<&'a PreSegmentOwned>,
    mds_cumsum: Vec<f64>,
}

impl<'a> PreSegmentView<'a> {
    pub fn new(
        pvalue: Option<f64>,
        rel_start: usize,
        rel_end: usize,
        parent: Arc<&'a PreSegmentOwned>,
    ) -> PreSegmentView<'a> {
        let iter = &parent.mds_cumsum[rel_start..rel_end];
        let first = iter.first().cloned().unwrap_or(0.0);
        let mds_cumsum = iter.into_iter().map(|x| *x - first).collect::<Vec<_>>();

        PreSegmentView {
            pvalue,
            rel_start,
            rel_end,
            mds_cumsum,
            parent,
        }
    }

    pub fn init_pvalue(&mut self) {
        if self.pvalue.is_none() {
            let (_, prob) = self.ks_test();
            self.pvalue = Some(prob)
        }
    }

    pub fn mds_cumsum(&self) -> &[f64] {
        &self.mds_cumsum
    }
    pub fn mds_orig(&self) -> &[f64] {
        &self.parent.mds_orig[self.rel_start..self.rel_end]
    }
    pub fn group_a(&self) -> &[f64] {
        &self.parent.group_a[self.rel_start..self.rel_end]
    }
    pub fn group_b(&self) -> &[f64] {
        &self.parent.group_b[self.rel_start..self.rel_end]
    }
    pub fn mean_diff(&self) -> f64 {
        self.group_a().mean() - self.group_b().mean()
    }
    pub fn start_pos(&self) -> u64 {
        self.parent.positions[self.rel_start]
    }
    pub fn end_pos(&self) -> u64 {
        self.parent.positions[self.rel_end - 1]
    }
    pub fn positions(&self) -> &[u64] {
        &self.parent.positions[self.rel_start..self.rel_end]
    }
    pub fn size(&self) -> usize {
        self.rel_end - self.rel_start
    }

    pub fn slice(&self, start: usize, end: usize) -> PreSegmentView<'a> {
        PreSegmentView::new(
            None,
            self.rel_start + start,
            if self.rel_start + end <= self.rel_end {
                self.rel_start + end
            } else {
                self.rel_end
            },
            self.parent.clone(),
        )
    }

    pub fn merge(&self, other: &Self) -> PreSegmentView<'a> {
        debug_assert!(other.rel_start.saturating_sub(self.rel_end) < 2);
        PreSegmentView::new(None, self.rel_start, other.rel_end, self.parent.clone())
    }

    pub fn ks_test(&self) -> (f64, f64) {
        ks2d_2sample(
            self.positions(),
            self.group_a(),
            self.positions(),
            self.group_b(),
        )
    }

    pub fn to_owned(&self) -> PreSegmentOwned {
        PreSegmentOwned {
            positions: self.positions().to_vec(),
            group_a: self.group_a().to_vec(),
            group_b: self.group_b().to_vec(),
            mds_cumsum: self.mds_cumsum().to_vec(),
            mds_orig: self.mds_orig().to_vec(),
        }
    }
}

impl PartialEq for PreSegmentOwned {
    fn eq(&self, other: &Self) -> bool {
        self.size() == other.size()
            && self.positions.first() == other.positions.first()
            && self.positions.last() == other.positions.last()
    }
}

impl PartialEq for PreSegmentView<'_> {
    fn eq(&self, other: &PreSegmentView) -> bool {
        if self.parent == other.parent {
            self.rel_start == other.rel_start && self.rel_end == other.rel_end
        } else {
            false
        }
    }
}

impl PartialOrd for PreSegmentView<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.parent == other.parent {
            self.rel_start.partial_cmp(&other.rel_start)
        } else {
            None
        }
    }
}

impl Eq for PreSegmentView<'_> {}

impl Ord for PreSegmentView<'_> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other)
            .unwrap_or_else(|| panic!("Different parent segments"))
    }
}

#[derive(Debug)]
pub struct PreSegmentOwned {
    group_a: Vec<f64>,
    group_b: Vec<f64>,
    mds_cumsum: Vec<f64>,
    mds_orig: Vec<f64>,
    positions: Vec<u64>,
}

impl PreSegmentOwned {
    pub fn new(positions: Vec<u64>, group_a: Vec<f64>, group_b: Vec<f64>) -> Self {
        let mds_orig = group_a
            .iter()
            .zip(group_b.iter())
            .map(|(a, b)| a - b)
            .collect_vec();
        assert_eq!(positions.len(), mds_orig.len());
        let mut cum = Vec::with_capacity(mds_orig.len());
        for &val in mds_orig.iter() {
            cum.push(cum.last().unwrap_or(&0.0) + val);
        }

        Self {
            group_a,
            group_b,
            mds_cumsum: cum,
            mds_orig,
            positions,
        }
    }

    pub fn split_by_dist(mut self, max_dist: u64) -> Vec<Self> {
        let mut split_idxs = arg_split_segment(&self.positions, max_dist);
        split_idxs.reverse();
        if split_idxs.is_empty() {
            vec![self]
        } else {
            let mut res = split_idxs.into_iter().fold(Vec::new(), |mut acc, idx| {
                // TODO not sure here
                let first_cumsum = self.mds_cumsum.get(idx).cloned().unwrap_or(0.0);
                let positions = self.positions.drain(idx..).collect_vec();
                let group_a = self.group_a.drain(idx..).collect_vec();
                let group_b = self.group_b.drain(idx..).collect_vec();
                let mds_orig = self.mds_orig.drain(idx..).collect_vec();
                let mds_cumsum = self
                    .mds_cumsum
                    .drain(idx..)
                    .map(|v| v - first_cumsum)
                    .collect_vec();

                acc.push(Self {
                    positions,
                    group_a,
                    group_b,
                    mds_cumsum,
                    mds_orig,
                });
                acc
            });
            res.reverse();
            res
        }
    }

    pub fn concat(mut self, mut other: PreSegmentOwned) -> PreSegmentOwned {
        self.group_a.append(&mut other.group_a);
        self.group_b.append(&mut other.group_b);
        self.mds_orig.append(&mut other.mds_orig);
        self.positions.append(&mut other.positions);
        self.mds_cumsum.append(&mut other.mds_cumsum);
        self
    }

    pub fn to_view(&self) -> PreSegmentView {
        PreSegmentView::new(None, 0, self.mds_orig.len(), Arc::new(self))
    }

    pub fn slice(&self, start: usize, end: usize) -> PreSegmentView {
        PreSegmentView::new(None, start, end, Arc::new(self))
    }

    pub fn length(&self) -> u64 {
        self.positions
            .first()
            .cloned()
            .unwrap_or(0)
            .saturating_sub(self.positions.last().cloned().unwrap_or(0))
    }

    pub fn size(&self) -> usize {
        self.positions.len()
    }

    pub fn group_a(&self) -> &Vec<f64> {
        &self.group_a
    }

    pub fn group_b(&self) -> &Vec<f64> {
        &self.group_b
    }

    pub fn mds_cumsum(&self) -> &Vec<f64> {
        &self.mds_cumsum
    }

    pub fn mds_orig(&self) -> &Vec<f64> {
        &self.mds_orig
    }

    pub fn positions(&self) -> &Vec<u64> {
        &self.positions
    }
}

pub fn arg_split_segment(positions: &[u64], max_dist: u64) -> Vec<usize> {
    (1..positions.len())
        .into_iter()
        .filter_map(|i| {
            if positions[i] - positions[i - 1] >= max_dist {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

fn tv_recurse_segment(
    segment: PreSegmentView,
    mut initial_l: f64,
    l_min: f64,
    l_coef: f64,
    min_cpg: usize,
    diff_threshold: f64,
    seg_tolerance: f64,
    merge_pvalue: f64,
) -> BTreeSet<PreSegmentView> {
    let mut result: Vec<PreSegmentView> = Vec::new();
    let mut segment_queue = vec![segment.clone()];

    loop {
        if segment_queue.is_empty() {
            break;
        }

        let mut cur_seg = segment_queue.pop().unwrap();
        cur_seg.init_pvalue();

        let diff_smoothed = izip!(
            condat(cur_seg.group_a(), initial_l),
            condat(cur_seg.group_b(), initial_l)
        )
        .map(|(a, b)| {
            let diff = a - b;
            if diff > diff_threshold {
                diff
            } else {
                0.0
            }
        })
        .collect_vec();

        let new_segments = extract_segments(&diff_smoothed, seg_tolerance)
            .into_iter()
            .map(|(start, end, _)| cur_seg.slice(start, end))
            .map(|mut s| {
                s.init_pvalue();
                s
            })
            .collect_vec();

        let (mut short_segments, mut better_segments) =
            new_segments
                .iter()
                .fold((Vec::new(), Vec::new()), |(mut short, mut better), s| {
                    if s.size() <= min_cpg {
                        short.push(s.clone())
                    } else if s.pvalue < cur_seg.pvalue {
                        better.push(s.clone())
                    }
                    (short, better)
                });

        result.append(&mut short_segments);
        if better_segments.len() == 0 {
            result.push(cur_seg);
        } else {
            segment_queue.append(&mut better_segments);
        }

        initial_l /= l_coef;
        if initial_l < l_min {
            break;
        }
    }

    let segmented = result
        .iter()
        .cloned()
        .sorted_by_key(|s| s.rel_start)
        .collect_vec();
    let mut last_end = 0;
    for s in segmented.iter() {
        if s.rel_start != last_end && s.rel_start > last_end {
            result.push(segment.slice(last_end, s.rel_start));
        }
        last_end = s.rel_end;
    }
    if last_end != segment.size() {
        result.push(segment.slice(last_end, segment.rel_end));
    }

    let merged = merge_segments(result, merge_pvalue);

    BTreeSet::from_iter(merged)
}

fn extract_segments(x: &[f64], tolerance: f64) -> Vec<(usize, usize, f64)> {
    let mut segments = Vec::new();
    if x.is_empty() {
        return segments;
    }
    let mut start = 0;
    let mut current_val = x[0];
    for (i, &val) in x.iter().enumerate() {
        if (val - current_val).abs() > tolerance {
            segments.push((start, i - 1, current_val));
            start = i;
            current_val = val;
        }
    }
    segments.push((start, x.len() - 1, current_val));
    segments
}

/// Merge adjacent segments if the t–test on their underlying data does not show a significant difference.
/// The t–test is performed on the original noisy data from each segment.
/// If the p–value exceeds p_threshold, the segments are merged.
fn merge_adjacent_segments(segments: Vec<PreSegmentView>, p_threshold: f64) -> Vec<PreSegmentView> {
    if segments.is_empty() {
        return vec![];
    }
    let mut merged_segments = Vec::new();
    let mut prev_seg = segments[0].clone();

    for cur_seg in segments.iter().skip(1) {
        if cur_seg.rel_start.saturating_sub(prev_seg.rel_end) < 2 {
            let (_ustat, p_value) = mann_whitney_u(prev_seg.mds_orig(), cur_seg.mds_orig());
            // If the difference is not significant, merge the segments.
            if p_value > p_threshold {
                prev_seg = prev_seg.merge(cur_seg);
            } else {
                merged_segments.push(prev_seg);
                prev_seg = cur_seg.clone();
            }
        } else {
            merged_segments.push(prev_seg);
            prev_seg = cur_seg.clone();
        }
    }
    merged_segments.push(prev_seg);
    merged_segments
}

/// Iteratively merge segments until no further merging occurs.
pub fn merge_segments(segments: Vec<PreSegmentView>, p_threshold: f64) -> Vec<PreSegmentView> {
    let mut current_segments = segments.to_vec();
    loop {
        current_segments.sort_by_key(|s| s.rel_start);
        let merged = merge_adjacent_segments(current_segments.clone(), p_threshold);
        if merged.len() == current_segments.len() {
            break;
        }
        current_segments = merged;
    }
    current_segments
}
