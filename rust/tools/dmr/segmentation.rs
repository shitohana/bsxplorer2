use crate::tools::dmr::data_structs::SegmentView;
use crate::tools::dmr::tv1d_clone::condat;
use crate::utils::mann_whitney_u;
use crate::utils::types::Context;
use itertools::{izip, Itertools};
use std::collections::BTreeSet;

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

pub fn tv_recurse_segment(
    segment: SegmentView,
    mut initial_l: f64,
    l_min: f64,
    l_coef: f64,
    min_cpg: usize,
    diff_threshold: f64,
    seg_tolerance: f64,
    merge_pvalue: f64,
) -> BTreeSet<SegmentView> {
    let mut result: Vec<SegmentView> = Vec::new();
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
fn merge_adjacent_segments(segments: Vec<SegmentView>, p_threshold: f64) -> Vec<SegmentView> {
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
pub fn merge_segments(segments: Vec<SegmentView>, p_threshold: f64) -> Vec<SegmentView> {
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

pub struct FilterConfig {
    pub context: Context,
    pub n_missing: usize,
    pub min_coverage: i16,
}
