use std::collections::BTreeSet;

use itertools::{izip, Itertools};

use crate::data_structs::enums::Context;
use crate::tools::dmr::data_structs::SegmentView;
use crate::tools::dmr::tv1d_clone::condat;
use crate::utils::mann_whitney_u;

/// Splits a segment based on the distance between positions.
///
/// Returns a vector of indices where the distance between consecutive positions
/// is greater than or equal to `max_dist`.
///
/// # Arguments
///
/// * `positions` - A slice of unsigned 64-bit integers representing positions.
/// * `max_dist` - The maximum allowed distance between consecutive positions.
///
/// # Returns
///
/// A vector of indices where the segment should be split.
/// # Panics
///
/// * If `positions` is empty.  This should be handled gracefully.
pub fn arg_split_segment(
    positions: &[u64],
    max_dist: u64,
) -> Vec<usize> {
    if positions.is_empty() {
        return Vec::new(); // Return an empty vector instead of panicking.
    }
    (1..positions.len())
        .filter(|&i| positions[i] - positions[i - 1] >= max_dist)
        .collect()
}

/// Recursively segments a [SegmentView] using Total Variation (TV)
/// regularization.
///
/// This function applies the Condat algorithm to smooth the data_structs within
/// the segment and identifies sub-segments based on the smoothed differences
/// between the two groups. It recursively processes these sub-segments until a
/// stopping criterion is met.
///
/// # Arguments
///
/// * `segment` - The `SegmentView` to segment.
/// * `initial_l` - The initial regularization parameter for the Condat
///   algorithm.
/// * `l_min` - The minimum value for the regularization parameter.
/// * `l_coef` - The coefficient by which `initial_l` is divided in each
///   iteration.
/// * `min_cpg` - The minimum number of CpGs required for a segment to be
///   considered valid.
/// * `diff_threshold` - The threshold for the smoothed difference between
///   groups.
/// * `seg_tolerance` -  The tolerance for merging adjacent segments with
///   similar values after the TV-filtering step of Condat.
/// * `merge_pvalue` - The p-value threshold for merging adjacent segments.
///
/// # Returns
///
/// A `BTreeSet` of `SegmentView`s representing the segmented regions. The
/// `BTreeSet` ensures the Segments are sorted and non-overlapping.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::mutable_key_type)]
pub fn tv_recurse_segment(
    segment: SegmentView,
    mut initial_l: f64,
    l_min: f64,
    l_coef: f64,
    min_cpg: usize,
    diff_threshold: f32,
    seg_tolerance: f32,
    merge_pvalue: f64,
) -> BTreeSet<SegmentView> {
    assert!(
        l_coef > 1.0,
        "l division coefficient must be bigger, than 1"
    );

    let mut result: Vec<SegmentView> = Vec::new();
    let mut segment_queue = vec![segment.clone()];
    let mut iterations = 0;
    const MAX_ITERATIONS: usize = 1000;

    loop {
        if segment_queue.is_empty() || iterations > MAX_ITERATIONS {
            break;
        }
        iterations += 1;

        #[allow(unused_mut)]
        let mut cur_seg = segment_queue.pop().unwrap();
        cur_seg.get_pvalue();

        let diff_smoothed = izip!(
            condat(cur_seg.group_a(), initial_l as f32),
            condat(cur_seg.group_b(), initial_l as f32)
        )
        .map(|(a, b)| {
            let diff = a - b;
            if diff.abs() > diff_threshold {
                diff
            }
            else {
                0.0
            }
        })
        .collect_vec();

        let new_segments = extract_segments(&diff_smoothed, seg_tolerance)
            .iter()
            .map(|(start, end, _)| cur_seg.slice(*start, *end))
            .collect_vec();

        let (mut short_segments, mut better_segments) = new_segments
            .iter()
            .fold((Vec::new(), Vec::new()), |(mut short, mut better), s| {
                if s.size() <= min_cpg {
                    short.push(s.clone())
                }
                else if s.get_pvalue() < cur_seg.get_pvalue() {
                    better.push(s.clone())
                }
                (short, better)
            });

        result.append(&mut short_segments);
        if better_segments.is_empty() {
            result.push(cur_seg);
        }
        else {
            segment_queue.append(&mut better_segments);
        }

        if l_coef <= 0.0 {
            // prevent infinite loop if l_coef is negative or 0.
            break;
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

/// Extracts segments from a slice of floats based on a tolerance.
///
/// This function identifies contiguous regions in `x` where the difference
/// between consecutive values does not exceed `tolerance`.
///
/// # Arguments
///
/// * `x` - A slice of floats.
/// * `tolerance` - The maximum allowed difference between consecutive values.
///
/// # Returns
///
/// A vector of tuples, where each tuple contains:
///   - The start index of the segment.
///   - The end index of the segment.
///   - The value of the segment (the value at the start index).
fn extract_segments(
    x: &[f32],
    tolerance: f32,
) -> Vec<(usize, usize, f32)> {
    let mut segments = Vec::new();
    if x.is_empty() {
        return segments;
    }
    if x.len() == 1 {
        return vec![(0, 0, x[0])];
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

/// Merge adjacent segments if the Mann-Whitney U test on their underlying
/// data_structs does not show a significant difference.
///
/// The test is performed on the original, unsmoothed methylation data_structs
/// from each segment. If the p-value exceeds `p_threshold`, the segments are
/// merged.
///
/// # Arguments
///
/// * `segments` - A vector of `SegmentView`s representing the segments to
///   merge.
/// * `p_threshold` - The p-value threshold for merging.
///
/// # Returns
///
/// A vector of `SegmentView`s representing the merged segments.
fn merge_adjacent_segments(
    segments: Vec<SegmentView>,
    p_threshold: f64,
) -> Vec<SegmentView> {
    if segments.is_empty() {
        return vec![];
    }
    let mut merged_segments = Vec::new();
    let mut sorted_segments = segments;
    sorted_segments.sort_by_key(|s| s.rel_start);
    let mut prev_seg = sorted_segments[0].clone();

    for cur_seg in sorted_segments.iter().skip(1) {
        // Check for overlap or adjacency.  If rel_start <= prev_seg.rel_end,
        // they are adjacent or overlapping.
        if cur_seg.rel_start <= prev_seg.rel_end {
            let (_ustat, p_value) =
                mann_whitney_u(prev_seg.mds_orig(), cur_seg.mds_orig());
            // If the difference is not significant, merge the segments.
            if p_value > p_threshold {
                prev_seg = prev_seg.merge(cur_seg);
            }
            else {
                merged_segments.push(prev_seg);
                prev_seg = cur_seg.clone();
            }
        }
        else {
            merged_segments.push(prev_seg);
            prev_seg = cur_seg.clone();
        }
    }
    merged_segments.push(prev_seg);
    merged_segments
}

/// Iteratively merge segments until no further merging occurs.
///
/// This function repeatedly calls `merge_adjacent_segments` until the number
/// of segments remains constant, indicating that no further merging is
/// possible.
///
/// # Arguments
///
/// * `segments` - A vector of `SegmentView`s representing the segments to
///   merge.
/// * `p_threshold` - The p-value threshold for merging.
///
/// # Returns
///
/// A vector of `SegmentView`s representing the final merged segments.
pub fn merge_segments(
    segments: Vec<SegmentView>,
    p_threshold: f64,
) -> Vec<SegmentView> {
    let mut current_segments = segments.to_vec();
    loop {
        current_segments.sort_by_key(|s| s.rel_start);
        let merged =
            merge_adjacent_segments(current_segments.clone(), p_threshold);
        if merged.len() == current_segments.len() {
            break;
        }
        current_segments = merged;
    }
    current_segments
}

/// Configuration for filtering data_structs.
pub struct FilterConfig {
    /// The context to consider (e.g., CpG, CHG, CHH).
    pub context:      Context,
    /// The maximum number of missing values allowed.
    pub n_missing:    usize,
    /// The minimum coverage required.
    pub min_coverage: i16,
}
