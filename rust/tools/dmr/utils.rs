use statrs::distribution::{ContinuousCDF, StudentsT};
use std::cmp::{max, min};

/// Filters intersecting indices between two sets of indices based on a threshold.
///
/// # Arguments
///
/// * `indices1` - A reference to a vector of tuples representing the first set of indices.
/// * `indices2` - A reference to a vector of tuples representing the second set of indices.
/// * `positions` - A reference to a vector of positions corresponding to the indices.
/// * `threshold` - A float value representing the intersection threshold.
///
/// # Returns
///
/// * `Vec<(usize, usize)>` - A vector of tuples representing the filtered intersecting indices.
pub fn filter_intersecting_indices(
    indices1: &[(usize, usize)],
    indices2: &[(usize, usize)],
    positions: &[u32],
    threshold: f64,
) -> Vec<(usize, usize)> {
    let mut result = Vec::new();
    let mut j = 0;

    for &(start1, end1) in indices1 {
        let pos1_start = positions[start1];
        let pos1_end = positions[end1];

        while j < indices2.len() && positions[indices2[j].1] < pos1_start {
            j += 1;
        }

        for &(start2, end2) in &indices2[j..] {
            let pos2_start = positions[start2];
            let pos2_end = positions[end2];

            if pos2_start > pos1_end {
                break;
            }

            let intersection_length = min(pos1_end, pos2_end) - max(pos1_start, pos2_start);
            let total_length = max(pos1_end, pos2_end) - min(pos1_start, pos2_start);
            let ratio = intersection_length as f64 / total_length as f64;

            if ratio > threshold {
                result.push((start1, end2));
            }
        }
    }

    result
}

fn mean_variance_welford(data: &[f64]) -> (f64, f64) {
    let mut n = 0;
    let mut mean = 0.0;
    let mut m2 = 0.0;
    for &x in data {
        n += 1;
        let delta = x - mean;
        mean += delta / (n as f64);
        let delta2 = x - mean;
        m2 += delta * delta2;
    }
    if n > 1 {
        (mean, m2 / (n as f64 - 1.0))
    } else {
        (mean, 0.0)
    }
}

/// Perform Welch’s t–test on two samples and return the two–tailed p–value.
/// (If there isn’t enough data, we return p = 1.0 so that no merge occurs.)
pub fn welch_t_test(sample1: &[f64], sample2: &[f64]) -> f64 {
    let n1 = sample1.len();
    let n2 = sample2.len();
    if n1 < 2 || n2 < 2 {
        return 1.0;
    }
    let (mean1, var1) = mean_variance_welford(sample1);
    let (mean2, var2) = mean_variance_welford(sample2);
    let se = (var1 / (n1 as f64) + var2 / (n2 as f64)).sqrt();
    if se == 0.0 {
        return 1.0;
    }
    let t_stat = (mean1 - mean2).abs() / se;
    let df = (var1 / (n1 as f64) + var2 / (n2 as f64)).powi(2)
        / ((var1 / (n1 as f64)).powi(2) / ((n1 - 1) as f64)
            + (var2 / (n2 as f64)).powi(2) / ((n2 - 1) as f64));
    let t_dist = StudentsT::new(0.0, 1.0, df).unwrap();
    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat));
    p_value
}

/// Extract segments from a (piecewise–constant) signal.
/// Two adjacent values are considered equal if their difference is below tol.
pub fn get_segments(x: &[f64], tol: f64) -> Vec<(usize, usize)> {
    let mut segments = Vec::new();
    if x.is_empty() {
        return segments;
    }
    let mut start = 0;
    let mut current_val = x[0];
    for (i, &val) in x.iter().enumerate() {
        if (val - current_val).abs() > tol {
            segments.push((start, i - 1));
            start = i;
            current_val = val;
        }
    }
    segments.push((start, x.len() - 1));
    segments
}

/// Merge adjacent segments if the t–test on their underlying data does not show a significant difference.
/// The t–test is performed on the original noisy data from each segment.
/// If the p–value exceeds p_threshold, the segments are merged.
fn merge_adjacent_segments(
    y: &[f64],
    segments: &[(usize, usize)],
    p_threshold: f64,
) -> Vec<(usize, usize)> {
    if segments.is_empty() {
        return vec![];
    }
    let mut merged_segments = Vec::new();
    let mut current_seg = segments[0];

    for seg in segments.iter().skip(1) {
        let sample1 = &y[current_seg.0..=current_seg.1];
        let sample2 = &y[seg.0..=seg.1];
        let p_value = welch_t_test(sample1, sample2);

        // If the difference is not significant, merge the segments.
        if p_value > p_threshold {
            let new_start = current_seg.0;
            let new_end = seg.1;
            current_seg = (new_start, new_end);
        } else {
            merged_segments.push(current_seg);
            current_seg = *seg;
        }
    }
    merged_segments.push(current_seg);
    merged_segments
}

/// Iteratively merge segments until no further merging occurs.
pub fn merge_segments(
    y: &[f64],
    segments: &[(usize, usize)],
    p_threshold: f64,
) -> Vec<(usize, usize)> {
    let mut current_segments = segments.to_vec();
    loop {
        let merged = merge_adjacent_segments(y, &current_segments, p_threshold);
        if merged.len() == current_segments.len() {
            break;
        }
        current_segments = merged;
    }
    current_segments
}
