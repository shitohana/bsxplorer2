use crate::tools::dmr::{tv1d_clone, utils};
use serde::Serialize;
use statrs::statistics::Statistics;
use std::{iter, ops};

#[derive(Debug, Clone, Serialize)]
pub struct PenaltySegmentModel {
    pub min_seg_length: usize,
    pub num_lambdas: usize,
    pub lambda_low: f64,
    pub penalty_weight: f64,
    pub seg_count_weight: f64,
    pub merge_pvalue: f64,
    pub segmentation_tol: f64,
    pub max_dist: u32,
    pub union_threshold: f64,
}

impl Default for PenaltySegmentModel {
    fn default() -> Self {
        Self {
            max_dist: 100,
            min_seg_length: 10,
            num_lambdas: 100,
            lambda_low: 1e-3,
            penalty_weight: 1e5,
            seg_count_weight: 1e4,
            merge_pvalue: 1e-2,
            segmentation_tol: 1e-6,
            union_threshold: 0.5,
        }
    }
}

impl PenaltySegmentModel {
    /// Optimize λ by grid search (in log–space) over a specified range.
    /// Returns the optimal λ and the corresponding denoised signal.
    ///
    /// You can tune:
    ///   - num_lambdas: grid resolution.
    ///   - lambda_low: lower bound (should be > 0).
    ///   - lambda_high: an upper bound (for example, the overall data range).
    pub fn denoise(&self, signal: &[f64]) -> (f64, Vec<f64>) {
        let mut best_lambda = self.lambda_low;
        let mut best_score = f64::INFINITY;
        let mut best_x = Vec::new();

        let data_min = signal.iter().cloned().fold(f64::INFINITY, f64::min);
        let data_max = signal.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let lambda_high = data_max - data_min; // roughly the overall data range.

        // Grid search in log-space.
        for i in 0..self.num_lambdas {
            let t = i as f64 / (self.num_lambdas - 1) as f64;
            let lambda = self.lambda_low * (lambda_high / self.lambda_low).powf(t);
            let denoised = tv1d_clone::condat(signal, lambda);
            let segments = self.get_segments(&denoised);
            let score = self.objective_function(signal, &denoised, &segments);

            // Uncomment the next line to see the grid search details.
            // println!("λ = {:.6} => {} segments, score = {:.3}", lambda, segments.len(), score);

            if score < best_score {
                best_score = score;
                best_lambda = lambda;
                best_x = denoised;
            }
        }
        (best_lambda, best_x)
    }

    pub fn map_segments(positions: &[u32], segments: Vec<(usize, usize, f64)>) -> Vec<u32> {
        Vec::from_iter(segments.into_iter().map(|(_, end, _)| positions[end]))
    }

    /// Extract segments from a (piecewise–constant) signal.
    /// Two adjacent values are considered equal if their difference is below tol.
    pub fn get_segments(&self, x: &[f64]) -> Vec<(usize, usize, f64)> {
        let mut segments = Vec::new();
        if x.is_empty() {
            return segments;
        }
        let mut start = 0;
        let mut current_val = x[0];
        for (i, &val) in x.iter().enumerate() {
            if (val - current_val).abs() > self.segmentation_tol {
                segments.push((start, i - 1, current_val));
                start = i;
                current_val = val;
            }
        }
        segments.push((start, x.len() - 1, current_val));
        segments
    }

    /// A modified objective function that penalizes both the residual error
    /// and (heavily) any segments shorter than a desired minimum length, as well as
    /// the overall number of segments. (Lower is better.)
    ///
    /// - rss: sum of squared residuals.
    /// - seg_penalty: for each segment shorter than min_seg_length, we add
    ///       penalty_weight * (min_seg_length - seg_length)².
    /// - count_penalty: seg_count_weight * (# segments)
    fn objective_function(&self, y: &[f64], x: &[f64], segments: &[(usize, usize, f64)]) -> f64 {
        let rss: f64 = y
            .iter()
            .zip(x.iter())
            .map(|(yi, xi)| (yi - xi).powi(2))
            .sum();
        let mut seg_penalty = 0.0;
        for &(start, end, _) in segments {
            let seg_len = end - start + 1;
            if seg_len < self.min_seg_length {
                let diff = self.min_seg_length as f64 - seg_len as f64;
                seg_penalty += self.penalty_weight * diff * diff;
            }
        }
        let count_penalty = self.seg_count_weight * segments.len() as f64;
        rss + seg_penalty + count_penalty
    }

    /// Merge adjacent segments if the t–test on their underlying data does not show a significant difference.
    /// The t–test is performed on the original noisy data from each segment.
    /// If the p–value exceeds p_threshold, the segments are merged.
    fn merge_adjacent_segments(
        &self,
        y: &[f64],
        segments: &[(usize, usize, f64)],
    ) -> Vec<(usize, usize, f64)> {
        if segments.is_empty() {
            return vec![];
        }
        let mut merged_segments = Vec::new();
        let mut current_seg = segments[0];

        for seg in segments.iter().skip(1) {
            let sample1 = &y[current_seg.0..=current_seg.1];
            let sample2 = &y[seg.0..=seg.1];
            let p_value = utils::welch_t_test(sample1, sample2);

            // If the difference is not significant, merge the segments.
            if p_value > self.merge_pvalue {
                let new_start = current_seg.0;
                let new_end = seg.1;
                let merged_data = &y[new_start..=new_end];
                let new_mean = merged_data.mean();
                current_seg = (new_start, new_end, new_mean);
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
        &self,
        y: &[f64],
        segments: &[(usize, usize, f64)],
    ) -> Vec<(usize, usize, f64)> {
        let mut current_segments = segments.to_vec();
        loop {
            let merged = self.merge_adjacent_segments(y, &current_segments);
            if merged.len() == current_segments.len() {
                break;
            }
            current_segments = merged;
        }
        current_segments
    }
}
