use crate::tools::dmr::{tv1d_clone, utils};
use log::debug;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, ParallelIterator,
};

#[derive(Clone, Debug)]
pub struct SureSegmentModelConfig {
    pub block_size: usize,
    pub seg_tolerance: f64,
    pub p_threshold: f64,
    pub coarse_steps: usize,
    pub refined_steps: usize,
    pub lambda_min: f64,
    pub lambda_max: f64,
    pub max_dist: u32,
    pub union_threshold: f64,
    pub base_null: (f64, f64),
    pub base_alt: (f64, f64, f64),
    pub epsilon: f64,
    pub max_iters: u64,
}

impl Default for SureSegmentModelConfig {
    fn default() -> Self {
        Self {
            coarse_steps: 20,
            refined_steps: 20,
            lambda_min: 0.001,
            lambda_max: 2.0,
            block_size: 20,
            seg_tolerance: 1e-6,
            p_threshold: 0.01,
            max_dist: 100,
            union_threshold: 0.7,
            base_null: (0.5, 0.1),
            base_alt: (0.55, 0.45, 0.2),
            epsilon: 1e-2,
            max_iters: 5000,
        }
    }
}

/// Computes a robust noise variance estimator using the median absolute deviation (MAD)
/// of the first differences of the signal.
/// Returns an estimate of σ².
fn robust_noise_variance(y: &[f64]) -> f64 {
    // Compute absolute differences between adjacent points.
    let mut diffs: Vec<f64> = y.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    if diffs.len() < 1 {
        return 0.0;
    }
    // Sort the differences to compute the median.
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if diffs.len() % 2 == 0 {
        let mid = diffs.len() / 2;
        (diffs[mid.saturating_sub(1)] + diffs[mid]) / 2.0
    } else {
        diffs[diffs.len() / 2]
    };
    // For Gaussian noise, the MAD is related to the standard deviation by 0.6745.
    let sigma = median / 0.6745;
    sigma * sigma
}

/// Computes an estimate of SURE (Stein’s Unbiased Risk Estimate)
/// for TV-denoising given:
/// - the original data `y`
/// - the denoised signal `f`
/// - the noise variance sigma² (`sigma2`)
///
/// Here we use the formula:
///     SURE(λ) = ||y - f||² + 2σ²·df - nσ²,
/// where the degrees of freedom (df) are approximated as the number
/// of constant segments in `f`.
fn compute_sure(y: &[f64], f: &[f64], sigma2: f64, config: &SureSegmentModelConfig) -> f64 {
    let n = y.len() as f64;
    let residual: f64 = y
        .par_iter()
        .zip(f.par_iter())
        .map(|(yi, fi)| (yi - fi).powi(2))
        .sum();

    // Count segments as changes in f.
    let tol = config.seg_tolerance;
    let mut segments = 1;
    for i in 1..f.len() {
        if (f[i] - f[i - 1]).abs() > tol {
            segments += 1;
        }
    }
    residual + 2.0 * sigma2 * (segments as f64) - n * sigma2
}

/// Computes a block-resampled SURE for a given lambda by partitioning the signal into blocks.
/// This helps account for spatial dependence in the methylation data.
fn block_resampled_sure(
    y: &[f64],
    lambda: f64,
    sigma2: f64,
    config: &SureSegmentModelConfig,
) -> f64 {
    let n = y.len();
    let num_blocks = (n + config.block_size - 1) / config.block_size; // Ceiling division

    let sure_sum: f64 = (0..num_blocks)
        .into_par_iter()
        .map(|i| {
            let start = i * config.block_size;
            let end = ((i + 1) * config.block_size).min(n);
            let block = &y[start..end];
            // Use the available tv_denoise and compute_sure functions on the block.
            let denoised_block = tv1d_clone::condat(block, lambda);
            let sure_block = compute_sure(block, &denoised_block, sigma2, &config);
            sure_block
        })
        .sum();

    sure_sum / (num_blocks as f64)
}

/// Performs a refined candidate grid search for the optimal λ.
/// First, a coarse grid is evaluated using block-resampled SURE, and then a refined grid
/// is searched around the best coarse candidate.
fn refined_candidate_grid(
    y: &[f64],
    sigma2: f64,
    config: &SureSegmentModelConfig,
) -> (f64, Vec<f64>) {
    // --- Coarse Grid Search ---
    let lambda_min = config.lambda_min;
    let lambda_max = config.lambda_max;
    let coarse_steps = config.coarse_steps;

    let coarse_candidates: Vec<(f64, f64, Vec<f64>)> = (0..coarse_steps)
        .into_par_iter()
        .map(|i| {
            let lambda =
                lambda_min + (lambda_max - lambda_min) * (i as f64) / ((coarse_steps - 1) as f64);
            let sure = block_resampled_sure(y, lambda, sigma2, config);
            let denoised = tv1d_clone::condat(y, lambda);
            (lambda, sure, denoised)
        })
        .collect();
    let (mut best_lambda, mut best_sure, mut best_denoised) = coarse_candidates
        .into_iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();

    // Updated refined grid search:
    let refinement_range = best_lambda * 0.5;
    let refined_lambda_min = if best_lambda > refinement_range {
        best_lambda - refinement_range
    } else {
        0.0
    };
    let refined_lambda_max = best_lambda + refinement_range;

    let refined_candidates: Vec<(f64, f64, Vec<f64>)> = (0..config.refined_steps)
        .into_par_iter()
        .map(|i| {
            let lambda = refined_lambda_min
                + (refined_lambda_max - refined_lambda_min) * (i as f64)
                    / ((config.refined_steps - 1) as f64);
            let sure = block_resampled_sure(y, lambda, sigma2, config);
            let denoised = tv1d_clone::condat(y, lambda);
            (lambda, sure, denoised)
        })
        .collect();
    if let Some((lambda, sure, denoised)) = refined_candidates
        .into_iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
    {
        if sure < best_sure {
            best_lambda = lambda;
            best_sure = sure;
            best_denoised = denoised;
        }
    }

    (best_lambda, best_denoised)
}

pub fn segment_signal_sure(y: &[f64], config: &SureSegmentModelConfig) -> Vec<(usize, usize)> {
    let sigma2 = robust_noise_variance(&y);
    let (opt_lambda, denoised_signal) = refined_candidate_grid(&y, sigma2, config);
    debug!("Optimal lambda: {:.6}", opt_lambda);

    let initial_segments = utils::get_segments(&denoised_signal, config.seg_tolerance);
    debug!("Initial segmentation: {} segments", initial_segments.len());
    let merged_segments = utils::merge_segments(&y, &initial_segments, config.p_threshold);
    debug!("Merged segmentation: {} segments", merged_segments.len());
    merged_segments
}
