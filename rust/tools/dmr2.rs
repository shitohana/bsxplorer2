use crate::data_structs::bsx_batch_group::EncodedBsxBatchGroup;
use crate::io::bsx::multiple_reader::MultiBsxFileReader;
use crate::io::bsx::region_read::BsxFileReader;
use crate::utils::types::Context;
use anyhow::anyhow;
use argmin::core::{CostFunction, Executor};
use argmin::solver::neldermead::NelderMead;
use bio_types::annot::refids::RefIDSet;
use itertools::Itertools;
use log::{debug, error};
use rayon::prelude::*;
use rust_lapper::Interval;
use serde::{Deserialize, Serialize};
use statrs::distribution::{ContinuousCDF, StudentsT};
use statrs::function::gamma::ln_gamma;
use statrs::statistics::Statistics;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::hash::Hash;
use std::io::{Read, Seek};
use std::sync::Arc;

/// An owned methylation region. This struct owns its data.
#[derive(Clone, Debug)]
pub struct MethylatedRegionOwned {
    pub positions: Vec<u32>,
    pub density: Vec<f64>,
    pub count_m: Vec<u32>,
    pub count_total: Vec<u32>,
}

impl MethylatedRegionOwned {
    /// Creates a new region, validating that all vectors have the same nonzero length.
    pub fn new(
        positions: Vec<u32>,
        density: Vec<f64>,
        count_m: Vec<u32>,
        count_total: Vec<u32>,
    ) -> Self {
        assert!(!positions.is_empty(), "Region cannot be empty");
        assert_eq!(positions.len(), density.len());
        assert_eq!(positions.len(), count_m.len());
        assert_eq!(positions.len(), count_total.len());
        Self {
            positions,
            density,
            count_m,
            count_total,
        }
    }

    /// Returns a read-only view of the entire region.
    pub fn as_view(&self) -> MethylatedRegionView {
        MethylatedRegionView {
            positions: &self.positions,
            density: &self.density,
            count_m: &self.count_m,
            count_total: &self.count_total,
        }
    }

    /// Returns a read-only view for the slice [start, end) without copying.
    pub fn slice(&self, start: usize, end: usize) -> MethylatedRegionView {
        assert!(start <= end, "Start must be <= end");
        assert!(
            end <= self.positions.len(),
            "Slice bounds exceed region size"
        );
        MethylatedRegionView {
            positions: &self.positions[start..end],
            density: &self.density[start..end],
            count_m: &self.count_m[start..end],
            count_total: &self.count_total[start..end],
        }
    }

    pub fn slice_multiple(&self, indices: &[(usize, usize)]) -> Vec<MethylatedRegionView> {
        indices
            .into_par_iter()
            .map(|(start, end)| self.slice(*start, end + 1))
            .collect()
    }
}

/// A lightweight view into a methylation region.
/// This struct borrows from an underlying owned region without extra allocations.
#[derive(Clone, Debug)]
pub struct MethylatedRegionView<'a> {
    pub positions: &'a [u32],
    pub density: &'a [f64],
    pub count_m: &'a [u32],
    pub count_total: &'a [u32],
}

impl<'a> MethylatedRegionView<'a> {
    pub fn to_owned(self) -> MethylatedRegionOwned {
        MethylatedRegionOwned {
            positions: self.positions.to_vec(),
            density: self.density.to_vec(),
            count_m: self.count_m.to_vec(),
            count_total: self.count_total.to_vec(),
        }
    }
    pub fn size(&self) -> usize {
        self.positions.len()
    }

    /// Returns the number of sites in this view.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn segment_indices(&self, config: &SegmentModelConfig) -> Vec<(usize, usize)> {
        let segment_indices = segment_signal(&self.density, config);
        segment_indices
            .into_par_iter()
            .map(|(start, end)| {
                let mut result = Vec::new();
                let mut current = start;
                for i in start..end {
                    if self.positions[i + 1] - self.positions[i] > config.max_dist {
                        result.push((current, i));
                        current = i + 1
                    }
                }

                result.push((current, end));
                result
            })
            .flatten()
            .collect()
    }

    /// Returns a read-only view for the slice [start, end) without copying.
    pub fn slice(&self, start: usize, end: usize) -> MethylatedRegionView {
        assert!(start <= end, "Start must be <= end");
        assert!(
            end <= self.positions.len(),
            "Slice bounds exceed region size"
        );
        MethylatedRegionView {
            positions: &self.positions[start..end],
            density: &self.density[start..end],
            count_m: &self.count_m[start..end],
            count_total: &self.count_total[start..end],
        }
    }

    pub fn slice_multiple(&self, indices: &[(usize, usize)]) -> Vec<MethylatedRegionView> {
        indices
            .into_par_iter()
            .map(|(start, end)| self.slice(*start, end + 1))
            .collect()
    }

    /// Splits the view at the provided index into two sub-views.
    pub fn split_at(&self, index: usize) -> (MethylatedRegionView<'a>, MethylatedRegionView<'a>) {
        let (pos_left, pos_right) = self.positions.split_at(index);
        let (dens_left, dens_right) = self.density.split_at(index);
        let (cm_left, cm_right) = self.count_m.split_at(index);
        let (ct_left, ct_right) = self.count_total.split_at(index);
        (
            MethylatedRegionView {
                positions: pos_left,
                density: dens_left,
                count_m: cm_left,
                count_total: ct_left,
            },
            MethylatedRegionView {
                positions: pos_right,
                density: dens_right,
                count_m: cm_right,
                count_total: ct_right,
            },
        )
    }

    /// Concatenates two views into a new owned region.
    /// Note: This operation does incur an allocation.
    pub fn concat(
        left: MethylatedRegionView<'a>,
        right: MethylatedRegionView<'a>,
    ) -> MethylatedRegionOwned {
        let mut positions = Vec::with_capacity(left.len() + right.len());
        positions.extend_from_slice(left.positions);
        positions.extend_from_slice(right.positions);

        let mut density = Vec::with_capacity(left.len() + right.len());
        density.extend_from_slice(left.density);
        density.extend_from_slice(right.density);

        let mut count_m = Vec::with_capacity(left.len() + right.len());
        count_m.extend_from_slice(left.count_m);
        count_m.extend_from_slice(right.count_m);

        let mut count_total = Vec::with_capacity(left.len() + right.len());
        count_total.extend_from_slice(left.count_total);
        count_total.extend_from_slice(right.count_total);

        MethylatedRegionOwned {
            positions,
            density,
            count_m,
            count_total,
        }
    }

    pub fn to_beta_binom_observations(&self) -> Vec<BetaBinomObservation> {
        self.count_m
            .iter()
            .zip(self.count_total.iter())
            .map(|(&count_m, &count_total)| {
                BetaBinomObservation::new(count_m as f64, count_total as f64).unwrap()
            })
            .collect()
    }
}

#[derive(Clone, Debug)]
pub struct SegmentModelConfig {
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

impl Default for SegmentModelConfig {
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
            base_alt: (0.5, 0.5, 0.1),
            epsilon: 1e-3,
            max_iters: 1000,
        }
    }
}

/// Fits the null model where all samples share the same p and phi values.
///
/// # Arguments
///
/// * `left` - A vector of `BetaBinomObservation` for the left group.
/// * `right` - A vector of `BetaBinomObservation` for the right group.
/// * `base_null` - A tuple containing the initial values for p and phi.
/// * `epsilon` - A small perturbation value for generating the initial simplex.
/// * `max_iters` - The maximum number of iterations for the Nelder-Mead solver.
///
/// # Returns
///
/// * `Result<f64, argmin::core::Error>` - The log-likelihood of the null model.
fn fit_null_model(
    left: Vec<BetaBinomObservation>,
    right: Vec<BetaBinomObservation>,
    base_null: (f64, f64),
    epsilon: f64,
    max_iters: u64,
) -> Result<f64, argmin::core::Error> {
    // Fit the NULL model (common p and phi for all samples)
    let mut all_observations = left;
    all_observations.extend(right);
    let null_model = NullBetaBinomCost {
        observations: all_observations,
    };
    // For the null model, we have 2 parameters.
    let base_null = vec![base_null.0, base_null.1];
    let simplex_null = generate_initial_simplex(base_null, epsilon);

    let solver = NelderMead::new(simplex_null);
    let res_null = Executor::new(null_model, solver)
        .configure(|state| state.max_iters(max_iters))
        .run()?;
    let best_null_params = res_null.state().best_param.clone().unwrap();
    let best_null_cost = res_null.state().best_cost;
    let log_likelihood_null = -best_null_cost;

    debug!("Null model best parameters: {:?}", best_null_params);
    debug!("Null model log-likelihood: {:.6}", log_likelihood_null);

    Ok(log_likelihood_null)
}

/// Fits the alternative model where group1 and group2 have different p values (p1 and p2), but share phi.
///
/// # Arguments
///
/// * `left` - A vector of `BetaBinomObservation` for group1.
/// * `right` - A vector of `BetaBinomObservation` for group2.
/// * `base_alt` - A tuple containing the initial values for p1, p2, and phi.
/// * `epsilon` - A small perturbation value for generating the initial simplex.
/// * `max_iters` - The maximum number of iterations for the Nelder-Mead solver.
///
/// # Returns
///
/// * `Result<f64, argmin::core::Error>` - The log-likelihood of the alternative model.
fn fit_alt_model(
    left: Vec<BetaBinomObservation>,
    right: Vec<BetaBinomObservation>,
    base_alt: (f64, f64, f64),
    epsilon: f64,
    max_iters: u64,
) -> Result<f64, argmin::core::Error> {
    // Fit the ALTERNATIVE model (different p for group1 and group2, shared phi)
    let alt_model = AltBetaBinomCost {
        observations_group1: left,
        observations_group2: right,
    };
    // For the alternative model, we have 3 parameters.
    let base_alt = vec![base_alt.0, base_alt.1, base_alt.2];
    let simplex_alt = generate_initial_simplex(base_alt, epsilon);

    let solver_alt = NelderMead::new(simplex_alt);
    let res_alt = Executor::new(alt_model, solver_alt)
        .configure(|state| state.max_iters(max_iters))
        .run()?;
    let best_alt_params = res_alt.state().best_param.clone().unwrap();
    let best_alt_cost = res_alt.state().best_cost;
    let log_likelihood_alt = -best_alt_cost;

    debug!("Alt model best parameters: {:?}", best_alt_params);
    debug!("Alt model log-likelihood: {:.6}", log_likelihood_alt);

    Ok(log_likelihood_alt)
}

/// Performs the Likelihood Ratio Test (LRT) given the log-likelihoods of the null and alternative models.
///
/// The LRT statistic is calculated as:
///     LRT = 2 * (log\_likelihood\_alt - log\_likelihood\_null)
///
/// The p-value is computed using the chi-square distribution with the specified degrees of freedom.
///
/// # Arguments
///
/// * `log_likelihood_null` - Log-likelihood of the null model.
/// * `log_likelihood_alt` - Log-likelihood of the alternative model.
/// * `df` - Degrees of freedom for the chi-square distribution.
///
/// # Returns
///
/// * `f64` - The p-value of the LRT.
fn likelihood_ratio_test(log_likelihood_null: f64, log_likelihood_alt: f64, df: usize) -> f64 {
    let test_statistic = 2.0 * (log_likelihood_alt - log_likelihood_null);
    let chi2_dist = statrs::distribution::ChiSquared::new(df as f64).unwrap();
    let p_value = 1.0 - chi2_dist.cdf(test_statistic);
    p_value
}

/// Helper function to generate a non-degenerate initial simplex for a given base vector.
/// For a parameter vector of dimension n, returns n+1 vertices with a small perturbation (epsilon).
fn generate_initial_simplex(base: Vec<f64>, epsilon: f64) -> Vec<Vec<f64>> {
    let n = base.len();
    let mut simplex = Vec::with_capacity(n + 1);
    simplex.push(base.clone());
    for i in 0..n {
        let mut vertex = base.clone();
        vertex[i] += epsilon; // add a small perturbation
        simplex.push(vertex);
    }
    simplex
}

/// A single observation: k methylated reads out of n total reads.
#[derive(Clone)]
struct BetaBinomObservation {
    k: f64,
    n: f64,
    ln_binom: f64,
}

impl BetaBinomObservation {
    /// Creates a new BetaBinomObservation after validating inputs.
    fn new(k: f64, n: f64) -> Result<Self, &'static str> {
        if k < 0.0 || n <= 0.0 || k > n {
            return Err("Invalid values: Ensure 0 <= k <= n and n > 0.");
        }
        let ln_binom = ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0);
        Ok(Self { k, n, ln_binom })
    }
}

/// Null model: all samples share the same p and phi.
struct NullBetaBinomCost {
    observations: Vec<BetaBinomObservation>,
}

impl CostFunction for NullBetaBinomCost {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let p = param[0];
        let phi = param[1];
        // Constrain: p in (0,1) and phi in (0,1)
        if !(0.0..1.0).contains(&p) || !(0.0..1.0).contains(&phi) {
            return Ok(f64::INFINITY);
        }
        let neg_log_likelihood: f64 = self
            .observations
            .par_iter()
            .map(|obs| -log_likelihood_single(obs, p, phi))
            .sum();
        Ok(neg_log_likelihood)
    }
}

/// Alternative model: group1 and group2 have different p values (p1 and p2), but share phi.
struct AltBetaBinomCost {
    observations_group1: Vec<BetaBinomObservation>,
    observations_group2: Vec<BetaBinomObservation>,
}

impl CostFunction for AltBetaBinomCost {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        let p1 = param[0];
        let p2 = param[1];
        let phi = param[2];
        if !(0.0..1.0).contains(&p1) || !(0.0..1.0).contains(&p2) || !(0.0..1.0).contains(&phi) {
            return Ok(f64::INFINITY);
        }
        let neg_log_likelihood1: f64 = self
            .observations_group1
            .par_iter()
            .map(|obs| -log_likelihood_single(obs, p1, phi))
            .sum();
        let neg_log_likelihood2: f64 = self
            .observations_group2
            .par_iter()
            .map(|obs| -log_likelihood_single(obs, p2, phi))
            .sum();
        Ok(neg_log_likelihood1 + neg_log_likelihood2)
    }
}

/// Compute the log-likelihood for one observation given parameters p and phi.
///
/// We use the reparameterization:
///   alpha = p * ((1 - phi) / phi)
///   beta  = (1 - p) * ((1 - phi) / phi)
fn log_likelihood_single(obs: &BetaBinomObservation, p: f64, phi: f64) -> f64 {
    let alpha = p * ((1.0 - phi) / phi);
    let beta = (1.0 - p) * ((1.0 - phi) / phi);
    let ln_b1 =
        ln_gamma(obs.k + alpha) + ln_gamma(obs.n - obs.k + beta) - ln_gamma(obs.n + alpha + beta);
    let ln_b0 = ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(alpha + beta);
    obs.ln_binom + ln_b1 - ln_b0
}

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

/// Computes a robust noise variance estimator using the median absolute deviation (MAD)
/// of the first differences of the signal.
/// Returns an estimate of σ².
fn robust_noise_variance(y: &[f64]) -> f64 {
    // Compute absolute differences between adjacent points.
    let mut diffs: Vec<f64> = y.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    // Sort the differences to compute the median.
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if diffs.len() % 2 == 0 {
        let mid = diffs.len() / 2;
        (diffs[mid - 1] + diffs[mid]) / 2.0
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
fn compute_sure(y: &[f64], f: &[f64], sigma2: f64, config: &SegmentModelConfig) -> f64 {
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
fn block_resampled_sure(y: &[f64], lambda: f64, sigma2: f64, config: &SegmentModelConfig) -> f64 {
    let n = y.len();
    let num_blocks = (n + config.block_size - 1) / config.block_size; // Ceiling division

    let sure_sum: f64 = (0..num_blocks)
        .into_par_iter()
        .map(|i| {
            let start = i * config.block_size;
            let end = ((i + 1) * config.block_size).min(n);
            let block = &y[start..end];
            // Use the available tv_denoise and compute_sure functions on the block.
            let denoised_block = tv1d::condat(block, lambda);
            let sure_block = compute_sure(block, &denoised_block, sigma2, &config);
            sure_block
        })
        .sum();

    sure_sum / (num_blocks as f64)
}

/// Performs a refined candidate grid search for the optimal λ.
/// First, a coarse grid is evaluated using block-resampled SURE, and then a refined grid
/// is searched around the best coarse candidate.
fn refined_candidate_grid(y: &[f64], sigma2: f64, config: &SegmentModelConfig) -> (f64, Vec<f64>) {
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
            let denoised = tv1d::condat(y, lambda);
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
            let denoised = tv1d::condat(y, lambda);
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
fn welch_t_test(sample1: &[f64], sample2: &[f64]) -> f64 {
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
fn get_segments(x: &[f64], tol: f64) -> Vec<(usize, usize)> {
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
fn merge_segments(y: &[f64], segments: &[(usize, usize)], p_threshold: f64) -> Vec<(usize, usize)> {
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

fn segment_signal(y: &[f64], config: &SegmentModelConfig) -> Vec<(usize, usize)> {
    let sigma2 = robust_noise_variance(&y);
    let (opt_lambda, denoised_signal) = refined_candidate_grid(&y, sigma2, config);
    debug!("Optimal lambda: {:.6}", opt_lambda);

    let initial_segments = get_segments(&denoised_signal, config.seg_tolerance);
    debug!("Initial segmentation: {} segments", initial_segments.len());
    let merged_segments = merge_segments(&y, &initial_segments, config.p_threshold);
    debug!("Merged segmentation: {} segments", merged_segments.len());
    merged_segments
}

#[derive(Debug, Clone)]
pub struct DmrConfig {
    pub context: Context,
    pub n_missing: usize,
    pub min_coverage: i16,
    pub diff_threshold: f64,
    pub min_cpgs: usize,
    pub segment_model: SegmentModelConfig,
}

impl DmrConfig {
    pub fn new(
        context: Context,
        n_missing: usize,
        min_coverage: i16,
        diff_threshold: f64,
        min_cpgs: usize,
        segment_model: SegmentModelConfig,
    ) -> Self {
        Self {
            context,
            n_missing,
            min_coverage,
            diff_threshold,
            min_cpgs,
            segment_model,
        }
    }

    pub fn try_finish<F, R>(&self, readers: Vec<(R, F)>) -> anyhow::Result<DmrIterator<R>>
    where
        F: Read + Seek + Send + Sync + 'static,
        R: Display + Eq + Hash + Clone + Default + std::fmt::Debug + Send + 'static,
    {
        let sample_mapping: HashMap<uuid::Uuid, (R, F)> = HashMap::from_iter(
            readers
                .into_iter()
                .map(|(group, handle)| (uuid::Uuid::new_v4(), (group, handle))),
        );

        let group_mapping: HashMap<uuid::Uuid, R> = HashMap::from_iter(
            sample_mapping
                .iter()
                .map(|(id, (group, _))| (id.clone(), group.clone())),
        );
        let group_order = {
            let group_set: HashSet<R> =
                HashSet::from_iter(group_mapping.values().map(|x| x.clone()));
            if group_set.len() != 2 {
                return Err(anyhow!(
                    "There should be only two groups of samples! ({:?})",
                    group_set
                ));
            }
            let groups_vec = group_set.into_iter().collect_vec();
            (groups_vec[0].clone(), groups_vec[1].clone())
        };
        let readers_mapping: HashMap<uuid::Uuid, BsxFileReader<F>> = HashMap::from_iter(
            sample_mapping
                .into_iter()
                .map(|(id, (_, handle))| (id, BsxFileReader::new(handle))),
        );

        let config_copy = self.clone();

        let (sender, receiver) = std::sync::mpsc::sync_channel(10);
        let last_chr = Arc::new(String::new());
        let ref_idset = RefIDSet::new();

        let group_mapping_clone = group_mapping.clone();
        let mut multi_reader =
            MultiBsxFileReader::try_new(readers_mapping).map_err(|e| anyhow!("{:?}", e))?;
        let reader_stat = ReaderMetadata::new(multi_reader.blocks_total());

        let join_handle = std::thread::spawn(move || {
            let group_mapping = group_mapping_clone;

            while let Some(batches) = multi_reader.next() {
                let labels = batches
                    .iter()
                    .map(|(id, _)| group_mapping.get(id).unwrap().clone())
                    .collect();
                let data = batches.into_iter().map(|(_, batch)| batch).collect();
                let group = EncodedBsxBatchGroup::try_new(data, Some(labels)).unwrap();
                sender.send(Some(group)).unwrap();
            }
            sender.send(None).unwrap();
        });

        let out = DmrIterator {
            config: config_copy,
            ref_idset,
            group_pair: group_order,
            leftover: None,
            regions_cache: Vec::new(),
            last_chr,
            receiver,
            reader_stat,
            _join_handle: join_handle,
        };

        Ok(out)
    }
}

impl Default for DmrConfig {
    fn default() -> Self {
        Self {
            context: Context::CG,
            n_missing: 0,
            min_coverage: 5,
            diff_threshold: 0.1,
            min_cpgs: 5,
            segment_model: SegmentModelConfig::default(),
        }
    }
}

struct ReaderMetadata {
    blocks_total: usize,
    current_block: usize,
}

impl ReaderMetadata {
    fn new(blocks_total: usize) -> Self {
        Self {
            blocks_total,
            current_block: 0,
        }
    }
}

pub struct DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default,
{
    config: DmrConfig,
    group_pair: (R, R),
    ref_idset: RefIDSet<Arc<String>>,
    leftover: Option<(MethylatedRegionOwned, MethylatedRegionOwned)>,
    regions_cache: Vec<DMRegion>,
    last_chr: Arc<String>,
    receiver: std::sync::mpsc::Receiver<Option<EncodedBsxBatchGroup<R>>>,
    reader_stat: ReaderMetadata,
    _join_handle: std::thread::JoinHandle<()>,
}

impl<R> DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + std::fmt::Debug + Sync,
{
    pub fn blocks_total(&self) -> usize {
        self.reader_stat.blocks_total
    }

    pub fn current_block(&self) -> usize {
        self.reader_stat.current_block
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
        let new_chr = self.ref_idset.intern(group.get_chr()?.as_str());

        // Divide the group into left and right groups
        let (group_left, group_right) = self.divide_groups(group)?;

        // Create methylated regions for the left and right groups
        let (mut region_left, mut region_right) = self.create_regions(&group_left, &group_right)?;

        // Handle leftover regions from the previous group
        if let Some((leftover_left, leftover_right)) = self.leftover.take() {
            if new_chr != self.last_chr {
                self.regions_cache.push(Self::bbinom_test(
                    self.last_chr.clone(),
                    leftover_left.as_view(),
                    leftover_right.as_view(),
                    Arc::new(self.config.segment_model.clone()),
                )?);
            } else {
                region_left =
                    MethylatedRegionView::concat(leftover_left.as_view(), region_left.as_view());
                region_right =
                    MethylatedRegionView::concat(leftover_right.as_view(), region_right.as_view());
            }
        }

        self.last_chr = new_chr.clone();

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
            .filter_context(self.config.context)?
            .mark_low_counts(self.config.min_coverage)?
            .filter_n_missing(self.config.n_missing)
    }

    fn divide_groups(
        &self,
        group: EncodedBsxBatchGroup<R>,
    ) -> anyhow::Result<(EncodedBsxBatchGroup<R>, EncodedBsxBatchGroup<R>)> {
        let mut individual_groups = group.split_groups().into_iter().collect_vec();
        if individual_groups.len() != 2 {
            return Err(anyhow!("Too many groups"));
        }
        let (group_left, group_right) = if individual_groups[0].0 == self.group_pair.0 {
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
    ) -> anyhow::Result<Vec<(usize, usize)>> {
        let segments_left = region_left.segment_indices(&self.config.segment_model);
        let segments_right = region_right.segment_indices(&self.config.segment_model);
        let positions = region_left.positions;
        Ok(filter_intersecting_indices(
            &segments_left,
            &segments_right,
            positions,
            self.config.segment_model.union_threshold,
        ))
    }

    fn update_cache(
        &mut self,
        chr: Arc<String>,
        region_left: MethylatedRegionOwned,
        region_right: MethylatedRegionOwned,
        intersecting_indices: Vec<(usize, usize)>,
    ) -> anyhow::Result<()> {
        let mut left = region_left.slice_multiple(&intersecting_indices);
        let mut right = region_right.slice_multiple(&intersecting_indices);

        self.leftover = Some((
            left.pop().unwrap().to_owned(),
            right.pop().unwrap().to_owned(),
        ));
        let config = Arc::new(self.config.segment_model.clone());

        let mut new_cache = left
            .clone()
            .into_par_iter()
            .zip(right.clone().into_par_iter())
            .map(|(left, right)| Self::bbinom_test(chr.clone(), left, right, config.clone()))
            .collect::<anyhow::Result<Vec<DMRegion>>>()?;

        self.regions_cache.append(&mut new_cache);
        Ok(())
    }

    fn bbinom_test(
        chr: Arc<String>,
        left: MethylatedRegionView,
        right: MethylatedRegionView,
        config: Arc<SegmentModelConfig>,
    ) -> anyhow::Result<DMRegion> {
        let left_obs = left.to_beta_binom_observations();
        let right_obs = right.to_beta_binom_observations();

        let null = fit_null_model(
            left_obs.clone(),
            right_obs.clone(),
            config.base_null,
            config.epsilon,
            config.max_iters,
        );
        let alt = fit_alt_model(
            left_obs,
            right_obs,
            config.base_alt,
            config.epsilon,
            config.max_iters,
        );

        if let (Ok(null), Ok(alt)) = (null, alt) {
            let p_value = likelihood_ratio_test(null, alt, 1);
            Ok(DMRegion::new(
                chr,
                left.positions.first().cloned().unwrap_or(0),
                left.positions.last().cloned().unwrap_or(0),
                p_value,
                left.density.mean(),
                right.density.mean(),
                left.size(),
            ))
        } else {
            Err(anyhow!("Error fitting models"))
        }
    }
}

impl<R> Iterator for DmrIterator<R>
where
    R: Display + Eq + Hash + Clone + Default + std::fmt::Debug + Sync,
{
    type Item = (usize, DMRegion);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(region) = self.regions_cache.pop() {
            Some((self.current_block(), region))
        } else {
            if let Ok(new_group) = self.receiver.recv() {
                match new_group {
                    Some(group) => {
                        self.reader_stat.current_block += 1;
                        if let Err(e) = self.process_group(group) {
                            error!("Error processing group: {}", e);
                        }
                    }
                    None => {}
                }
            } else {
                if let Some(leftover) = self.leftover.take() {
                    let p_value = Self::bbinom_test(
                        self.last_chr.clone(),
                        leftover.0.as_view(),
                        leftover.1.as_view(),
                        Arc::new(self.config.segment_model.clone()),
                    )
                    .expect("Error fitting models");
                    self.regions_cache.push(p_value);
                }
            }
            self.next()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DMRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub p_value: f64,
    pub meth_left: f64,
    pub meth_right: f64,
    pub n_cytosines: usize,
}

impl DMRegion {
    fn new(
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
