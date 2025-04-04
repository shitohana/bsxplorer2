use itertools::Itertools;
use ndarray::{arr0, Array1, ArrayView1};
use num::Float;
use std::ops::{Div, Mul, Sub};

use crate::tools::dimred::to_num;

const PRECISION_LIMIT: f64 = 1e-6;
// const VARIANCE_MIN_RATIO: f64 = 1e-2;
const MEAN_SCALE: u64 = 1000;

// TODO add pseudocounts
// meth_pseudo = sum(methylated_reads) + α₀
// total_pseudo = sum(total_reads) + α₀ + β₀
// p_pseudo = meth_pseudo / total_pseudo
/// BetaBinomial distribution parameters (alpha, beta)
pub struct BetaBinomParams<N: Float>(pub N, pub N);

/// Functions for estimating alpha and beta parameters
fn estimate_alpha_beta_from_mean_var<N: Float>(
    mean: N,
    variance: N,
    mean_trials: N,
) -> (N, N) {
    let binomial_var = mean * (to_num::<_, N>(1) - mean) / mean_trials;
    // Extra-binomial variance component (overdispersion)
    let extra_var: N = to_num(
        (variance - binomial_var)
            .to_f64()
            .unwrap()
            .max(0.0),
    );

    if extra_var > to_num(0) {
        // Overdispersion factor
        let phi = to_num::<_, N>(1) + extra_var / binomial_var;

        // Estimate alpha and beta
        let temp = {
            let unbounded = if mean_trials > to_num(1) {
                (phi - to_num(1)) / (mean_trials - to_num(1))
            }
            else {
                to_num(PRECISION_LIMIT)
            };
            to_num(bound_prec(unbounded.to_f64().unwrap()))
        };

        let alpha = mean * (to_num::<_, N>(1) - temp) / temp;
        let beta =
            (to_num::<_, N>(1) - mean) * (to_num::<_, N>(1) - temp) / temp;
        (alpha, beta)
    }
    else {
        // Use fixed scale for low variance
        (
            mean * to_num(MEAN_SCALE),
            (to_num::<_, N>(1) - mean) * to_num(MEAN_SCALE),
        )
    }
}

/// Calculate the empirical sample size needed
fn empiric_sample_size<N: Float>(
    alpha_est: N,
    beta_est: N,
    confidence: N,
    precision: N,
) -> usize {
    todo!("Switch from using z_table");
    // let z: N =
    //     to_num(z_table::lookup((1.0 + confidence.to_f32().unwrap()) / 2.0));
    // let factor = to_num::<_, N>(1) / alpha_est
    //     + to_num::<_, N>(1) / beta_est
    //     + to_num::<_, N>(1) / (alpha_est + beta_est);
    // let var_alpha = (alpha_est.powi(2)) * factor;
    // let var_beta = (beta_est.powi(2)) * factor;
    //
    // let denominator = precision.powi(2) * alpha_est.powi(2);
    // let n_alpha: N = (to_num::<_, N>(z.powi(2)) * var_alpha) / denominator;
    // let n_beta: N = (to_num::<_, N>(z.powi(2)) * var_beta) / denominator;
    //
    // n_alpha
    //     .max(n_beta)
    //     .ceil()
    //     .to_usize()
    //     .unwrap_or(usize::MAX)
}

/// Calculate the theoretical sample size needed
fn theoretical_sample_size<N: Float>(
    alpha_est: N,
    beta_est: N,
    n_trials: ArrayView1<N>,
    precision: N,
) -> usize {
    let mu = alpha_est / (alpha_est + beta_est);
    let rho: N =
        to_num::<i32, N>(1) / (alpha_est + beta_est + to_num::<i32, N>(1));
    let q: N = to_num::<i32, N>(1) - mu;

    // OPTIMIZATION: Array operations are vectorized and more efficient than
    // elementwise calculation
    let var_sum = (n_trials.mul(arr0(mu)).mul(arr0(q))
        * (n_trials
            .sub(arr0(to_num::<_, N>(1)))
            .mul(arr0(rho)))
        + arr0(to_num::<_, N>(1)))
    .sum();
    let avg_var = var_sum / to_num(n_trials.len());

    let n_alpha = alpha_est.powi(2) / (precision.powi(2) * avg_var);
    let n_beta = beta_est.powi(2) / (precision.powi(2) * avg_var);

    n_alpha
        .max(n_beta)
        .ceil()
        .to_usize()
        .unwrap_or(usize::MAX)
}

fn bound_prec(value: f64) -> f64 {
    value
        .min(1.0 - PRECISION_LIMIT)
        .max(PRECISION_LIMIT)
}

impl<N: Float + ndarray::ScalarOperand> BetaBinomParams<N> {
    pub fn new(
        alpha: N,
        beta: N,
    ) -> Self {
        assert!(alpha > to_num(0), "Alpha must be positive");
        assert!(beta > to_num(0), "Beta must be positive");
        Self(alpha, beta)
    }

    /// Create a new BetaBinomial distribution from counts of successes and
    /// total trials.
    ///
    /// # Arguments
    ///
    /// * `count_success` - A slice of counts of successes.
    /// * `count_total` - A slice of counts of total trials.
    /// * `confidence` - The desired confidence level.
    /// * `precision` - The desired precision.
    ///
    /// # Returns
    ///
    /// A tuple (Self, empiric_n, theoretical_n) containing the BetaBinomial
    /// distribution parameters and the sample sizes needed for empirical
    /// and theoretical estimations.
    pub fn from_counts<M: num::PrimInt>(
        count_success: &[M],
        count_total: &[M],
        confidence: N,
        precision: N,
    ) -> (Self, usize, usize) {
        Self::from_data_n(
            &count_success
                .iter()
                .cloned()
                .map(to_num::<M, N>)
                .collect_vec(),
            &count_total
                .iter()
                .cloned()
                .map(to_num::<M, N>)
                .collect_vec(),
            confidence,
            precision,
        )
    }

    pub fn from_variable_counts<M: num::PrimInt>(
        count_success: &[M],
        count_total: &[M],
        confidence: N,
        precision: N,
    ) -> (Self, usize, usize) {
        // OPTIMIZATION: Pre-allocate Array1 for better performance
        let count_success_adj = count_success
            .iter()
            .map(|val| to_num::<M, N>(*val))
            .map(|val| val + to_num(1))
            .collect::<Array1<N>>();

        let count_total_adj = count_total
            .iter()
            .map(|val| to_num::<M, N>(*val))
            .map(|val| val + to_num(1))
            .collect::<Array1<N>>();

        let initial_weights = count_success_adj.clone();
        let success_ratio = &count_success_adj / &count_total_adj;

        let (mu, r) = Self::get_params_weights_ratio(
            initial_weights.view(),
            success_ratio.view(),
            count_total_adj.view(),
        );

        let adj_weights = &count_total_adj
            / ((&count_total_adj - to_num::<_, N>(1)) * r + to_num::<_, N>(1));

        let common_var_part = (&success_ratio - mu).map_mut(|x| x.powi(2));

        let initial_variance = {
            let len = count_total_adj.len();
            let coef = len / (len - 1);
            (initial_weights
                .clone()
                .map_mut(|x| x.powi(2))
                * common_var_part.view())
            .sum()
                / initial_weights.sum().powi(2)
                * to_num::<_, N>(coef)
        };
        let adj_variance = {
            let sq_adj_weights = adj_weights
                .clone()
                .map_mut(|x| x.powi(2));
            let numerator = (&sq_adj_weights * common_var_part).sum();
            let denominator = adj_weights.sum().powi(2) - sq_adj_weights.sum();
            numerator / denominator
        };

        let (mu, r) = if initial_variance < adj_variance {
            (mu, r)
        }
        else {
            Self::get_params_weights_ratio(
                adj_weights.view(),
                success_ratio.view(),
                count_total_adj.view(),
            )
        };

        let alpha = mu * (to_num::<_, N>(1) - r) / r;
        let beta = (to_num::<_, N>(1) - mu) * (to_num::<_, N>(1) - r) / r;

        let params = Self::new(alpha, beta);
        let count_total_array = Array1::from_vec(count_total_adj.to_vec());

        let empiric_n = empiric_sample_size(alpha, beta, confidence, precision);
        let theoretical_n = theoretical_sample_size(
            alpha,
            beta,
            count_total_array.view(),
            precision,
        );

        (params, empiric_n, theoretical_n)
    }

    fn get_params_weights_ratio(
        weights: ArrayView1<N>,
        success_ratio: ArrayView1<N>,
        n_trials: ArrayView1<N>,
    ) -> (N, N) {
        let weights_sum = weights.sum();
        let mean_ratio = (&weights * &success_ratio).sum() / weights_sum;
        let var_ratio = (&weights
            * (&success_ratio - mean_ratio).map_mut(|x| x.powi(2)))
        .sum();

        let r = {
            let p = mean_ratio;
            let q = to_num::<_, N>(1) - mean_ratio;
            let one_m_wmean = weights.div(weights_sum) * to_num::<_, N>(-1)
                + to_num::<_, N>(1);
            let sum_wn_onemwmean = (&weights / &n_trials * &one_m_wmean).sum();
            let numerator = var_ratio - p * q * sum_wn_onemwmean;
            let denominator =
                ((&weights * &one_m_wmean).sum() - sum_wn_onemwmean) * p * q;
            numerator / denominator
        };

        (mean_ratio, r)
    }

    /// Create a new BetaBinomial distribution from counts of successes and
    /// total trials.
    ///
    /// # Arguments
    ///
    /// * `data` - A slice of counts of successes.
    /// * `n_trials` - A slice of counts of total trials.
    /// * `confidence` - The desired confidence level.
    /// * `precision` - The desired precision.
    ///
    /// # Returns
    ///
    /// A tuple (Self, empiric_n, theoretical_n) containing the BetaBinomial
    /// distribution parameters and the sample sizes needed for empirical
    /// and theoretical estimations.
    pub fn from_data_n(
        data: &[N],
        n_trials: &[N],
        confidence: N,
        precision: N,
    ) -> (Self, usize, usize) {
        // OPTIMIZATION: Calculate mean and variance in one pass through the
        // data
        let n_elements = to_num(data.len());
        let (sum, sum_sq, _props): (N, N, Vec<N>) = data
            .iter()
            .zip(n_trials)
            .map(|(s, t)| *s / *t)
            .fold(
                (to_num(0), to_num(0), Vec::with_capacity(data.len())),
                |(sum_acc, sum_sq_acc, mut props), p| {
                    props.push(p);
                    (sum_acc + p, sum_sq_acc + p.powi(2), props)
                },
            );

        let mean = sum / n_elements;
        let variance = sum_sq / n_elements - mean.powi(2);

        Self::from_mean_var(mean, variance, n_trials, confidence, precision)
    }

    /// Create a new BetaBinomial distribution from mean and variance.
    ///
    /// # Arguments
    ///
    /// * `mean` - The mean of the distribution.
    /// * `variance` - The variance of the distribution.
    /// * `trials` - A slice of counts of total trials.
    /// * `confidence` - The desired confidence level.
    /// * `precision` - The desired precision.
    ///
    /// # Returns
    ///
    /// A tuple (Self, empiric_n, theoretical_n) containing the BetaBinomial
    /// distribution parameters and the sample sizes needed for empirical
    /// and theoretical estimations.
    pub fn from_mean_var(
        mean: N,
        variance: N,
        trials: &[N],
        confidence: N,
        precision: N,
    ) -> (Self, usize, usize) {
        let mean_trials = trials
            .iter()
            .fold(to_num(0), |acc: N, x| acc + *x)
            / to_num(trials.len());

        let (alpha_est, beta_est) =
            estimate_alpha_beta_from_mean_var(mean, variance, mean_trials);

        let params = BetaBinomParams::new(alpha_est, beta_est);
        let n_trials_array = Array1::from_vec(trials.to_vec());

        let empiric_n =
            empiric_sample_size(alpha_est, beta_est, confidence, precision);
        let theoretical_n = theoretical_sample_size(
            alpha_est,
            beta_est,
            n_trials_array.view(),
            precision,
        );

        (params, empiric_n, theoretical_n)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use assert_approx_eq::assert_approx_eq;

    use crate::tools::dimred::mom::{bound_prec,
                                    BetaBinomParams,
                                    PRECISION_LIMIT};
    #[test]
    fn test_new() {
        // Test valid parameters
        let params = BetaBinomParams::<f64>::new(2.0, 3.0);
        assert_eq!(params.0, 2.0);
        assert_eq!(params.1, 3.0);

        // Test with different numeric type
        let params = BetaBinomParams::<f32>::new(2.0, 3.0);
        assert_eq!(params.0, 2.0);
        assert_eq!(params.1, 3.0);
    }

    #[test]
    #[should_panic(expected = "Alpha must be positive")]
    fn test_new_invalid_alpha() {
        // This should panic with "Alpha must be positive"
        BetaBinomParams::<f64>::new(0.0, 3.0);
    }

    #[test]
    #[should_panic(expected = "Beta must be positive")]
    fn test_new_invalid_beta() {
        // This should panic with "Beta must be positive"
        BetaBinomParams::<f64>::new(2.0, 0.0);
    }
    /// Helper function to generate beta-binomial distributed samples
    /// alpha and beta are parameters of the beta distribution
    /// n_trials is the number of trials for each sample
    /// n_samples is the number of samples to generate
    pub fn generate_beta_binomial_samples(
        alpha: f64,
        beta: f64,
        n_trials: u32,
        n_samples: usize,
    ) -> (Vec<u32>, Vec<u32>) {
        use rand::prelude::*;
        use rand_distr::{Beta, Binomial};

        let mut rng = thread_rng();
        let beta_dist = Beta::new(alpha, beta).unwrap();

        let mut successes = Vec::with_capacity(n_samples);
        let mut totals = Vec::with_capacity(n_samples);

        for _ in 0..n_samples {
            // Sample p from Beta(alpha, beta)
            let p = beta_dist.sample(&mut rng);

            // Sample number of successes from Binomial(n_trials, p)
            let binom = Binomial::new(n_trials as u64, p).unwrap();
            let success_count = binom.sample(&mut rng) as u32;

            successes.push(success_count);
            totals.push(n_trials);
        }
        (successes, totals)
    }

    #[test]
    fn test_from_counts() {
        // Generate sample from beta-binomial with known parameters
        let true_alpha = 25.0;
        let true_beta = 15.0;
        let n_trials = 30u32;
        let n_samples = 100; // Large enough for estimation

        let (successes, totals) = generate_beta_binomial_samples(
            true_alpha, true_beta, n_trials, n_samples,
        );

        // Test with sufficient data
        let (params, empiric_n, theoretical_n) =
            BetaBinomParams::<f64>::from_counts(&successes, &totals, 0.95, 0.1);

        // We should get a valid result with reasonable parameters
        assert!(empiric_n <= n_samples); // Should be sufficient
        assert!(params.0 > 0.0); // Alpha should be positive
        assert!(params.1 > 0.0); // Beta should be positive
        assert!(theoretical_n > 0); // Should be a positive number

        // Alpha should be larger than beta since true_alpha > true_beta
        assert!(params.0 > params.1);

        // Test with insufficient data
        let small_size = 5;
        let (small_successes, small_totals) = generate_beta_binomial_samples(
            true_alpha, true_beta, n_trials, small_size,
        );

        let (params, empiric_n, _) = BetaBinomParams::<f64>::from_counts(
            &small_successes,
            &small_totals,
            0.99, // Higher confidence
            0.01, // Higher precision
        );

        // This should indicate insufficient samples
        assert!(empiric_n > small_size);
    }

    #[test]
    fn test_from_data_n() {
        // Generate sample with known parameters
        let true_alpha = 8.0;
        let true_beta = 2.0;
        let n_trials = 15u32;
        let n_samples = 100;

        let (successes, totals) = generate_beta_binomial_samples(
            true_alpha, true_beta, n_trials, n_samples,
        );

        // Convert to floating point
        let successes_f64: Vec<f64> = successes
            .iter()
            .map(|&x| x as f64)
            .collect();
        let totals_f64: Vec<f64> = totals
            .iter()
            .map(|&x| x as f64)
            .collect();

        let (params, empiric_n, theoretical_n) =
            BetaBinomParams::<f64>::from_data_n(
                &successes_f64,
                &totals_f64,
                0.95,
                0.1,
            );

        // Check that we got a valid result
        assert!(empiric_n <= n_samples);
        assert!(params.0 > 0.0);
        assert!(params.1 > 0.0);
        assert!(theoretical_n > 0);

        // For this distribution (alpha > beta), the estimated parameters ratio
        // should reflect that
        let ratio = params.0 / params.1;
        assert!(ratio > 2.0 && ratio < 6.0);
    }

    #[test]
    fn test_from_mean_var() {
        // Test case with moderate overdispersion
        let mean = 0.5f64; // 50% success rate
        let variance = 0.03f64; // Variance higher than binomial
        let trials = vec![20.0f64; 100]; // 100 samples with 20 trials each

        let (params, empiric_n, theoretical_n) =
            BetaBinomParams::<f64>::from_mean_var(
                mean, variance, &trials, 0.95, 0.1,
            );

        // Check that we got a valid result
        assert!(empiric_n > 0);
        assert!(params.0 > 0.0);
        assert!(params.1 > 0.0);
        assert!(theoretical_n > 0);
    }

    #[test]
    fn test_bound_prec() {
        // Test bounding of values
        assert_approx_eq!(bound_prec(0.5), 0.5); // Within bounds
        assert_approx_eq!(bound_prec(0.0), PRECISION_LIMIT); // Lower bound
        assert_approx_eq!(bound_prec(1.0), 1.0 - PRECISION_LIMIT); // Upper bound
        assert_approx_eq!(bound_prec(-1.0), PRECISION_LIMIT); // Below lower bound
        assert_approx_eq!(bound_prec(2.0), 1.0 - PRECISION_LIMIT); // Above upper bound
    }
}
