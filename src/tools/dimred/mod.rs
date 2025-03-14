/// ****************************************************************************
/// * Copyright (c) 2025
/// The Prosperity Public License 3.0.0
///
/// Contributor: [shitohana](https://github.com/shitohana)
///
/// Source Code: https://github.com/shitohana/BSXplorer
/// ***************************************************************************

/// ****************************************************************************
/// * Copyright (c) 2025
/// ***************************************************************************
pub mod mle;

use itertools::Itertools;
use num::NumCast;

enum DistanceMetric {
    Eucledian,
    Cosine,
    Canberra,
    Minkowsky(usize),
}

trait MetricAccumulator<N: num::Float>: Clone {
    fn finalize(&self) -> N;
    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    );
}

#[derive(Clone)]
struct MinkowskiAccum<N: num::Float> {
    power:  usize,
    powsum: N,
}

impl<N: num::Float> MinkowskiAccum<N> {
    fn new(power: usize) -> Self {
        Self {
            power,
            powsum: N::from(0.0).unwrap(),
        }
    }
}

impl<N: num::Float> MetricAccumulator<N> for MinkowskiAccum<N> {
    fn finalize(&self) -> N {
        self.powsum
            .powf(N::from(1f64 / self.power as f64).unwrap())
    }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.powsum = self.powsum + (xi - yi).abs().powi(self.power as i32)
    }
}

#[derive(Clone)]
struct CanberraAccum<N: num::Float> {
    sum: N,
}

impl<N: num::Float> CanberraAccum<N> {
    fn new() -> Self {
        Self {
            sum: N::from(0.0).unwrap(),
        }
    }
}

impl<N: num::Float> MetricAccumulator<N> for CanberraAccum<N> {
    fn finalize(&self) -> N { self.sum }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.sum = self.sum + (xi - yi).abs() / (xi.abs() + yi.abs());
    }
}

#[derive(Clone)]
struct EucledianAccum<N: num::Float> {
    sqsum: N,
}

impl<N: num::Float> EucledianAccum<N> {
    fn new() -> Self {
        Self {
            sqsum: N::from(0.0).unwrap(),
        }
    }
}
impl<N: num::Float> MetricAccumulator<N> for EucledianAccum<N> {
    fn finalize(&self) -> N { self.sqsum.sqrt() }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.sqsum = self.sqsum + (xi - yi).powi(2)
    }
}

#[derive(Clone)]
struct CosineAccum<N: num::Float> {
    dot_prod:    N,
    right_sqsum: N,
    left_sqsum:  N,
}

impl<N: num::Float> CosineAccum<N> {
    fn new() -> Self {
        Self {
            dot_prod:    N::from(0.0).unwrap(),
            right_sqsum: N::from(0.0).unwrap(),
            left_sqsum:  N::from(0.0).unwrap(),
        }
    }
}

impl<N: num::Float> MetricAccumulator<N> for CosineAccum<N> {
    fn finalize(&self) -> N {
        N::from(1.0).unwrap()
            - self.dot_prod / (self.right_sqsum.sqrt() * self.left_sqsum.sqrt())
    }

    fn add_ppoints(
        &mut self,
        xi: N,
        yi: N,
    ) {
        self.dot_prod = self.dot_prod + xi * yi;
        self.right_sqsum = self.right_sqsum + xi.powi(2);
        self.left_sqsum = self.left_sqsum + yi.powi(2);
    }
}

struct MethCountsDataOwned {
    positions:      Vec<u32>,
    count_m:        Vec<u16>,
    count_total:    Vec<u16>,
    density:        Vec<f64>,
    density_cumsum: Vec<f64>,
}
// TODO add pseudocounts
// meth_pseudo = sum(methylated_reads) + α₀
// total_pseudo = sum(total_reads) + α₀ + β₀
// p_pseudo = meth_pseudo / total_pseudo
/// BetaBinomial distribution parameters (alpha, beta)
pub struct BetaBinomParams<N: num::Float>(pub N, pub N);

impl<N: num::Float> BetaBinomParams<N> {
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
    /// A tuple (Self, bool) containing the BetaBinomial distribution parameters
    /// and a boolean indicating if the sample size is adequate for
    /// effective approximation.
    pub fn from_counts<M: num::PrimInt>(
        count_success: &[M],
        count_total: &[M],
        confidence: N,
        precision: N,
    ) -> (Self, bool) {
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

    /// Create a new BetaBinomial distribution from counts of successes and
    /// total trials.
    ///
    /// # Arguments
    ///
    /// * `data` - A slice of counts of successes.
    /// * `confidence` - The desired confidence level.
    /// * `precision` - The desired precision.
    ///
    /// # Returns
    ///
    /// A tuple (Self, bool) containing the BetaBinomial distribution parameters
    /// and a boolean indicating if the sample size is adequate for
    /// effective approximation.
    pub fn from_data_n(
        data: &[N],
        n_trials: &[N],
        confidence: N,
        precision: N,
    ) -> (Self, bool) {
        let n_elements = to_num(data.len());
        let prop = data
            .iter()
            .zip(n_trials)
            .map(|(s, t)| *s / *t)
            .collect_vec();
        let mean = prop
            .iter()
            .fold(to_num(0), |acc: N, p| acc + *p)
            / n_elements;
        let variance = prop
            .iter()
            .map(|val| val.powi(2))
            .fold(to_num(0), |acc: N, x| acc + x)
            / n_elements
            - mean.powi(2);

        Self::from_mean_var(mean, variance, n_trials, confidence, precision)
    }

    /// Create a new BetaBinomial distribution from counts of successes and
    /// total trials.
    ///
    /// # Arguments
    ///
    /// * `mean` - The mean of the distribution.
    /// * `variance` - The variance of the distribution.
    /// * `confidence` - The desired confidence level.
    /// * `precision` - The desired precision.
    ///
    /// # Returns
    ///
    /// A tuple (Self, bool) containing the BetaBinomial distribution parameters
    /// and a boolean indicating if the sample size is adequate for
    /// effective approximation.
    pub fn from_mean_var(
        mean: N,
        variance: N,
        trials: &[N],
        confidence: N,
        precision: N,
    ) -> (Self, bool) {
        let mean_trials = trials
            .iter()
            .fold(to_num(0), |acc: N, x| acc + *x)
            / to_num(trials.len());
        let binomial_var = mean * (to_num::<_, N>(1) - mean) / mean_trials;
        // Extra-binomial variance component (overdispersion)
        let extra_var: N = to_num(
            (variance - binomial_var)
                .to_f64()
                .unwrap()
                .max(0.0),
        );

        let (alpha_est, beta_est) = if extra_var > to_num(0) {
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

            let a = mean * (to_num::<_, N>(1) - temp) / temp;
            let b =
                (to_num::<_, N>(1) - mean) * (to_num::<_, N>(1) - temp) / temp;
            (a, b)
        }
        else {
            (
                mean * to_num(MEAN_SCALE),
                (to_num::<_, N>(1) - mean) * to_num(MEAN_SCALE),
            )
        };

        let z: N =
            to_num(z_table::lookup((1.0 + confidence.to_f32().unwrap()) / 2.0));
        let factor = to_num::<_, N>(1) / alpha_est
            + to_num::<_, N>(1) / beta_est
            + to_num::<_, N>(1) / (alpha_est + beta_est);
        let var_alpha = (alpha_est.powi(2)) * factor;
        let var_beta = (beta_est.powi(2)) * factor;

        let denominator = precision.powi(2) * alpha_est.powi(2);
        let n_alpha: N = (to_num::<_, N>(z.powi(2)) * var_alpha) / denominator;
        let n_beta: N = (to_num::<_, N>(z.powi(2)) * var_beta) / denominator;

        let effective_n = n_alpha.max(n_beta).ceil();
        let enough_samples = effective_n <= to_num(trials.len());

        (BetaBinomParams::new(alpha_est, beta_est), enough_samples)
    }
}

const PRECISION_LIMIT: f64 = 1e-6;
const VARIANCE_MIN_RATIO: f64 = 1e-2;
const MEAN_SCALE: u64 = 1000;

fn to_num<F, T>(x: F) -> T
where
    T: NumCast,
    F: num::ToPrimitive, {
    T::from(x).unwrap()
}

enum Trials<'a, N: num::Float> {
    Single(N),
    Multiple(&'a [N]),
}

fn bound_prec(value: f64) -> f64 {
    value
        .min(1.0 - PRECISION_LIMIT)
        .max(PRECISION_LIMIT)
}

#[cfg(test)]
pub(crate) mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn test_minkowski_accum() {
        // Test with p=1 (Manhattan distance)
        let mut accum = MinkowskiAccum::<f64>::new(1);
        accum.add_ppoints(1.0, 4.0);
        accum.add_ppoints(2.0, 6.0);
        assert_approx_eq!(accum.finalize(), 7.0); // |1-4| + |2-6| = 3 + 4 = 7

        // Test with p=2 (Euclidean distance)
        let mut accum = MinkowskiAccum::<f64>::new(2);
        accum.add_ppoints(0.0, 3.0);
        accum.add_ppoints(0.0, 4.0);
        assert_approx_eq!(accum.finalize(), 5.0); // sqrt((0-3)^2 + (0-4)^2) =
                                                  // sqrt(9 + 16) = 5
    }

    #[test]
    fn test_canberra_accum() {
        let mut accum = CanberraAccum::<f64>::new();
        accum.add_ppoints(1.0, 3.0);
        accum.add_ppoints(2.0, 5.0);
        // |1-3|/(|1|+|3|) + |2-5|/(|2|+|5|) = 2/4 + 3/7 = 0.5 + 0.428... = 0.928...
        assert_approx_eq!(accum.finalize(), 0.5 + 3.0 / 7.0, 1e-10);
    }

    #[test]
    fn test_eucledian_accum() {
        let mut accum = EucledianAccum::<f64>::new();
        accum.add_ppoints(1.0, 4.0);
        accum.add_ppoints(2.0, 6.0);
        // sqrt((1-4)^2 + (2-6)^2) = sqrt(9 + 16) = sqrt(25) = 5
        assert_approx_eq!(accum.finalize(), 5.0);
    }

    #[test]
    fn test_cosine_accum() {
        let mut accum = CosineAccum::<f64>::new();
        accum.add_ppoints(1.0, 0.0);
        accum.add_ppoints(0.0, 1.0);
        // dot product = 0, |a| = 1, |b| = 1, cosine similarity = 0, distance =
        // 1-0 = 1
        assert_approx_eq!(accum.finalize(), 1.0);

        // Test case for identical vectors (should give distance 0)
        let mut accum = CosineAccum::<f64>::new();
        accum.add_ppoints(2.0, 2.0);
        accum.add_ppoints(3.0, 3.0);
        assert_approx_eq!(accum.finalize(), 0.0, 1e-10);
    }

    #[test]
    fn test_with_different_numeric_types() {
        // Test with f32
        let mut accum = EucledianAccum::<f32>::new();
        accum.add_ppoints(1.0, 4.0);
        assert_approx_eq!(accum.finalize(), 3.0);

        // Test with f64
        let mut accum = EucledianAccum::<f64>::new();
        accum.add_ppoints(1.0, 4.0);
        assert_approx_eq!(accum.finalize(), 3.0);
    }

    #[test]
    fn test_clone_functionality() {
        let mut accum = EucledianAccum::<f64>::new();
        accum.add_ppoints(1.0, 4.0);

        let cloned = accum.clone();
        assert_approx_eq!(cloned.finalize(), 3.0);

        // Ensure original is unchanged
        assert_approx_eq!(accum.finalize(), 3.0);

        // Modify original after cloning
        accum.add_ppoints(2.0, 6.0);
        assert_approx_eq!(accum.finalize(), 5.0);

        // Ensure clone is still the same
        assert_approx_eq!(cloned.finalize(), 3.0);
    }

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
        let (params, is_sufficient) =
            BetaBinomParams::<f64>::from_counts(&successes, &totals, 0.95, 0.1);

        // We should get a valid result with reasonable parameters
        assert!(is_sufficient);
        assert!(params.0 > 0.0); // Alpha should be positive
        assert!(params.1 > 0.0); // Beta should be positive

        // Alpha should be larger than beta since true_alpha > true_beta
        assert!(params.0 > params.1);

        // Test with insufficient data (should return false for is_sufficient)
        let small_size = 5;
        let (small_successes, small_totals) = generate_beta_binomial_samples(
            true_alpha, true_beta, n_trials, small_size,
        );

        let (params, is_sufficient) = BetaBinomParams::<f64>::from_counts(
            &small_successes,
            &small_totals,
            0.99, // Higher confidence
            0.01, // Higher precision
        );

        // This should return false for is_sufficient
        assert!(!is_sufficient);
    }

    #[test]
    fn test_from_data_n() {
        // Generate sample with known parameters
        let true_alpha = 5.0;
        let true_beta = 5.0;
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

        let (params, is_sufficient) = BetaBinomParams::<f64>::from_data_n(
            &successes_f64,
            &totals_f64,
            0.95,
            0.1,
        );

        // Check that we got a valid result
        assert!(is_sufficient);
        assert!(params.0 > 0.0);
        assert!(params.1 > 0.0);

        // For symmetric distribution (alpha = beta), the estimated parameters
        // should be roughly equal
        let ratio = params.0 / params.1;
        assert!(
            ratio > 0.5 && ratio < 2.0,
            "Estimated alpha/beta ratio should be close to 1 for symmetric \
             distribution"
        );
    }

    #[test]
    fn test_from_mean_var() {
        // Test case with moderate overdispersion
        let mean = 0.5f64; // 70% success rate
        let variance = 0.03f64; // Variance higher than binomial
        let trials = vec![20.0f64; 100]; // 20 samples with 10 trials each

        let (params, is_sufficient) = BetaBinomParams::<f64>::from_mean_var(
            mean, variance, &trials, 0.95, 0.1,
        );

        // Check that we got a valid result
        assert!(is_sufficient);
        assert!(params.0 > 0.0);
        assert!(params.1 > 0.0);
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
