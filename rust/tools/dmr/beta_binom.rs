use argmin::core::{CostFunction, Executor};
use argmin::solver::neldermead::NelderMead;
use log::debug;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use statrs::distribution::ContinuousCDF;
use statrs::function::gamma::ln_gamma;

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
pub fn fit_null_model(
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
pub fn fit_alt_model(
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
pub fn likelihood_ratio_test(log_likelihood_null: f64, log_likelihood_alt: f64, df: usize) -> f64 {
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
pub struct BetaBinomObservation {
    k: f64,
    n: f64,
    ln_binom: f64,
}

impl BetaBinomObservation {
    /// Creates a new BetaBinomObservation after validating inputs.
    pub(crate) fn new(k: f64, n: f64) -> Result<Self, &'static str> {
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
