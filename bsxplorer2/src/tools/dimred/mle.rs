// Copyright (c) 2025.
// The Prosperity Public License 3.0.0
// Contributor: [shitohana](https://github.com/shitohana)
// Source Code: https://github.com/shitohana/BSXplorer

use argmin::core::{CostFunction, Error, Gradient};
use statrs::function::factorial::ln_factorial;
use statrs::function::gamma::{digamma, ln_gamma};

/// This is not a good estimator, and currently MoM nearly allways
/// performs better
use crate::tools::dimred::mom::BetaBinomParams;

#[derive(Clone)]
pub struct MethylationData {
    methylated_reads: Vec<u32>,
    total_reads: Vec<u32>,
}

#[derive(Clone)]
pub struct BetaBinomialMLE {
    data: MethylationData,
}

impl BetaBinomialMLE {
    fn new(
        methylated_reads: Vec<u32>,
        total_reads: Vec<u32>,
    ) -> Self {
        BetaBinomialMLE {
            data: MethylationData {
                methylated_reads,
                total_reads,
            },
        }
    }

    fn estimate_initial_params(
        &self,
        confidence: f64,
        error: f64,
    ) -> Vec<f64> {
        let (params, ..) = BetaBinomParams::from_counts(
            &self.data.methylated_reads,
            &self.data.total_reads,
            confidence,
            error,
        );
        vec![
            params.0 / (params.0 + params.1),
            1.0 / (params.0 + params.1 + 1.0),
        ]
    }
}

impl CostFunction for BetaBinomialMLE {
    // Instead of (α, β), use (μ, φ) where:
    // μ = α/(α+β) = mean methylation
    // φ = 1/(α+β+1) = overdispersion parameter
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(
        &self,
        params: &Self::Param,
    ) -> Result<Self::Output, Error> {
        let mu = params[0]
            .min(1.0 - f64::EPSILON.sqrt())
            .max(0.0 + f64::EPSILON.sqrt()); // Mean between 0 and 1
        let phi = params[1]
            .min(1.0 - f64::EPSILON.sqrt())
            .max(0.0 + f64::EPSILON.sqrt()); // Overdispersion between 0 and 1

        // Convert to alpha and beta
        let alpha = mu * (1.0 - phi) / phi;
        let beta = (1.0 - mu) * (1.0 - phi) / phi;

        let ln_beta_ab =
            ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(alpha + beta);

        let log_likelihood = self
            .data
            .methylated_reads
            .iter()
            .zip(self.data.total_reads.iter())
            .filter(|(_, &n)| n > 0)
            .map(|(&k, &n)| {
                let comb_numerator = ln_factorial(n as u64);
                let comb_denominator =
                    ln_factorial(k as u64) + ln_factorial((n - k) as u64);
                let k = k as f64;
                let n = n as f64;

                comb_numerator - comb_denominator
                    + ln_gamma(k + alpha)
                    + ln_gamma(n - k + beta)
                    - ln_gamma(n + alpha + beta)
                    - ln_beta_ab
            })
            .sum::<f64>();

        Ok(-log_likelihood)
    }
}

impl Gradient for BetaBinomialMLE {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    fn gradient(
        &self,
        params: &Self::Param,
    ) -> Result<Self::Gradient, Error> {
        let mu = params[0];
        let phi = params[1];

        // Convert to alpha and beta
        let alpha = mu * (1.0 - phi) / phi;
        let beta = (1.0 - mu) * (1.0 - phi) / phi;

        let dig_ab = digamma(alpha + beta);
        let dig_a = digamma(alpha);
        let dig_b = digamma(beta);

        let (d_mu_nscaled, d_phi_nscaled) = self
            .data
            .methylated_reads
            .iter()
            .zip(self.data.total_reads.iter())
            .filter(|(_, &n)| n > 0)
            .map(|(&k, &n)| {
                let k = k as f64;
                let n = n as f64;
                let dig_k_alpha = digamma(k + alpha);
                let dig_n_k_beta = digamma(n - k + beta);
                let dig_n_alpha_beta = digamma(n + alpha + beta);

                // Accumulate partial derivative with respect to mu
                let dmu = dig_k_alpha - dig_n_k_beta - dig_a + dig_b;
                let dphi = mu * (dig_a - dig_k_alpha)
                    + (1.0 - mu) * (dig_b - dig_n_k_beta)
                    + dig_n_alpha_beta
                    - dig_ab;

                (dmu, dphi)
            })
            .fold((0.0, 0.0), |(acc_a, acc_b), (da, db)| {
                (acc_a + da, acc_b + db)
            });
        let res_mu = -d_mu_nscaled * (1.0 - phi) / phi;
        let res_phi = -d_phi_nscaled / phi.powi(2);
        Ok(vec![res_mu, res_phi])
    }
}

// Estimates the parameters of a Beta-Binomial distribution using LBFGS
// optimization.
//
// Returns a tuple of (mu, phi) parameters where:
// - mu: Mean methylation rate (between 0 and 1)
// - phi: Overdispersion parameter (between 0 and 1)
// pub fn estimate_beta_binomial_params(
//     methylated_reads: Vec<u32>,
//     total_reads: Vec<u32>,
//     confidence: f64,
//     error: f64,
// ) -> Result<(f64, f64), Box<dyn std::error::Error>> {
//     let problem = BetaBinomialMLE::new(methylated_reads, total_reads);
//     let initial_params = problem.estimate_initial_params(confidence, error);
//
//     let linesearch = MoreThuenteLineSearch::new();
//     let solver = LBFGS::new(linesearch, 10)
//         .with_tolerance_grad(1e-12)?
//         .with_tolerance_cost(1e-12)?;
//
//     let res = Executor::new(problem, solver)
//         .configure(|state| {
//             state
//                 .param(initial_params)
//                 .max_iters(1000)
//         })
//         .run()?;
//
//     let best_params = res
//         .state()
//         .best_param
//         .clone()
//         .ok_or("No parameters found")?;
//
//     Ok((best_params[0], best_params[1]))
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::tools::dimred::tests::generate_beta_binomial_samples;
//
//     #[test]
//     fn test_parameter_estimation() -> Result<(), Box<dyn std::error::Error>>
// {         // True parameters for simulation
//         let true_alpha = 1.5;
//         let true_beta = 4.0;
//         let n_trials = 20;
//         let n_samples = 10;
//
//         let (successes, totals) = generate_beta_binomial_samples(
//             true_alpha, true_beta, n_trials, n_samples,
//         );
//
//         let (mu, phi) =
//             estimate_beta_binomial_params(successes, totals, 0.95, 0.05)?;
//
//         // Convert back to alpha and beta for comparison
//         let alpha = mu * (1.0 - phi) / phi;
//         let beta = (1.0 - mu) * (1.0 - phi) / phi;
//
//         println!("\nParameter Estimation Results:");
//         println!(
//             "True parameters: α = {:.4}, β = {:.4}",
//             true_alpha, true_beta
//         );
//         println!("Estimated parameters: α = {:.4}, β = {:.4}", alpha, beta);
//         println!("In transformed space: μ = {:.4}, φ = {:.4}", mu, phi);
//
//         // Basic assertion to ensure estimates are reasonably close
//         assert!(
//             (alpha - true_alpha).abs() < 1.0,
//             "Alpha estimate too far from true value"
//         );
//         assert!(
//             (beta - true_beta).abs() < 1.0,
//             "Beta estimate too far from true value"
//         );
//
//         Ok(())
//     }
// }
