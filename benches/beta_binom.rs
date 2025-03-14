/// ****************************************************************************
/// * Copyright (c) 2025
/// The Prosperity Public License 3.0.0
///
/// Contributor: [shitohana](https://github.com/shitohana)
///
/// Source Code: https://github.com/shitohana/BSXplorer
/// ***************************************************************************

/// ***********************************************************************
/// *****
/// * Copyright (c) 2025
/// ***********************************************************************
/// ****
use core::f64;

use bsxplorer2::tools::dimred::mle::estimate_beta_binomial_params;
use bsxplorer2::tools::dimred::BetaBinomParams;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use num::integer::sqrt;
use num::FromPrimitive;
use statrs::statistics::Statistics;

struct ValueCollector<T> {
    values: Vec<T>,
}

impl<T: Clone + num::FromPrimitive + num::ToPrimitive> ValueCollector<T> {
    fn new() -> Self { ValueCollector { values: Vec::new() } }

    fn bench_function<F>(
        &mut self,
        f: F,
    ) where
        F: FnMut() -> T, {
        let mut func = f;
        let result = func();
        self.values.push(result);
    }

    fn mean(&self) -> f64 {
        self.values
            .iter()
            .map(|val| val.to_f64().unwrap())
            .mean()
    }
}

fn generate_beta_binomial_samples(
    alpha: f64,
    beta: f64,
    n_trials: u32,
    n_samples: usize,
) -> (Vec<u32>, Vec<u32>) {
    use rand::prelude::*;
    use rand_distr::{Beta, Binomial};

    let mut rng = rand::thread_rng();
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

pub fn benchmark_from_counts(c: &mut Criterion) {
    // Setup multiple dataset sizes to measure performance scaling
    let sizes = [100, 1000, 10000];

    for size in sizes.iter() {
        let group_name = format!("from_counts_size_{}", size);
        let mut group = c.benchmark_group(group_name);

        // Generate dataset of appropriate size
        let true_alpha = 20.0;
        let true_beta = 10.0;
        let n_trials = 50u32;
        let (successes, totals) = generate_beta_binomial_samples(
            true_alpha, true_beta, n_trials, *size,
        );

        // Benchmark with different confidence/precision combinations
        let confidence_precision_pairs = [
            (0.90, 0.1),  // Lower confidence, lower precision
            (0.95, 0.05), // Medium confidence, medium precision
            (0.99, 0.01), // High confidence, high precision
        ];

        for (confidence, precision) in confidence_precision_pairs.iter() {
            let bench_name = format!("conf_{}_prec_{}", confidence, precision);

            group.bench_function(bench_name, |b| {
                b.iter(|| {
                    BetaBinomParams::<f64>::from_counts(
                        black_box(&successes),
                        black_box(&totals),
                        black_box(*confidence),
                        black_box(*precision),
                    )
                });
            });
        }

        group.finish();
    }
}
pub fn benchmark_parameter_estimation(c: &mut Criterion) {
    // Define sample sizes from 10 to 1000
    let sample_sizes = [10, 20, 50, 100, 200, 500, 1000];
    let true_alpha = 2.0;
    let true_beta = 5.0;
    let n_trials = 30;

    let mut group = c.benchmark_group("parameter_estimation");

    for &size in &sample_sizes {
        group.bench_function(format!("sample_size_{}", size), |b| {
            b.iter(|| {
                let (successes, totals) = generate_beta_binomial_samples(
                    true_alpha, true_beta, n_trials, size,
                );

                let result = estimate_beta_binomial_params(
                    black_box(successes),
                    black_box(totals),
                    black_box(0.95),
                    black_box(0.05),
                )
                .unwrap();

                result
            });
        });

        // Calculate distance from true CDF to estimated CDF
        group.bench_function(format!("cdf_distance_size_{}", size), |b| {
            b.iter(|| {
                let (successes, totals) = generate_beta_binomial_samples(
                    true_alpha, true_beta, n_trials, size,
                );

                let (mu, phi) = estimate_beta_binomial_params(
                    black_box(successes),
                    black_box(totals),
                    black_box(0.95),
                    black_box(0.05),
                )
                .unwrap();

                let est_alpha =
                    (mu * (1.0 - phi) / phi).max(f64::EPSILON.sqrt());
                let est_beta =
                    ((1.0 - mu) * (1.0 - phi) / phi).max(f64::EPSILON.sqrt());

                // Calculate CDF distance (L1 norm over all possible outcomes)
                let mut distance = 0.0;
                for k in 0..=n_trials {
                    // Calculate true vs estimated probability mass at k
                    let true_pmf = true_beta_binomial_pmf(
                        k, n_trials, true_alpha, true_beta,
                    );
                    let est_pmf = true_beta_binomial_pmf(
                        k, n_trials, est_alpha, est_beta,
                    );
                    distance += (true_pmf - est_pmf).abs();
                }

                distance
            });
        });
    }
    group.finish();
}

// Helper function to calculate Beta-Binomial PMF
fn true_beta_binomial_pmf(
    k: u32,
    n: u32,
    alpha: f64,
    beta: f64,
) -> f64 {
    use statrs::function::beta;

    let nck = choose(n, k) as f64;
    let beta_ratio = beta::beta(
        (k as f64 + alpha).max(f64::EPSILON.sqrt()),
        (n as f64 - k as f64 + beta).max(f64::EPSILON.sqrt()),
    ) / beta::beta(alpha, beta);

    nck * beta_ratio
}

// Binomial coefficient calculation
fn choose(
    n: u32,
    k: u32,
) -> u64 {
    if k > n {
        return 0;
    }

    let k = k.min(n - k); // Take advantage of symmetry
    let mut result = 1u64;

    for i in 0..k {
        result = result * (n - i) as u64 / (i + 1) as u64;
    }

    result
}

criterion_group!(
    benches,
    benchmark_parameter_estimation,
    benchmark_from_counts
);
criterion_main!(benches);
