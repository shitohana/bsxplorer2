/// ****************************************************************************
/// * Copyright (c) 2025
/// ***************************************************************************

/// ***********************************************************************
/// *****
/// * Copyright (c) 2025
/// ***********************************************************************
/// ****
use core::f64;
use std::cell::RefCell;
use std::path::PathBuf;
use std::time::Instant;

use bsxplorer2::tools::dimred::mle::estimate_beta_binomial_params;
use bsxplorer2::tools::dimred::BetaBinomParams;
use itertools::{multiunzip, Itertools};
use plotters::prelude::*;
use plotters::style::Color;
use quickcheck::TestResult;
use rand::prelude::*;

thread_local! {
    static RESULTS: RefCell<Vec<(f64, f64, f64, f64, u128, u128)>> = RefCell::new(Vec::new());
    static RNG: RefCell<ThreadRng> = RefCell::new(rand::thread_rng());
}
const N_TRIALS: u32 = 30;
const SAMPLE_SIZE: u32 = 100;

fn test_with_timing<F, I, O>(
    mut function: F,
    input: I,
) -> (O, u128)
where
    F: FnMut(I) -> O, {
    let start = Instant::now();
    let result = function(input);
    let elapsed = start.elapsed().as_nanos();
    (result, elapsed)
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

fn test_fun(input: (f64, f64)) -> (f64, f64, u128, u128) {
    let (mu, phi) = input;
    let true_alpha = mu * (1.0 - phi) / phi;
    let true_beta = (1.0 - mu) * (1.0 - phi) / phi;
    let (successes, totals) = generate_beta_binomial_samples(
        true_alpha,
        true_beta,
        N_TRIALS,
        SAMPLE_SIZE as usize,
    );

    let start = Instant::now();
    let (params, _) =
        BetaBinomParams::from_counts(&successes, &totals, 0.95, 0.05);
    let elapsed_mom = start.elapsed().as_nanos();
    let start = Instant::now();
    let (mu, phi) =
        estimate_beta_binomial_params(successes, totals, 0.95, 0.05).unwrap();
    let elapsed_mle = start.elapsed().as_nanos();

    let est_alpha = (mu * (1.0 - phi) / phi).max(f64::EPSILON.sqrt());
    let est_beta = ((1.0 - mu) * (1.0 - phi) / phi).max(f64::EPSILON.sqrt());

    // Calculate CDF distance (L1 norm over all possible outcomes)
    let mut mle_distance = 0.0;
    for k in 0..=N_TRIALS {
        // Calculate true vs estimated probability mass at k
        let true_pmf =
            true_beta_binomial_pmf(k, N_TRIALS, true_alpha, true_beta);
        let est_pmf = true_beta_binomial_pmf(k, N_TRIALS, est_alpha, est_beta);
        mle_distance += (true_pmf - est_pmf).abs();
    }
    let est_alpha = params.0;
    let est_beta = params.1;

    let mut mom_distance = 0.0;
    for k in 0..=N_TRIALS {
        // Calculate true vs estimated probability mass at k
        let true_pmf =
            true_beta_binomial_pmf(k, N_TRIALS, true_alpha, true_beta);
        let est_pmf = true_beta_binomial_pmf(k, N_TRIALS, est_alpha, est_beta);
        mom_distance += (true_pmf - est_pmf).abs();
    }
    (mom_distance, mle_distance, elapsed_mom, elapsed_mle)
}

fn test_property() -> TestResult {
    let mut mu = 0.0;
    let mut alpha = 0.0;

    RNG.with(|rng| {
        let mut rng = rng.borrow_mut();
        mu = rng.gen_range(f64::EPSILON.sqrt()..(1.0 - f64::EPSILON.sqrt()));
        alpha = rng.gen_range((1.0 - mu - f64::EPSILON.sqrt())..200.0);
    });

    let true_alpha = alpha;
    let true_beta = alpha / mu - alpha;
    let phi = 1.0 / (true_alpha + true_beta + 1.0);

    let (result, _) = test_with_timing(test_fun, (mu, phi));
    let (mom, mle, elapsed_mom, elapsed_mle) = result;

    // Store results in thread-local storage
    RESULTS.with(|results| {
        results.borrow_mut().push((
            mu,
            phi,
            mom,
            mle,
            elapsed_mom,
            elapsed_mle,
        ));
    });

    TestResult::from_bool(true)
}

const MAX_DIGITS: usize = 10;
const N_TESTS: usize = 101;
#[test]
fn main() {
    for i in 0..N_TESTS {
        test_property();
    }
    let output_dir = PathBuf::from(
        &(env!("CARGO_MANIFEST_DIR").to_owned() + "/target/tests"),
    );
    std::fs::create_dir(output_dir.clone()).unwrap_or(());
    let output_file = output_dir.join(format!(
        "{}.svg",
        std::env::current_exe()
            .ok()
            .and_then(|p| {
                p.file_name()
                    .map(|n| n.to_string_lossy().to_string())
            })
            .unwrap_or_else(|| "unknown".to_string())
    ));

    RESULTS.with(|results| {
        let results = results.borrow().clone();
        let (_, _, mom, mle, elapsed_mom, elapsed_mle): (
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
            Vec<_>,
        ) = multiunzip(results.into_iter());

        let root =
            SVGBackend::new(&output_file, (1024, 768)).into_drawing_area();

        // Convert elapsed to ms for better readability
        let elapsed_mom_ms: Vec<f64> = elapsed_mom
            .iter()
            .map(|&t| t as f64 / 1_000_000.0)
            .collect();
        let elapsed_mle_ms: Vec<f64> = elapsed_mle
            .iter()
            .map(|&t| t as f64 / 1_000_000.0)
            .collect();

        // Create a layout with 2 plots
        let areas = root.split_evenly((1, 2));

        // First plot: Accuracy comparison (mom vs mle)
        {
            let mut chart = ChartBuilder::on(&areas[0])
                .caption("Accuracy Comparison", ("sans-serif", 20).into_font())
                .margin(5)
                .x_label_area_size(40)
                .y_label_area_size(60)
                .build_cartesian_2d(
                    0.0..mom.len() as f64,
                    0.0..mom
                        .iter()
                        .chain(mle.iter())
                        .fold(0.0f64, |a, &b| a.max(b))
                        * 1.1,
                )
                .unwrap();

            chart
                .configure_mesh()
                .x_desc("Sample Index")
                .y_desc("Distance (L1 norm)")
                .draw()
                .unwrap();

            chart
                .draw_series(LineSeries::new(
                    mom.iter()
                        .enumerate()
                        .map(|(i, &v)| (i as f64, v)),
                    &RED,
                ))
                .unwrap()
                .label("Method of Moments")
                .legend(|(x, y)| {
                    PathElement::new(vec![(x, y), (x + 20, y)], &RED)
                });

            chart
                .draw_series(LineSeries::new(
                    mle.iter()
                        .enumerate()
                        .map(|(i, &v)| (i as f64, v)),
                    &BLUE,
                ))
                .unwrap()
                .label("Maximum Likelihood")
                .legend(|(x, y)| {
                    PathElement::new(vec![(x, y), (x + 20, y)], &BLUE)
                });

            chart
                .configure_series_labels()
                .background_style(&WHITE.mix(0.8))
                .border_style(&BLACK)
                .draw()
                .unwrap();
        }

        // Second plot: Performance comparison (elapsed_mom vs elapsed_mle)
        {
            let mut chart = ChartBuilder::on(&areas[1])
                .caption(
                    "Performance Comparison",
                    ("sans-serif", 20).into_font(),
                )
                .margin(5)
                .x_label_area_size(40)
                .y_label_area_size(60)
                .build_cartesian_2d(
                    0.0..elapsed_mom_ms.len() as f64,
                    0.0..elapsed_mom_ms
                        .iter()
                        .chain(elapsed_mle_ms.iter())
                        .fold(0.0f64, |a, &b| a.max(b))
                        * 1.1,
                )
                .unwrap();

            chart
                .configure_mesh()
                .x_desc("Sample Index")
                .y_desc("Time (ms)")
                .draw()
                .unwrap();

            chart
                .draw_series(LineSeries::new(
                    elapsed_mom_ms
                        .iter()
                        .enumerate()
                        .map(|(i, &v)| (i as f64, v)),
                    &RED,
                ))
                .unwrap()
                .label("Method of Moments")
                .legend(|(x, y)| {
                    PathElement::new(vec![(x, y), (x + 20, y)], &RED)
                });

            chart
                .draw_series(LineSeries::new(
                    elapsed_mle_ms
                        .iter()
                        .enumerate()
                        .map(|(i, &v)| (i as f64, v)),
                    &BLUE,
                ))
                .unwrap()
                .label("Maximum Likelihood")
                .legend(|(x, y)| {
                    PathElement::new(vec![(x, y), (x + 20, y)], &BLUE)
                });

            chart
                .configure_series_labels()
                .background_style(&WHITE.mix(0.8))
                .border_style(&BLACK)
                .draw()
                .unwrap();
        }

        println!("Generated plot at {}", output_file.display());
        // root.fill(&WHITE).unwrap();
    })
}
