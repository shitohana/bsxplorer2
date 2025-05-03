use std::fmt;

use log::*;
use num::Float;
use statrs::distribution::{ContinuousCDF, Normal};
use statrs::statistics::Statistics;

/// Calculates Pearson correlation coefficient between two variables.
pub fn pearson_r<X, Y>(
    x: &[X],
    y: &[Y],
) -> f64
where
    X: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps
        + fmt::Debug,
    Y: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps
        + fmt::Debug, {
    if x.len() != y.len() {
        warn!(
            "Cannot calculate Pearson's r: x length ({}) doesn't match y \
             length ({})",
            x.len(),
            y.len()
        );
        return 0.0;
    }

    if x.is_empty() {
        warn!("Cannot calculate Pearson's r: empty arrays");
        return 0.0;
    }

    let x_f64 = x
        .iter()
        .map(|x| x.to_f64().unwrap())
        .collect::<Vec<_>>();
    let y_f64 = y
        .iter()
        .map(|y| y.to_f64().unwrap())
        .collect::<Vec<_>>();

    let x_mean = x_f64.iter().mean();
    let y_mean = y_f64.iter().mean();

    // Calculate numerator (covariance)
    let numerator = x_f64
        .iter()
        .zip(y_f64.iter())
        .map(|(valx, valy)| (valx - x_mean) * (valy - y_mean))
        .sum::<f64>();

    // Calculate denominator (product of standard deviations)
    let denominator = {
        let x_dev: f64 = x_f64
            .iter()
            .map(|valx| (valx - x_mean).powi(2))
            .sum();
        let y_dev: f64 = y_f64
            .iter()
            .map(|valy| (valy - y_mean).powi(2))
            .sum();
        (x_dev * y_dev).sqrt()
    };

    if denominator == 0.0 {
        debug!("Denominator is zero, returning r=0");
        return 0.0;
    }

    let r = numerator / denominator;
    debug!("Pearson's r = {:.4}", r);
    r
}

/// Represents an observation in the Mann-Whitney U test
#[derive(Debug)]
struct Observation<F: Float> {
    /// The observed value
    value: F,
    /// Which group the observation belongs to (0 for group1, 1 for group2)
    group: usize,
    /// The assigned rank of this observation
    rank:  f64,
}

/// Performs Mann-Whitney U test.
/// A non-parametric test for distribution differences.
pub fn mann_whitney_u<F: Float>(
    group1: &[F],
    group2: &[F],
) -> (f64, f64) {
    info!(
        "Performing Mann-Whitney U test: group1={}, group2={}",
        group1.len(),
        group2.len()
    );

    let n1 = F::from(group1.len()).unwrap();
    let n2 = F::from(group2.len()).unwrap();
    let n_total = n1 + n2;

    if group1.is_empty() || group2.is_empty() {
        warn!("Mann-Whitney U test: one or both groups are empty");
        return (0.0, 1.0);
    }

    // Combine observations from both groups
    let mut observations: Vec<Observation<F>> = group1
        .iter()
        .map(|&v| {
            Observation {
                value: v,
                group: 0,
                rank:  0.0,
            }
        })
        .chain(group2.iter().map(|&v| {
            Observation {
                value: v,
                group: 1,
                rank:  0.0,
            }
        }))
        .collect();

    // Sort observations by value
    observations.sort_by(|a, b| {
        a.value
            .partial_cmp(&b.value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Assign ranks, handling ties by averaging
    let mut tie_groups: Vec<usize> = Vec::new();
    let mut i = 0;
    while i < observations.len() {
        let start = i;
        let mut end = i + 1;

        // Find all tied values
        while end < observations.len()
            && (observations[end].value - observations[start].value).abs()
                < F::from(1e-6).unwrap()
        {
            end += 1;
        }

        let count = end - start;
        let avg_rank = (start as f64 + 1.0 + end as f64) / 2.0; // ranks are 1-indexed

        // Assign average rank to all tied values
        for val in observations
            .iter_mut()
            .take(end)
            .skip(start)
        {
            val.rank = avg_rank;
        }

        // Record ties for variance correction
        if count > 1 {
            tie_groups.push(count);
        }

        i = end;
    }

    // Compute the sum of ranks for group1
    let r1: F = F::from(
        observations
            .iter()
            .filter(|obs| obs.group == 0)
            .map(|obs| obs.rank)
            .sum::<f64>(),
    )
    .unwrap();

    // Calculate U statistics
    let u1 = r1 - n1 * (n1 + F::from(1.0).unwrap()) / F::from(2.0).unwrap();
    let u2 = n1 * n2 - u1;
    // Use the smaller U value for the test
    let u_stat = u1.min(u2);

    // Mean of U under H0
    let mean_u = n1 * n2 / F::from(2.0).unwrap();

    // Variance with tie correction
    let tie_sum = F::from(
        tie_groups
            .iter()
            .map(|&t| (t * t * t - t) as f64)
            .sum::<f64>(),
    )
    .unwrap();
    let variance_u = n1 * n2 / F::from(12.0).unwrap()
        * ((n_total + F::from(1.0).unwrap())
            - tie_sum / (n_total * (n_total - F::from(1.0).unwrap())));

    // Apply continuity correction and compute Z-score
    let z = if variance_u > F::from(0.0).unwrap() {
        (u_stat - mean_u + F::from(0.5).unwrap()) / variance_u.sqrt()
    }
    else {
        warn!("Variance is zero in Mann-Whitney U test");
        F::from(0.0).unwrap()
    };
    let z = z.to_f64().unwrap();

    // Calculate two-tailed p-value
    let normal = Normal::new(0.0, 1.0).unwrap();
    let p_value = 2.0 * normal.cdf(-z.abs());
    let p_value = p_value.min(1.0); // ensure p-value doesn't exceed 1

    (u_stat.to_f64().unwrap(), p_value)
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use crate::utils::{mann_whitney_u, pearson_r};

    #[test]
    fn test_utest() {
        let group1 = vec![1.5, 2.3, 3.1, 4.8, 5.7, 5.6];
        let group2 = vec![2.0, 3.5, 3.8, 4.0, 6.2, 3.5];

        let (u, p) = mann_whitney_u(&group1, &group2);
        let (_uleft, pleft) = mann_whitney_u(&group2, &group1);
        assert_approx_eq!(pleft, p);
        println!("Mann–Whitney U statistic: {:.3}", u);
        println!("Two-tailed p–value: {:.5}", p);
    }

    #[test]
    fn pearson_r_test() {
        let x = vec![1, 2, 3, 4, 5, 6];
        let y = vec![6, 7, 8, 9, 10, 11];
        assert_eq!(pearson_r(&x, &y), 1f64);
    }
}
