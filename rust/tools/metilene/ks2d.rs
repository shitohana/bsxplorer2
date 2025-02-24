use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::distribution::{ContinuousCDF, Normal};
use statrs::statistics::Statistics;

/// ks stands for Kolmogorov-Smirnov, 2d for 2 dimensional,
/// 2s for 2 samples.
///
/// KS test for goodness-of-fit on two 2D samples. Tests the hypothesis that
/// the two samples are from the same distribution.
///
/// ### Returns
/// A tuple of two floats. First, the two-sample K-S statistic.
///  If this value is higher than the significance level of the hypothesis,
///  it is rejected. Second, the significance level of *d*. Small values of
///  prob show that the two samples are significantly different.
pub fn ks2d_2sample<X, Y>(x1: &[X], y1: &[Y], x2: &[X], y2: &[Y]) -> (f64, f64)
where
    X: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps
        + Sync
        + Send,
    Y: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps
        + Sync
        + Send,
{
    assert_eq!(x1.len(), y1.len());
    assert_eq!(x2.len(), y2.len());
    let d1 = (0..x1.len())
        .into_par_iter()
        .map(|i| {
            let (fpp1, fmp1, fpm1, fmm1) = count_quads_generic(x1, y1, x1[i], y1[i]);
            let (fpp2, fmp2, fpm2, fmm2) = count_quads_generic(x2, y2, x1[i], y1[i]);
            vec![
                (fpp1 - fpp2).abs(),
                (fmp1 - fmp2).abs(),
                (fpm1 - fpm2).abs(),
                (fmm1 - fmm2).abs(),
            ]
            .max()
        })
        .max_by(|a, b| a.total_cmp(&b))
        .unwrap_or(0f64);

    let d2 = (0..x2.len())
        .into_par_iter()
        .map(|i| {
            let (fpp1, fmp1, fpm1, fmm1) = count_quads_generic(x1, y1, x2[i], y2[i]);
            let (fpp2, fmp2, fpm2, fmm2) = count_quads_generic(x2, y2, x2[i], y2[i]);
            vec![
                (fpp1 - fpp2).abs(),
                (fmp1 - fmp2).abs(),
                (fpm1 - fpm2).abs(),
                (fmm1 - fmm2).abs(),
            ]
            .max()
        })
        .max_by(|a, b| a.total_cmp(&b))
        .unwrap_or(0f64);

    let d = (d1 + d2) / 2.0;
    let sqen = ((x1.len() * x2.len()) as f64 / (x1.len() + x2.len()) as f64).sqrt();
    let r1 = pearson_r(x1, y1);
    let r2 = pearson_r(x2, y2);
    let rr = (1.0 - (r1 * r1 + r2 * r2) / 2.0).sqrt();
    let prob = ks_prob(
        d * sqen / (1.0 + rr * (0.25 - 0.75 / sqen)),
        Default::default(),
        Default::default(),
    );
    (d, prob)
}

/// Calculates pearson correlation coefficient for set of points
pub fn pearson_r<X, Y>(x: &[X], y: &[Y]) -> f64
where
    X: num::ToPrimitive + num::FromPrimitive + PartialOrd + Copy + Clone + num::traits::NumOps,
    Y: num::ToPrimitive + num::FromPrimitive + PartialOrd + Copy + Clone + num::traits::NumOps,
{
    assert_eq!(x.len(), y.len());

    let x_f64 = x.iter().map(|x| x.to_f64().unwrap()).collect::<Vec<_>>();
    let y_f64 = y.iter().map(|x| x.to_f64().unwrap()).collect::<Vec<_>>();

    let x_mean = x_f64.iter().mean();
    let y_mean = y_f64.iter().mean();

    let numerator = x_f64
        .iter()
        .zip(y_f64.iter())
        .map(|(valx, valy)| (valx - x_mean) * (valy - y_mean))
        .sum::<f64>();
    let denominator = {
        let x_dev: f64 = x_f64.iter().map(|valx| (valx - x_mean).powi(2)).sum();
        let y_dev: f64 = y_f64.iter().map(|valy| (valy - y_mean).powi(2)).sum();
        (x_dev * y_dev).sqrt()
    };
    if denominator == 0.0 {
        return 0.0;
    }
    numerator / denominator
}

/// Computes the probabilities of finding points in each 4 quadrant
/// defined by a vertical and horizontal lines crossing the point, by counting
/// the proportion of points in Arr2D in each quadrant.
///
/// ### Returns
/// A tuple of 4 floats.  The probabilities of finding a point in
/// each quadrants, with point as the origin.  p stands for positive, n for
/// negative, with the first and second positions meaning the x and y
/// directions respectively.
fn count_quads_generic<X, Y>(x: &[X], y: &[Y], x0: X, y0: Y) -> (f64, f64, f64, f64)
where
    X: num::ToPrimitive + num::FromPrimitive + PartialOrd + Copy + Clone + num::traits::NumOps,
    Y: num::ToPrimitive + num::FromPrimitive + PartialOrd + Copy + Clone + num::traits::NumOps,
{
    assert_eq!(x.len(), y.len());
    let points_len = x.len();
    let ff = 1.0 / points_len as f64;
    let fpp = (0..points_len)
        .into_iter()
        .filter(|cur_idx| (x[*cur_idx] > x0) && (y[*cur_idx] > y0))
        .count() as f64
        * ff;
    let fnp = (0..points_len)
        .into_iter()
        .filter(|cur_idx| (x[*cur_idx] < x0) && (y[*cur_idx] > y0))
        .count() as f64
        * ff;
    let fpn = (0..points_len)
        .into_iter()
        .filter(|cur_idx| (x[*cur_idx] > x0) && (y[*cur_idx] < y0))
        .count() as f64
        * ff;
    let fnn = (0..points_len)
        .into_iter()
        .filter(|cur_idx| (x[*cur_idx] < x0) && (y[*cur_idx] < y0))
        .count() as f64
        * ff;
    (fpp, fnp, fpn, fnn)
}

/// Computes the value of the KS probability function, as a function of
/// alam, the D statistic. From *Numerical recipes in C* page 623:
/// the K–S statistic useful is that its distribution in the case of the null
/// hypothesis (data sets drawn from the same distribution) can be calculated,
/// at least to useful approximation, thus giving the significance of any
/// observed nonzero value of D.' (D being the KS statistic).
///
/// - alam: D statistic.
/// - iter: Number of iterations to be perfomed. On non-convergence, returns 1.0.
/// - prec: Convergence criteria of the qks. Stops converging if that precision
///         is attained.
///
/// ### Returns
///
/// A float. The significance level of the observed D statistic.
fn ks_prob(alam: f64, iter: Option<usize>, prec: Option<f64>) -> f64 {
    let mut toadd: Vec<f64> = vec![1.0];
    let mut qks = 0f64;
    let mut j = 1;
    let iter = iter.unwrap_or(100);
    let prec = prec.unwrap_or(1e-17);

    while (j < iter) && ((toadd.last().unwrap_or(&f64::MAX)).abs() > prec * 2.0) {
        let new_elem = 2.0
            * (-1.0f64).powi((j as i32 - 1) as i32)
            * (-2.0 * j.pow(2) as f64 * alam.powi(2)).exp();
        qks += new_elem;
        toadd.push(new_elem);
        j += 1
    }
    if (j == iter) || (qks > 1.0) {
        // If no convergence after j iter, return 1.0
        return 1.0;
    }
    if qks < prec {
        0.0
    } else {
        qks
    }
}

#[derive(Debug)]
struct Observation {
    value: f64,
    group: usize, // 0 for group1, 1 for group2
    rank: f64,
}

/// Implements the Mann–Whitney U test for two independent samples of f64 values.
/// Returns a tuple: (U statistic, two–tailed p–value).
pub fn mann_whitney_u(group1: &[f64], group2: &[f64]) -> (f64, f64) {
    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;
    let n_total = (group1.len() + group2.len()) as f64;

    // Combine observations from both groups.
    let mut observations: Vec<Observation> = group1
        .iter()
        .map(|&v| Observation {
            value: v,
            group: 0,
            rank: 0.0,
        })
        .chain(group2.iter().map(|&v| Observation {
            value: v,
            group: 1,
            rank: 0.0,
        }))
        .collect();

    // Sort observations by value.
    observations.sort_by(|a, b| a.value.partial_cmp(&b.value).unwrap());

    // Assign ranks, averaging in the case of ties.
    let mut tie_groups: Vec<usize> = Vec::new();
    let mut i = 0;
    while i < observations.len() {
        let start = i;
        let mut end = i + 1;
        while end < observations.len()
            && (observations[end].value - observations[start].value).abs() < 1e-6
        {
            end += 1;
        }
        let count = end - start;
        let avg_rank = (start as f64 + 1.0 + end as f64) / 2.0; // ranks are 1-indexed
        for j in start..end {
            observations[j].rank = avg_rank;
        }
        if count > 1 {
            tie_groups.push(count);
        }
        i = end;
    }

    // Compute the sum of ranks for group1.
    let r1: f64 = observations
        .iter()
        .filter(|obs| obs.group == 0)
        .map(|obs| obs.rank)
        .sum();

    // Calculate U statistics.
    let u1 = r1 - n1 * (n1 + 1.0) / 2.0;
    let u2 = n1 * n2 - u1;
    // Use the smaller U value for the test.
    let u_stat = u1.min(u2);

    // Mean of U under H0.
    let mean_u = n1 * n2 / 2.0;

    // Variance with tie correction.
    let tie_sum: f64 = tie_groups.iter().map(|&t| (t * t * t - t) as f64).sum();
    let variance_u = n1 * n2 / 12.0 * ((n_total + 1.0) - tie_sum / (n_total * (n_total - 1.0)));

    // Apply continuity correction.
    // Since we take the minimum U, it is always <= mean_u.
    let z = if variance_u > 0.0 {
        (u_stat - mean_u + 0.5) / variance_u.sqrt()
    } else {
        0.0
    };

    // Use statrs to calculate the two-tailed p–value.
    let normal = Normal::new(0.0, 1.0).unwrap();
    let p_value = 2.0 * normal.cdf(-z.abs());

    (u_stat, p_value.min(1.0)) // ensure the p-value does not exceed 1.0
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_utest() {
        let group1 = vec![1.5, 2.3, 3.1, 4.8, 5.7, 5.6];
        let group2 = vec![2.0, 3.5, 3.8, 4.0, 6.2, 3.5];

        let (u, p) = mann_whitney_u(&group1, &group2);

        println!("Mann–Whitney U statistic: {:.3}", u);
        println!("Two-tailed p–value: {:.5}", p);
    }
}
