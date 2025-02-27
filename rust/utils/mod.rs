use crate::data_structs::region::GenomicPosition;
use itertools::Itertools;
use polars::datatypes::PlIndexMap;
use polars::prelude::*;
use std::sync::Arc;

pub mod types;

#[cfg(feature = "python")]
mod python {
    #[macro_export]
    macro_rules! wrap_polars_result {
        ($expression: expr) => {{
            match $expression {
                Ok(v) => Ok(v),
                Err(e) => Err(PyPolarsErr::from(e).into()),
            }
        }};
    }

    #[macro_export]
    macro_rules! wrap_box_result {
        ($error: ty, $expression: expr) => {{
            match $expression {
                Ok(v) => Ok(v),
                Err(e) => Err(<$error>::new_err(e.to_string())),
            }
        }};
    }
    pub use {wrap_box_result, wrap_polars_result};
}
#[cfg(feature = "python")]
pub(crate) use python::{wrap_box_result, wrap_polars_result};

macro_rules! polars_schema {
    ( $($name: expr => $dtype: expr),* ) => {
        {
            let mut fields: Vec<(PlSmallStr, DataType)> = Vec::new();
            $(
                fields.push(($name.into(), $dtype));
            )*
            Schema::from_iter(fields)
        }
    };
}
use num::ToPrimitive;
pub(crate) use polars_schema;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::distribution::{ContinuousCDF, Normal};
use statrs::statistics::Statistics;

pub fn array_to_schema(array: &[(&str, DataType)]) -> Schema {
    Schema::from(PlIndexMap::from_iter(
        array.iter().cloned().map(|(k, v)| (PlSmallStr::from(k), v)),
    ))
}

pub fn get_categorical_dtype(categories: Vec<String>) -> DataType {
    let categories = polars::export::arrow::array::Utf8ViewArray::from_vec(
        categories.iter().map(String::as_str).collect_vec(),
        ArrowDataType::Utf8View,
    );
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Enum(Some(rev_mapping), CategoricalOrdering::Physical)
}

pub(crate) fn schema_from_arrays(names: &[&str], dtypes: &[DataType]) -> Schema {
    Schema::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
}

pub(crate) fn hashmap_from_arrays<'a>(
    names: &[&'a str],
    dtypes: &[DataType],
) -> PlHashMap<&'a str, DataType> {
    PlHashMap::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
}

pub(crate) fn first_position(
    data: &DataFrame,
    chr_col: &str,
    pos_col: &str,
) -> PolarsResult<GenomicPosition<u64>> {
    let chr = data
        .column(chr_col)?
        .as_series()
        .unwrap()
        .first()
        .as_any_value()
        .cast(&DataType::String)
        .str_value()
        .to_string();
    let pos = data
        .column(pos_col)?
        .cast(&DataType::UInt64)?
        .u64()?
        .first()
        .unwrap();
    Ok(GenomicPosition::new(chr, pos))
}

pub(crate) fn last_position(
    data: &DataFrame,
    chr_col: &str,
    pos_col: &str,
) -> PolarsResult<GenomicPosition<u64>> {
    let chr = data
        .column(chr_col)?
        .as_series()
        .unwrap()
        .last()
        .as_any_value()
        .cast(&DataType::String)
        .str_value()
        .to_string();
    let pos = data
        .column(pos_col)?
        .cast(&DataType::UInt64)?
        .u64()?
        .last()
        .unwrap_or(0);
    Ok(GenomicPosition::new(chr, pos))
}

pub fn encode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(strand_col).eq(lit("+")))
            .then(lit(true))
            .when(col(strand_col).eq(lit("-")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias("strand"),
    )
}

pub fn encode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(context_col).eq(lit("CG")))
            .then(lit(true))
            .when(col(context_col).eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias(context_col),
    )
}

pub fn decode_strand(lazy_frame: LazyFrame, strand_col: &str, result_name: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(strand_col).eq(lit(true)))
            .then(lit("+"))
            .when(col(strand_col).eq(lit(false)))
            .then(lit("-"))
            .otherwise(lit("."))
            .cast(DataType::String)
            .alias(result_name),
    )
}

pub fn decode_context(lazy_frame: LazyFrame, context_col: &str, result_name: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(context_col).eq(lit(true)))
            .then(lit("CG"))
            .when(col(context_col).eq(lit(false)))
            .then(lit("CHG"))
            .otherwise(lit("CHH"))
            .cast(DataType::String)
            .alias(result_name),
    )
}

pub(crate) fn f64_to_u64_scaled(value: f64) -> u64 {
    (value * u64::MAX as f64) as u64
}

pub(crate) fn u64_to_f64_scaled(value: u64) -> f64 {
    value as f64 / u64::MAX as f64
}
// < 0.5% error
pub(crate) const GROUPING_POWER: u8 = 8;

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
    use crate::utils::{ks2d_2sample, mann_whitney_u, pearson_r};
    use assert_approx_eq::assert_approx_eq;
    use itertools::Itertools;
    use rand::distributions::Distribution;
    use rand::thread_rng;

    #[test]
    fn test_utest() {
        let group1 = vec![1.5, 2.3, 3.1, 4.8, 5.7, 5.6];
        let group2 = vec![2.0, 3.5, 3.8, 4.0, 6.2, 3.5];

        let (u, p) = mann_whitney_u(&group1, &group2);
        let (uleft, pleft) = mann_whitney_u(&group2, &group1);
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

    #[test]
    fn ks_test() {
        let x1 = vec![1, 2, 3, 4, 5, 6];
        let y1 = vec![1, 4, 9, 16, 25, 36];

        let x2 = vec![1, 2, 3, 4, 5, 6];
        let y2 = vec![1, 4, 9, 16, 25, 36];
        let (prob, _) = ks2d_2sample(&x1, &y1, &x2, &y2);
        assert_approx_eq!(prob, 0f64);

        let x1 = vec![1, 2, 3, 4, 5, 6];
        let y1 = statrs::distribution::Normal::standard()
            .sample_iter(thread_rng())
            .take(x1.len())
            .collect_vec();

        let x2 = vec![1, 2, 3, 4, 5, 6];
        let y2 = statrs::distribution::Beta::new(1f64, 3f64)
            .unwrap()
            .sample_iter(thread_rng())
            .take(x2.len())
            .collect_vec();
        let (prob, d) = ks2d_2sample(&x1, &y1, &x2, &y2);
        let (prob_reverse, d) = ks2d_2sample(&x2, &y2, &x1, &y1);
        assert_approx_eq!(prob, prob_reverse)
    }
}
