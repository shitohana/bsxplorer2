use std::fmt;
use std::sync::Arc;

use anyhow::Context;
use itertools::Itertools;
use log::{debug, info, warn};
use num::Float;
use polars::datatypes::PlIndexMap;
use polars::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::distribution::{ContinuousCDF, Normal};
use statrs::statistics::Statistics;

use crate::data_structs::region::GenomicPosition;

pub mod types;


/// Converts an array of name-datatype pairs to a Polars Schema
///
/// # Arguments
/// * `array` - Slice of tuples containing column name and data_structs type
///
/// # Returns
/// A Polars Schema constructed from the provided array
pub fn array_to_schema(array: &[(&str, DataType)]) -> Schema {
    debug!("Creating schema from array with {} elements", array.len());
    Schema::from(PlIndexMap::from_iter(
        array
            .iter()
            .cloned()
            .map(|(k, v)| (PlSmallStr::from(k), v)),
    ))
}

/// Creates a categorical data_structs type from a list of categories
///
/// # Arguments
/// * `categories` - Vector of strings representing the categories
///
/// # Returns
/// A DataType::Enum with the provided categories
pub fn get_categorical_dtype(categories: Vec<String>) -> DataType {
    debug!(
        "Creating categorical data_structs type with {} categories",
        categories.len()
    );
    let categories = polars::export::arrow::array::Utf8ViewArray::from_vec(
        categories
            .iter()
            .map(String::as_str)
            .collect_vec(),
        ArrowDataType::Utf8View,
    );
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Enum(Some(rev_mapping), CategoricalOrdering::Physical)
}

/// Creates a schema from separate arrays of names and data_structs types
///
/// # Arguments
/// * `names` - Array of column names
/// * `dtypes` - Array of data_structs types corresponding to the names
///
/// # Returns
/// A Polars Schema with the provided column names and data_structs types
pub(crate) fn schema_from_arrays(
    names: &[&str],
    dtypes: &[DataType],
) -> Schema {
    if names.len() != dtypes.len() {
        warn!(
            "Mismatch between names and dtypes array lengths: {} vs {}",
            names.len(),
            dtypes.len()
        );
    }
    Schema::from_iter(
        names
            .iter()
            .cloned()
            .map_into()
            .zip(dtypes.iter().cloned()),
    )
}

/// Creates a hashmap from separate arrays of names and data_structs types
///
/// # Arguments
/// * `names` - Array of column names
/// * `dtypes` - Array of data_structs types corresponding to the names
///
/// # Returns
/// A PlHashMap with column names as keys and data_structs types as values
pub(crate) fn hashmap_from_arrays<'a>(
    names: &[&'a str],
    dtypes: &[DataType],
) -> PlHashMap<&'a str, DataType> {
    if names.len() != dtypes.len() {
        warn!(
            "Mismatch between names and dtypes array lengths: {} vs {}",
            names.len(),
            dtypes.len()
        );
    }
    PlHashMap::from_iter(
        names
            .iter()
            .cloned()
            .map_into()
            .zip(dtypes.iter().cloned()),
    )
}

// ======================= Genomic Position Utilities =======================
/// Extracts the first genomic position from a DataFrame
///
/// # Arguments
/// * `data_structs` - The DataFrame containing position data_structs
/// * `chr_col` - Name of the chromosome column
/// * `pos_col` - Name of the position column
///
/// # Returns
/// A Result containing the first GenomicPosition or a Polars error
pub(crate) fn first_position(
    data: &DataFrame,
    chr_col: &str,
    pos_col: &str,
) -> anyhow::Result<GenomicPosition<u64>> {
    debug!(
        "Extracting first position from DataFrame using columns: {}, {}",
        chr_col, pos_col
    );

    let chr = data
        .column(chr_col)
        .with_context(|| {
            format!("Failed to find chromosome column '{}'", chr_col)
        })?
        .as_series()
        .unwrap()
        .first()
        .as_any_value()
        .cast(&DataType::String)
        .str_value()
        .to_string();

    let pos = data
        .column(pos_col)
        .with_context(|| {
            format!("Failed to find position column '{}'", pos_col)
        })?
        .cast(&DataType::UInt64)?
        .u64()?
        .first()
        .unwrap();

    Ok(GenomicPosition::new(chr, pos))
}

/// Extracts the last genomic position from a DataFrame
///
/// # Arguments
/// * `data_structs` - The DataFrame containing position data_structs
/// * `chr_col` - Name of the chromosome column
/// * `pos_col` - Name of the position column
///
/// # Returns
/// A Result containing the last GenomicPosition or a Polars error
pub(crate) fn last_position(
    data: &DataFrame,
    chr_col: &str,
    pos_col: &str,
) -> anyhow::Result<GenomicPosition<u64>> {
    debug!(
        "Extracting last position from DataFrame using columns: {}, {}",
        chr_col, pos_col
    );

    let chr = data
        .column(chr_col)
        .with_context(|| {
            format!("Failed to find chromosome column '{}'", chr_col)
        })?
        .as_series()
        .unwrap()
        .last()
        .as_any_value()
        .cast(&DataType::String)
        .str_value()
        .to_string();

    let pos = data
        .column(pos_col)
        .with_context(|| {
            format!("Failed to find position column '{}'", pos_col)
        })?
        .cast(&DataType::UInt64)?
        .u64()?
        .last()
        .unwrap_or(0);

    Ok(GenomicPosition::new(chr, pos))
}

// ======================= Data Encoding/Decoding Functions
// =======================
/// Encodes strand information as boolean values
///
/// Converts "+" to true, "-" to false, and everything else to null
///
/// # Arguments
/// * `lazy_frame` - The input LazyFrame
/// * `strand_col` - Name of the strand column to encode
///
/// # Returns
/// A LazyFrame with an additional "strand" column containing encoded values
pub fn encode_strand(
    lazy_frame: LazyFrame,
    strand_col: &str,
) -> LazyFrame {
    info!("Encoding strand column '{}' to boolean values", strand_col);
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

/// Encodes methylation context information as boolean values
///
/// Converts "CG" to true, "CHG" to false, and everything else to null
///
/// # Arguments
/// * `lazy_frame` - The input LazyFrame
/// * `context_col` - Name of the context column to encode
///
/// # Returns
/// A LazyFrame with the context column containing encoded values
pub fn encode_context(
    lazy_frame: LazyFrame,
    context_col: &str,
) -> LazyFrame {
    info!(
        "Encoding context column '{}' to boolean values",
        context_col
    );
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

/// Decodes boolean strand information back to string representation
///
/// Converts true to "+", false to "-", and null to "."
///
/// # Arguments
/// * `lazy_frame` - The input LazyFrame
/// * `strand_col` - Name of the boolean strand column to decode
/// * `result_name` - Name for the resulting string column
///
/// # Returns
/// A LazyFrame with an additional column containing decoded strand values
pub fn decode_strand(
    lazy_frame: LazyFrame,
    strand_col: &str,
    result_name: &str,
) -> LazyFrame {
    info!(
        "Decoding strand column '{}' to string column '{}'",
        strand_col, result_name
    );
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

/// Decodes boolean context information back to string representation
///
/// Converts true to "CG", false to "CHG", and null to "CHH"
///
/// # Arguments
/// * `lazy_frame` - The input LazyFrame
/// * `context_col` - Name of the boolean context column to decode
/// * `result_name` - Name for the resulting string column
///
/// # Returns
/// A LazyFrame with an additional column containing decoded context values
pub fn decode_context(
    lazy_frame: LazyFrame,
    context_col: &str,
    result_name: &str,
) -> LazyFrame {
    info!(
        "Decoding context column '{}' to string column '{}'",
        context_col, result_name
    );
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

// ======================= Statistical Tests =======================
/// Two-dimensional, two-sample Kolmogorov-Smirnov test
///
/// Tests the hypothesis that two 2D samples are from the same distribution.
///
/// # Arguments
/// * `x1`, `y1` - First sample's x and y coordinates
/// * `x2`, `y2` - Second sample's x and y coordinates
///
/// # Returns
/// A tuple of (KS statistic, p-value).
/// Small p-values indicate the samples are significantly different.
pub fn ks2d_2sample<X, Y>(
    x1: &[X],
    y1: &[Y],
    x2: &[X],
    y2: &[Y],
) -> (f64, f64)
where
    X: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps
        + Sync
        + Send
        + fmt::Debug,
    Y: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps
        + Sync
        + Send
        + fmt::Debug, {
    info!(
        "Performing 2D two-sample KS test: sample1={}, sample2={}",
        x1.len(),
        x2.len()
    );

    // Validate input dimensions
    if x1.len() != y1.len() {
        warn!(
            "Input error: x1 length ({}) doesn't match y1 length ({})",
            x1.len(),
            y1.len()
        );
        return (0.0, 1.0);
    }
    if x2.len() != y2.len() {
        warn!(
            "Input error: x2 length ({}) doesn't match y2 length ({})",
            x2.len(),
            y2.len()
        );
        return (0.0, 1.0);
    }

    // Find maximum difference using first sample's points as references
    let d1 = (0..x1.len())
        .into_par_iter()
        .map(|i| {
            let (fpp1, fmp1, fpm1, fmm1) =
                count_quads_generic(x1, y1, x1[i], y1[i]);
            let (fpp2, fmp2, fpm2, fmm2) =
                count_quads_generic(x2, y2, x1[i], y1[i]);
            vec![
                (fpp1 - fpp2).abs(),
                (fmp1 - fmp2).abs(),
                (fpm1 - fpm2).abs(),
                (fmm1 - fmm2).abs(),
            ]
            .into_iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap_or(0.0)
        })
        .max_by(|a, b| a.total_cmp(b))
        .unwrap_or(0.0);

    // Find maximum difference using second sample's points as references
    let d2 = (0..x2.len())
        .into_par_iter()
        .map(|i| {
            let (fpp1, fmp1, fpm1, fmm1) =
                count_quads_generic(x1, y1, x2[i], y2[i]);
            let (fpp2, fmp2, fpm2, fmm2) =
                count_quads_generic(x2, y2, x2[i], y2[i]);
            vec![
                (fpp1 - fpp2).abs(),
                (fmp1 - fmp2).abs(),
                (fpm1 - fpm2).abs(),
                (fmm1 - fmm2).abs(),
            ]
            .into_iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap_or(0.0)
        })
        .max_by(|a, b| a.total_cmp(b))
        .unwrap_or(0.0);

    // Calculate test statistic and p-value
    let d = (d1 + d2) / 2.0;
    let sqen =
        ((x1.len() * x2.len()) as f64 / (x1.len() + x2.len()) as f64).sqrt();
    let r1 = pearson_r(x1, y1);
    let r2 = pearson_r(x2, y2);
    let rr = (1.0 - (r1 * r1 + r2 * r2) / 2.0).sqrt();

    let prob = ks_prob(
        d * sqen / (1.0 + rr * (0.25 - 0.75 / sqen)),
        Some(100),   // default iterations
        Some(1e-17), // default precision
    );

    debug!("KS2D test results: d={:.4}, p={:.6}", d, prob);
    (d, prob)
}

/// Calculates Pearson correlation coefficient between two variables
///
/// # Arguments
/// * `x` - Array of x values
/// * `y` - Array of y values
///
/// # Returns
/// Pearson's r coefficient (-1 to 1)
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

    debug!(
        "Calculating Pearson correlation for {} data_structs points",
        x.len()
    );

    // Convert to f64 for calculations
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

/// Computes quadrant probabilities for points relative to a reference point
///
/// # Arguments
/// * `x`, `y` - Arrays of x and y coordinates
/// * `x0`, `y0` - Reference point coordinates
///
/// # Returns
/// Probabilities for each quadrant (++, -+, +-, --) relative to the reference
/// point
fn count_quads_generic<X, Y>(
    x: &[X],
    y: &[Y],
    x0: X,
    y0: Y,
) -> (f64, f64, f64, f64)
where
    X: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps,
    Y: num::ToPrimitive
        + num::FromPrimitive
        + PartialOrd
        + Copy
        + Clone
        + num::traits::NumOps, {
    if x.len() != y.len() {
        warn!(
            "Cannot count quadrants: x length ({}) doesn't match y length ({})",
            x.len(),
            y.len()
        );
        return (0.0, 0.0, 0.0, 0.0);
    }

    let points_len = x.len();
    if points_len == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let ff = 1.0 / points_len as f64;

    // Count points in each quadrant
    let fpp = (0..points_len)
        .filter(|cur_idx| (x[*cur_idx] > x0) && (y[*cur_idx] > y0))
        .count() as f64
        * ff;

    let fnp = (0..points_len)
        .filter(|cur_idx| (x[*cur_idx] < x0) && (y[*cur_idx] > y0))
        .count() as f64
        * ff;

    let fpn = (0..points_len)
        .filter(|cur_idx| (x[*cur_idx] > x0) && (y[*cur_idx] < y0))
        .count() as f64
        * ff;

    let fnn = (0..points_len)
        .filter(|cur_idx| (x[*cur_idx] < x0) && (y[*cur_idx] < y0))
        .count() as f64
        * ff;

    (fpp, fnp, fpn, fnn)
}

/// Computes the KS probability function value
///
/// Based on Numerical Recipes algorithm for the significance level
/// of an observed Kolmogorov-Smirnov D statistic.
///
/// # Arguments
/// * `alam` - D statistic value
/// * `iter` - Maximum number of iterations (default: 100)
/// * `prec` - Convergence precision (default: 1e-17)
///
/// # Returns
/// The significance level (p-value)
fn ks_prob(
    alam: f64,
    iter: Option<usize>,
    prec: Option<f64>,
) -> f64 {
    let mut toadd: Vec<f64> = vec![1.0];
    let mut qks = 0f64;
    let mut j = 1;
    let iter = iter.unwrap_or(100);
    let prec = prec.unwrap_or(1e-17);

    debug!(
        "Computing KS probability for D={}, max_iter={}, precision={}",
        alam, iter, prec
    );

    // Iteratively compute probability until convergence
    while (j < iter) && (toadd.last().unwrap_or(&f64::MAX).abs() > prec * 2.0) {
        let new_elem = 2.0
            * (-1.0f64).powi(j as i32 - 1)
            * (-2.0 * j.pow(2) as f64 * alam.powi(2)).exp();
        qks += new_elem;
        toadd.push(new_elem);
        j += 1;
    }

    if j == iter {
        warn!(
            "KS probability calculation did not converge after {} iterations",
            iter
        );
        return 1.0;
    }

    if qks > 1.0 {
        warn!(
            "KS probability calculation resulted in value > 1.0: {}",
            qks
        );
        return 1.0;
    }

    if qks < prec {
        debug!("KS probability close to zero (< {})", prec);
        0.0
    }
    else {
        qks
    }
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

/// Performs Mann-Whitney U test for two independent samples
///
/// A non-parametric test for determining whether two samples come from the same
/// distribution.
///
/// # Arguments
/// * `group1` - First sample values
/// * `group2` - Second sample values
///
/// # Returns
/// A tuple containing (U statistic, two-tailed p-value)
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
            .unwrap_or_else(|| {
                std::cmp::Ordering::Equal
            })
    });

    // Assign ranks, handling ties by averaging
    debug!("Assigning ranks to {} observations", observations.len());
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
        for j in start..end {
            observations[j].rank = avg_rank;
        }

        // Record ties for variance correction
        if count > 1 {
            tie_groups.push(count);
            debug!("Found tie group of size {}", count);
        }

        i = end;
    }

    // Compute the sum of ranks for group1
    let r1: F = F::from(observations
        .iter()
        .filter(|obs| obs.group == 0)
        .map(|obs| obs.rank)
        .sum::<f64>()).unwrap();

    // Calculate U statistics
    let u1 = r1 - n1 * (n1 + F::from(1.0).unwrap()) / F::from(2.0).unwrap();
    let u2 = n1 * n2 - u1;
    // Use the smaller U value for the test
    let u_stat = u1.min(u2);

    // Mean of U under H0
    let mean_u = n1 * n2 / F::from(2.0).unwrap();

    // Variance with tie correction
    let tie_sum = F::from(tie_groups
        .iter()
        .map(|&t| (t * t * t - t) as f64)
        .sum::<f64>()).unwrap();
    let variance_u = n1 * n2 / F::from(12.0).unwrap()
        * ((n_total + F::from(1.0).unwrap()) - tie_sum / (n_total * (n_total - F::from(1.0).unwrap())));

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
