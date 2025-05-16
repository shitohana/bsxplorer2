use itertools::Itertools;
use polars::frame::DataFrame;
use polars::prelude::*;

use super::{BsxBatch, BsxColumns};
use crate::plsmallstr;

pub fn merge_replicates(
    mut batches: Vec<BsxBatch>,
    count_agg: fn(Vec<&Column>) -> Column,
    density_agg: fn(Vec<&Column>) -> Column,
) -> anyhow::Result<BsxBatch> {
    if batches.is_empty() {
        anyhow::bail!("batches cannot be empty");
    }
    else if batches.len() == 1 {
        return Ok(batches.pop().unwrap());
    }
    else {
        // Assume all batches have the same length and identical chr, pos,
        // strand, context This should be guaranteed by the caller or
        // previous steps (e.g., alignment)
        let first_batch = unsafe { batches.get_unchecked(0) };
        let chr_col = first_batch.data().column(BsxColumns::Chr.as_str())?;
        let pos_col = first_batch.data().column(BsxColumns::Position.as_str())?;
        let strand_col = first_batch.data().column(BsxColumns::Strand.as_str())?;
        let context_col = first_batch.data().column(BsxColumns::Context.as_str())?;

        let count_m_col = count_agg(
            batches
                .iter()
                .map(|b| b.data().column(BsxColumns::CountM.as_str()).unwrap())
                .collect_vec(),
        )
        .with_name(plsmallstr!(BsxColumns::CountM.as_str())); // Ensure correct name

        let count_total_col = count_agg(
            batches
                .iter()
                .map(|b| b.data().column(BsxColumns::CountTotal.as_str()).unwrap())
                .collect_vec(),
        )
        .with_name(plsmallstr!(BsxColumns::CountTotal.as_str())); // Ensure correct name

        let density_col = density_agg(
            batches
                .iter()
                .map(|b| b.data().column(BsxColumns::Density.as_str()).unwrap())
                .collect_vec(),
        )
        .with_name(plsmallstr!(BsxColumns::Density.as_str())); // Ensure correct name

        // Ensure columns have the correct data type for the target batch type B
        let count_m_col = count_m_col.cast(&BsxColumns::CountM.dtype())?;
        let count_total_col = count_total_col.cast(&BsxColumns::CountTotal.dtype())?;
        let density_col = density_col.cast(&BsxColumns::Density.dtype())?;

        let df = DataFrame::from_iter([
            chr_col.to_owned(),
            pos_col.to_owned(),
            strand_col.to_owned(),
            context_col.to_owned(),
            count_m_col,
            count_total_col,
            density_col,
        ]);

        // Create the batch without checks, assuming input batches were valid
        // and aligned

        Ok(unsafe { BsxBatch::new_unchecked(df) })
    }
}

pub fn create_caregorical_dtype<S, P>(chr_values: P) -> DataType
where
    S: AsRef<str>,
    P: AsRef<[Option<S>]>, {
    use polars::export::arrow::array::Utf8ViewArray;
    let categories = Utf8ViewArray::from_slice(chr_values);
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Categorical(Some(rev_mapping), CategoricalOrdering::Physical)
}

pub const fn create_empty_categorical_dtype() -> DataType {
    DataType::Categorical(None, CategoricalOrdering::Physical)
}

#[macro_export]
macro_rules! name_dtype_tuple {
    ($enum_var: expr) => {
        ($enum_var.as_str().into(), $enum_var.dtype())
    };
}
#[macro_export]
macro_rules! create_empty_series {
    ($col: ident) => {
        Series::new_empty(plsmallstr!(BsxCol::$col.as_str()), &BsxCol::$col.dtype())
    };
}
#[macro_export]
macro_rules! get_col_fn {
    ($name: ident, $col: expr, $col_fn: ident, $rettype: ty) => {
        pub fn $name(&self) -> &$rettype {
            self.data().column($col).unwrap().$col_fn().unwrap()
        }
    };
}

pub(crate) use {create_empty_series, get_col_fn, name_dtype_tuple};


#[cfg(test)]
mod tests {
    use polars::df;
    use polars::prelude::*;

    use super::*;

    // --- Helper Functions for Aggregation ---

    fn sum_agg(cols: Vec<&Column>) -> Column {
        if cols.is_empty() {
            panic!("sum_agg received empty input");
        }

        let series_vec: Vec<Series> = cols
            .into_iter()
            .map(|c| c.to_owned().as_materialized_series().clone())
            .collect();

        let initial_series = series_vec[0].clone();

        let sum_series = series_vec
            .into_iter()
            .skip(1) // Skip the first one as it's the initial value
            .fold(initial_series, |acc, s| {
                (&acc + &s).expect("Failed to sum Series in sum_agg") // Add error handling if needed
            });

        sum_series.into() // Convert the final Series back to a Column
    }

    fn mean_agg(cols: Vec<&Column>) -> Column {
        if cols.is_empty() {
            return Series::new_empty("mean_agg_empty".into(), &DataType::Float64)
                .into();
        }

        let series_vec: PolarsResult<Vec<Series>> = cols
            .into_iter()
            .map(|c| {
                c.to_owned()
                    .as_materialized_series()
                    .cast(&DataType::Float64)
            })
            .collect();

        let series_vec = match series_vec {
            Ok(vec) => vec,
            Err(e) => {
                panic!("Failed to cast column to Float64 in mean_agg: {}", e);
            },
        };

        let n = series_vec.len() as f64; // Use f64 for division

        let initial_series = series_vec[0].clone();

        let total_sum = series_vec
                .into_iter()
                .skip(1) // Skip the first one as it's the initial value
                .fold(initial_series, |acc, s| {
                    // Perform addition between Series
                    // This should work now as all series are Float64
                    (&acc + &s).expect("Failed to sum Series in mean_agg")
                });

        // Divide the sum by the number of columns
        let mean_series = &total_sum / n;
        mean_series.into() // Convert the final Series back to a Column
    }

    // --- Helper Functions to Create Test Batches ---

    #[allow(unsafe_code)]
    fn create_decoded_batch_1() -> BsxBatch {
        let df = df!(
            BsxColumns::Chr.as_str() => &["chr1", "chr1"],
            BsxColumns::Position.as_str() => &[10u64, 20],
            BsxColumns::Strand.as_str() => &["+", "-"],
            BsxColumns::Context.as_str() => &["CG", "CHG"],
            BsxColumns::CountM.as_str() => &[2u32, 4],
            BsxColumns::CountTotal.as_str() => &[10u32, 8],
            BsxColumns::Density.as_str() => &[0.2f64, 0.5],
        )
        .unwrap();
        unsafe { BsxBatch::new_unchecked(df) }
    }

    #[allow(unsafe_code)]
    fn create_decoded_batch_2() -> BsxBatch {
        let df = df!(
            BsxColumns::Chr.as_str() => &["chr1", "chr1"],
            BsxColumns::Position.as_str() => &[10u64, 20],
            BsxColumns::Strand.as_str() => &["+", "-"],
            BsxColumns::Context.as_str() => &["CG", "CHG"],
            BsxColumns::CountM.as_str() => &[3u32, 6],
            BsxColumns::CountTotal.as_str() => &[10u32, 12],
            BsxColumns::Density.as_str() => &[0.3f64, 0.5],
        )
        .unwrap();
        unsafe { BsxBatch::new_unchecked(df) }
    }

    // --- Test Cases ---

    #[test]
    fn test_merge_replicates_decoded() -> anyhow::Result<()> {
        let batch1 = create_decoded_batch_1();
        let batch2 = create_decoded_batch_2();
        let batches = vec![batch1, batch2];

        let merged_batch: BsxBatch = merge_replicates(batches, sum_agg, mean_agg)?;

        // Expected results
        let expected_df = df!(
            BsxColumns::Chr.as_str() => &["chr1", "chr1"],
            BsxColumns::Position.as_str() => &[10u64, 20],
            BsxColumns::Strand.as_str() => &["+", "-"],
            BsxColumns::Context.as_str() => &["CG", "CHG"],
            BsxColumns::CountM.as_str() => &[5u32, 10],
            BsxColumns::CountTotal.as_str() => &[20u32, 20],
            BsxColumns::Density.as_str() => &[0.25f64, 0.5],
        )?;

        // Compare DataFrames
        assert_eq!(merged_batch.data(), &expected_df);
        // Ensure the type is correct
        assert_eq!(merged_batch.count_m().dtype(), &BsxColumns::CountM.dtype());
        assert_eq!(merged_batch.density().dtype(), &BsxColumns::Density.dtype());

        Ok(())
    }

    #[test]
    fn test_merge_replicates_empty() {
        let batches: Vec<BsxBatch> = vec![];
        let result = merge_replicates(batches, sum_agg, mean_agg);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("batches cannot be empty"));
    }

    #[test]
    fn test_merge_replicates_single() -> anyhow::Result<()> {
        let batch1 = create_decoded_batch_1();
        let batches = vec![batch1.clone()]; // Clone batch1 for comparison
        let result_batch: BsxBatch = merge_replicates(batches, sum_agg, mean_agg)?;
        assert_eq!(result_batch, batch1); // Should return the original batch
        Ok(())
    }
}
