use itertools::Itertools;
use polars::frame::DataFrame;
use polars::prelude::Column;

use crate::data_structs::batch::traits::{colnames, BsxBatchMethods};

#[allow(unsafe_code)]
pub fn merge_replicates<B: BsxBatchMethods>(
    mut batches: Vec<B>,
    count_agg: fn(Vec<&Column>) -> Column,
    density_agg: fn(Vec<&Column>) -> Column,
) -> anyhow::Result<B> {
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
        let first_batch_data = unsafe { batches.get_unchecked(0).data() };
        let chr_col = first_batch_data.column(colnames::CHR_NAME)?;
        let pos_col = first_batch_data.column(colnames::POS_NAME)?;
        let strand_col = first_batch_data.column(colnames::STRAND_NAME)?;
        let context_col = first_batch_data.column(colnames::CONTEXT_NAME)?;

        let count_m_col = count_agg(
            batches
                .iter()
                .map(|b| {
                    b.data()
                        .column(colnames::COUNT_M_NAME)
                        .unwrap()
                })
                .collect_vec(),
        )
        .with_name(colnames::COUNT_M_NAME.into()); // Ensure correct name

        let count_total_col = count_agg(
            batches
                .iter()
                .map(|b| {
                    b.data()
                        .column(colnames::COUNT_TOTAL_NAME)
                        .unwrap()
                })
                .collect_vec(),
        )
        .with_name(colnames::COUNT_TOTAL_NAME.into()); // Ensure correct name

        let density_col = density_agg(
            batches
                .iter()
                .map(|b| {
                    b.data()
                        .column(colnames::DENSITY_NAME)
                        .unwrap()
                })
                .collect_vec(),
        )
        .with_name(colnames::DENSITY_NAME.into()); // Ensure correct name

        // Ensure columns have the correct data type for the target batch type B
        let count_m_col = count_m_col.cast(&B::count_type())?;
        let count_total_col = count_total_col.cast(&B::count_type())?;
        let density_col = density_col.cast(&B::density_type())?;

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
        let batch = unsafe { B::new_unchecked(df) };
        Ok(batch)
    }
}

#[cfg(test)]
mod tests {
    use polars::df;
    use polars::prelude::*;

    use super::*;
    use crate::data_structs::batch::decoded::BsxBatch;
    use crate::data_structs::batch::encoded::EncodedBsxBatch;
    use crate::data_structs::batch::traits::colnames::*;
    use crate::data_structs::batch::traits::BsxBatchMethods;

    // --- Helper Functions for Aggregation ---

    fn sum_agg(cols: Vec<&Column>) -> Column {
        if cols.is_empty() {
            // Handle empty input, perhaps return an empty column or error
            // For now, let's assume at least one column based on
            // merge_replicates logic
            panic!("sum_agg received empty input");
        }

        // Convert columns to Series
        let series_vec: Vec<Series> = cols
            .into_iter()
            .map(|c| {
                c.to_owned()
                    .as_materialized_series()
                    .clone()
            })
            .collect();

        // Use fold to sum all series together
        // Start with the first series as the initial accumulator
        let initial_series = series_vec[0].clone();

        let sum_series = series_vec
            .into_iter()
            .skip(1) // Skip the first one as it's the initial value
            .fold(initial_series, |acc, s| {
                // Perform addition between Series
                // This might panic if types are incompatible, but merge_replicates
                // should provide columns of the same type (e.g., all CountType)
                (&acc + &s).expect("Failed to sum Series in sum_agg") // Add error handling if needed
            });

        sum_series.into() // Convert the final Series back to a Column
    }

    fn mean_agg(cols: Vec<&Column>) -> Column {
        if cols.is_empty() {
            // Handle empty input
            // Return an empty Float64 column as a reasonable default
            return Series::new_empty(
                "mean_agg_empty".into(),
                &DataType::Float64,
            )
            .into();
        }

        // Convert columns to Series and cast to Float64
        let series_vec: PolarsResult<Vec<Series>> = cols
            .into_iter()
            .map(|c| {
                c.to_owned()
                    .as_materialized_series()
                    .cast(&DataType::Float64)
            })
            .collect();

        // Check for casting errors
        let series_vec = match series_vec {
            Ok(vec) => vec,
            Err(e) => {
                // Handle casting error, perhaps panic or return an error Column
                // For now, panic as it indicates unexpected input type
                panic!("Failed to cast column to Float64 in mean_agg: {}", e);
            },
        };

        let n = series_vec.len() as f64; // Use f64 for division

        // Use fold to sum all series together
        // Start with the first series as the initial accumulator
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
            CHR_NAME => &["chr1", "chr1"],
            POS_NAME => &[10u64, 20],
            STRAND_NAME => &["+", "-"],
            CONTEXT_NAME => &["CG", "CHG"],
            COUNT_M_NAME => &[2u32, 4],
            COUNT_TOTAL_NAME => &[10u32, 8],
            DENSITY_NAME => &[0.2f64, 0.5],
        )
        .unwrap();
        unsafe { BsxBatch::new_unchecked(df) }
    }

    #[allow(unsafe_code)]
    fn create_decoded_batch_2() -> BsxBatch {
        let df = df!(
            CHR_NAME => &["chr1", "chr1"], // Identical metadata
            POS_NAME => &[10u64, 20],
            STRAND_NAME => &["+", "-"],
            CONTEXT_NAME => &["CG", "CHG"],
            COUNT_M_NAME => &[3u32, 6], // Different counts/density
            COUNT_TOTAL_NAME => &[10u32, 12],
            DENSITY_NAME => &[0.3f64, 0.5],
        )
        .unwrap();
        unsafe { BsxBatch::new_unchecked(df) }
    }

    #[allow(unsafe_code)]
    fn create_encoded_batch_1() -> EncodedBsxBatch {
        let chr_series = Series::new(CHR_NAME.into(), &["chr1", "chr1"])
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
            .unwrap();
        let df = df!(
            CHR_NAME => chr_series,
            POS_NAME => &[10u32, 20],
            STRAND_NAME => &[true, false], // "+" -> true, "-" -> false
            CONTEXT_NAME => &[true, false], // "CG" -> true, "CHG" -> false
            COUNT_M_NAME => &[2i16, 4],
            COUNT_TOTAL_NAME => &[10i16, 8],
            DENSITY_NAME => &[0.2f32, 0.5],
        )
        .unwrap();
        unsafe { EncodedBsxBatch::new_unchecked(df) }
    }

    #[allow(unsafe_code)]
    fn create_encoded_batch_2() -> EncodedBsxBatch {
        let chr_series = Series::new(CHR_NAME.into(), &["chr1", "chr1"])
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
            .unwrap();
        let df = df!(
            CHR_NAME => chr_series, // Identical metadata
            POS_NAME => &[10u32, 20],
            STRAND_NAME => &[true, false],
            CONTEXT_NAME => &[true, false],
            COUNT_M_NAME => &[3i16, 6], // Different counts/density
            COUNT_TOTAL_NAME => &[10i16, 12],
            DENSITY_NAME => &[0.3f32, 0.5],
        )
        .unwrap();
        unsafe { EncodedBsxBatch::new_unchecked(df) }
    }

    // --- Test Cases ---

    #[test]
    fn test_merge_replicates_decoded() -> anyhow::Result<()> {
        let batch1 = create_decoded_batch_1();
        let batch2 = create_decoded_batch_2();
        let batches = vec![batch1, batch2];

        let merged_batch: BsxBatch =
            merge_replicates(batches, sum_agg, mean_agg)?;

        // Expected results
        let expected_df = df!(
            CHR_NAME => &["chr1", "chr1"],
            POS_NAME => &[10u64, 20],
            STRAND_NAME => &["+", "-"],
            CONTEXT_NAME => &["CG", "CHG"],
            COUNT_M_NAME => &[5u32, 10], // Sum: 2+3, 4+6
            COUNT_TOTAL_NAME => &[20u32, 20], // Sum: 10+10, 8+12
            DENSITY_NAME => &[0.25f64, 0.5], // Mean: (0.2+0.3)/2, (0.5+0.5)/2
        )?;

        // Compare DataFrames
        assert_eq!(merged_batch.data(), &expected_df);
        // Ensure the type is correct
        assert_eq!(merged_batch.count_m().dtype(), &BsxBatch::count_type());
        assert_eq!(merged_batch.density().dtype(), &BsxBatch::density_type());

        Ok(())
    }

    #[test]
    fn test_merge_replicates_encoded() -> anyhow::Result<()> {
        let batch1 = create_encoded_batch_1();
        let batch2 = create_encoded_batch_2();
        let batches = vec![batch1, batch2];

        let merged_batch: EncodedBsxBatch =
            merge_replicates(batches, sum_agg, mean_agg)?;

        // Expected results
        let expected_chr_series =
            Series::new(CHR_NAME.into(), &["chr1", "chr1"]).cast(
                &DataType::Categorical(None, CategoricalOrdering::Physical),
            )?;
        let expected_df = df!(
            CHR_NAME => expected_chr_series,
            POS_NAME => &[10u32, 20],
            STRAND_NAME => &[true, false],
            CONTEXT_NAME => &[true, false],
            COUNT_M_NAME => &[5i16, 10], // Sum: 2+3, 4+6
            COUNT_TOTAL_NAME => &[20i16, 20], // Sum: 10+10, 8+12
            DENSITY_NAME => &[0.25f32, 0.5f32], // Mean: (0.2+0.3)/2, (0.5+0.5)/2
        )?;

        // Compare DataFrames
        assert_eq!(merged_batch.data(), &expected_df);
        // Ensure the type is correct
        assert_eq!(
            merged_batch.count_m().dtype(),
            &EncodedBsxBatch::count_type()
        );
        assert_eq!(
            merged_batch.density().dtype(),
            &EncodedBsxBatch::density_type()
        );

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
        let result_batch: BsxBatch =
            merge_replicates(batches, sum_agg, mean_agg)?;
        assert_eq!(result_batch, batch1); // Should return the original batch
        Ok(())
    }
}
