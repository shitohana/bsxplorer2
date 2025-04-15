use super::{
    builder::BsxBatchBuilder,
    decoded::BsxBatch,
    encoded::EncodedBsxBatch,
    traits::{BsxBatchMethods, BsxTypeTag},
};
use crate::data_structs::batch::builder::ContextData;
use crate::data_structs::batch::traits::{colnames, BatchType};
use crate::utils::types::{Context, IPCEncodedEnum, Strand};
use polars::prelude::*;
use polars::series::IsSorted;
use std::marker::PhantomData;

/// A lazy representation of a BSX batch for efficient query operations
pub struct LazyBsxBatch<T: BsxTypeTag + BsxBatchMethods> {
    data: LazyFrame,
    _phantom: PhantomData<T>,
}

impl<T: BsxTypeTag + BsxBatchMethods> LazyBsxBatch<T> {
    pub fn collect(self) -> PolarsResult<T> {
        let data = self.data.collect()?;
        Ok(unsafe {T::new_unchecked(data)})
    }
    /// Applies a filter expression to the batch
    fn filter(
        self,
        predicate: Expr,
    ) -> Self {
        Self {
            data: self.data.filter(predicate),
            _phantom: PhantomData::<T>,
        }
    }

    /// Creates a LazyBsxBatch from an existing LazyFrame
    fn from_lazy(lazy: LazyFrame) -> Self {
        Self {
            data: lazy,
            _phantom: PhantomData::<T>,
        }
    }

    /// Filters positions less than the specified value
    pub fn filter_pos_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(colnames::POS_NAME).lt(lit(value)))
    }

    /// Filters positions greater than the specified value
    pub fn filter_pos_gt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(colnames::POS_NAME).gt(lit(value)))
    }

    /// Filters entries with coverage less than the specified value
    pub fn filter_coverage_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(colnames::COUNT_TOTAL_NAME).lt(lit(value)))
    }

    /// Filters entries by strand value
    pub fn filter_strand(
        self,
        value: Strand,
    ) -> Self {
        let expr = match T::type_name() {
            "decoded" => col(colnames::STRAND_NAME).eq(lit(value.to_string())),
            "encoded" => col(colnames::STRAND_NAME).eq(value
                .to_bool()
                .map(|v| lit(v))
                .unwrap_or(lit(NULL))),
            other => unimplemented!("Type {}", other),
        };
        self.filter(expr)
    }

    /// Filters entries by context value
    pub fn filter_context(
        self,
        value: Context,
    ) -> Self {
        let expr = match T::type_name() {
            "decoded" => col(colnames::CONTEXT_NAME).eq(lit(value.to_string())),
            "encoded" => col(colnames::CONTEXT_NAME).eq(value
                .to_bool()
                .map(|v| lit(v))
                .unwrap_or(lit(NULL))),
            other => unimplemented!("Type {}", other),
        };
        self.filter(expr)
    }

    /// Marks entries with coverage below threshold as zero and adjusts related values
    pub fn mark_low_coverage(self, threshold: u32) -> Self {
        let res = self.data.with_column(
            when(col("count_total").lt(lit(threshold)))
                .then(lit(0))
                .otherwise(col("count_total"))
                .alias("count_total"),
        )
            .with_columns([
                when(col("count_total").eq(lit(0)))
                    .then(lit(0))
                    .otherwise(col("count_m"))
                    .alias("count_m"),
                when(col("count_total").eq(lit(0)))
                    .then(lit(f64::NAN))
                    .otherwise(col("density"))
                    .alias("density"),
            ]);
        Self::from_lazy(res)
    }

    pub fn align_with_contexts(mut self, context_data: ContextData, chr_val: &str) -> Self {
        let (pos_col, strand_col, context_col) = context_data.take();
        let mut pos_col = Series::from_vec(colnames::POS_NAME.into(), pos_col);
        pos_col.set_sorted_flag(IsSorted::Ascending);

        let strand_col = match T::type_enum() {
            BatchType::Decoded => {
                Series::from_iter(
                    strand_col.into_iter().map(|v| v.to_string())
                )
            },
            BatchType::Encoded => {
                Series::from_iter(
                    strand_col.into_iter().map(|v| v.to_bool())
                )
            }
        }.with_name(colnames::STRAND_NAME.into());

        let context_col = match T::type_enum() {
            BatchType::Decoded => {
                Series::from_iter(
                    context_col.into_iter().map(|v| v.to_string())
                )
            },
            BatchType::Encoded => {
                Series::from_iter(
                    context_col.into_iter().map(|v| v.to_bool())
                )
            }
        }.with_name(colnames::CONTEXT_NAME.into());

        let context_df = DataFrame::from_iter([pos_col, strand_col, context_col])
            .lazy()
            .cast(
                PlHashMap::from_iter([
                    (colnames::POS_NAME, T::pos_type()),
                    (colnames::STRAND_NAME, T::strand_type()),
                    (colnames::CONTEXT_NAME, T::context_type()),
                ]),
                true
            );
        let self_selected = self.data.drop(
            [col(colnames::CONTEXT_NAME), col(colnames::STRAND_NAME)]
        );
        let joined = context_df
            .left_join(
                self_selected,
                col(colnames::POS_NAME),
                col(colnames::POS_NAME)
            )
            .with_columns([
                col(colnames::CHR_NAME).fill_null(lit(chr_val)),
                col(colnames::COUNT_M_NAME).fill_null(lit(0)).alias(colnames::COUNT_M_NAME),
                col(colnames::COUNT_TOTAL_NAME).fill_null(lit(0)).alias(colnames::COUNT_TOTAL_NAME),
                col(colnames::DENSITY_NAME).fill_null(lit(f32::NAN))
            ]);
        Self::from_lazy(joined)
    }
}

impl<T: BsxTypeTag + BsxBatchMethods> TryFrom<LazyBsxBatch<T>> for BsxBatch {
    type Error = anyhow::Error;

    /// Attempts to convert a lazy BSX batch to a decoded BSX batch
    fn try_from(lazy_batch: LazyBsxBatch<T>) -> Result<Self, Self::Error> {
        let df = lazy_batch.data.collect()?;
        BsxBatchBuilder::no_checks().build(df)
    }
}

impl<T: BsxTypeTag + BsxBatchMethods> TryFrom<LazyBsxBatch<T>>
    for EncodedBsxBatch
{
    type Error = anyhow::Error;

    /// Attempts to convert a lazy BSX batch to an encoded BSX batch
    fn try_from(lazy_batch: LazyBsxBatch<T>) -> Result<Self, Self::Error> {
        let df = lazy_batch.data.collect()?;
        unsafe { Ok(EncodedBsxBatch::new_unchecked(df)) }
    }
}

impl<B: BsxBatchMethods> From<B> for LazyBsxBatch<B> {
    fn from(batch: B) -> Self {
        LazyBsxBatch::from_lazy(batch.take().lazy())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::batch::{
        builder::BsxBatchBuilder, decoded::BsxBatch, encoded::EncodedBsxBatch,
        traits::colnames::*, traits::BsxBatchMethods,
    };
    use crate::utils::get_categorical_dtype;
    use crate::utils::types::{Context, Strand};
    use polars::chunked_array::cast::CastOptions;
    use polars::df;
    use polars::prelude::*;
    use std::ops::Not;

    // --- Helper Functions ---

    fn create_test_decoded_df() -> DataFrame {
        df!(
            CHR_NAME => &["chr1", "chr1", "chr1", "chr1", "chr1"],
            POS_NAME => &[10u64, 20, 30, 40, 50],
            STRAND_NAME => &["+", "+", "-", "+", "-"],
            CONTEXT_NAME => &["CG", "CHG", "CHH", "CG", "CHG"],
            COUNT_M_NAME => &[1u32, 2, 3, 4, 5],
            COUNT_TOTAL_NAME => &[5u32, 8, 6, 8, 10],
            DENSITY_NAME => &[0.2f64, 0.25, 0.5, 0.5, 0.5],
        )
        .unwrap()
    }

    fn create_test_encoded_df() -> DataFrame {
        create_test_encoded_batch().take()
    }

    fn create_test_decoded_batch() -> BsxBatch {
        BsxBatchBuilder::no_checks()
            .build(create_test_decoded_df())
            .expect("Failed to build decoded batch")
    }

     fn create_test_encoded_batch() -> EncodedBsxBatch {
        // We build from a decoded DF, the builder handles encoding
         BsxBatchBuilder::encode_batch(create_test_decoded_batch(), get_categorical_dtype(vec!["chr1".to_string()]))
            .expect("Failed to build encoded batch")
    }


    // --- Test Cases ---

    #[test]
    fn test_from_decoded_and_collect() {
        let batch = create_test_decoded_batch();
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch.clone());
        let collected = lazy_batch.collect().unwrap();
        assert_eq!(collected, batch);
    }

    #[test]
    fn test_from_encoded_and_collect() {
        let batch = create_test_encoded_batch();
        let lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch.clone());
        let collected = lazy_batch.collect().unwrap();
        // Note: Encoded batch might have different physical representation for chr,
        // but the data should be logically equivalent. Direct equality works here
        // because we built it deterministically.
        assert_eq!(collected, batch);
    }

    #[test]
    fn test_filter_pos_lt_decoded() {
        let batch = create_test_decoded_batch();
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch.clone());
        let filtered_lazy = lazy_batch.filter_pos_lt(30u64); // Keep pos 10, 20
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_decoded_df().slice(0, 2);
        let expected_batch = BsxBatch::try_from(expected_df).unwrap();
        assert_eq!(collected, expected_batch);
    }

    #[test]
    fn test_filter_pos_lt_encoded() {
        let batch = create_test_encoded_batch();
        let lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch.clone());
        let filtered_lazy = lazy_batch.filter_pos_lt(30u32); // Keep pos 10, 20
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_encoded_df().slice(0, 2);
        let expected_batch = unsafe { EncodedBsxBatch::new_unchecked(expected_df) }; // Assume valid after slice
        assert_eq!(collected.data(), expected_batch.data());
    }


    #[test]
    fn test_filter_pos_gt_decoded() {
        let batch = create_test_decoded_batch();
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch.clone());
        let filtered_lazy = lazy_batch.filter_pos_gt(30u64); // Keep pos 40, 50
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_decoded_df().slice(3, 2);
         let expected_batch = BsxBatch::try_from(expected_df).unwrap();
        assert_eq!(collected, expected_batch);
    }

     #[test]
    fn test_filter_pos_gt_encoded() {
        let batch = create_test_encoded_batch();
        let lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch.clone());
        let filtered_lazy = lazy_batch.filter_pos_gt(30u32); // Keep pos 40, 50
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_encoded_df().slice(3, 2);
        let expected_batch = unsafe { EncodedBsxBatch::new_unchecked(expected_df) };
        assert_eq!(collected.data(), expected_batch.data());
    }

    #[test]
    fn test_filter_coverage_lt_decoded() {
        let batch = create_test_decoded_batch();
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch.clone());
        // Keep total_cov < 8 (rows with 5, 6) -> indices 0, 2
        let filtered_lazy = lazy_batch.filter_coverage_lt(8u32);
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_decoded_df().filter(&create_test_decoded_df().column(COUNT_TOTAL_NAME).unwrap().lt(&Scalar::new(DataType::UInt32, AnyValue::UInt32(8u32)).into_column("test".into())).unwrap()).unwrap();
        let expected_batch = BsxBatch::try_from(expected_df).unwrap();
        assert_eq!(collected, expected_batch);
    }

    #[test]
    fn test_filter_coverage_lt_encoded() {
        let batch = create_test_encoded_batch();
        let lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch.clone());
         // Keep total_cov < 8 (rows with 5, 6) -> indices 0, 2
        let filtered_lazy = lazy_batch.filter_coverage_lt(8i16);
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_encoded_df().filter(&create_test_encoded_df().column(COUNT_TOTAL_NAME).unwrap().lt(&Scalar::new(DataType::Int16, AnyValue::Int16(8i16)).into_column("test".into())).unwrap()).unwrap();
        let expected_batch = unsafe { EncodedBsxBatch::new_unchecked(expected_df) };
        assert_eq!(collected.data(), expected_batch.data());
    }


    #[test]
    fn test_filter_strand_decoded() {
        let batch = create_test_decoded_batch();
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch.clone());
        // Keep strand == "-" -> indices 2, 4
        let filtered_lazy = lazy_batch.filter_strand(Strand::Reverse);
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_decoded_df().filter(&create_test_decoded_df().column(STRAND_NAME).unwrap().equal(&Scalar::new(DataType::String, AnyValue::String("-")).into_column("test".into())).unwrap()).unwrap();
        let expected_batch = BsxBatch::try_from(expected_df).unwrap();
        assert_eq!(collected, expected_batch);
    }

    #[test]
    fn test_filter_strand_encoded() {
        let batch = create_test_encoded_batch();
        let lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch.clone());
        // Keep strand == false (Minus) -> indices 2, 4
        let filtered_lazy = lazy_batch.filter_strand(Strand::Reverse);
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_encoded_df().filter(&create_test_encoded_df().column(STRAND_NAME).unwrap().equal(&Scalar::new(DataType::Boolean, AnyValue::Boolean(false)).into_column("test".into())).unwrap()).unwrap();
         let expected_batch = unsafe { EncodedBsxBatch::new_unchecked(expected_df) };
        assert_eq!(collected.data(), expected_batch.data());
    }

     #[test]
    fn test_filter_context_decoded() {
        let batch = create_test_decoded_batch();
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch.clone());
        // Keep context == "CG" -> indices 0, 3
        let filtered_lazy = lazy_batch.filter_context(Context::CG);
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_decoded_df().filter(&create_test_decoded_df().column(CONTEXT_NAME).unwrap().equal(&Scalar::new(DataType::String, AnyValue::String("CG")).into_column("test".into())).unwrap()).unwrap();
        let expected_batch = BsxBatch::try_from(expected_df).unwrap();
        assert_eq!(collected, expected_batch);
    }

     #[test]
    fn test_filter_context_encoded() {
        let batch = create_test_encoded_batch();
        let lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch.clone());
        // Keep context == true (CG) -> indices 0, 3
        let filtered_lazy = lazy_batch.filter_context(Context::CG);
        let collected = filtered_lazy.collect().unwrap();

        let expected_df = create_test_encoded_df().filter(&create_test_encoded_df().column(CONTEXT_NAME).unwrap().equal(&Scalar::new(DataType::Boolean, AnyValue::Boolean(true)).into_column("test".into())).unwrap()).unwrap();
        let expected_batch = unsafe { EncodedBsxBatch::new_unchecked(expected_df) };
        assert_eq!(collected.data(), expected_batch.data());
    }

    #[test]
    fn test_mark_low_coverage_decoded() {
        let batch = create_test_decoded_batch();
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch);
        let marked_lazy = lazy_batch.mark_low_coverage(7); // Mark rows with total_cov 5, 6 (indices 0, 2)
        let collected = marked_lazy.collect().unwrap();

        let expected_df = df!(
            CHR_NAME => &["chr1", "chr1", "chr1", "chr1", "chr1"],
            POS_NAME => &[10u64, 20, 30, 40, 50],
            STRAND_NAME => &["+", "+", "-", "+", "-"],
            CONTEXT_NAME => &["CG", "CHG", "CHH", "CG", "CHG"],
            // Index 0, 2: count_m becomes 0
            COUNT_M_NAME => &[0u32, 2, 0, 4, 5],
             // Index 0, 2: count_total becomes 0
            COUNT_TOTAL_NAME => &[0u32, 8, 0, 8, 10],
            // Index 0, 2: density becomes NaN
            DENSITY_NAME => &[f64::NAN, 0.25, f64::NAN, 0.5, 0.5],
        )
        .unwrap();
         let expected_batch = BsxBatch::try_from(expected_df).unwrap();

        // Compare DataFrames directly due to NaN
        assert!(collected.data().equals_missing(expected_batch.data()));
    }

     #[test]
    fn test_mark_low_coverage_encoded() {
        let batch = create_test_encoded_batch();
        let lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch);
        let marked_lazy = lazy_batch.mark_low_coverage(7); // Mark rows with total_cov 5, 6 (indices 0, 2)
        let collected = marked_lazy.collect().unwrap();

        let chr_series = Series::new(CHR_NAME.into(), &["chr1", "chr1", "chr1", "chr1", "chr1"])
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
            .unwrap();
        let expected_df = df!(
            CHR_NAME => chr_series,
            POS_NAME => &[10u32, 20, 30, 40, 50],
            STRAND_NAME => &[Some(true), Some(true), Some(false), Some(true), Some(false)],
            CONTEXT_NAME => &[Some(true), Some(false), None, Some(true), Some(false)],
            COUNT_M_NAME => &[0i16, 2, 0, 4, 5],       // Indices 0, 2 are 0
            COUNT_TOTAL_NAME => &[0i16, 8, 0, 8, 10], // Indices 0, 2 are 0
            DENSITY_NAME => &[f32::NAN, 0.25f32, f32::NAN, 0.5f32, 0.5f32], // Indices 0, 2 are NaN
        )
        .unwrap();
         let expected_batch = unsafe{ EncodedBsxBatch::new_unchecked(expected_df) };

        // Compare DataFrames directly due to NaN and potential categorical diffs
        assert!(collected.data().equals_missing(expected_batch.data()));
    }

    #[test]
    fn test_align_with_contexts_decoded() {
        let batch = create_test_decoded_batch(); // pos: 10, 20, 30, 40, 50
        let lazy_batch = LazyBsxBatch::<BsxBatch>::from(batch);

        // Context data: include existing pos (20, 40), new pos (15, 55), different context/strand for existing pos
        let context_positions = vec![15u32, 20, 40, 55];
        let context_strands = vec![Strand::Forward, Strand::Reverse, Strand::Forward, Strand::Reverse]; // Original 20 was +, 40 was +
        let context_contexts = vec![Context::CHH, Context::CG, Context::CHG, Context::CHH]; // Original 20 was CHG, 40 was CG
        let context_data = ContextData::new(context_positions, context_strands, context_contexts);

        let aligned_lazy = lazy_batch.align_with_contexts(context_data, "chr1");
        let collected = aligned_lazy.collect().unwrap();

        // Expected DataFrame:
        // - Rows for pos 15, 55 added with new context/strand, counts 0, density NaN, chr Null.
        // - Rows for pos 20, 40 keep original counts/density/chr, but update context/strand.
        // - Rows for pos 10, 30, 50 are dropped.
        let expected_df = df!(
            // Note: Chr is Null for rows originating only from context_data
            // Note: Original chr column was String
            CHR_NAME => &[None::<String>, Some("chr1".to_string()), Some("chr1".to_string()), None::<String>],
            POS_NAME => &[15u64, 20, 40, 55], // Type matches BsxBatch::pos_type()
            STRAND_NAME => &["+".to_string(), "-".to_string(), "+".to_string(), "-".to_string()], // Type matches BsxBatch::strand_type()
            CONTEXT_NAME => &["CHH".to_string(), "CG".to_string(), "CHG".to_string(), "CHH".to_string()], // Type matches BsxBatch::context_type()
            // Original counts for pos 20 (idx 1), 40 (idx 3)
            COUNT_M_NAME => &[0u32, 2, 4, 0], // Filled nulls with 0
            COUNT_TOTAL_NAME => &[0u32, 8, 8, 0], // Filled nulls with 0
            // Original density for pos 20, 40
            DENSITY_NAME => &[f64::NAN, 0.25, 0.5, f64::NAN], // Filled nulls with NaN
        )
        .unwrap();
        // Reorder collected columns to match expected df for comparison
        let collected_df = collected.data().select([
             CHR_NAME, POS_NAME, STRAND_NAME, CONTEXT_NAME, COUNT_M_NAME, COUNT_TOTAL_NAME, DENSITY_NAME
        ]).unwrap();

        // Need equals_missing because of NaN
        assert!(collected_df.equals_missing(&expected_df));
    }

    #[test]
    fn test_align_with_contexts_encoded() {
        let batch = create_test_encoded_batch(); // pos: 10, 20, 30, 40, 50
        let chr_dtype = batch.chr().dtype().clone();
        let mut lazy_batch = LazyBsxBatch::<EncodedBsxBatch>::from(batch.clone());

        // Context data: include existing pos (20, 40), new pos (15, 55), different context/strand for existing pos
        let context_positions = vec![15u32, 20, 40, 55];
        let context_strands = vec![Strand::Forward, Strand::Reverse, Strand::Forward, Strand::Reverse]; // Original 20 was +, 40 was + -> true, true
        let context_contexts = vec![Context::CHH, Context::CG, Context::CHG, Context::CHH]; // Original 20 was CHG, 40 was CG -> false, true

        let context_data = ContextData::new(context_positions, context_strands, context_contexts);

        let aligned_lazy = lazy_batch.align_with_contexts(context_data, "chr1");
        let collected = aligned_lazy.collect().unwrap();

        // Expected DataFrame:
        // - Rows for pos 15, 55 added with new context/strand, counts 0, density NaN, chr Null.
        // - Rows for pos 20, 40 keep original counts/density/chr, but update context/strand.
        // - Rows for pos 10, 30, 50 are dropped.
        let expected_df = df!(
            // Note: Chr is Null for rows originating only from context_data
            // Note: Original chr column was String
            CHR_NAME => &[Some("chr1".to_string()), Some("chr1".to_string()), Some("chr1".to_string()), Some("chr1".to_string())],
            POS_NAME => &[15u64, 20, 40, 55], // Type matches BsxBatch::pos_type()
            STRAND_NAME => &["+".to_string(), "-".to_string(), "+".to_string(), "-".to_string()], // Type matches BsxBatch::strand_type()
            CONTEXT_NAME => &["CHH".to_string(), "CG".to_string(), "CHG".to_string(), "CHH".to_string()], // Type matches BsxBatch::context_type()
            // Original counts for pos 20 (idx 1), 40 (idx 3)
            COUNT_M_NAME => &[0u32, 2, 4, 0], // Filled nulls with 0
            COUNT_TOTAL_NAME => &[0u32, 8, 8, 0], // Filled nulls with 0
            // Original density for pos 20, 40
            DENSITY_NAME => &[f64::NAN, 0.25, 0.5, f64::NAN], // Filled nulls with NaN
        )
        .unwrap();

        let expected_batch = BsxBatchBuilder::encode_batch(
            unsafe { BsxBatch::new_unchecked(expected_df) },
            chr_dtype
        ).unwrap();

         // Reorder collected columns to match expected df for comparison
        let collected_df = collected.data().select([
             CHR_NAME, POS_NAME, STRAND_NAME, CONTEXT_NAME, COUNT_M_NAME, COUNT_TOTAL_NAME, DENSITY_NAME
        ]).unwrap();

        // Need equals_missing because of NaN and potential categorical differences
        assert_eq!(collected_df, expected_batch.into());
    }
}