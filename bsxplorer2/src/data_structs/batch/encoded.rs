use super::traits::{colnames, BatchType, BsxTypeTag};
use crate::data_structs::batch::traits::BsxBatchMethods;
use itertools::Itertools;
use once_cell::sync::OnceCell;
use polars::prelude::*;

/// Encoded version of [BsxBatch]
///
/// Encodes
/// 1. `context` as [bool]
/// 2. `strand` as [bool]
///
/// Other types are also downcasted
#[derive(Debug, Clone)]
pub struct EncodedBsxBatch {
    data: DataFrame,
    chr: OnceCell<String>,
}

impl Eq for EncodedBsxBatch {}
impl PartialEq for EncodedBsxBatch {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.data == other.data
    }
}

impl EncodedBsxBatch {
    /// Returns a reference to the categorical chromosome column
    pub fn chr(&self) -> &CategoricalChunked {
        self.data
            .column(colnames::CHR_NAME)
            .unwrap()
            .categorical()
            .unwrap()
    }

    /// Creates a new EncodedBsxBatch with pre-initialized fields
    pub(crate) fn new_with_fields(
        data: DataFrame,
        chr: String,
    ) -> EncodedBsxBatch {
        let chr = OnceCell::with_value(chr);
        Self { data, chr }
    }
}

impl From<EncodedBsxBatch> for DataFrame {
    fn from(batch: EncodedBsxBatch) -> Self {
        batch.data
    }
}

impl BsxBatchMethods for EncodedBsxBatch {
    type ChrType = UInt32Type;
    type PosType = UInt32Type;
    type StrandType = BooleanType;
    type ContextType = BooleanType;
    type CountType = Int16Type;
    type DensityType = Float32Type;

    fn chr(&self) -> &ChunkedArray<Self::ChrType> {
        self.data
            .column(colnames::CHR_NAME)
            .unwrap()
            .categorical()
            .unwrap()
            .physical()
    }
    fn position(&self) -> &ChunkedArray<Self::PosType> {
        self.data
            .column(colnames::POS_NAME)
            .unwrap()
            .u32()
            .unwrap()
    }
    fn strand(&self) -> &ChunkedArray<Self::StrandType> {
        self.data
            .column(colnames::STRAND_NAME)
            .unwrap()
            .bool()
            .unwrap()
    }

    fn context(&self) -> &ChunkedArray<Self::ContextType> {
        self.data
            .column(colnames::CONTEXT_NAME)
            .unwrap()
            .bool()
            .unwrap()
    }

    fn count_m(&self) -> &ChunkedArray<Self::CountType> {
        self.data
            .column(colnames::COUNT_M_NAME)
            .unwrap()
            .i16()
            .unwrap()
    }

    fn count_total(&self) -> &ChunkedArray<Self::CountType> {
        self.data
            .column(colnames::COUNT_TOTAL_NAME)
            .unwrap()
            .i16()
            .unwrap()
    }

    fn density(&self) -> &ChunkedArray<Self::DensityType> {
        self.data
            .column(colnames::DENSITY_NAME)
            .unwrap()
            .f32()
            .unwrap()
    }

    unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        EncodedBsxBatch {
            data: data_frame,
            chr: OnceCell::new(),
        }
    }

    fn data(&self) -> &DataFrame {
        &self.data
    }

    fn data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }

    fn take(self) -> DataFrame
    where
        Self: Sized,
    {
        self.data
    }

    fn chr_val(&self) -> anyhow::Result<&str> {
        let val = self
            .chr
            .get_or_try_init(|| {
                if self.data.is_empty() {
                    return Err(anyhow::anyhow!("no data"));
                }
                let rev_map = self.chr().get_rev_map();
                let first_value = self.chr().physical().first();
                first_value
                    .map(|v| rev_map.get(v).to_string())
                    .ok_or(anyhow::anyhow!("No data!"))
            })?
            .as_str();
        Ok(val)
    }
}

impl BsxTypeTag for EncodedBsxBatch {
    fn type_name() -> &'static str {
        "encoded"
    }
    fn type_enum() -> BatchType {
        BatchType::Encoded
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::batch::traits::colnames::*;
    use polars::df;
    use polars::prelude::*;

    // Helper function to create a sample DataFrame suitable for EncodedBsxBatch
    fn create_test_df() -> DataFrame {
        // Raw data similar to decoded, but will be converted
        let chr_series =
            Series::new(CHR_NAME.into(), &["chr1", "chr1", "chr1"])
                .cast(&DataType::Categorical(
                    None,
                    CategoricalOrdering::Physical,
                ))
                .unwrap();
        let pos_series = Series::new(POS_NAME.into(), &[10u32, 20, 30]);
        // Encode strand: "+" -> true, "-" -> false
        let strand_series =
            Series::new(STRAND_NAME.into(), &[true, true, false]);
        // Encode context: "CG" -> true, others -> false
        let context_series =
            Series::new(CONTEXT_NAME.into(), &[true, false, false]);
        let count_m_series = Series::new(COUNT_M_NAME.into(), &[5i16, 10, 15]);
        let count_total_series =
            Series::new(COUNT_TOTAL_NAME.into(), &[10i16, 20, 30]);
        let density_series =
            Series::new(DENSITY_NAME.into(), &[0.5f32, 0.5, 0.5]);

        df!(
            colnames::CHR_NAME => chr_series,
            colnames::POS_NAME => pos_series,
            colnames::STRAND_NAME => strand_series,
            colnames::CONTEXT_NAME => context_series,
            colnames::COUNT_M_NAME => count_m_series,
            colnames::COUNT_TOTAL_NAME => count_total_series,
            colnames::DENSITY_NAME => density_series,
        )
        .unwrap()
    }

    // Helper function to create an EncodedBsxBatch (unchecked for simplicity in tests)
    fn create_test_batch() -> EncodedBsxBatch {
        unsafe { EncodedBsxBatch::new_unchecked(create_test_df()) }
    }

    // Test the specific `chr` getter returning CategoricalChunked
    #[test]
    fn test_chr_categorical() {
        let batch = create_test_batch();
        let expected_dtype =
            DataType::Categorical(None, CategoricalOrdering::Physical);
        assert_eq!(batch.chr().dtype(), &expected_dtype);
        assert_eq!(batch.chr().len(), 3);
    }

    // Test BsxBatchMethods::chr returning physical representation
    #[test]
    fn test_bsx_chr() {
        let batch = create_test_batch();
        // Assuming "chr1" is the first category, its physical representation is 0
        let expected = Series::new(CHR_NAME.into(), &[0u32, 0, 0]);
        assert_eq!(
            <EncodedBsxBatch as BsxBatchMethods>::chr(&batch)
                .clone()
                .into_series(),
            expected
        );
    }

    #[test]
    fn test_position() {
        let batch = create_test_batch();
        let expected = Series::new(POS_NAME.into(), &[10u32, 20, 30]);
        assert_eq!(batch.position().clone().into_series(), expected);
    }

    #[test]
    fn test_strand() {
        let batch = create_test_batch();
        let expected =
            Series::new(STRAND_NAME.into(), &[true, true, false]);
        assert_eq!(batch.strand().clone().into_series(), expected);
    }

    #[test]
    fn test_context() {
        let batch = create_test_batch();
        let expected =
            Series::new(CONTEXT_NAME.into(), &[true, false, false]);
        assert_eq!(batch.context().clone().into_series(), expected);
    }

    #[test]
    fn test_count_m() {
        let batch = create_test_batch();
        let expected =
            Series::new(COUNT_M_NAME.into(), &[5i16, 10, 15]);
        assert_eq!(batch.count_m().clone().into_series(), expected);
    }

    #[test]
    fn test_count_total() {
        let batch = create_test_batch();
        let expected =
            Series::new(COUNT_TOTAL_NAME.into(), &[10i16, 20, 30]);
        assert_eq!(
            batch
                .count_total()
                .clone()
                .into_series(),
            expected
        );
    }

    #[test]
    fn test_density() {
        let batch = create_test_batch();
        let expected =
            Series::new(DENSITY_NAME.into(), &[0.5f32, 0.5, 0.5]);
        assert_eq!(batch.density().clone().into_series(), expected);
    }

    #[test]
    fn test_new_unchecked() {
        let df = create_test_df();
        let batch = unsafe { EncodedBsxBatch::new_unchecked(df.clone()) };
        assert_eq!(batch.data(), &df);
        assert!(batch.chr.get().is_none()); // Chr cache should not be initialized
    }

    #[test]
    fn test_new_with_fields() {
        let df = create_test_df();
        let chr_str = "chr1".to_string();
        let batch =
            EncodedBsxBatch::new_with_fields(df.clone(), chr_str.clone());
        assert_eq!(batch.data(), &df);
        assert!(batch.chr.get().is_some());
        assert_eq!(batch.chr.get().unwrap(), &chr_str);
    }

    #[test]
    fn test_data() {
        let df = create_test_df();
        let batch = unsafe { EncodedBsxBatch::new_unchecked(df.clone()) };
        assert_eq!(batch.data(), &df);
    }

    #[test]
    fn test_data_mut() {
        let df = create_test_df();
        let mut batch = unsafe { EncodedBsxBatch::new_unchecked(df.clone()) };
        let df_mut = batch.data_mut();
        // Example modification
        *df_mut = df!(CHR_NAME => Series::new(CHR_NAME.into(), &["chrX"])
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap())
        .unwrap();
        assert_eq!(batch.data().shape(), (1, 1));
        // Check if the modification took place (may need more specific checks)
    }

    #[test]
    fn test_take() {
        let df = create_test_df();
        let batch = unsafe { EncodedBsxBatch::new_unchecked(df.clone()) };
        let taken_df = batch.take();
        assert_eq!(taken_df, df);
    }

    #[test]
    fn test_chr_val_ok_uncached() {
        let batch = create_test_batch(); // Created with new_unchecked, so cache is empty
        assert!(batch.chr.get().is_none());
        // First call - calculates from rev_map
        assert_eq!(batch.chr_val().unwrap(), "chr1");
        // Second call - should now use cached value
        assert!(batch.chr.get().is_some());
        assert_eq!(batch.chr_val().unwrap(), "chr1");
    }

    #[test]
    fn test_chr_val_ok_cached() {
        let df = create_test_df();
        let chr_str = "chr1".to_string();
        let batch = EncodedBsxBatch::new_with_fields(df, chr_str.clone());
        assert!(batch.chr.get().is_some());
        assert_eq!(batch.chr_val().unwrap(), "chr1");
        // Ensure it's still the same cached value
        assert_eq!(batch.chr.get().unwrap(), &chr_str);
    }

    #[test]
    fn test_chr_val_empty() {
        let empty_df = df!(
             colnames::CHR_NAME => Series::new_empty(colnames::CHR_NAME.into(), &DataType::Categorical(None, CategoricalOrdering::Physical)),
            colnames::POS_NAME => Series::new_empty(colnames::POS_NAME.into(), &DataType::UInt32),
            colnames::STRAND_NAME => Series::new_empty(colnames::STRAND_NAME.into(), &DataType::Boolean),
            colnames::CONTEXT_NAME => Series::new_empty(colnames::CONTEXT_NAME.into(), &DataType::Boolean),
            colnames::COUNT_M_NAME => Series::new_empty(colnames::COUNT_M_NAME.into(), &DataType::Int16),
            colnames::COUNT_TOTAL_NAME => Series::new_empty(colnames::COUNT_TOTAL_NAME.into(), &DataType::Int16),
            colnames::DENSITY_NAME => Series::new_empty(colnames::DENSITY_NAME.into(), &DataType::Float32),
        )
        .unwrap();
        let batch = unsafe { EncodedBsxBatch::new_unchecked(empty_df) };
        assert!(batch.chr_val().is_err());
        assert!(batch
            .chr_val()
            .unwrap_err()
            .to_string()
            .contains("no data"));
    }

    #[test]
    fn test_from_encoded_bsxbatch() {
        let df = create_test_df();
        let batch = unsafe { EncodedBsxBatch::new_unchecked(df.clone()) };
        let converted_df: DataFrame = batch.into();
        assert_eq!(converted_df, df);
    }

    #[test]
    fn test_type_tag() {
        assert_eq!(EncodedBsxBatch::type_name(), "encoded");
        assert_eq!(EncodedBsxBatch::type_enum(), BatchType::Encoded);
    }
}
