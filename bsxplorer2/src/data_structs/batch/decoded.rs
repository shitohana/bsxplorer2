use crate::data_structs::batch::builder::BsxBatchBuilder;
use crate::data_structs::batch::traits::{colnames, BsxTypeTag};
use crate::data_structs::batch::traits::{BatchType, BsxBatchMethods};
use anyhow::anyhow;
use once_cell::sync::OnceCell;
use polars::prelude::*;

/// A batch of BSX data stored in a DataFrame with the following guarantees:
/// 1. Non-null chr and position columns
/// 2. Single unique chromosome value
/// 3. Ascending sorted positions
#[derive(Clone, Debug)]
pub struct BsxBatch {
    data: DataFrame,
    chr: OnceCell<String>,
}

impl Eq for BsxBatch {}
impl PartialEq for BsxBatch {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.data == other.data
    }
}

impl BsxBatchMethods for BsxBatch {
    type ChrType = StringType;
    type PosType = UInt64Type;
    type StrandType = StringType;
    type ContextType = StringType;
    type CountType = UInt32Type;
    type DensityType = Float64Type;

    fn chr(&self) -> &ChunkedArray<Self::ChrType> {
        self.data
            .column(colnames::CHR_NAME)
            .unwrap()
            .str()
            .unwrap()
    }
    fn position(&self) -> &ChunkedArray<Self::PosType> {
        self.data
            .column(colnames::POS_NAME)
            .unwrap()
            .u64()
            .unwrap()
    }
    fn strand(&self) -> &ChunkedArray<Self::StrandType> {
        self.data
            .column(colnames::STRAND_NAME)
            .unwrap()
            .str()
            .unwrap()
    }
    fn context(&self) -> &ChunkedArray<Self::ContextType> {
        self.data
            .column(colnames::CONTEXT_NAME)
            .unwrap()
            .str()
            .unwrap()
    }
    fn count_m(&self) -> &ChunkedArray<Self::CountType> {
        self.data
            .column(colnames::COUNT_M_NAME)
            .unwrap()
            .u32()
            .unwrap()
    }
    fn count_total(&self) -> &ChunkedArray<Self::CountType> {
        self.data
            .column(colnames::COUNT_TOTAL_NAME)
            .unwrap()
            .u32()
            .unwrap()
    }

    fn density(&self) -> &ChunkedArray<Self::DensityType> {
        self.data
            .column(colnames::DENSITY_NAME)
            .unwrap()
            .f64()
            .unwrap()
    }

    unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        BsxBatch {
            data: data_frame,
            chr: Default::default(),
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
        self.chr
            .get_or_try_init(|| {
                if self.data.is_empty() {
                    return Err(anyhow!("no data"));
                }
                let first = self
                    .chr()
                    .first()
                    .ok_or_else(|| anyhow!("no data"))?;
                Ok(first.to_string())
            })
            .map(String::as_str)
    }
}

/// Implementation to create a BsxBatch from a DataFrame with validation
impl TryFrom<DataFrame> for BsxBatch {
    type Error = anyhow::Error;

    fn try_from(value: DataFrame) -> Result<Self, Self::Error> {
        BsxBatchBuilder::default()
            .with_check_duplicates(true)
            .with_rechunk(true)
            .with_check_nulls(true)
            .with_check_sorted(true)
            .with_check_single_chr(true)
            .build(value)
    }
}

/// Implementation to convert a BsxBatch back to its inner DataFrame
impl From<BsxBatch> for DataFrame {
    fn from(b: BsxBatch) -> Self {
        b.data
    }
}

impl BsxTypeTag for BsxBatch {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn type_name() -> &'static str {
        "decoded"
    }
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn type_enum() -> BatchType {
        BatchType::Decoded
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::batch::traits::colnames;
    use polars::df;

    // Helper function to create a sample DataFrame
    fn create_test_df() -> DataFrame {
        df!(
            colnames::CHR_NAME => &["chr1", "chr1", "chr1"],
            colnames::POS_NAME => &[10u64, 20, 30],
            colnames::STRAND_NAME => &["+", "+", "-"],
            colnames::CONTEXT_NAME => &["CG", "CHG", "CHH"],
            colnames::COUNT_M_NAME => &[5u32, 10, 15],
            colnames::COUNT_TOTAL_NAME => &[10u32, 20, 30],
            colnames::DENSITY_NAME => &[0.5f64, 0.5, 0.5],
        )
        .unwrap()
    }

    // Helper function to create a BsxBatch (unchecked for simplicity in tests)
    fn create_test_batch() -> BsxBatch {
        unsafe { BsxBatch::new_unchecked(create_test_df()) }
    }

    #[test]
    fn test_chr() {
        let batch = create_test_batch();
        let expected =
            Series::new(colnames::CHR_NAME.into(), &["chr1", "chr1", "chr1"]);
        assert_eq!(batch.chr().clone().into_series(), expected);
    }

    #[test]
    fn test_position() {
        let batch = create_test_batch();
        let expected = Series::new(colnames::POS_NAME.into(), &[10u64, 20, 30]);
        assert_eq!(batch.position().clone().into_series(), expected);
    }

    #[test]
    fn test_strand() {
        let batch = create_test_batch();
        let expected =
            Series::new(colnames::STRAND_NAME.into(), &["+", "+", "-"]);
        assert_eq!(batch.strand().clone().into_series(), expected);
    }

    #[test]
    fn test_context() {
        let batch = create_test_batch();
        let expected =
            Series::new(colnames::CONTEXT_NAME.into(), &["CG", "CHG", "CHH"]);
        assert_eq!(batch.context().clone().into_series(), expected);
    }

    #[test]
    fn test_count_m() {
        let batch = create_test_batch();
        let expected =
            Series::new(colnames::COUNT_M_NAME.into(), &[5u32, 10, 15]);
        assert_eq!(batch.count_m().clone().into_series(), expected);
    }

    #[test]
    fn test_count_total() {
        let batch = create_test_batch();
        let expected =
            Series::new(colnames::COUNT_TOTAL_NAME.into(), &[10u32, 20, 30]);
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
            Series::new(colnames::DENSITY_NAME.into(), &[0.5f64, 0.5, 0.5]);
        assert_eq!(batch.density().clone().into_series(), expected);
    }

    #[test]
    fn test_take() {
        let df = create_test_df();
        let batch = unsafe { BsxBatch::new_unchecked(df.clone()) };
        let taken_df = batch.take();
        assert_eq!(taken_df, df);
    }

    #[test]
    fn test_chr_val_ok() {
        let batch = create_test_batch();
        // First call
        assert_eq!(batch.chr_val().unwrap(), "chr1");
        // Second call (should use cached value)
        assert_eq!(batch.chr_val().unwrap(), "chr1");
    }

    #[test]
    fn test_chr_val_empty() {
        let empty_df = df!(
            colnames::CHR_NAME => Series::new_empty(colnames::CHR_NAME.into(), &DataType::String),
            colnames::POS_NAME => Series::new_empty(colnames::POS_NAME.into(), &DataType::UInt64),
            colnames::STRAND_NAME => Series::new_empty(colnames::STRAND_NAME.into(), &DataType::String),
            colnames::CONTEXT_NAME => Series::new_empty(colnames::CONTEXT_NAME.into(), &DataType::String),
            colnames::COUNT_M_NAME => Series::new_empty(colnames::COUNT_M_NAME.into(), &DataType::UInt32),
            colnames::COUNT_TOTAL_NAME => Series::new_empty(colnames::COUNT_TOTAL_NAME.into(), &DataType::UInt32),
            colnames::DENSITY_NAME => Series::new_empty(colnames::DENSITY_NAME.into(), &DataType::Float64),
        ).unwrap();
        let batch = unsafe { BsxBatch::new_unchecked(empty_df) };
        assert!(batch.chr_val().is_err());
        // Check error kind if possible/needed, e.g., contains("no data")
        assert!(batch
            .chr_val()
            .unwrap_err()
            .to_string()
            .contains("no data"));
    }

    #[test]
    fn test_try_from_valid() {
        let df = create_test_df();
        let batch = BsxBatch::try_from(df);
        assert!(batch.is_ok());
        assert_eq!(batch.unwrap().chr_val().unwrap(), "chr1");
    }

    #[test]
    fn test_try_from_invalid_chromosome() {
        let df = df!(
            colnames::CHR_NAME => &["chr1", "chr2", "chr1"], // Multiple chromosomes
            colnames::POS_NAME => &[10u64, 20, 30],
            colnames::STRAND_NAME => &["+", "+", "-"],
            colnames::CONTEXT_NAME => &["CG", "CHG", "CHH"],
            colnames::COUNT_M_NAME => &[5u32, 10, 15],
            colnames::COUNT_TOTAL_NAME => &[10u32, 20, 30],
            colnames::DENSITY_NAME => &[0.5f64, 0.5, 0.5],
        )
        .unwrap();
        let batch = BsxBatch::try_from(df);
        assert!(batch.is_err());
        assert!(batch
            .unwrap_err()
            .to_string()
            .contains("Multiple chromosomes"));
    }
}
