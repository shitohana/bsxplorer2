//! Module contains BsxBatches which represent single chromosome methylation
//! data_structs
//!
//! BsxBatch stores data_structs in [DataFrame] with
//! * Non-null chr and position
//! * Single unique chromosome
//! * Ascending sorted positions
//!
//! Module contains two main structs:
//! 1. [BsxBatch] struct for storing non-encoded data_structs
//! 2. [EncodedBsxBatch] struct for storing encoded data_structs (context and
//! strand encoded as boolean)
//!
//! and [BsxBatchMethods] trait
//!
//! Both BsxBatch and EncodedBsxBatch can be filtered using context and/or
//! strand
use crate::data_structs::batch::builder::BsxBatchBuilder;
use crate::data_structs::batch::traits::BsxTypeTag;
use crate::data_structs::batch::traits::{BsxBatchMethods, BsxColNames};
use crate::utils::types::{IPCEncodedEnum, Strand};
use anyhow::anyhow;
use once_cell::sync::OnceCell;
use polars::prelude::*;

/// DataFrame with
/// 1. Non-null chr and position
/// 2. Single unique chr
/// 3. Ascending sorted positions
#[derive(Clone, PartialEq, Debug)]
pub struct BsxBatch {
    data: DataFrame,
    chr: OnceCell<String>,
    start: OnceCell<u32>,
    end: OnceCell<u32>,
}

impl BsxBatch {
    pub const CHR_DTYPE: DataType = DataType::String;
    pub const POS_DTYPE: DataType = DataType::UInt64;
    pub const STRAND_DTYPE: DataType = DataType::String;
    pub const CONTEXT_DTYPE: DataType = DataType::String;
    pub const COUNT_DTYPE: DataType = DataType::UInt32;
    pub const DENSITY_DTYPE: DataType = DataType::Float64;
}

impl BsxBatch {
    /// Returns expected types of columns
    pub(crate) const fn col_types() -> &'static [DataType] {
        &[
            Self::CHR_DTYPE,
            Self::POS_DTYPE,
            Self::STRAND_DTYPE,
            Self::CONTEXT_DTYPE,
            Self::COUNT_DTYPE,
            Self::COUNT_DTYPE,
            Self::DENSITY_DTYPE,
        ]
    }

    pub fn empty() -> Self {
        unsafe { Self::new_unchecked(DataFrame::empty_with_schema(&BsxBatch::schema())) }
    }

    pub fn split_at(
        self,
        index: usize,
    ) -> (Self, Self) {
        let (a, b) = self.data.split_at(index as i64);
        unsafe {(
            Self::new_unchecked(a),
            Self::new_unchecked(b),
        )}
    }

    /// Returns expected schema of [BsxBatch]
    pub fn schema() -> Schema {
        use crate::data_structs::batch::traits::BsxColNames;
        use crate::utils::schema_from_arrays;
        schema_from_arrays(&Self::col_names(), Self::col_types())
    }

    /// Returns expected schema of [BsxBatch] as [PlHashMap]
    pub fn hashmap() -> PlHashMap<&'static str, DataType> {
        use crate::data_structs::batch::traits::BsxColNames;
        use crate::utils::hashmap_from_arrays;
        hashmap_from_arrays(&Self::col_names(), Self::col_types())
    }

    pub fn position(&self) -> &ChunkedArray<UInt64Type> {
        self.data
            .column(Self::POS_NAME)
            .unwrap()
            .u64()
            .unwrap()
    }

    pub fn chr(&self) -> &ChunkedArray<StringType> {
        self.data
            .column(Self::CHR_NAME)
            .unwrap()
            .str()
            .unwrap()
    }

    pub fn strand(&self) -> &ChunkedArray<StringType> {
        self.data
            .column(Self::STRAND_NAME)
            .unwrap()
            .str()
            .unwrap()
    }

    pub fn context(&self) -> &ChunkedArray<StringType> {
        self.data
            .column(Self::CONTEXT_NAME)
            .unwrap()
            .str()
            .unwrap()
    }

    pub fn count_m(&self) -> &ChunkedArray<UInt32Type> {
        self.data
            .column(Self::COUNT_M_NAME)
            .unwrap()
            .u32()
            .unwrap()
    }

    pub fn count_total(&self) -> &ChunkedArray<UInt32Type> {
        self.data
            .column(Self::COUNT_TOTAL_NAME)
            .unwrap()
            .u32()
            .unwrap()
    }

    pub fn density(&self) -> &ChunkedArray<Float64Type> {
        self.data
            .column(Self::DENSITY_NAME)
            .unwrap()
            .f64()
            .unwrap()
    }
}

impl BsxColNames for BsxBatch {}

impl BsxBatchMethods for BsxBatch {
    unsafe fn new_unchecked(data_frame: DataFrame) -> Self {
        BsxBatch {
            data: data_frame,
            chr: Default::default(),
            start: Default::default(),
            end: Default::default(),
        }
    }
    /// Filters [BsxBatch] by `context` and `strand`
    ///
    /// If both `context` and `strand` are [None], returns [BsxBatch]
    fn filter(
        self,
        expr: Expr,
    ) -> anyhow::Result<Self>
    where
        Self: Sized,
    {
        let df = self.data.lazy().filter(expr).collect()?;
        BsxBatchBuilder::no_checks().build_decoded(df)
    }

    /// Returns reference to inner [DataFrame]
    #[inline]
    fn data(&self) -> &DataFrame {
        &self.data
    }

    fn chr_val(&self) -> anyhow::Result<&str> {
        self.chr.get_or_try_init(|| {
                let first = self.chr().first().ok_or_else(|| anyhow!("no data"))?;
                Ok(first.to_string())
            }).map(String::as_str)

    }

    /// Returns mutable reference to inner [DataFrame]
    #[inline]
    fn data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }

    fn start_pos(&self) -> Option<u32> {
        self.start.get_or_try_init(|| self.position()
            .first()
            .map(|v| v as u32)
            .ok_or(anyhow!("no data"))
        ).ok().cloned()

    }
    fn end_pos(&self) -> Option<u32> {
        self.end.get_or_try_init(|| self.position()
            .last()
            .map(|v| v as u32)
            .ok_or(anyhow!("no data"))
        ).ok().cloned()
    }
}

impl TryFrom<DataFrame> for BsxBatch {
    type Error = anyhow::Error;

    fn try_from(value: DataFrame) -> Result<Self, Self::Error> {
        BsxBatchBuilder::default()
            .with_check_duplicates(true)
            .with_rechunk(true)
            .with_check_nulls(true)
            .with_check_sorted(true)
            .build_decoded(value)
    }
}

impl From<BsxBatch> for DataFrame {
    fn from(b: BsxBatch) -> Self {
        b.data
    }
}

impl BsxTypeTag for BsxBatch {
    fn type_name() -> &'static str {
        "decoded"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_structs::batch::traits::BsxBatchMethods;
    use crate::data_structs::region::RegionCoordinates;
    use crate::utils::get_categorical_dtype;

    fn dummy_batch() -> BsxBatch {
        BsxBatch::try_from(
            df![
                "chr" => ["1", "1", "1"],
                "position" => [1, 2, 3],
                "strand" => [".", ".", "."],
                "context" => ["CG", "CHG", "CHH"],
                "count_m" => [0,0,0],
                "count_total" => [0,0,0],
                "density" => [0, 1, 0]
            ]
            .unwrap()
            .lazy()
            .cast(BsxBatch::hashmap(), true)
            .collect()
            .unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn test_position() {
        let batch = dummy_batch();
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded =
            BsxBatchBuilder::encode_batch(batch.clone(), chr_dtype).unwrap();

        assert_eq!(batch.start_pos().unwrap(), encoded.start_pos().unwrap());
        assert_eq!(batch.end_pos().unwrap(), encoded.end_pos().unwrap());
    }

    #[test]
    fn test_encode_decode() {
        let batch = dummy_batch();
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded =
            BsxBatchBuilder::encode_batch(batch.clone(), chr_dtype).unwrap();
        let decoded = BsxBatchBuilder::decode_batch(encoded).unwrap();
        assert_eq!(batch, decoded);
    }

    #[test]
    fn test_trim() {
        let batch = BsxBatch::try_from(
            df![
                "chr" => ["1", "1", "1", "1", "1", "1", "1"],
                "position" => [2, 3, 4, 5, 7, 9, 10],
                "strand" => [".", ".", ".", ".", ".", ".", "."],
                "context" => ["CG", "CHG", "CHH", "CG", "CHG", "CHH", "CG"],
                "count_m" => [0,0,0,0,0,0,0],
                "count_total" => [0,0,0,0,0,0,0],
                "density" => [0, 1, 0, 0, 1, 0, 0]
            ]
            .unwrap()
            .lazy()
            .cast(BsxBatch::hashmap(), true)
            .collect()
            .unwrap(),
        )
        .unwrap();
        let chr_dtype = get_categorical_dtype(vec!["1".into()]);
        let encoded =
            BsxBatchBuilder::encode_batch(batch.clone(), chr_dtype).unwrap();

        let trimmed = encoded
            .trim_region(&RegionCoordinates::new(
                "1".to_string(),
                4,
                9,
                Strand::None,
            ))
            .unwrap();
        assert_eq!(trimmed.start_pos().unwrap(), 4);
        assert_eq!(trimmed.end_pos().unwrap(), 9);

        let trimmed = encoded
            .trim_region(&RegionCoordinates::new(
                "1".to_string(),
                1,
                9,
                Strand::None,
            ))
            .unwrap();
        assert_eq!(trimmed.start_pos().unwrap(), 2);
        assert_eq!(trimmed.end_pos().unwrap(), 9);

        let trimmed = encoded
            .trim_region(&RegionCoordinates::new(
                "1".to_string(),
                4,
                12,
                Strand::None,
            ))
            .unwrap();
        assert_eq!(trimmed.start_pos().unwrap(), 4);
        assert_eq!(trimmed.end_pos().unwrap(), 10);
    }
}
