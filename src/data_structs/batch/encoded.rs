use super::traits::{colnames, BsxTypeTag};
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
#[derive(Debug, Clone, PartialEq)]
pub struct EncodedBsxBatch {
    data: DataFrame,
    chr: OnceCell<String>,
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
        Self {
            data,
            chr,
        }
    }
}

impl From<EncodedBsxBatch> for DataFrame {
    fn from(batch: EncodedBsxBatch) -> Self {
        batch.data
    }
}

impl BsxBatchMethods for EncodedBsxBatch where {
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

    fn chr_val(&self) -> anyhow::Result<&str> {
        let val = self
            .chr
            .get_or_try_init(|| {
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
}
