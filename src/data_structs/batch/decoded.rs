use crate::data_structs::batch::builder::BsxBatchBuilder;
use crate::data_structs::batch::traits::BsxBatchMethods;
use crate::data_structs::batch::traits::{colnames, BsxTypeTag};
use anyhow::anyhow;
use once_cell::sync::OnceCell;
use polars::prelude::*;

/// A batch of BSX data stored in a DataFrame with the following guarantees:
/// 1. Non-null chr and position columns
/// 2. Single unique chromosome value
/// 3. Ascending sorted positions
#[derive(Clone, PartialEq, Debug)]
pub struct BsxBatch {
    data: DataFrame,
    chr: OnceCell<String>
}

impl BsxBatchMethods for BsxBatch where {
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
            chr: Default::default()
        }
    }

    fn data(&self) -> &DataFrame {
        &self.data
    }

    fn chr_val(&self) -> anyhow::Result<&str> {
        self.chr.get_or_try_init(|| {
                let first = self.chr().first().ok_or_else(|| anyhow!("no data"))?;
                Ok(first.to_string())
            }).map(String::as_str)

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
            .build_decoded(value)
    }
}

/// Implementation to convert a BsxBatch back to its inner DataFrame
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
