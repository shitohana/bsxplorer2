use crate::data_structs::coords::{Contig, GenomicPosition};
use crate::data_structs::enums::Strand;

use super::builder::BsxBatchBuilder;
use anyhow::anyhow;
use polars::datatypes::BooleanChunked;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::prelude::*;

pub mod colnames {
    pub const CHR_NAME: &str = "chr";
    /// Position column name
    pub const POS_NAME: &str = "position";
    /// Strand column name
    pub const STRAND_NAME: &str = "strand";
    /// Context column name
    pub const CONTEXT_NAME: &str = "context";
    /// Methylated count column name
    pub const COUNT_M_NAME: &str = "count_m";
    /// Total count column name
    pub const COUNT_TOTAL_NAME: &str = "count_total";
    /// Density column name
    pub const DENSITY_NAME: &str = "density";

    pub const fn col_names() -> [&'static str; 7] {
        [
            CHR_NAME,
            POS_NAME,
            STRAND_NAME,
            CONTEXT_NAME,
            COUNT_M_NAME,
            COUNT_TOTAL_NAME,
            DENSITY_NAME,
        ]
    }
}
use colnames::*;
use crate::data_structs::batch::LazyBsxBatch;

/// Trait for common methods for [BsxBatch] and [EncodedBsxBatch]
pub trait BsxBatchMethods: BsxTypeTag + Eq + PartialEq {
    /// Type for chromosome data
    type ChrType: PolarsDataType;
    /// Type for position data
    type PosType: PolarsIntegerType;
    /// Type for strand data
    type StrandType: PolarsDataType;
    /// Type for context data
    type ContextType: PolarsDataType;
    /// Type for count data
    type CountType: PolarsIntegerType;
    /// Type for density data
    type DensityType: PolarsFloatType;

    /// Access chromosome column
    fn chr(&self) -> &ChunkedArray<Self::ChrType>;
    fn position(&self) -> &ChunkedArray<Self::PosType>;
    fn strand(&self) -> &ChunkedArray<Self::StrandType>;
    fn context(&self) -> &ChunkedArray<Self::ContextType>;
    fn count_m(&self) -> &ChunkedArray<Self::CountType>;
    fn count_total(&self) -> &ChunkedArray<Self::CountType>;

    fn density(&self) -> &ChunkedArray<Self::DensityType>;

    /// Get chromosome data type
    fn chr_type() -> DataType
    where
        Self: Sized,
    {
        Self::ChrType::get_dtype()
    }
    /// Get position data type
    fn pos_type() -> DataType
    where
        Self: Sized,
    {
        Self::PosType::get_dtype()
    }
    /// Get strand data type
    fn strand_type() -> DataType
    where
        Self: Sized,
    {
        Self::StrandType::get_dtype()
    }
    /// Get context data type
    fn context_type() -> DataType
    where
        Self: Sized,
    {
        Self::ContextType::get_dtype()
    }
    /// Get count data type
    fn count_type() -> DataType
    where
        Self: Sized,
    {
        Self::CountType::get_dtype()
    }
    /// Get density data type
    fn density_type() -> DataType
    where
        Self: Sized,
    {
        Self::DensityType::get_dtype()
    }

    /// Create schema for the batch data
    fn schema() -> Schema
    where
        Self: Sized,
    {
        Schema::from_iter([
            (CHR_NAME.into(), Self::ChrType::get_dtype()),
            (POS_NAME.into(), Self::PosType::get_dtype()),
            (STRAND_NAME.into(), Self::StrandType::get_dtype()),
            (CONTEXT_NAME.into(), Self::ContextType::get_dtype()),
            (COUNT_M_NAME.into(), Self::CountType::get_dtype()),
            (COUNT_TOTAL_NAME.into(), Self::CountType::get_dtype()),
            (DENSITY_NAME.into(), Self::DensityType::get_dtype()),
        ])
    }

    /// Create hashmap of column names to data types
    fn hashmap() -> PlHashMap<&'static str, DataType>
    where
        Self: Sized,
    {
        PlHashMap::from_iter([
            (CHR_NAME, Self::ChrType::get_dtype()),
            (POS_NAME, Self::PosType::get_dtype()),
            (STRAND_NAME, Self::StrandType::get_dtype()),
            (CONTEXT_NAME, Self::ContextType::get_dtype()),
            (COUNT_M_NAME, Self::CountType::get_dtype()),
            (COUNT_TOTAL_NAME, Self::CountType::get_dtype()),
            (DENSITY_NAME, Self::DensityType::get_dtype()),
        ])
    }

    /// Create a new batch from a DataFrame without checks
    ///
    /// # Safety
    ///
    /// This function assumes that the DataFrame is valid and that the columns are of the correct type.
    unsafe fn new_unchecked(data_frame: DataFrame) -> Self
    where
        Self: Sized;

    /// Create an empty batch
    fn empty() -> Self
    where
        Self: Sized,
    {
        unsafe {
            Self::new_unchecked(DataFrame::empty_with_schema(&Self::schema()))
        }
    }
    fn is_empty(&self) -> bool {
        self.data().is_empty()
    }

    /// Split batch at specified index
    fn split_at(
        self,
        index: usize,
    ) -> (Self, Self)
    where
        Self: Sized,
    {
        let (a, b) = self.data().split_at(index as i64);
        unsafe { (Self::new_unchecked(a), Self::new_unchecked(b)) }
    }

    /// Returns reference to inner [DataFrame]
    fn data(&self) -> &DataFrame;
    fn data_mut(&mut self) -> &mut DataFrame;

    fn take(self) -> DataFrame
    where
        Self: Sized;

    fn lazy(self) -> LazyBsxBatch<Self> where Self: Sized {
        LazyBsxBatch::from(self)
    }
    /// Get chromosome value as string
    fn chr_val(&self) -> anyhow::Result<&str>;

    /// Get start position
    fn start_pos(&self) -> Option<u32> {
        self.data()
            .column(POS_NAME)
            .unwrap()
            .get(self.data().height() - 1)
            .map(|v| {
                v.cast(&DataType::UInt32)
                    .try_extract()
                    .unwrap()
            })
            .ok()
    }
    /// Get end position
    fn end_pos(&self) -> Option<u32> {
        self.data()
            .column(POS_NAME)
            .unwrap()
            .get(self.height() - 1)
            .map(|v| {
                v.cast(&DataType::UInt32)
                    .try_extract()
                    .unwrap()
            })
            .ok()
    }

    fn start_gpos(&self) -> anyhow::Result<GenomicPosition<&str, u32>> {
        let pos = self
            .start_pos()
            .ok_or(anyhow!("no data"))?;
        let chr = self.chr_val()?;
        let gpos = GenomicPosition::new(chr, pos);
        Ok(gpos)
    }

    fn end_gpos(&self) -> anyhow::Result<GenomicPosition<&str, u32>> {
        let pos = self
            .end_pos()
            .ok_or(anyhow!("no data"))?;
        let chr = self.chr_val()?;
        let gpos = GenomicPosition::new(chr, pos);
        Ok(gpos)
    }

    /// Vertically stack with another batch
    fn vstack(
        &self,
        other: &Self,
    ) -> anyhow::Result<Self>
    where
        Self: Sized,
    {
        let new_data = self.data().vstack(other.data())?;
        let res = BsxBatchBuilder::all_checks().check_modify(new_data)?;
        Ok(unsafe { Self::new_unchecked(res) })
    }

    fn extend(
        &mut self,
        other: &Self,
    ) -> anyhow::Result<()> {
        self.data_mut().extend(other.data())?;
        BsxBatchBuilder::all_checks().check(self.data())?;
        Ok(())
    }

    /// Filter batch based on boolean mask
    fn filter_mask(
        &self,
        mask: &BooleanChunked,
    ) -> PolarsResult<Self>
    where
        Self: Sized,
    {
        Ok(unsafe { Self::new_unchecked(self.data().filter(mask)?) })
    }

    /// Returns number of rows
    fn height(&self) -> usize {
        self.data().height()
    }

    /// Convert batch to genomic contig
    fn as_contig(&self) -> anyhow::Result<Contig<String, u32>> {
        let start = self
            .start_pos()
            .ok_or(anyhow!("no data"))?;
        let end = self
            .end_pos()
            .ok_or(anyhow!("no data"))?;
        let chr = self.chr_val()?;

        Ok(Contig::new(chr.to_owned(), start, end + 1, Strand::None))
    }
}

/// Trait for BSX type identification
pub trait BsxTypeTag {
    /// Get type name as string
    fn type_name() -> &'static str;
    fn type_enum() -> BatchType;
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum BatchType {
    Decoded,
    Encoded,
}
