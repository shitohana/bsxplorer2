use super::builder::BsxBatchBuilder;
use crate::data_structs::batch::decoded::BsxBatch;
use anyhow::anyhow;
use bio_types::annot::contig::Contig;
use bio_types::strand::NoStrand;
use polars::datatypes::BooleanChunked;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::prelude::*;

/// Defines column names for BSX data structures
pub trait BsxColNames {
    /// Chromosome column name
    const CHR_NAME: &'static str = "chr";
    /// Position column name
    const POS_NAME: &'static str = "position";
    /// Strand column name
    const STRAND_NAME: &'static str = "strand";
    /// Context column name
    const CONTEXT_NAME: &'static str = "context";
    /// Methylated count column name
    const COUNT_M_NAME: &'static str = "count_m";
    /// Total count column name
    const COUNT_TOTAL_NAME: &'static str = "count_total";
    /// Density column name
    const DENSITY_NAME: &'static str = "density";

    /// Returns an array of all column names
    fn col_names() -> [&'static str; 7] {
        [
            Self::CHR_NAME,
            Self::POS_NAME,
            Self::STRAND_NAME,
            Self::CONTEXT_NAME,
            Self::COUNT_M_NAME,
            Self::COUNT_TOTAL_NAME,
            Self::DENSITY_NAME,
        ]
    }

}

/// Trait for common methods for [BsxBatch] and [EncodedBsxBatch]
pub trait BsxBatchMethods
{
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
    /// Access position column
    fn position(&self) -> &ChunkedArray<Self::PosType>;
    /// Access strand column
    fn strand(&self) -> &ChunkedArray<Self::StrandType>;
    /// Access context column
    fn context(&self) -> &ChunkedArray<Self::ContextType>;
    /// Access methylated count column
    fn count_m(&self) -> &ChunkedArray<Self::CountType>;
    /// Access total count column
    fn count_total(&self) -> &ChunkedArray<Self::CountType>;
    /// Access density column
    fn density(&self) -> &ChunkedArray<Self::DensityType>;

    /// Get chromosome data type
    fn chr_type() -> DataType where Self: Sized {
        Self::ChrType::get_dtype()
    }
    /// Get position data type
    fn pos_type() -> DataType where Self: Sized {
        Self::PosType::get_dtype()
    }
    /// Get strand data type
    fn strand_type() -> DataType where Self: Sized {
        Self::StrandType::get_dtype()
    }
    /// Get context data type
    fn context_type() -> DataType where Self: Sized {
        Self::ContextType::get_dtype()
    }
    /// Get count data type
    fn count_type() -> DataType where Self: Sized {
        Self::CountType::get_dtype()
    }
    /// Get density data type
    fn density_type() -> DataType where Self: Sized {
        Self::DensityType::get_dtype()
    }

    /// Create schema for the batch data
    fn schema() -> Schema where Self: Sized + BsxColNames {
        Schema::from_iter([
            (Self::CHR_NAME.into(), Self::ChrType::get_dtype()),
            (Self::POS_NAME.into(), Self::PosType::get_dtype()),
            (Self::STRAND_NAME.into(), Self::StrandType::get_dtype()),
            (Self::CONTEXT_NAME.into(), Self::ContextType::get_dtype()),
            (Self::COUNT_M_NAME.into(), Self::CountType::get_dtype()),
            (Self::COUNT_TOTAL_NAME.into(), Self::CountType::get_dtype()),
            (Self::DENSITY_NAME.into(), Self::DensityType::get_dtype()),
        ])
    }

    /// Create hashmap of column names to data types
    fn hashmap() -> PlHashMap<&'static str, DataType> where Self: Sized + BsxColNames {
        PlHashMap::from_iter([
            (Self::CHR_NAME, Self::ChrType::get_dtype()),
            (Self::POS_NAME, Self::PosType::get_dtype()),
            (Self::STRAND_NAME, Self::StrandType::get_dtype()),
            (Self::CONTEXT_NAME, Self::ContextType::get_dtype()),
            (Self::COUNT_M_NAME, Self::CountType::get_dtype()),
            (Self::COUNT_TOTAL_NAME, Self::CountType::get_dtype()),
            (Self::DENSITY_NAME, Self::DensityType::get_dtype()),
        ])
    }

    /// Create a new batch from a DataFrame without checks
    unsafe fn new_unchecked(data_frame: DataFrame) -> Self where Self: Sized;

    /// Create an empty batch
    fn empty() -> Self where Self: Sized + BsxColNames {
        unsafe { Self::new_unchecked(DataFrame::empty_with_schema(&Self::schema())) }
    }

    /// Split batch at specified index
    fn split_at(
        self,
        index: usize,
    ) -> (Self, Self) where Self: Sized {
        let (a, b) = self.data().split_at(index as i64);
        unsafe {(
            Self::new_unchecked(a),
            Self::new_unchecked(b),
        )}
    }

    /// Returns reference to inner [DataFrame]
    fn data(&self) -> &DataFrame;

    /// Get chromosome value as string
    fn chr_val(&self) -> anyhow::Result<&str>;

    /// Get start position
    fn start_pos(&self) -> Option<u32>
    where
        u32: TryFrom<<Self::PosType as PolarsDataType>::OwnedPhysical>
    {
        self.position()
            .first()
            .map(|v| u32::try_from(v).unwrap_or_else(|_| panic!("Failed to cast to u32")))
    }
    /// Get end position
    fn end_pos(&self) -> Option<u32>
                      where
                          u32: TryFrom<<Self::PosType as PolarsDataType>::OwnedPhysical>
    {
        self.position()
            .last()
            .map(|v| u32::try_from(v).unwrap_or_else(|_| panic!("Failed to cast to u32")))
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
        let res = BsxBatchBuilder::all_checks().check(new_data)?;
        Ok(unsafe { Self::new_unchecked(res) })
    }

    /// Filter batch based on boolean mask
    fn filter_mask(
        &self,
        mask: &BooleanChunked,
    ) -> PolarsResult<Self>
    where Self: Sized
    {
        Ok(unsafe { Self::new_unchecked(self.data().filter(mask)?) })
    }

    /// Returns number of rows
    fn height(&self) -> usize {
        self.data().height()
    }

    /// Convert batch to genomic contig
    fn as_contig(&self) -> anyhow::Result<Contig<&str, NoStrand>>
    where
        u32: TryFrom<<Self::PosType as PolarsDataType>::OwnedPhysical>
    {
        let start = self
            .start_pos()
            .ok_or(anyhow!("no data"))?;
        let end = self
            .end_pos()
            .ok_or(anyhow!("no data"))?;
        let chr = self.chr_val()?;

        Ok(Contig::new(
            chr,
            start as isize,
            (end + 1 - start) as usize,
            NoStrand::Unknown,
        ))
    }
}

/// Trait for BSX type identification
pub trait BsxTypeTag {
    /// Get type name as string
    fn type_name() -> &'static str;
}
