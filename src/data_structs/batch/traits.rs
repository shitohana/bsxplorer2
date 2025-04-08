use super::builder::BsxBatchBuilder;
use crate::data_structs::batch::decoded::BsxBatch;
use anyhow::anyhow;
use bio_types::annot::contig::Contig;
use bio_types::strand::NoStrand;
use polars::datatypes::BooleanChunked;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::prelude::*;

pub trait BsxColNames {
    const CHR_NAME: &'static str = "chr";
    const POS_NAME: &'static str = "position";
    const STRAND_NAME: &'static str = "strand";
    const CONTEXT_NAME: &'static str = "context";
    const COUNT_M_NAME: &'static str = "count_m";
    const COUNT_TOTAL_NAME: &'static str = "count_total";
    const DENSITY_NAME: &'static str = "density";

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
pub trait BsxBatchMethods {
    type ChrType: PolarsDataType;
    type PosType: PolarsDataType;
    type StrandType: PolarsDataType;
    type ContextType: PolarsDataType;
    type CountType: PolarsDataType;
    type DensityType: PolarsDataType;

    fn chr(&self) -> &ChunkedArray<Self::ChrType>;
    fn position(&self) -> &ChunkedArray<Self::PosType>;
    fn strand(&self) -> &ChunkedArray<Self::StrandType>;
    fn context(&self) -> &ChunkedArray<Self::ContextType>;
    fn count_m(&self) -> &ChunkedArray<Self::CountType>;
    fn count_total(&self) -> &ChunkedArray<Self::CountType>;
    fn density(&self) -> &ChunkedArray<Self::DensityType>;

    fn chr_type() -> DataType where Self: Sized {
        Self::ChrType::get_dtype()
    }
    fn pos_type() -> DataType where Self: Sized {
        Self::PosType::get_dtype()
    }
    fn strand_type() -> DataType where Self: Sized {
        Self::StrandType::get_dtype()
    }
    fn context_type() -> DataType where Self: Sized {
        Self::ContextType::get_dtype()
    }
    fn count_type() -> DataType where Self: Sized {
        Self::CountType::get_dtype()
    }
    fn density_type() -> DataType where Self: Sized {
        Self::DensityType::get_dtype()
    }

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

    unsafe fn new_unchecked(data_frame: DataFrame) -> Self where Self: Sized;

    fn empty() -> Self where Self: Sized + BsxColNames {
        unsafe { Self::new_unchecked(DataFrame::empty_with_schema(&Self::schema())) }
    }

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

    /// Returns mutable reference to inner [DataFrame]
    fn data_mut(&mut self) -> &mut DataFrame;

    fn chr_val(&self) -> anyhow::Result<&str>;

    fn start_pos(&self) -> Option<u32>;
    fn end_pos(&self) -> Option<u32>;

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

    fn as_contig(&self) -> anyhow::Result<Contig<&str, NoStrand>> {
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

pub trait BsxTypeTag {
    fn type_name() -> &'static str;
}
