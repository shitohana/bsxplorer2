use anyhow::anyhow;
use bio_types::annot::contig::Contig;
use bio_types::strand::NoStrand;
use polars::datatypes::BooleanChunked;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::prelude::Expr;

use super::builder::BsxBatchBuilder;

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
    unsafe fn new_unchecked(data_frame: DataFrame) -> Self
    where
        Self: Sized;
    /// Filters [BsxBatch] by `context` and `strand`
    ///
    /// If both `context` and `strand` are [None], returns [BsxBatch]
    fn filter(
        self,
        expr: Expr,
    ) -> anyhow::Result<Self>
    where
        Self: Sized;

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
