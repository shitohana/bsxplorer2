use crate::data_structs::region::GenomicPosition;
use crate::utils::types::{Context, Strand};
use log::warn;
use polars::error::PolarsError;
use polars::frame::DataFrame;
use polars::prelude::{Expr, SortMultipleOptions};

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

    /// Extends [BsxBatch] by other
    ///
    /// Prints warning if first position of other is not less than last of self,
    /// but still sorts data_structs by position.
    ///
    /// Returns error if
    /// 1. Chromosome names do not match
    /// 2. Chromosome columns non-unique
    /// 3. Resulting data_structs still contains duplicates
    fn extend(
        &mut self,
        other: &Self,
    ) -> anyhow::Result<()>
    where
        Self: Sized 
    {
        if self.chr_val()? != other.chr_val()? {
            return Err(PolarsError::ComputeError(
                "Chromosomes in batches differ".into(),
            )
            .into());
        }
        if self.end_pos() >= other.start_pos() {
            warn!(
                "First position in other batch ({:?}) must be less than last \
                 position in the current ({:?})",
                self.end_pos(), other.start_pos()
            );
            self.data_mut().extend(other.data())?;
            self.data_mut().sort_in_place(
                ["position"],
                SortMultipleOptions::default()
                    .with_order_descending(false)
                    .with_multithreaded(true),
            )?;

            if self
                .data()
                .column("position")?
                .n_unique()?
                != self.data().height()
            {
                println!(
                    "{} != {}",
                    self.data()
                        .column("position")?
                        .n_unique()?,
                    self.data().height()
                );
                println!("{}", self.data());
                return Err(PolarsError::ComputeError(
                    "Position values are duplicated".into(),
                )
                .into());
            }
        }

        Ok(self.data_mut().extend(other.data())?)
    }

    /// Returns number of rows
    fn height(&self) -> usize { self.data().height() }
}