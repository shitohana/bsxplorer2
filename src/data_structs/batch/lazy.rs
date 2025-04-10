use crate::data_structs::batch::traits::colnames;
use crate::utils::types::{Context, IPCEncodedEnum, Strand};
use polars::prelude::*;
use std::marker::PhantomData;

use super::{
    builder::BsxBatchBuilder,
    decoded::BsxBatch,
    encoded::EncodedBsxBatch,
    traits::{BsxBatchMethods, BsxTypeTag},
};

/// A lazy representation of a BSX batch for efficient query operations
pub struct LazyBsxBatch<T: BsxTypeTag + BsxBatchMethods> {
    data: LazyFrame,
    _phantom: PhantomData<T>,
}

impl<T: BsxTypeTag + BsxBatchMethods> LazyBsxBatch<T> {
    /// Applies a filter expression to the batch
    fn filter(
        self,
        predicate: Expr,
    ) -> Self {
        Self {
            data: self.data.filter(predicate),
            _phantom: PhantomData::<T>,
        }
    }

    /// Creates a LazyBsxBatch from an existing LazyFrame
    fn from_lazy(lazy: LazyFrame) -> Self {
        Self {
            data: lazy,
            _phantom: PhantomData::<T>,
        }
    }

    /// Filters positions less than the specified value
    pub fn filter_pos_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(colnames::POS_NAME).lt(lit(value)))
    }

    /// Filters positions greater than the specified value
    pub fn filter_pos_gt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(colnames::POS_NAME).gt(lit(value)))
    }

    /// Filters entries with coverage less than the specified value
    pub fn filter_coverage_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(colnames::COUNT_TOTAL_NAME).lt(lit(value)))
    }

    /// Filters entries by strand value
    pub fn filter_strand(
        self,
        value: Strand,
    ) -> Self {
        let expr = match T::type_name() {
            "decoded" => col(colnames::STRAND_NAME).eq(lit(value.to_string())),
            "encoded" => col(colnames::STRAND_NAME).eq(value
                .to_bool()
                .map(|v| lit(v))
                .unwrap_or(lit(NULL))),
            other => unimplemented!("Type {}", other),
        };
        self.filter(expr)
    }

    /// Filters entries by context value
    pub fn filter_context(
        self,
        value: Context,
    ) -> Self {
        let expr = match T::type_name() {
            "decoded" => col(colnames::CONTEXT_NAME).eq(lit(value.to_string())),
            "encoded" => col(colnames::CONTEXT_NAME).eq(value
                .to_bool()
                .map(|v| lit(v))
                .unwrap_or(lit(NULL))),
            other => unimplemented!("Type {}", other),
        };
        self.filter(expr)
    }

    /// Marks entries with coverage below threshold as zero and adjusts related values
    pub fn mark_low_coverage(self, threshold: u32) -> Self {
        let res = self.data.with_column(
            when(col("count_total").lt(lit(threshold)))
                .then(lit(0))
                .otherwise(col("count_total"))
                .alias("count_total"),
        )
            .with_columns([
                when(col("count_total").eq(lit(0)))
                    .then(lit(0))
                    .otherwise(col("count_m"))
                    .alias("count_m"),
                when(col("count_total").eq(lit(0)))
                    .then(lit(f64::NAN))
                    .otherwise(col("density"))
                    .alias("density"),
            ]);
        Self::from_lazy(res)
    }
}

impl From<BsxBatch> for LazyBsxBatch<BsxBatch> {
    /// Converts a decoded BSX batch to its lazy representation
    fn from(batch: BsxBatch) -> Self {
        Self {
            data: DataFrame::from(batch).lazy(),
            _phantom: PhantomData::<BsxBatch>,
        }
    }
}

impl From<EncodedBsxBatch> for LazyBsxBatch<EncodedBsxBatch> {
    /// Converts an encoded BSX batch to its lazy representation
    fn from(batch: EncodedBsxBatch) -> Self {
        Self {
            data: DataFrame::from(batch).lazy(),
            _phantom: PhantomData::<EncodedBsxBatch>,
        }
    }
}

impl<T: BsxTypeTag + BsxBatchMethods> TryFrom<LazyBsxBatch<T>> for BsxBatch {
    type Error = anyhow::Error;

    /// Attempts to convert a lazy BSX batch to a decoded BSX batch
    fn try_from(lazy_batch: LazyBsxBatch<T>) -> Result<Self, Self::Error> {
        let df = lazy_batch.data.collect()?;
        BsxBatchBuilder::no_checks().build_decoded(df)
    }
}

impl<T: BsxTypeTag + BsxBatchMethods> TryFrom<LazyBsxBatch<T>>
    for EncodedBsxBatch
{
    type Error = anyhow::Error;

    /// Attempts to convert a lazy BSX batch to an encoded BSX batch
    fn try_from(lazy_batch: LazyBsxBatch<T>) -> Result<Self, Self::Error> {
        let df = lazy_batch.data.collect()?;
        unsafe { Ok(EncodedBsxBatch::new_unchecked(df)) }
    }
}
