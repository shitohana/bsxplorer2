use polars::prelude::*;
use std::marker::PhantomData;

use crate::utils::types::{Context, IPCEncodedEnum, Strand};

use super::{
    builder::BsxBatchBuilder,
    decoded::BsxBatch,
    encoded::EncodedBsxBatch,
    traits::{BsxBatchMethods, BsxColNames, BsxTypeTag},
};

pub struct LazyBsxBatch<T: BsxTypeTag + BsxBatchMethods> {
    data: LazyFrame,
    _phantom: PhantomData<T>,
}

impl<T: BsxTypeTag + BsxBatchMethods> LazyBsxBatch<T> {
    fn filter(
        self,
        predicate: Expr,
    ) -> Self {
        Self {
            data: self.data.filter(predicate),
            _phantom: PhantomData::<T>,
        }
    }
    pub fn filter_pos_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(BsxBatch::POS_NAME).lt(lit(value)))
    }
    pub fn filter_pos_gt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(BsxBatch::POS_NAME).gt(lit(value)))
    }
    pub fn filter_coverage_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(col(BsxBatch::COUNT_TOTAL_NAME).lt(lit(value)))
    }
    pub fn filter_strand(
        self,
        value: Strand,
    ) -> Self {
        let expr = match T::type_name() {
            "decoded" => col(BsxBatch::STRAND_NAME).eq(lit(value.to_string())),
            "encoded" => col(EncodedBsxBatch::STRAND_NAME).eq(value
                .to_bool()
                .map(|v| lit(v))
                .unwrap_or(lit(NULL))),
            other => unimplemented!("Type {}", other),
        };
        self.filter(expr)
    }
    pub fn filter_context(
        self,
        value: Context,
    ) -> Self {
        let expr = match T::type_name() {
            "decoded" => col(BsxBatch::CONTEXT_NAME).eq(lit(value.to_string())),
            "encoded" => col(EncodedBsxBatch::CONTEXT_NAME).eq(value
                .to_bool()
                .map(|v| lit(v))
                .unwrap_or(lit(NULL))),
            other => unimplemented!("Type {}", other),
        };
        self.filter(expr)
    }
}

impl From<BsxBatch> for LazyBsxBatch<BsxBatch> {
    fn from(batch: BsxBatch) -> Self {
        Self {
            data: DataFrame::from(batch).lazy(),
            _phantom: PhantomData::<BsxBatch>,
        }
    }
}

impl From<EncodedBsxBatch> for LazyBsxBatch<EncodedBsxBatch> {
    fn from(batch: EncodedBsxBatch) -> Self {
        Self {
            data: DataFrame::from(batch).lazy(),
            _phantom: PhantomData::<EncodedBsxBatch>,
        }
    }
}

impl<T: BsxTypeTag + BsxBatchMethods> TryFrom<LazyBsxBatch<T>> for BsxBatch {
    type Error = anyhow::Error;

    fn try_from(lazy_batch: LazyBsxBatch<T>) -> Result<Self, Self::Error> {
        let df = lazy_batch.data.collect()?;
        BsxBatchBuilder::no_checks().build_decoded(df)
    }
}

impl<T: BsxTypeTag + BsxBatchMethods> TryFrom<LazyBsxBatch<T>>
    for EncodedBsxBatch
{
    type Error = anyhow::Error;

    fn try_from(lazy_batch: LazyBsxBatch<T>) -> Result<Self, Self::Error> {
        let df = lazy_batch.data.collect()?;
        unsafe { Ok(EncodedBsxBatch::new_unchecked(df)) }
    }
}
