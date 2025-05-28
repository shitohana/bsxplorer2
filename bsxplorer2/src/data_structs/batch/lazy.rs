use itertools::Itertools;
use polars::prelude::*;

use super::{
    BsxBatch,
    BsxColumns as BsxCol,
};
use crate::data_structs::enums::{
    Context,
    Strand,
};

#[derive(Clone)]
pub struct LazyBsxBatch {
    data: LazyFrame,
}

impl From<BsxBatch> for LazyBsxBatch {
    fn from(batch: BsxBatch) -> Self {
        Self {
            data: batch.into_inner().lazy(),
        }
    }
}

impl LazyBsxBatch {
    pub fn collect(self) -> PolarsResult<BsxBatch> {
        let data = self
            .data
            .collect()?
            .select(BsxCol::schema().iter_names_cloned().collect_vec())?;
        Ok(unsafe { BsxBatch::new_unchecked(data) })
    }

    fn filter(
        self,
        predicate: Expr,
    ) -> Self {
        Self {
            data: self.data.filter(predicate),
        }
    }

    /// Filters positions less than the specified value.
    pub fn filter_pos_lt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(BsxCol::Position.col().lt(lit(value)))
    }

    /// Filters positions greater than the specified value.
    pub fn filter_pos_gt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(BsxCol::Position.col().gt(lit(value)))
    }

    /// Filters entries with coverage greater than the specified value.
    pub fn filter_coverage_gt<N: Literal>(
        self,
        value: N,
    ) -> Self {
        self.filter(BsxCol::CountTotal.col().gt(lit(value)))
    }

    /// Filters entries by strand value.
    pub fn filter_strand(
        self,
        value: Strand,
    ) -> Self {
        self.filter(match value {
            Strand::Forward | Strand::Reverse => {
                BsxCol::Strand
                    .col()
                    .eq(lit(Option::<bool>::from(value).unwrap()))
            },
            _ => BsxCol::Context.col().is_null(),
        })
    }

    /// Filters entries by context value.
    pub fn filter_context(
        self,
        value: Context,
    ) -> Self {
        self.filter(match value {
            Context::CG | Context::CHG => {
                BsxCol::Context
                    .col()
                    .eq(lit(Option::<bool>::from(value).unwrap()))
            },
            _ => BsxCol::Context.col().is_null(),
        })
    }
}
