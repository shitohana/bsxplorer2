mod bsx_batch;
mod data_batch;
mod data_schema;
mod reader;
mod report_schema;

mod utils {
    use crate::region::RegionCoordinates;
    use log::warn;
    use polars::prelude::*;
    use std::cmp::Ordering;
    use std::fmt::{Display, Formatter};
    use std::ops::{Range, RangeInclusive, Shr, Sub};

    pub(crate) fn encode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(strand_col).eq(lit("+")))
                .then(lit(true))
                .when(col(strand_col).eq(lit("-")))
                .then(lit(false))
                .otherwise(lit(NULL))
                .cast(DataType::Boolean),
        )
    }

    pub(crate) fn encode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(context_col).eq(lit("CG")))
                .then(lit(true))
                .when(col(context_col).eq(lit("CHG")))
                .then(lit(false))
                .otherwise(lit(NULL))
                .cast(DataType::Boolean),
        )
    }

    pub(crate) fn decode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(strand_col).eq(lit(true)))
                .then(lit("+"))
                .when(col(strand_col).eq(lit(false)))
                .then(lit("-"))
                .otherwise(".")
                .cast(DataType::String),
        )
    }

    pub(crate) fn decode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
        lazy_frame.with_column(
            when(col(context_col).eq(lit(true)))
                .then(lit("CG"))
                .when(col(context_col).eq(lit(false)))
                .then(lit("CHG"))
                .otherwise("CHH")
                .cast(DataType::String),
        )
    }

    #[derive(Debug, Clone)]
    pub struct GenomicPosition {
        chr: String,
        position: u64,
    }

    impl Display for GenomicPosition {
        fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}:{}", self.chr, self.position)
        }
    }

    impl GenomicPosition {
        pub fn new(chr: String, position: u64) -> GenomicPosition {
            GenomicPosition { chr, position }
        }
    }

    impl PartialEq for GenomicPosition {
        fn eq(&self, other: &Self) -> bool {
            self.chr == other.chr && self.position == other.position
        }
    }

    impl PartialOrd for GenomicPosition {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            if self.chr == other.chr {
                Some(self.position.cmp(&other.position))
            } else {
                None
            }
        }
    }

    impl Shr for GenomicPosition {
        type Output = Option<RegionCoordinates>;

        fn shr(self, rhs: Self) -> Self::Output {
            let chr_cmp = self.partial_cmp(&rhs);
            if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Greater {
                Some(RegionCoordinates::new(
                    self.chr,
                    rhs.position as u32,
                    self.position as u32,
                ))
            } else if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Less {
                warn!(
                    "Bad range operation for {} and {} regions. Lhs must be less than rhs.",
                    self, rhs
                );
                Some(RegionCoordinates::new(
                    self.chr,
                    self.position as u32,
                    rhs.position as u32,
                ))
            } else {
                None
            }
        }
    }
}
