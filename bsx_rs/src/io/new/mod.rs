mod report_schema;
mod bsx_batch;
mod reader;
mod report_batch;

mod utils {
    use std::cmp::Ordering;
    use std::ops::{Sub};
    use polars::prelude::*;
    use crate::region::RegionCoordinates;

    pub(crate) fn encode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
        lazy_frame
            .with_column(
                when(col(strand_col).eq(lit("+")))
                    .then(lit(true))
                    .when(col(strand_col).eq(lit("-")))
                    .then(lit(false))
                    .otherwise(lit(NULL))
                    .cast(DataType::Boolean)
            )
    }

    pub(crate) fn encode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
        lazy_frame
            .with_column(
                when(col(context_col).eq(lit("CG")))
                    .then(lit(true))
                    .when(col(context_col).eq(lit("CHG")))
                    .then(lit(false))
                    .otherwise(lit(NULL))
                    .cast(DataType::Boolean)
            )
    }

    pub(crate) fn decode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
        lazy_frame
            .with_column(
                when(col(strand_col).eq(lit(true)))
                    .then(lit("+"))
                    .when(col(strand_col).eq(lit(false)))
                    .then(lit("-"))
                    .otherwise(".")
                    .cast(DataType::String)
            )
    }

    pub(crate) fn decode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
        lazy_frame
            .with_column(
                when(col(context_col).eq(lit(true)))
                    .then(lit("CG"))
                    .when(col(context_col).eq(lit(false)))
                    .then(lit("CHG"))
                    .otherwise("CHH")
                    .cast(DataType::String)
            )
    }

    pub struct GenomicPosition {
        chr: String,
        position: usize
    }

    impl GenomicPosition {
        pub fn new(chr: String, position: usize) -> GenomicPosition {
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
    
    impl Sub for GenomicPosition {
        type Output = Option<RegionCoordinates>;

        fn sub(self, rhs: Self) -> Self::Output {
            let chr_cmp = self.partial_cmp(&rhs);
            if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Greater {
                Some(RegionCoordinates::new(self.chr, rhs.position as u32, self.position as u32))
            } else if chr_cmp.is_some() && chr_cmp.unwrap() == Ordering::Less {
                Some(RegionCoordinates::new(self.chr, self.position as u32, rhs.position as u32))
            } else {
                None
            }
        }
    }
}