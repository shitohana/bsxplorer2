use crate::region::GenomicPosition;
use polars::prelude::*;

#[derive(Debug)]
pub struct BsxBatch(DataFrame);

impl BsxBatch {
    pub fn new(data_frame: DataFrame) -> Self {
        BsxBatch(data_frame)
    }
    
    pub fn data(&self) -> &DataFrame {
        &self.0
    }
    
    pub(crate) const fn col_names() -> &'static [&'static str] {
        &[
            "chr",
            "position",
            "strand",
            "context",
            "count_m",
            "count_total",
            "density",
        ]
    }

    pub(crate) const fn col_types() -> &'static [DataType] {
        &[
            DataType::String,
            DataType::UInt64,
            DataType::String,
            DataType::String,
            DataType::UInt32,
            DataType::UInt32,
            DataType::Float64,
        ]
    }

    pub(crate) const fn chr_col() -> &'static str {
        "chr"
    }

    pub(crate) const fn pos_col() -> &'static str {
        "position"
    }

    pub fn schema() -> Schema {
        use crate::io::report::report_batch_utils::schema_from_arrays;
        schema_from_arrays(Self::col_names(), Self::col_types())
    }

    pub fn hashmap() -> PlHashMap<&'static str, DataType> {
        use crate::io::report::report_batch_utils::hashmap_from_arrays;
        hashmap_from_arrays(Self::col_names(), Self::col_types())
    }

    pub fn first_position(&self) -> PolarsResult<GenomicPosition> {
        use crate::io::report::report_batch_utils::first_position;
        first_position(&self.0, Self::chr_col(), Self::pos_col())
    }
    pub fn last_position(&self) -> PolarsResult<GenomicPosition> {
        use crate::io::report::report_batch_utils::last_position;
        last_position(&self.0, Self::chr_col(), Self::pos_col())
    }
}