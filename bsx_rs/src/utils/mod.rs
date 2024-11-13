use polars::datatypes::PlIndexMap;
use polars::prelude::{DataType, PlSmallStr, Schema};

pub mod types;

pub fn array_to_schema(array: &[(&str, DataType)]) -> Schema {
    Schema::from(PlIndexMap::from_iter(
        array
            .iter()
            .cloned()
            .map(|(k, v)| (PlSmallStr::from(k), v)),
    ))
}
