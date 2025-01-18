use std::sync::Arc;
use itertools::Itertools;
use polars::datatypes::PlIndexMap;
use polars::prelude::*;

pub mod types;

pub fn array_to_schema(array: &[(&str, DataType)]) -> Schema {
    Schema::from(PlIndexMap::from_iter(
        array.iter().cloned().map(|(k, v)| (PlSmallStr::from(k), v)),
    ))
}

pub fn get_categorical_dtype(categories: Vec<String>) -> DataType {
    let categories = polars::export::arrow::array::Utf8ViewArray::from_vec(
        categories.iter().map(String::as_str).collect_vec(),
        ArrowDataType::Utf8View
    );
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Enum(Some(rev_mapping), CategoricalOrdering::Physical)
}