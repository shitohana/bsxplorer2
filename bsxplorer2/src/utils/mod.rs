use std::sync::Arc;

use itertools::Itertools;
use log::warn;
use polars::datatypes::PlIndexMap;
use polars::prelude::*;

pub mod types;
mod stats;
pub use stats::*;

/// Converts an array of name-datatype pairs to a Polars Schema.
pub fn array_to_schema(array: &[(&str, DataType)]) -> Schema {
    Schema::from(PlIndexMap::from_iter(
        array
            .iter()
            .cloned()
            .map(|(k, v)| (PlSmallStr::from(k), v)),
    ))
}

/// Creates a categorical data_structs type from a list of categories.
pub fn get_categorical_dtype(categories: Vec<String>) -> DataType {
    let categories = polars::export::arrow::array::Utf8ViewArray::from_vec(
        categories
            .iter()
            .map(String::as_str)
            .collect_vec(),
        ArrowDataType::Utf8View,
    );
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Enum(Some(rev_mapping), CategoricalOrdering::Physical)
}

/// Creates a schema from separate arrays of names and data_structs types.
pub(crate) fn schema_from_arrays(
    names: &[&str],
    dtypes: &[DataType],
) -> Schema {
    if names.len() != dtypes.len() {
        warn!(
            "Mismatch between names and dtypes array lengths: {} vs {}",
            names.len(),
            dtypes.len()
        );
    }
    Schema::from_iter(
        names
            .iter()
            .cloned()
            .map_into()
            .zip(dtypes.iter().cloned()),
    )
}

/// Creates a hashmap from separate arrays of names and data_structs types.
pub(crate) fn hashmap_from_arrays<'a>(
    names: &[&'a str],
    dtypes: &[DataType],
) -> PlHashMap<&'a str, DataType> {
    if names.len() != dtypes.len() {
        warn!(
            "Mismatch between names and dtypes array lengths: {} vs {}",
            names.len(),
            dtypes.len()
        );
    }
    PlHashMap::from_iter(
        names
            .iter()
            .cloned()
            .map_into()
            .zip(dtypes.iter().cloned()),
    )
}
