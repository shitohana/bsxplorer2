use crate::region::GenomicPosition;
use itertools::Itertools;
use polars::datatypes::PlIndexMap;
use polars::prelude::*;
use std::sync::Arc;

pub mod types;

pub fn array_to_schema(array: &[(&str, DataType)]) -> Schema {
    Schema::from(PlIndexMap::from_iter(
        array.iter().cloned().map(|(k, v)| (PlSmallStr::from(k), v)),
    ))
}

pub fn get_categorical_dtype(categories: Vec<String>) -> DataType {
    let categories = polars::export::arrow::array::Utf8ViewArray::from_vec(
        categories.iter().map(String::as_str).collect_vec(),
        ArrowDataType::Utf8View,
    );
    let rev_mapping = Arc::new(RevMapping::build_local(categories));
    DataType::Enum(Some(rev_mapping), CategoricalOrdering::Physical)
}

pub(crate) fn schema_from_arrays(names: &[&str], dtypes: &[DataType]) -> Schema {
    Schema::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
}

pub(crate) fn hashmap_from_arrays<'a>(
    names: &[&'a str],
    dtypes: &[DataType],
) -> PlHashMap<&'a str, DataType> {
    PlHashMap::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
}

pub(crate) fn first_position(
    data: &DataFrame,
    chr_col: &str,
    pos_col: &str,
) -> PolarsResult<GenomicPosition<u64>> {
    let chr = data
        .column(chr_col)?
        .as_series()
        .unwrap()
        .first()
        .as_any_value()
        .cast(&DataType::String)
        .str_value()
        .to_string();
    let pos = data
        .column(pos_col)?
        .cast(&DataType::UInt64)?
        .u64()?
        .first()
        .unwrap();
    Ok(GenomicPosition::new(chr, pos))
}

pub(crate) fn last_position(
    data: &DataFrame,
    chr_col: &str,
    pos_col: &str,
) -> PolarsResult<GenomicPosition<u64>> {
    let chr = data
        .column(chr_col)?
        .as_series()
        .unwrap()
        .last()
        .as_any_value()
        .cast(&DataType::String)
        .str_value()
        .to_string();
    let pos = data
        .column(pos_col)?
        .cast(&DataType::UInt64)?
        .u64()?
        .last()
        .unwrap();
    Ok(GenomicPosition::new(chr, pos))
}

pub fn encode_strand(lazy_frame: LazyFrame, strand_col: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(strand_col).eq(lit("+")))
            .then(lit(true))
            .when(col(strand_col).eq(lit("-")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias("strand"),
    )
}

pub fn encode_context(lazy_frame: LazyFrame, context_col: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(context_col).eq(lit("CG")))
            .then(lit(true))
            .when(col(context_col).eq(lit("CHG")))
            .then(lit(false))
            .otherwise(lit(NULL))
            .cast(DataType::Boolean)
            .alias(context_col),
    )
}

pub fn decode_strand(lazy_frame: LazyFrame, strand_col: &str, result_name: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(strand_col).eq(lit(true)))
            .then(lit("+"))
            .when(col(strand_col).eq(lit(false)))
            .then(lit("-"))
            .otherwise(lit("."))
            .cast(DataType::String)
            .alias(result_name),
    )
}

pub fn decode_context(lazy_frame: LazyFrame, context_col: &str, result_name: &str) -> LazyFrame {
    lazy_frame.with_column(
        when(col(context_col).eq(lit(true)))
            .then(lit("CG"))
            .when(col(context_col).eq(lit(false)))
            .then(lit("CHG"))
            .otherwise(lit("CHH"))
            .cast(DataType::String)
            .alias(result_name),
    )
}
