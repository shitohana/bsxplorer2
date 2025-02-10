use crate::data_structs::region::GenomicPosition;
use itertools::Itertools;
use polars::datatypes::PlIndexMap;
use polars::prelude::*;
use std::sync::Arc;

pub mod types;

#[cfg(feature = "python")]
mod python {
    #[macro_export]
    macro_rules! wrap_polars_result {
        ($expression: expr) => {{
            match $expression {
                Ok(v) => Ok(v),
                Err(e) => Err(PyPolarsErr::from(e).into()),
            }
        }};
    }

    #[macro_export]
    macro_rules! wrap_box_result {
        ($error: ty, $expression: expr) => {{
            match $expression {
                Ok(v) => Ok(v),
                Err(e) => Err(<$error>::new_err(e.to_string())),
            }
        }};
    }
    pub use {wrap_box_result, wrap_polars_result};
}
#[cfg(feature = "python")]
pub(crate) use python::{wrap_box_result, wrap_polars_result};

macro_rules! polars_schema {
    ( $($name: expr => $dtype: expr),* ) => {
        {
            let mut fields: Vec<(PlSmallStr, DataType)> = Vec::new();
            $(
                fields.push(($name.into(), $dtype));
            )*
            Schema::from_iter(fields)
        }
    };
}
pub(crate) use polars_schema;

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
        .unwrap_or(0);
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

pub(crate) fn f64_to_u64_scaled(value: f64) -> u64 {
    (value * u64::MAX as f64) as u64
}

pub(crate) fn u64_to_f64_scaled(value: u64) -> f64 {
    value as f64 / u64::MAX as f64
}
// < 0.5% error
pub(crate) const GROUPING_POWER: u8 = 8;
