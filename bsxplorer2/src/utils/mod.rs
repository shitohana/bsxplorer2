//! This module contains various utility functions and helper macros used
//! throughout the bsxplorer2 crate.
//!
//! It provides common functionalities that are not specific to a particular
//! module but are required by multiple components, promoting code reuse
//! and maintainability.
//!
//! Key functionalities include:
//!
//! - Statistical functions, such as Pearson correlation and Mann-Whitney U
//!   test.
//! - Handling of Polars data types, including creation of categorical types
//!   from chromosome lists and schema/hashmap generation from arrays.
//! - Macros for common struct operations (e.g., getter functions, builder-style
//!   `with_*` methods).
//! - Reading chromosome names from FASTA (`.fa`) and FASTA index (`.fai`)
//!   files.
//! - Utility functions for converting between float and integer representations
//!   of methylation density values.

use std::io::{
    BufReader,
    Read,
};
use std::sync::Arc;

use itertools::Itertools;
use log::warn;
use noodles::fasta::io::Indexer;
use once_cell::sync::Lazy;
use polars::prelude::*;

mod stats;
use rayon::{
    ThreadPool,
    ThreadPoolBuilder,
};
pub use stats::*;

pub static THREAD_POOL: Lazy<ThreadPool> = Lazy::new(|| {
    let num_threads: Option<usize> = std::env::var("BSX_NUM_THREADS")
        .ok()
        .and_then(|str| str.parse::<usize>().ok());
    ThreadPoolBuilder::new()
        .num_threads(num_threads.unwrap_or(0))
        .build()
        .expect("Failed to create thread pool")
});

pub fn n_threads() -> usize {
    THREAD_POOL.current_num_threads()
}

/// Creates a categorical data_structs type from a list of categories.
pub fn get_categorical_dtype(categories: Vec<String>) -> DataType {
    let categories = polars::export::arrow::array::Utf8ViewArray::from_vec(
        categories.iter().map(String::as_str).collect_vec(),
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
    Schema::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
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
    PlHashMap::from_iter(names.iter().cloned().map_into().zip(dtypes.iter().cloned()))
}

#[macro_export]
macro_rules! plsmallstr {
    ($string: expr) => {
        PlSmallStr::from($string)
    };
    () => {
        PlSmallStr::from("")
    };
    () => {
        PlSmallStr::from("")
    };
}

#[macro_export]
macro_rules! getter_fn {
    ($field_name: ident, $field_type: ty) => {
        #[cfg_attr(coverage_nightly, coverage(off))]
        pub fn $field_name(&self) -> &$field_type {
            &self.$field_name
        }
    };
}
pub use getter_fn;

#[macro_export]
macro_rules! with_field_fn {
    ($field_name: ident, $field_type: ty) => {
        paste::paste! {
            #[cfg_attr(coverage_nightly, coverage(off))]
            pub fn [<with_$field_name>](mut self, value: $field_type) -> Self {
            self.$field_name = value;
            self
            }
        }
    };
}

pub fn read_chrs_from_fai<R: Read>(reader: R) -> anyhow::Result<Vec<String>> {
    let records: Vec<noodles::fasta::fai::Record> =
        noodles::fasta::fai::io::Reader::new(BufReader::new(reader))
            .read_index()?
            .into();
    Ok(records
        .into_iter()
        .map(|r| String::from_utf8_lossy(r.name()).to_string())
        .collect())
}

pub fn read_chrs_from_fa<R: Read>(reader: R) -> anyhow::Result<Vec<String>> {
    let mut indexer = Indexer::new(BufReader::new(reader));
    let mut records = Vec::new();

    while let Some(record) = indexer.index_record()? {
        records.push(record);
    }

    Ok(records
        .into_iter()
        .map(|r| String::from_utf8_lossy(r.name()).to_string())
        .collect())
}
