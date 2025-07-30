//! Contains the data structures and logic for representing and manipulating
//! batches of genomic methylation data.
//!
//! A `BsxBatch` is the primary data structure, holding methylation data for a
//! contiguous genomic region on a single chromosome (contig) in a Polars
//! `DataFrame`.
//!
//! The data within a `BsxBatch` is structured according to the `BsxColumns`
//! schema, which defines the standard column names and data types required
//! (e.g., chromosome, position, methylation counts, density, etc.).
//!
//! The module provides the following key components:
//!
//! - [`BsxBatch`]: The core struct wrapping a Polars `DataFrame` and providing
//!   domain-specific methods for accessing columns, performing genomic
//!   operations (like slicing and reporting), calculating methylation
//!   statistics, and transforming data (e.g., binomial test).
//! - [`BsxBatchBuilder`]: A builder pattern for creating validated `BsxBatch`
//!   instances, especially useful when loading data from various report formats
//!   (Bismark, CGmap, etc.) or concatenating batches. It includes configurable
//!   data validation checks (e.g., for sorted positions, duplicates, single
//!   chromosome).
//! - [`LazyBsxBatch`]: A wrapper around a Polars `LazyFrame` for `BsxBatch`,
//!   allowing for lazy evaluation of operations like filtering, which can be
//!   more efficient for large datasets.
//! - [`BsxColumns`]: An enum defining the canonical column names and data types
//!   used within a `BsxBatch` DataFrame.
//! - `utils`: Contains utility functions, such as `merge_replicates` for
//!   combining data from multiple samples covering the same genomic regions.

mod base;
mod builder;
mod lazy;
mod schema;
mod utils;

pub use base::*;
pub use builder::*;
pub use lazy::*;
pub use schema::*;
pub use utils::*;

#[cfg(test)]
mod tests;
