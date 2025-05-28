//! This module contains the core data structures used throughout the
//! `bsxplorer2` crate for representing and manipulating biological data,
//! particularly genomic annotations and methylation information.
//!
//! It provides specialized structures built on top of Rust's standard library
//! and external crates like `polars`, `bio-rs`, and `id_tree`, designed for
//! efficiency and ease of use in bioinformatics contexts.
//!
//! Key components of this module include:
//!
//! - [`annotation`]: Structures for handling genomic features and annotations,
//!   such as genes, transcripts, and regulatory elements. Includes types like
//!   [`GffEntry`], [`GffEntryAttributes`], and a hierarchical storage mechanism
//!   [`HcAnnotStore`].
//! - [`batch`]: Structures for managing batches of genomic methylation data,
//!   typically for a contiguous region on a single chromosome. The central type
//!   is [`BsxBatch`], which wraps a `polars::DataFrame`, along with builders
//!   and lazy evaluation wrappers.
//! - [`coords`]: Basic structures for representing genomic locations:
//!   [`GenomicPosition`] for single points and [`Contig`] for genomic regions.
//! - [`ContextData`], which stores genomic context (like CG, CHG, CHH sites)
//!   derived from sequence information.
//! - Common enumerations used across the application, such as [`Context`] for
//!   methylation context and [`Strand`] for genomic strand.
//! - Structures ([`MethylationStats`], [`MethylationStatFlat`]) for holding and
//!   summarizing methylation statistics for genomic regions or entire samples.
//! - [`typedef`]: Defines type aliases for common data types used for
//!   positions, counts, densities, and sequence names to improve code
//!   readability and maintainability.


pub mod annotation;
pub mod batch;
mod context_data;
pub mod coords;
mod enums;
mod methstats;
pub mod typedef;

#[cfg(test)]
mod tests;

pub use context_data::ContextData;
pub use enums::{
    Context,
    Strand,
};
pub use methstats::{
    MethylationStatFlat,
    MethylationStats,
};
