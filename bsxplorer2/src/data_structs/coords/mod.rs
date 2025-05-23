//! This module defines data structures for representing genomic coordinates.
//!
//! It provides two main structures:
//!
//! - [`GenomicPosition`]: Represents a single point on a genome,
//!   specified by a sequence name (e.g., chromosome) and a position.
//! - [`Contig`]: Represents a genomic region, defined by a sequence name,
//!   a start position, an end position, and an optional strand.
//!
//! These structures are fundamental for handling location-based data within the
//! crate.


mod contig;
mod gpos;

pub use contig::Contig;
pub use gpos::GenomicPosition;

#[cfg(test)]
mod tests;
