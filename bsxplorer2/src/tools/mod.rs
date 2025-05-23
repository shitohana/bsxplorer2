//! This module provides various tools and algorithms for analyzing genomic and
//! methylation data within the bsxplorer2 crate.
//!
//! It includes functionality for identifying differentially methylated regions (DMRs),
//! performing segmentation on data streams, and gathering genome-wide statistics.
//!
//! Key submodules:
//!
//! - [`dimred`]: Contains algorithms for dimensionality reduction and segmentation,
//!   such as the PELT algorithm and Total Variation based methods, useful for identifying
//!   segments or changepoints in sequential genomic data.
//! - [`dmr`]: Provides tools specifically designed for the analysis of Differential
//!   Methylated Regions, including configuration, data structures for representing
//!   potential DMRs, and an iterator-based approach for processing data streams
//!   and identifying DMRs.
pub mod dimred;
pub mod dmr;
mod stat_gather;
