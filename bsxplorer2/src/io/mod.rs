//! This module provides comprehensive functionality for handling file input
//! and output within the bsxplorer2 crate. It focuses on reading and
//! writing genomic data, particularly in the custom `.bsx` format based on
//! Apache Arrow IPC, and various common methylation report formats.
//!
//! Key submodules include:
//!
//! - [`bsx`]: Manages the reading and writing of `.bsx` files, offering
//!   optimized sequential and regional access through structures like
//!   [`BsxFileReader`], [`BsxFileWriter`], [`BatchIndex`], and
//!   [`RegionReader`].
//! - [`compression`]: Provides utilities for handling compressed data streams,
//!   used in conjunction with file readers and writers (feature-gated).
//! - [`report`]: Handles the parsing and writing of different methylation
//!   report formats (e.g., Bismark, BedGraph, Coverage), including features for
//!   aligning report data with reference genome contexts using FASTA/FAI.

pub mod bsx;
#[cfg(feature = "compression")]
pub mod compression;
pub mod report;
