//! This module provides functionality for reading and writing various
//! methylation report file formats, such as Bismark, CgMap, BedGraph, and
//! Coverage reports.
//!
//! It offers tools to parse these formats efficiently, handle different data
//! structures, and potentially align report data with a reference genome
//! (FASTA) for context information.
//!
//! Key components:
//!
//! - [`ReportType`]: An enum defining the supported report formats, their
//!   schemas, and parsing options.
//! - [`ReportReader`]: A reader for processing report files, capable of
//!   handling large datasets through chunking, buffering, and optional context
//!   alignment from a FASTA file. It leverages multithreading for performance.
//! - [`ReportReaderBuilder`]: A builder pattern for configuring the
//!   [`ReportReader`] with various options like chunk size, FASTA/FAI paths,
//!   batch size, number of threads, low memory mode, and compression.
//! - [`ReportWriter`]: A writer for converting internal data structures (like
//!   [`BsxBatch`]) into the specified report formats and writing them to a
//!   sink, using a batched CSV writer.
//!
//! This module integrates with `polars` for efficient data handling and
//! processing, `rayon` for parallelization, and optionally `bio` and `noodles`
//! for FASTA indexing and sequence reading.

mod read;
mod schema;
mod write;

pub use read::{
    ReportReader,
    ReportReaderBuilder,
};
pub use schema::ReportType;
pub use write::ReportWriter;
