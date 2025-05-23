//! This module provides functionality for reading and writing files in the
//! custom `.bsx` format, which is based on the Apache Arrow IPC file format.
//!
//! The `.bsx` format is designed to store genomic data in a structured,
//! batch-oriented manner, allowing for efficient sequential and random access
//! to data by genomic region.
//!
//! Key components of this module include:
//!
//! - [`BsxFileReader`]: A reader optimized for parallel access to data
//!   batches within a `.bsx` file, leveraging memory mapping and multithreading.
//! - [`BsxFileWriter`]: A writer for creating `.bsx` files from data batches,
//!   supporting different data sources for schema information (e.g., FASTA index)
//!   and optional compression.
//! - [`BatchIndex`]: An in-memory index mapping genomic regions (contigs) to
//!   the indices of the data batches that contain data for those regions. This
//!   index is associated with the `.bsx` file and is used by
//!   the [`RegionReader`].
//! - [`RegionReader`]: A higher-level reader that uses the [`BatchIndex`] and
//!   the [`BsxFileReader`] to efficiently retrieve data for specific genomic
//!   regions, managing caching of batches to minimize I/O.
//!
//! The module utilizes `polars-arrow` for the core Arrow IPC serialization/deserialization,
//! `memmap2` for efficient file mapping, and `rayon` for parallel processing.
//!
//! ## Basic reading
//! ```
//! # use std::path::PathBuf;
//!
//! use std::fs::File;
//! use bsxplorer2::io::bsx::BsxFileReader;
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! # let filepath = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx");
//! 
//! let handle = File::open(filepath)?;
//! let mut reader = BsxFileReader::try_new(handle)?;
//!
//! for batch_res in reader.iter() {
//!     let batch = batch_res?;
//!     assert!(!batch.is_empty())
//! }
//! # Ok(())
//! # }
//! ```
//!
//! ## Reading regions
//! ```
//! # use std::path::PathBuf;
//! use std::fs::File;
//! use bsxplorer2::data_structs::annotation::HcAnnotStore;
//! use bsxplorer2::io::bsx::{BsxFileReader, RegionReader};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! # let filepath = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/report.bsx");
//! # let bed_filepath = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/genes.bed");
//!
//! let bsx_reader = BsxFileReader::try_new(File::open(filepath)?)?;
//!
//! let mut region_reader = RegionReader::from_reader(bsx_reader)?;
//!
//! let annotation = HcAnnotStore::from_bed(File::open(bed_filepath)?)?;
//! let contigs = region_reader.index().sort(
//!     annotation.iter().map(|entry| entry.contig().clone())
//! ).collect::<Vec<_>>();
//! 
//! for contig_res in region_reader.iter_contigs(&contigs) {
//!     let contig_data = contig_res?;
//! }
//! # Ok(())
//! # }
//! ```
mod index;
mod read;
mod region;
mod write;

pub use index::BatchIndex;
// pub use read::{BsxFileIterator, BsxFileReader};
pub use read::{BsxFileIterator, BsxFileReader};
pub use region::RegionReader;
pub use write::BsxFileWriter;
