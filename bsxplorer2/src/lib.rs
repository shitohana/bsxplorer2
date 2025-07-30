//! # bsxplorer2
//!
//! `bsxplorer2` is a Rust library and command-line tool for efficient handling,
//! processing, and analysis of genomic methylation data. It is designed to
//! work with various common methylation report formats (like Bismark, CGmap,
//! BedGraph, Coverage) and introduces a custom, optimized `.bsx` file format
//! based on Apache Arrow IPC for fast random access and stream processing.
//!
//! The crate provides core data structures to represent genomic regions and
//! methylation data, efficient I/O mechanisms, and analytical tools,
//! particularly for identifying differentially methylated regions (DMRs) and
//! performing data segmentation.
//!
//! If you do not want to use bsxplorer as crate, check out
//! [**bsxplorer CLI tool**](https://crates.io/crates/bsxplorer-ci)!
//!
//! ## Key Features
//!
//! * **Efficient Data Structures**: Optimized structures for genomic
//!   coordinates ([`Contig`], [`GenomicPosition`]) and methylation data batches
//!   ([`BsxBatch`]) leveraging Polars DataFrames for in-memory performance.
//! * **Custom `.bsx` Format**: A high-performance, indexed binary format based
//!   on Apache Arrow IPC for storing methylation data, enabling fast sequential
//!   and random access to genomic regions.
//! * **Flexible I/O**: Readers and writers for `.bsx` files ([`BsxFileReader`],
//!   [`BsxFileWriter`], [`RegionReader`]) and common methylation report formats
//!   ([`ReportReader`], [`ReportWriter`]), with support for compression
//!   (feature-gated) and FASTA/FAI for context alignment.
//! * **Parallel Processing**: Leverages Rayon for parallel operations and
//!   Polars' native multi-threading for data processing tasks, improving
//!   performance on multi-core systems.
//! * **DMR Identification**: Tools for detecting Differentially Methylated
//!   Regions based on statistical methods and segmentation algorithms.
//! * **Data Segmentation**: Implementation of segmentation algorithms (e.g.,
//!   PELT) for identifying changepoints in genomic data streams.
//! * **Integration**: Uses established libraries like `bio-rs` for biological
//!   formats (GFF, BED, FASTA) and `polars` for data manipulation.
//!
//! Number of threads to be used can be configured with setting
//! `BSX_NUM_THREADS` environment variable.
//!
//! ## Structure
//!
//! The crate is organized into several modules:
//!
//! * [`data_structs`]: Defines the fundamental data types used throughout the
//!   crate, including genomic coordinates ([`Contig`], [`GenomicPosition`]),
//!   methylation data batches ([`BsxBatch`]), genomic annotations
//!   ([`GffEntry`], [`HcAnnotStore`]), and methylation statistics
//!   ([`MethylationStats`]).
//! * [`io`]: Handles file input and output, including the custom `.bsx` format
//!   and various methylation report formats.
//! * [`tools`]: Contains higher-level analytical tools, such as those for
//!   Differential Methylation Region (DMR) analysis (`dmr`) and data
//!   segmentation (`dimred`). (feature-gated by default)
//! * [`utils`]: Provides common utility functions, including statistical
//!   methods and helpers for Polars data types.
//!
//! ## Installation
//!
//! ```bash
//! # Add as a dependency to your Cargo.toml
//! cargo add bsxplorer2
//! ```
//!
//! To enable the `tools` module and compression features:
//!
//! ```bash
//! # Add with features
//! cargo add bsxplorer2 --features tools,compression
//! ```
//!
//! ## Usage
//!
//! ### Reading a `.bsx` file
//!
//! ```no_run
//! use std::fs::File;
//! use bsxplorer2::prelude::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let file = File::open("path/to/your/file.bsx")?;
//!     let mut reader = BsxFileReader::try_new(file)?;
//!
//!     for batch_result in reader.iter() {
//!         let batch = batch_result?;
//!         // Process the batch (BsxBatch)
//!         println!(
//!             "Read batch from {} with {} sites",
//!             batch.seqname().unwrap_or("unknown"),
//!             batch.len()
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Reading a specific genomic region from `.bsx`
//!
//! ```no_run
//! use std::fs::File;
//! use bsxplorer2::prelude::*;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let file = File::open("path/to/your/indexed_file.bsx")?;
//!     let bsx_reader = BsxFileReader::try_new(file)?;
//!     let mut region_reader = RegionReader::from_reader(bsx_reader)?;
//!
//!     // Define the region of interest
//!     let region = Contig::new(BsxSmallStr::from("chr1"), 1000, 5000, Strand::None);
//!
//!     if let Some(batch_result) = region_reader.query(region.clone(), None)? {
//!          // Process the batch for the requested region
//!          println!(
//!             "Read data for region {} with {} sites",
//!             region,
//!             batch_result.len()
//!         );
//!     } else {
//!         println!("No data found for region {}", region);
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ### Reading a methylation report file (e.g., Bismark)
//!
//! ```no_run
//! use std::fs::File;
//! use bsxplorer2::prelude::*;
//! use std::path::PathBuf;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let file_path = PathBuf::from("path/to/your/report.txt");
//!     let fasta_path = PathBuf::from("path/to/your/genome.fa");
//!
//!     // Use the builder to configure the reader, including FASTA for context
//!     let mut reader = ReportReaderBuilder::default()
//!         .with_report_type(ReportType::Bismark)
//!         .with_fasta_path(Some(fasta_path))
//!         .build(file_path)?;
//!
//!     for batch_result in reader {
//!         let batch = batch_result?;
//!         // Process the BsxBatch converted from the report
//!         println!(
//!             "Processed batch from {} with {} sites",
//!             batch.seqname().unwrap_or("unknown"),
//!             batch.len()
//!         );
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Writing a `.bsx` file
//!
//! ```no_run
//! use std::fs::File;
//! use bsxplorer2::prelude::*;
//! use std::path::PathBuf;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let file = File::create("output.bsx")?;
//!     let fasta_index = PathBuf::from("path/to/your/genome.fa.fai");
//!
//!     // Create writer, using FASTA index to get chromosome names
//!     let mut writer = BsxFileWriter::try_from_sink_and_fai(file, fasta_index, None, None)?;
//!
//!     // Create some dummy data (replace with your actual data processing)
//!     let batch1 = BsxBatch::try_from_columns(
//!         "chr1", None, // None uses default categorical dtype
//!         vec![100, 101, 150],
//!         vec![true, false, true],
//!         vec![Some(true), Some(false), Some(true)],
//!         vec![5, 10, 15],
//!         vec![10, 20, 30]
//!     )?;
//!
//!      let batch2 = BsxBatch::try_from_columns(
//!         "chr2", None, // None uses default categorical dtype
//!         vec![500, 550],
//!         vec![true, false],
//!         vec![Some(true), Some(false)],
//!         vec![20, 25],
//!         vec![40, 50]
//!     )?;
//!
//!     // Write batches
//!     writer.write_batch(batch1)?;
//!     writer.write_batch(batch2)?;
//!
//!     // Close the writer to finalize the file
//!     writer.close()?;
//!
//!     println!("Successfully wrote output.bsx");
//!
//!     Ok(())
//! }
//! ```
#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

#[ctor::ctor]
fn init() {
    if let Ok(n) = std::env::var("BSX_NUM_THREADS") {
        std::env::set_var("POLARS_MAX_THREADS", n)
    }
}

pub mod data_structs;
pub mod io;
pub mod prelude;
pub mod utils;

#[cfg(feature = "tools")]
#[cfg_attr(coverage_nightly, coverage(off))]
pub mod tools;

#[allow(unused_imports)]
use prelude::*;
