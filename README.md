# BSXplorer

A high-performance, Rust-based library for bisulfite sequencing data analysis and DNA methylation research.

![License](https://img.shields.io/badge/license-Prosperity-blue)
![Version](https://img.shields.io/badge/version-0.1.0-green)


<!-- mtoc-start -->

* [Overview](#overview)
* [Features](#features)
* [Installation](#installation)
* [Usage](#usage)
  * [Basic Example: Reading and Processing Methylation Data](#basic-example-reading-and-processing-methylation-data)
* [Console Application](#console-application)
* [BSX Format (IPC File Format)](#bsx-format-ipc-file-format)
  * [Performance Benefits](#performance-benefits)
  * [Compression Capabilities](#compression-capabilities)
  * [Data Organization](#data-organization)
  * [Integration Advantages](#integration-advantages)
* [DMR Identification Benchmark](#dmr-identification-benchmark)
* [Roadmap](#roadmap)
* [License](#license)
* [Acknowledgements](#acknowledgements)

<!-- mtoc-end -->

## Overview

BSXplorer is a comprehensive toolkit for analyzing bisulfite sequencing data, focusing on efficient processing, 
statistical analysis, and identification of differentially methylated regions (DMRs). Built with performance in mind, 
it leverages Rust's memory safety and concurrency features to handle large-scale methylation datasets effectively.

## Features

- **Efficient Data Structures**
    - Optimized storage and processing of methylation data using Polars DataFrames
    - Memory-efficient encoding of methylation contexts and strand information
    - Support for batch processing of large datasets

- **Versatile I/O Support**
    - Custom BSX file (Apache IPC File) format for efficient methylation data storage
    - Support for popular methylation report formats:
        - Bismark methylation extractor output
        - CG methylation map (CgMap)
        - BedGraph methylation density format
        - Coverage reports with methylated/unmethylated counts
    - FASTA sequence integration for genomic context analysis

- **Methylation Analysis Tools**
    - Context-specific methylation analysis (CG, CHG, CHH)
    - Strand-specific methylation patterns
    - Comprehensive methylation statistics calculation
    - Coverage distribution analysis

- **Differentially Methylated Region (DMR) Detection**
    - Advanced total variation segmentation algorithm
    - Mann-Whitney U statistical testing for DMR validation
    - Configurable DMR parameters (minimum coverage, p-value thresholds, etc.)
    - Region filtering and merging capabilities

- **Statistical Methods**
    - Beta-binomial distribution modeling for methylation data
    - Method of Moments (MoM) estimation for distribution parameters
    - Kolmogorov-Smirnov and Mann-Whitney U non-parametric tests
    - Dimensionality reduction techniques for methylation patterns

- **Performance Optimizations**
    - Parallel processing with Rayon for CPU-intensive operations
    - Memory-efficient data representations
    - Batch processing for large datasets
    - Optimized algorithms for DMR detection

## Installation

Add BSXplorer to your Rust project by including it in your `Cargo.toml`:

```toml
[dependencies]
bsxplorer = "0.1.0"
```

Documentation is available at [docs.rs](https://docs.rs/bsxplorer2)

## Usage

### Basic Example: Reading and Processing Methylation Data

```rust
use bsxplorer::io::bsx::read::BsxFileReader;
use bsxplorer::data_structs::bsx_batch::BsxBatchMethods;
use bsxplorer::utils::types::Context;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Open a BSX file
    let mut reader = BsxFileReader::new(std::fs::File::open("sample.bsx")?);
    
    // Process the first batch
    if let Some(batch_result) = reader.next() {
        let batch = batch_result?;
        
        // Filter for CG context only
        let cg_batch = batch.filter(Some(Context::CG), None);
        
        // Calculate methylation statistics
        let stats = cg_batch.get_methylation_stats()?;
        println!("Mean methylation: {}", stats.mean_methylation());
        
        // Access positions and methylation values
        let positions = cg_batch.get_position_vals()?;
        let methylation = cg_batch.get_density_vals()?;
        
        println!("Analyzed {} CpG sites", positions.len());
    }
    
    Ok(())
}
```

## Console Application
BSXplorer includes a powerful command-line interface for direct interaction with methylation data. The console 
application provides convenient access to the library's core functionality without requiring Rust programming knowledge.

[Detailed command descriptions](console/README.md).

## BSX Format (IPC File Format)
BSXplorer utilizes Arrow's Interprocess Communication (IPC) file format as the foundation for its custom BSX format, delivering significant advantages for methylation data processing:

### Performance Benefits

- Memory Efficiency: Column-oriented storage dramatically reduces memory footprint compared to traditional formats
- Zero-Copy Reading: Data can be accessed without redundant copying between memory regions
- Parallel Processing: Format supports concurrent access patterns for multi-threaded operations
- Vectorized Operations: Enables CPU-optimized SIMD instructions for faster data processing

### Compression Capabilities

- Multiple Compression Options: Supports both LZ4 (faster) and ZSTD (better compression ratio)
- Column-Level Compression: Each column is compressed independently, optimizing for data characteristics
- Minimal Decompression Overhead: Selective decompression of only required columns

### Data Organization

- Efficient Categorical Encoding: Methylation contexts and strands are stored as enumerated values, not strings
- Batched Storage: Data is organized in batches for efficient in-memory processing
- Type-Aware Storage: Numeric types are stored in their binary representation, not as text

### Integration Advantages

- Cross-Platform Compatibility: Works consistently across operating systems
- Language Interoperability: Can be read by any language with Arrow bindings (Python, R, etc.)
- Schema Enforcement: Strong typing prevents data corruption and format inconsistencies
- Metadata Support: Embedded metadata for tracking experimental conditions and processing steps

The BSX format combines these advantages into a specialized format optimized for methylation data, ensuring the best 
possible performance for complex analytical tasks.

## DMR Identification Benchmark

We've evaluated our DMR identification model F1-score, using benchmarking dataset from 
_C. Kreutz et al., ‘A blind and independent benchmark study for detecting differentially 
methylated regions in plants’, Bioinformatics, vol. 36, no. 11, pp. 3314–3321, Jun. 2020, 
doi: 10.1093/bioinformatics/btaa191._

![](https://private-user-images.githubusercontent.com/43905117/427238513-85507124-6347-4b51-930c-a153466c1646.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NDMwMjU4NDYsIm5iZiI6MTc0MzAyNTU0NiwicGF0aCI6Ii80MzkwNTExNy80MjcyMzg1MTMtODU1MDcxMjQtNjM0Ny00YjUxLTkzMGMtYTE1MzQ2NmMxNjQ2LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTAzMjYlMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwMzI2VDIxNDU0NlomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTU5ZWYzMGU5MmU0MTY5NDMyODk5NGJkNWZmODU3YTQ2NWY1MjVhMDE0NzUzZDQ4ZTIzNmFhY2Y3YmZmZjgyZWImWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.cut1umPWh_WNMdObsMgwXs5HIdRrHY-apxp_hp1a_VM)

## Roadmap

BSXplorer is under active development. Future plans include:

- Enhanced visualization capabilities for methylation patterns
- Integration with genome browser formats (BigWig, BigBed)
- Support for single-cell bisulfite sequencing analysis
- Integration with genomic annotation data (genes, regulatory elements)
- Machine learning models for methylation pattern prediction
- Web interface for interactive analysis
- Additional statistical methods for differential methylation analysis

## License

This project is licensed under the Prosperity Public License 3.0.0 - see the [LICENSE](LICENSE.md) file for details.

## Acknowledgements

- The total variation segmentation algorithm is based on work by Laurent Condat
- Statistical methods draw from established techniques in bioinformatics literature
- Parts of the codebase leverage community-developed libraries including bio-types, polars, and rayon

---

Created by [shitohana](https://github.com/shitohana) - Empowering methylation analysis through efficient computational 
methods.
