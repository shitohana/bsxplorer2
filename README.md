# BSXplorer

A high-performance, Rust-based library for bisulfite sequencing data analysis and DNA methylation research.

![License](https://img.shields.io/badge/license-Prosperity-blue)
![Version](https://img.shields.io/badge/version-0.1.0-green)

## Overview

BSXplorer is a comprehensive toolkit for analyzing bisulfite sequencing data, focusing on efficient processing, statistical analysis, and identification of differentially methylated regions (DMRs). Built with performance in mind, it leverages Rust's memory safety and concurrency features to handle large-scale methylation datasets effectively.

The library provides a complete pipeline for methylation analysis, from raw data processing to advanced statistical testing and visualization, supporting various input formats commonly used in bisulfite sequencing research.

## Features

- **Efficient Data Structures**
    - Optimized storage and processing of methylation data using Polars DataFrames
    - Memory-efficient encoding of methylation contexts and strand information
    - Support for batch processing of large datasets

- **Versatile I/O Support**
    - Custom BSX file format for efficient methylation data storage
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

## Roadmap

BSXplorer is under active development. Future plans include:

- Enhanced visualization capabilities for methylation patterns
- Integration with genome browser formats (BigWig, BigBed)
- Support for single-cell bisulfite sequencing analysis
- Integration with genomic annotation data (genes, regulatory elements)
- Machine learning models for methylation pattern prediction
- Web interface for interactive analysis
- Additional statistical methods for differential methylation analysis
- Support for hydroxymethylation (5hmC) data

## Contributing

Contributions to BSXplorer are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the Prosperity Public License 3.0.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

- The total variation segmentation algorithm is based on work by Laurent Condat
- Statistical methods draw from established techniques in bioinformatics literature
- Parts of the codebase leverage community-developed libraries including bio-types, polars, and rayon

---

Created by [shitohana](https://github.com/shitohana) - Empowering methylation analysis through efficient computational methods.