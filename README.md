# BSXplorer2: Accelerating DNA Methylation Analysis

[![Documentation](https://docs.rs/bsxplorer2/badge.svg)](https://docs.rs/bsxplorer2)
[![Version](https://img.shields.io/crates/v/bsxplorer2)](https://crates.io/crates/bsxplorer2)
[![codecov](https://codecov.io/github/shitohana/bsxplorer2/graph/badge.svg)](https://codecov.io/github/shitohana/bsxplorer2)

![License](https://img.shields.io/github/license/shitohana/bsxplorer2)
![Downloads](https://img.shields.io/crates/dr/bsxplorer2)

A cutting-edge, high-performance toolkit built in Rust for bisulfite sequencing 
data analysis and DNA methylation research.


<!-- mtoc-start -->

* [Overview](#overview)
* [Features](#features)
* [Components](#components)
  * [Core Rust Library](#core-rust-library)
  * [Python Wrapper (bsx2)](#python-wrapper-bsx2)
  * [Console Application (bsxplorer)](#console-application-bsxplorer)
* [Installation](#installation)
* [BSX Format (Arrow IPC File Format)](#bsx-format-arrow-ipc-file-format)
  * [Performance Benefits](#performance-benefits)
  * [Compression Capabilities](#compression-capabilities)
  * [Data Organization](#data-organization)
  * [Integration Advantages](#integration-advantages)
* [Roadmap](#roadmap)
* [License](#license)
* [Acknowledgements](#acknowledgements)

<!-- mtoc-end -->

## Overview

BSXplorer2 is designed from the ground up for speed and efficiency, enabling 
researchers and developers to process and analyze large-scale bisulfite sequencing 
datasets with unprecedented performance. By leveraging Rust's powerful features 
and integrating with modern data processing libraries like Polars and Arrow, 
BSXplorer2 provides a robust and scalable solution for identifying differentially 
methylated regions (DMRs), calculating methylation statistics, and handling various 
report formats.

Whether you prefer command-line tools for quick analyses or a programmatic interface 
for complex pipelines, BSXplorer2 offers flexible access through its console binary 
and Python bindings.

For detailed documentation and usage examples please refer to the 
[documentation](https://docs.rs/bsxplorer2).

## Features

✨ **Core Capabilities for High-Impact Research**

-   **Blazing Fast Data Handling:** Process massive datasets efficiently using 
memory-optimized data structures and native parallelization.
-   **Comprehensive Report Support:** Seamlessly work with Bismark, CGmap, BedGraph, 
and Coverage formats, plus our high-performance BSX format.
-   **Context-Aware Analysis:** Drill down into CG, CHG, and CHH methylation patterns.
-   **Advanced DMR Detection:** Pinpoint differentially methylated regions using 
cutting-edge segmentation and statistical methods.
-   **Robust Statistics:** Calculate detailed methylation statistics, coverage distributions, 
and apply sophisticated statistical tests.

⚡ **Engineered for Performance**

-   **Rust Native Speed:** Built on Rust for maximum performance and reliability.
-   **Polars & Arrow Integration:** Leverage column-oriented processing for speed and 
memory efficiency.
-   **Parallel Execution:** Utilize multi-core processors effectively with Rayon.

🤝 **User-Friendly & Accessible**

-   **Intuitive Console App:** Perform common tasks easily with the `bsxplorer` 
command-line tool.
-   **Flexible Python API:** Build custom analysis workflows using the `bsx2` Python library.
-   **Detailed Documentation:** Get started quickly with clear guides and examples.

## Components

BSXplorer2 is composed of three main parts:

### Core Rust Library
The heart of BSXplorer2, containing all the core data structures, algorithms, 
and file format implementations. Designed for high performance and low-level 
control.
Explore the Rust source code: [@src](@file:bsxplorer2_dev/bsxplorer2)

### Python Wrapper (bsx2)
Provides idiomatic Python bindings to the core Rust library using PyO3. Enables 
seamless integration with the Python data science ecosystem (Polars, NumPy, SciPy, 
Matplotlib, Plotly, Pydantic). Ideal for building complex analysis pipelines and 
interactive data exploration in Jupyter notebooks.
Find the Python package here: [@python](@file:bsxplorer2_dev/python)

### Console Application (bsxplorer)
A standalone command-line tool built on the Rust library. Offers convenient commands 
for file format conversion, DMR calling, validation, and more. Perfect for scripting 
and integrating into existing bioinformatics workflows without writing Rust or Python 
code.
Check out the console source and commands: [@console](@file:bsxplorer2_dev/console)

## Installation

### For the Console Application (`bsxplorer`)
Install the console binary directly using Cargo:

```bash
cargo install --locked bsxplorer-ci
```

Ensure your Cargo bin directory is in your system's PATH.

### For the Python Library (`bsx2`)

🚧 WIP! Python library is currently being actively developed! 🚧

## Usage

Dive into analyzing your methylation data using the `bsxplorer` console 
application or the `bsx2` Python library.

*   **Console Usage:** Get detailed help and command examples in the 
[Console Application README](console/README.md).
*   **Python Usage:** Explore the `bsx2` package documentation (coming soon!) and 
examples within the [@python](@file:bsxplorer2_dev/python) directory to use the Python API.

## BSX Format (Arrow IPC File Format)

BSXplorer2 introduces the BSX file format, leveraging the power of Apache Arrow's Interprocess 
Communication (IPC) format. This isn't just another file type; it's a foundation for 
highly efficient methylation data processing:

### Performance Benefits

-   **Memory Efficiency:** Column-oriented storage significantly reduces memory footprint 
compared to row-based formats.
-   **Zero-Copy Reading:** Data can be accessed in memory without expensive copying, 
boosting speed.
-   **Parallel Processing:** Designed for concurrent access, perfectly complementing 
multi-threaded operations.
-   **Vectorized Operations:** Enables leveraging modern CPU instructions (SIMD) for 
faster calculations.

### Compression Capabilities

-   **Flexible Compression:** Supports LZ4 (optimized for speed) and ZSTD (optimized 
for compression ratio).
-   **Column-Specific:** Compression is applied per column, adapting to different data types.
-   **Efficient Decompression:** Only necessary columns are decompressed, minimizing 
overhead.

### Data Organization

-   **Efficient Categorical Encoding:** Methylation contexts (CG, CHG, CHH) and strands 
are stored as efficient categorical types, not verbose strings.
-   **Batched Storage:** Data is chunked into logical batches for efficient processing 
in memory.
-   **Type-Aware:** Data types (integers, floats, booleans) are stored in optimized 
binary representations.

### Integration Advantages

-   **Cross-Platform:** Works consistently across various operating systems.
-   **Language Interoperability:** Accessible from any language with robust Arrow 
bindings (Python, R, Java, etc.).
-   **Schema Enforcement:** Strict schema ensures data integrity and prevents 
format ambiguities.
-   **Rich Metadata:** Supports embedding custom metadata for better data tracking 
and provenance.

The BSX format is purpose-built for methylation data, providing the optimal storage 
solution for BSXplorer2's high-performance analytical tasks.

## Roadmap

BSXplorer2 is under active development. Future plans include:

-   [x] High-performance file format support (BSX, Bismark, CGmap, BedGraph, Coverage) including 
reading, writing, conversion, validation, and sorting.
-   [x] Efficient indexing and region-based querying for BSX files.
-   [x] Core DMR identification algorithm implementation.
-   [x] Basic methylation statistics calculation.
-   [ ] Enhanced visualization tools within the Python library.
-   [ ] Tighter integration and utilities for genomic annotation data (genes, regulatory elements).
-   [ ] Exploration of a web-based interactive analysis interface.
-   [ ] Expansion of statistical methods for sophisticated differential methylation analysis.
-   [ ] Implement Metagene profile generation.

Contributions and feature requests are welcome!

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.

## Acknowledgements

-   The foundational work for the total variation segmentation algorithm is inspired by Laurent Condat.
-   Statistical implementations draw upon established techniques from bioinformatics literature.
-   We gratefully acknowledge the contributions of community-developed libraries, including `bio-types`, `polars`, `pyo3`, and `rayon`, which are integral to BSXplorer2.

---

Created by [shitohana](https://github.com/shitohana) - Empowering your DNA methylation research with speed and precision.
