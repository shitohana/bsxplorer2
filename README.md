# BSXplorer2: Accelerating DNA Methylation Analysis

[![Documentation](https://docs.rs/bsxplorer2/badge.svg)](https://docs.rs/bsxplorer2)
[![Version](https://img.shields.io/crates/v/bsxplorer2)](https://crates.io/crates/bsxplorer2)
[![codecov](https://codecov.io/github/shitohana/bsxplorer2/graph/badge.svg)](https://codecov.io/github/shitohana/bsxplorer2)

![License](https://img.shields.io/github/license/shitohana/bsxplorer2)
![Downloads](https://img.shields.io/crates/dr/bsxplorer2)

A high-performance toolkit for bisulfite sequencing data analysis and DNA methylation research.

<!-- mtoc-start -->

* [Overview](#overview)
* [Features](#features)
* [Components](#components)
  * [Core Rust Library](#core-rust-library)
  * [Python Wrapper (bsx2)](#python-wrapper-bsx2)
  * [Console Application (bsxplorer)](#console-application-bsxplorer)
* [Installation](#installation)
* [Usage](#usage)
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

BSXplorer2 is designed for fast, scalable methylation analysis. The project combines:

- a Rust core for storage, indexing, querying, and analysis primitives
- a Python package for visualization, clustering, and workflow integration
- a console application for command-line usage

For crate-level reference, see the [docs.rs documentation](https://docs.rs/bsxplorer2).

## Features

**Core Capabilities**

- High-performance support for BSX, Bismark, CGmap, BedGraph, and Coverage-like workflows
- Context-aware methylation analysis for CG, CHG, and CHH
- Efficient region-based querying through indexed BSX files
- Statistical and exploratory analysis building blocks for downstream pipelines

**Performance-Oriented Design**

- Rust-native implementation for core storage and compute paths
- Polars and Arrow integration for efficient columnar workflows
- Region-level and block-level access patterns designed for large datasets

**Python Visualization Layer**

- Metagene aggregation from arbitrary contigs
- Annotation-driven metagene aggregation from `RegionReader + HcAnnotStore`
- HoloViews renderers for line, heatmap, box, and violin metagene plots
- Clustering utilities with PCA, dendrogram, and cluster metagene views
- Chromosome methylation map support

## Components

BSXplorer2 is composed of three main parts:

### Core Rust Library

The heart of BSXplorer2, containing core data structures, algorithms, and file format
implementations. Designed for high performance and low-level control.

Source: [bsxplorer2](bsxplorer2)

### Python Wrapper (bsx2)

Idiomatic Python bindings and higher-level workflow helpers built on top of the Rust
core. This layer now includes the main visualization and clustering API.

Source and docs: [python](python), [python/README.md](python/README.md)

### Console Application (bsxplorer)

A standalone command-line tool built on the Rust library for conversion, validation,
and scripted analysis workflows.

Source and commands: [console](console), [console/README.md](console/README.md)

## Installation

### Console Application (`bsxplorer`)

Install the console binary with Cargo:

```bash
cargo install --locked bsxplorer-ci
```

### Python Library (`bsx2`)

Install the Python package with Poetry:

```bash
cd python
poetry install
```

## Usage

- Console usage: see [console/README.md](console/README.md)
- Python usage: see [python/README.md](python/README.md) for metagene, clustering,
  and chromosome-map workflows

## BSX Format (Arrow IPC File Format)

BSXplorer2 uses the BSX file format, built on Arrow IPC, as a storage layer for efficient
methylation data access and analysis.

### Performance Benefits

- Column-oriented storage for lower memory overhead
- Efficient indexed region access
- Good fit for vectorized and batched operations
- Cross-language interoperability through Arrow-compatible tooling

### Compression Capabilities

- Support for LZ4 and ZSTD-backed workflows
- Compression aligned with columnar storage patterns
- Efficient selective decompression for relevant data slices

### Data Organization

- Context and strand data stored in efficient typed representations
- Batched layout suited for large analytical workloads
- Explicit schema and metadata support

### Integration Advantages

- Cross-platform operation
- Interoperability with Python and other Arrow-aware ecosystems
- Clear schema enforcement for data integrity

## Roadmap

BSXplorer2 is under active development. Current state:

- [x] High-performance file format support and conversion workflows
- [x] Efficient indexing and region-based querying for BSX files
- [x] Core DMR identification building blocks
- [x] Basic methylation statistics calculation
- [x] Python visualization tools for metagene, clustering, and chromosome maps
- [x] Metagene profile generation
- [ ] Deeper utilities for richer genomic annotation workflows
- [ ] Expanded statistical methods for more advanced differential methylation analysis
- [ ] Broader interactive and web-facing analysis surfaces

## License

This project is licensed under the MIT License. See [LICENSE.md](LICENSE.md).

## Acknowledgements

- The total variation segmentation work draws on ideas from Laurent Condat
- Statistical implementation choices are informed by established bioinformatics literature
- The project relies on key ecosystem libraries including `bio-types`, `polars`, `pyo3`, and `rayon`
