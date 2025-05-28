This binary provides the command-line interface (`bsxplorer`) for the 
BSXplorer2 suite of tools.
It serves as the user-facing application that leverages the high-performance 
`bsxplorer2` core library
to perform various tasks related to bisulfite sequencing data analysis and 
management.

Key functionalities provided by the `bsxplorer` CLI include:

- **`convert`**: Convert methylation data between different file formats 
(e.g., Bismark, CgMap, BedGraph, Coverage, BSX).
- **`dmr`**: Identify differentially methylated regions between two groups 
of samples using rigorous statistical testing and segmentation.
- **`dimred`**: Perform dimensionality reduction analysis to identify regions 
with high methylation variability across samples using changepoint detection.
- **`validate`**: Ensure the structural integrity and consistency of multiple 
BSX files before joint analysis.
- **`sort`**: Reorder genomic data within a BSX file based on chromosome order.

## Installation

You can install `bsxplorer` CLI tool with:

```shell
cargo install --locked bsxplorer-ci
```

After installation, `bsxplorer` binary should become available.

Explore possible commands with

```shell
bsxplorer help
```