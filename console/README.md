# bsxplorer2-ci


<!-- mtoc-start -->

* [Installation](#installation)
* [bsxplorer convert](#bsxplorer-convert)
* [bsxplorer dmr](#bsxplorer-dmr)
* [bsxplorer stats](#bsxplorer-stats)

<!-- mtoc-end -->

## Installation

```commandline
cargo install --locked bsxplorer-ci
```

After installation, bsxplorer executable will be available in your PATH as `bsxplorer`

For more detailed information and benchmarks, please refer to [bsxplorer2](../README.md)

## bsxplorer convert

```text
bsxplorer convert --help
BSXplorer report type conversion tool

Usage: bsxplorer convert [OPTIONS] --output <OUTPUT> --from <FROM_TYPE> --into <INTO_TYPE> <INPUT>

Arguments:
  <INPUT>  Path of the input file.

Options:
  -o, --output <OUTPUT>                Path for the generated output file.
  -f, --from <FROM_TYPE>               [default: bismark] [possible values: bsx, bismark, cg-map, bed-graph, coverage]
  -i, --into <INTO_TYPE>               [default: bsx] [possible values: bsx, bismark, cg-map, bed-graph, coverage]
  -C, --compression <IPC_COMPRESSION>  [default: zstd] [possible values: lz4, zstd, none]
      --batch-size <BATCH_SIZE>        Size of raw batches. [default: 2097152]
      --progress                       Display progress bar (Disable if you need clean pipeline logs).
      --threads <THREADS>              Number of threads to use. [default: 1]
      --verbose                        Verbose output.
  -h, --help                           Print help

REPORT ARGS:
      --low-memory
          Use less RAM, but elongate computation.
  -c, --chunk <CHUNK_SIZE>
          Number of rows in the output batches (Important when converting to bsx format). [default: 10000]
      --fa <FASTA_PATH>
          Path to the reference sequence file. Obligatory when converting BedGraph or Coverage.
      --fai <FAI_PATH>
          Path to the fasta index. Obligatory when converting BedGraph or Coverage.
      --batch-per-read <BATCH_PER_READ>
          Number of batches to read simultaneously. Affects RAM usage. [default: 8]
```

Examples:

Convert from Bismark methylation report to BSX file format

```commandline
bsxplorer convert --from bismark --into bsx -o report.bsx --fai example.fa.fai -c 20000 bismark_report.CX_report.txt
```

Convert from Bismark to BedGraph

```commandline
bsxplorer convert --from bismark --into bed-graph -o report.bedGraph bismark_report.CX_report.txt
```

Convert from BSX file format to Bismark

```commandline
bsxplorer convert --from bsx --into Bismark -o report.CX_report.txt bsx_report.bsx
```

## bsxplorer dmr

```text
BSXplorer DMR identification algorithm.

Usage: bsxplorer dmr [OPTIONS] --group-a <GROUP_A> --group-b <GROUP_B> --output <OUTPUT>

Options:
      --progress           Display progress bar (Disable if you need clean pipeline logs).
      --threads <THREADS>  Number of threads to use. [default: 1]
      --verbose            Verbose output.
  -A, --group-a <GROUP_A>  Paths to BSX files of the first sample group.
  -B, --group-b <GROUP_B>  Paths to BSX files of the second sample group.
  -o, --output <OUTPUT>    Prefix for the generated output files.
  -f, --force              Automatically confirm selected paths.
  -h, --help               Print help

FILTER ARGS:
  -c, --context <CONTEXT>
          Select cytosine methylation context. Only cytosines in this context will be used for DMR calling. CG/CHG/CHH. [default: cg] [possible values: cg, chg, chh]
  -n, --n-missing <N_MISSING>
          Set missing values threshold. Cytosines with no data_structs in more than specified number of samples will be discarded. [default: 0]
  -v, --min-coverage <MIN_COVERAGE>
          Set coverage threshold. Cytosines with coverage below this value in any of the samples will be discarded. [default: 5]
  -m, --min-cytosines <MIN_CYTOSINES>
          Set minimum number of cytosines threshold. DMRs with cytosine count below this value will be discarded. [default: 10]
  -d, --diff-threshold <DIFF_THRESHOLD>
          Set minimum difference threshold. DMRs with an absolute difference in methylation proportion between the two groups smaller than this value will be discarded. [default: 0.05]
  -p, --padj <PADJ>
          Adjusted P-value threshold for DMR identification using 2D-Kolmogorov-Smirnov test. Segments with a p-value smaller than specified will be reported as DMRs. [default: 0.05]
      --pmethod <PMETHOD>
          [default: bh] [possible values: bonf, bh, by, none]

SEGMENTATION ARGS:
  -D, --max-dist <MAX_DIST>    Maximum distance between adjacent cytosines in a segment.  Cytosines further apart than this distance will be in separate segments. [default: 100]
  -L, --initial-l <INITIAL_L>  Initial regularization parameter for the Condat algorithm.  Larger values result in stronger smoothing. [default: 2]
  -l, --l-min <L_MIN>          Minimum value for the regularization parameter.  The regularization parameter is decreased during segmentation until it is smaller than this value. [default: 0.001]
      --coef <L_COEF>          Coefficient by which `initial_l` is divided in each iteration of the segmentation algorithm. Smaller values perform more segmentation iterations. [default: 1.5]
      --tolerance <TOLERANCE>  Tolerance for merging adjacent segments after the Total Variation denoising step (Condat's algorithm).  Smaller values result in more segments being merged. Should be very small to avoid over-segmentation after denoising. [default: 0.000001]
      --merge-p <MERGE_P>      Mann-Whitney U-test P-value threshold for merging adjacent segments during recursive segmentation. Smaller p-values result in more iterations and fewer falsely merged segments. [default: 0.01]
```

## bsxplorer stats

```text
bsxplorer stats --help
Compute methylation statistics.

Usage: bsxplorer stats [OPTIONS] --output <OUTPUT> <INPUT>

Arguments:
  <INPUT>  Path of the input file.

Options:
  -o, --output <OUTPUT>              Path for the generated output file.
  -m, --mode <MODE>                  Stats mode. [default: genomewide] [possible values: genomewide, regions]
  -f, --format <FORMAT>              Annotation format. [default: gff] [possible values: gff, bed]
  -a, --annot-path <ANNOT_PATH>      Path for the generated output file.
      --feature-type <FEATURE_TYPE>  Feature type to filter. [default: gene]
      --threads <THREADS>            Number of threads to use. [default: 1]
      --progress                     Display progress bar (Disable if you need clean pipeline logs).
      --threads <THREADS>            Number of threads to use. [default: 1]
      --verbose                      Verbose output.
  -h, --help                         Print help

If mode is set to `regions`, annotation file must be provided.

Output is JSON for genome-wide mode and TSV for regions mode.
```

Examples:

For genome-wide mode:
```commandline
bsxplorer stats --output stats.json report.ipc
```

For regions mode:
```commandline
bsxplorer stats --output stats.tsv --threads 12 --mode regions --format gff -a genomic.gff report.ipc
```
