# Rust Region Aggregation Backend

## Purpose

BSX2 now has an additive Rust-backed path for aggregating methylated and
unmethylated counts by arbitrary genomic regions. The target use cases are DMRs,
genes, promoters, BED regions, and external caller intervals.

The backend accelerates the heavy count aggregation step. It does not change the
DMR caller, Regional Evidence Model, p-value calculation, confidence intervals,
evidence scoring, QC interpretation, or visualization logic.

## Architecture

The public flow is:

```text
Python RegionSignal API
  -> Rust RegionCountAggregator binding
  -> BSX2 RegionReader / Contig / BsxBatch
  -> .bsx indexed methylation records
```

The Rust module is `bsxplorer2::tools::region_aggregation`. It reuses existing
BSX2 structures:

- `RegionReader` for indexed `.bsx` access.
- `Contig` for genomic intervals.
- `BsxBatch` columns for `position`, `strand`, `context`, `count_m`, and
  `count_total`.

The Python wrapper is `bsx2.analysis.region_signal_rust`. The high-level
fallback-aware interface remains `bsx2.analysis.region_signal`.

## What Moved To Rust

Rust now performs:

- region-by-region `.bsx` queries;
- context and strand filtering;
- per-region sums of `mC`, `uC`, and `total`;
- `n_cytosines`;
- `mean_methylation`;
- coverage QC labels.

The Rust side returns pandas-friendly rows with this schema:

```text
region_id, seqname, chrom, start, end, strand, context, sample_id,
mC, uC, total, n_cytosines, mean_methylation, coverage_qc
```

## What Remains In Python

Python still owns:

- RegionSignal table normalization;
- seqname alias harmonization;
- pandas fallback aggregation;
- Regional Evidence statistical testing;
- beta-binomial validation;
- evidence scoring and class assignment;
- QC/report summaries;
- plotting and thesis output tables.

## API Examples

```python
from bsx2.analysis import RegionSignalConfig, aggregate_region_signal

result = aggregate_region_signal(
    regions_df,
    "/path/to/sample.bsx",
    RegionSignalConfig(backend="auto", context="CG"),
    sample_id="sample_1",
)
```

```python
result = aggregate_region_signal(
    regions_df,
    "/path/to/sample.bsx",
    RegionSignalConfig(backend="rust", strand_policy="region_strand"),
)
```

```python
result = aggregate_region_signal(
    regions_df,
    counts_df,
    RegionSignalConfig(backend="pandas"),
)
```

Backend discovery:

```python
from bsx2.analysis import available_region_signal_backends

available_region_signal_backends()
```

## CLI Example

```bash
python scripts/aggregate_region_signal.py \
  --regions regions.tsv \
  --counts sample.bsx \
  --backend auto \
  --context CG \
  --out region_signal.tsv \
  --qc-out region_signal_qc.tsv \
  --summary-out region_signal_summary.md
```

For tabular count files, use:

```bash
python scripts/aggregate_region_signal.py \
  --regions regions.tsv \
  --counts counts.tsv \
  --backend pandas \
  --out region_signal.tsv \
  --qc-out region_signal_qc.tsv
```

## Limitations

- This is a first production-hardening step for count aggregation.
- Statistical testing is not implemented in Rust.
- Performance depends on existing `.bsx` indexing and `RegionReader` cache
  behavior.
- If the Rust extension is unavailable, `backend="auto"` falls back to pandas
  when a pandas-compatible count table is supplied.
- Full production-scale benchmarking on the B. napus six-sample run remains a
  separate future benchmark.

## Thesis Wording

The Rust/chunked backend accelerates region-level methylation count aggregation
for `.bsx` methylation files while preserving the same downstream Regional
Evidence statistical model. The methodical contribution is engineering
separation of high-throughput count aggregation from evidence scoring, not a new
DMR calling model.
