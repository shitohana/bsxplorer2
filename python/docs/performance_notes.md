# Performance Notes

This note summarizes the current scaling profile of the Python clustering and
visualization layer.

## What Was Optimized

The clustering metagene path now avoids per-bin `ProfileBin` allocation inside the
hot loop that materializes gene profiles.

Current fast path in `bsx2.clustering.gene_profile`:

- merge genes into chromosome query blocks
- query one BSX block at a time
- build cumulative sums once per block
- build metagene interval matrices per query block
- aggregate all bins for a block of genes with batched `searchsorted` lookups
- finalize retained matrices by reusing NaN masks and filling missing retained values in place
- build clustering tables from numpy-backed columns instead of eagerly converting most arrays to Python lists

Current render/composition fast path in `bsx2.viz`:

- heatmap row matrices are preallocated instead of materialized as row lists and stacked later
- line profiles use direct per-window accumulation for mean/min/max instead of concatenating all points first
- chromosome-wide profiles aggregate windows as arrays and reduce once per chromosome instead of merging Python dict entries for every window

This reduces Python overhead in the per-gene aggregation step, but it does not change
the dominant cost center for many real datasets: selective BSX reads and block
materialization.

## Practical Guidance

### Repeated runs

Enable persistent block caching for repeated clustering runs:

```bash
poetry run cluster-bsx ... --query-block-cache
```

This is usually the highest-impact optimization available without changing the data
pipeline.

If cache serialization starts to matter more than BSX reads, switch to
uncompressed cache payloads:

```bash
poetry run cluster-bsx ... --query-block-cache --query-block-cache-mode uncompressed
```

This increases cache size on disk but reduces CPU spent compressing and
decompressing `.npz` payloads.

### Query block planning

You can reduce `RegionReader.query(...)` fragmentation by merging nearby gene spans
into larger chromosome query blocks:

```bash
poetry run cluster-bsx ... --query-block-merge-gap-bp 1000
```

This is most useful when genes are dense and the default exact-span merge strategy
produces many small adjacent region queries.

### Hierarchical clustering

Hierarchical clustering is exact and quadratic in the number of retained genes.

Current default policy:

- hierarchical clustering is enabled
- it runs exactly up to `5,000` retained genes
- above that limit it falls back to a diagnostic subsample when configured
- if `cluster_source=hierarchical`, non-exact fallback is treated as an error

CLI knobs:

- `--no-hierarchical`
- `--hierarchical-mode`
- `--hierarchical-max-genes`
- `--hierarchical-subsample-genes`

Use `--hierarchical-mode exact` to force a full dendrogram when the extra quadratic
cost is acceptable.

### Silhouette score

Exact silhouette is also quadratic in retained gene count because it builds a full
pairwise distance matrix on the PCA embedding.

Current default policy:

- silhouette scoring is enabled
- it runs exactly up to `2,000` retained genes
- above that limit it switches to stratified sampled silhouette

CLI knobs:

- `--no-silhouette`
- `--silhouette-mode`
- `--silhouette-max-samples`
- `--silhouette-exact-max-genes`

Use `--silhouette-mode exact` only when validating small or medium matrices and the
extra distance work is acceptable.

### PCA

PCA now supports solver selection:

- `exact`
- `truncated`
- `auto`

Current default policy:

- `auto` uses exact SVD on smaller retained matrices
- above the configured matrix-size threshold it switches to truncated SVD
- if truncated SVD is not valid for the requested component count, it falls back to exact

CLI knobs:

- `--pca-solver`
- `--pca-exact-max-matrix-size`

If the matrix approaches millions of values, the practical options today are:

- reduce retained genes
- reduce profile bins
- skip hierarchical clustering first
- force `--pca-solver truncated` if exact SVD becomes the dominant backend cost

## Current Scaling Model

In rough order, the Python pipeline tends to become expensive in this sequence:

1. selective BSX block reads
2. exact hierarchical clustering
3. exact silhouette score
4. dense PCA on very large matrices

## Benchmark Scenarios

Synthetic perf-smoke fixtures now track three reference matrix shapes:

- `small`: `1,000 x 100`
- `medium`: `5,000 x 100`
- `large`: `10,000 x 100`

Default test runs execute only the `small` scenarios. To include the extended
`medium` and `large` perf smoke suite, run tests with:

```bash
BSX2_RUN_PERF=1 poetry run pytest tests/perf
```

The perf suite is intentionally a non-regression contract for policy selection
and successful execution, not a hard wall-clock benchmark.

## Reference Policy Contract

The reference perf scenarios are expected to exercise auto policies like this:

- `small (1,000 x 100)`: exact silhouette, exact hierarchical, exact PCA
- `medium (5,000 x 100)`: sampled silhouette, exact hierarchical, exact PCA
- `large (10,000 x 100)`: sampled silhouette, subsampled hierarchical, truncated PCA

These expectations are checked in `tests/perf` so that changes to scaling policy
become explicit review decisions instead of accidental regressions.

The metrics emitted by `cluster-bsx` now include:

- `performance_analysis`
- `scalability_guidance`
- `scalability_limits`
- `optimization_hotspots`
- `perf_non_regression_contract`

These fields are intended to make scaling behavior explicit in the generated
`metrics.json`.
