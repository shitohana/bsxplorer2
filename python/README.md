# bsx2

`bsx2` is the Python package for `bsxplorer2`: fast bisulfite sequencing data IO, normalized-region aggregation, clustering utilities, and HoloViews-based visualization.

## Status

Implemented in the current Python layer:

- report and `.bsx` IO
- scaled aggregation from arbitrary regions
- metagene aggregation from `RegionReader + HcAnnotStore`
- HoloViews plots for normalized region and annotation-layout profiles:
  - line plot
  - heatmap
  - box plot
  - violin plot
- clustering plots and helpers:
  - PCA embedding plot
  - hierarchical dendrogram
  - cluster metagene summary
- chromosome methylation map

Out of scope for the public plotting API:

- built-in HTML helper generation
- Plotly-specific wrapper modules from the old prototype

The plotting layer now returns HoloViews objects. Rendering, saving, and backend-specific export are left to the user.

Further docs:

- [Plotting support matrix](docs/plotting_support_matrix.md)
- [Upgrade notes](docs/upgrade_notes.md)
- [Performance notes](docs/performance_notes.md)
- [Release notes](docs/release_notes_aw25.md)
- [Code inventory summary](docs/code_inventory_summary.md)
- [Manual acceptance checklist](docs/manual_acceptance_checklist.csv)

## Installation

The project uses Poetry.

```bash
cd python
poetry install
```

## Development Commands

Run tests:

```bash
cd python
poetry run pytest tests
```

Run Ruff:

```bash
cd python
poetry run ruff check src tests
poetry run ruff format src tests
```

Run the clustering CLI:

```bash
cd python
poetry run cluster-bsx --help
```

For repeated clustering runs on the same BSX file, consider enabling the persistent
query-block cache:

```bash
cd python
poetry run cluster-bsx ... --query-block-cache
```

## Plotting API

### Support matrix

| Surface | Status | Notes |
| --- | --- | --- |
| `bsx2.viz` compute + HoloViews renderers | Supported | Main public plotting API |
| `bsx2.plots` | Removed | Legacy prototype surface |

See the full [plotting support matrix](docs/plotting_support_matrix.md) for the explicit contract.

### Profile data preparation

Use one of the data constructors first:

- `bsx2.viz.compute_discrete_regions(...)`
- `bsx2.viz.compute_from_annot(...)`

Both return `DiscreteRegionData`, which can be rendered independently of how it was computed.

`compute_discrete_regions(...)` scales each input region independently into the
relative interval `[0, 1]` and aggregates signal by that shared coordinate. In
other words, the ordinary arbitrary-region path is a scaled region profile, not
an implicit `upstream/body/downstream` metagene unless you explicitly compose a
segmented layout yourself.

### Annotation layouts

`compute_from_annot(...)` supports two modes:

- single-part aggregation via `feature_type=...`
- ordered multi-part aggregation via `layout=AnnotProfileLayout(...)`

The `layout=` mode is the preferred path for annotation-driven metagenes. Supported
`AnnotProfilePart.source` values are:

- `"feature"`: resolve a specific annotation feature type per gene
- `"gene"`: use the gene body itself
- `"flank5"`: synthesize a 5' flank from the gene anchor
- `"flank3"`: synthesize a 3' flank from the gene anchor

Parts are collected gene-by-gene and concatenated in declared order, so the API is not
restricted to `promoter/body/terminator`.

### Profile plots

Public plotting wrappers:

- `bsx2.viz.line_plot(...) -> hv.Curve | hv.Overlay`
- `bsx2.viz.heatmap(...) -> hv.Image | hv.Overlay`
- `bsx2.viz.box_plot(...) -> hv.BoxWhisker`
- `bsx2.viz.violin_plot(...) -> hv.Violin`

### Clustering plots

For `GeneClusterResult`:

- `bsx2.viz.build_gene_embedding_data(...)`
- `bsx2.viz.build_gene_dendrogram_data(...)`
- `bsx2.viz.build_cluster_metagene_data(...)`
- `bsx2.viz.cluster_metagene_plot(...)`
- `bsx2.viz.GeneEmbeddingPlotComposer`
- `bsx2.viz.GeneDendrogramPlotComposer`

### Clustering performance

The clustering pipeline is still usually dominated by selective BSX reads rather than
PCA or KMeans. For larger retained gene sets:

- enable `--query-block-cache` for repeated runs
- use `--query-block-cache-mode uncompressed` if cache CPU becomes noticeable
- use `--query-block-merge-gap-bp` to merge nearby gene spans into fewer region queries
- expect exact hierarchical clustering and exact silhouette to become the first
  quadratic backend steps
- note that silhouette switches to sampled mode by default above `2,000` retained genes
- note that hierarchical clustering switches away from full exact mode above `5,000`
  retained genes unless the policy is overridden
- note that PCA can switch from exact to truncated SVD in auto mode once the retained
  matrix exceeds the configured size threshold
- note that `tests/perf` runs only the `small` reference scenarios by default; set
  `BSX2_RUN_PERF=1` to include `medium` and `large` non-regression perf smoke runs

See [performance notes](docs/performance_notes.md) for the current scaling model and
CLI knobs.

### Chromosome methylation map

Public chromosome-wide API:

- `bsx2.viz.compute_chromosome_methylation_map_data(...)`
- `bsx2.viz.build_chromosome_methylation_map(...)`
- `bsx2.viz.chromosome_methylation_map(...)`

`chromosome_methylation_map(...)` is the one-step wrapper. It computes chromosome-window aggregates and returns a HoloViews plot.

### Export

The preferred flow is:

1. compute data with `bsx2.viz`
2. render HoloViews objects
3. save/export with the backend you actually need

## Usage Examples

### Scaled region profile from arbitrary contigs

```python
from bsx2 import Contig, RegionReader, Strand
from bsx2.viz import (
    compute_discrete_regions,
    heatmap,
    line_plot,
)

reader = RegionReader("/path/to/report.bsx")
contigs = [
    Contig("chr1", 100_000, 120_000, Strand.Forward),
    Contig("chr2", 50_000, 70_000, Strand.Reverse),
]

drd = compute_discrete_regions(
    reader,
    contigs,
)

curve = line_plot(drd, name="sample")
hm = heatmap(drd)
```

### Metagene from annotations

```python
from bsx2 import HcAnnotStore, RegionReader
from bsx2.viz import (
    AnnotProfileLayout,
    AnnotProfilePart,
    box_plot,
    compute_from_annot,
    violin_plot,
)

reader = RegionReader("/path/to/report.bsx")
annot = HcAnnotStore.from_gff("/path/to/annot.gff")
layout = AnnotProfileLayout(
    (
        AnnotProfilePart("promoter", 25, source="flank5", flank_bp=2000),
        AnnotProfilePart("gene", 50, source="gene"),
        AnnotProfilePart("terminator", 25, source="flank3", flank_bp=2000),
    )
)
segments = list(layout.segments)

drd = compute_from_annot(
    reader,
    annot,
    layout=layout,
    segments=segments,
)

box = box_plot(drd, segments=segments)
violin = violin_plot(drd, segments=segments)
```

### Cluster metagene

```python
from bsx2.viz import cluster_metagene_plot

plot = cluster_metagene_plot(
    cluster_result,
    cluster_ids=[0, 1],
)
```

If you need the intermediate data object first:

```python
from bsx2.viz import build_cluster_metagene_data, build_cluster_metagene_plot

data = build_cluster_metagene_data(cluster_result, cluster_ids=[0], collapse=False)
plot = build_cluster_metagene_plot(data)
```

### PCA and dendrogram plots

```python
from bsx2.viz import (
    GeneDendrogramPlotComposer,
    GeneEmbeddingPlotComposer,
    build_gene_dendrogram_data,
    build_gene_embedding_data,
)

embedding_data = build_gene_embedding_data(cluster_result)
pca_plot = GeneEmbeddingPlotComposer().add_data(embedding_data).finish()

dendrogram_data = build_gene_dendrogram_data(cluster_result)
if dendrogram_data is not None:
    dendrogram_plot = GeneDendrogramPlotComposer().add_data(dendrogram_data).finish()
```

### Chromosome methylation map

```python
from bsx2 import Context, RegionReader
from bsx2.viz import chromosome_methylation_map

reader = RegionReader("/path/to/report.bsx")

plot = chromosome_methylation_map(
    reader,
    context=Context.CG,
    bin_size_bp=50_000,
    chr_lengths={"chr1": 24_000_000, "chr2": 19_000_000},
    name="sample",
)
```

## Notes

- `compute_discrete_regions(...)` produces a scaled region profile on a shared relative coordinate.
- True multi-part profiles remain explicit: use `AnnotProfileLayout(...)` when you need promoter/gene/terminator style structure.
- Data preparation and visualization are intentionally separated.
- Smaller plotting components do not depend on parent structures beyond the data they receive.
- HoloViews objects can be customized further with `.opts(...)`.
- Built-in HTML helper generation is intentionally not part of the main `bsx2.viz` surface.

## Example Scripts

- [examples/metagene_from_annot.py](examples/metagene_from_annot.py) builds line, heatmap, box, and violin plots from `RegionReader + HcAnnotStore`.
- [examples/metagene_from_contigs.py](examples/metagene_from_contigs.py) shows arbitrary-region metagene computation from explicit contigs.
- [examples/clustering_plots.py](examples/clustering_plots.py) builds PCA, dendrogram, and cluster metagene plots from a full clustering run.
- [examples/chromosome_map.py](examples/chromosome_map.py) builds a chromosome methylation map.

## Style

- runtime validation: `beartype`
- tests: `pytest`
- lint/format: `ruff`
- docstrings: NumPy style
