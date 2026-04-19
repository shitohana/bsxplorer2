# Plotting Support Matrix

This document defines the intended public plotting surface for the `bsx2` Python
package as of the current `bsx2.viz` architecture.

## Public Status

| Surface | Status | Notes |
| --- | --- | --- |
| `bsx2.viz.compute_discrete_regions(...)` | Supported | Compute-first API for arbitrary contigs. Returns `DiscreteRegionData`. |
| `bsx2.viz.compute_from_annot(...)` | Supported | Main annotation-driven metagene API. Supports single-part `feature_type=...` and preferred multipart `layout=...`. |
| `bsx2.viz.AnnotProfileLayout` / `AnnotProfilePart` | Supported | First-class layout system for gene-centric annotation metagenes. |
| `bsx2.viz.line_plot(...)` | Supported | Returns HoloViews line output from `DiscreteRegionData`. |
| `bsx2.viz.heatmap(...)` | Supported | Returns HoloViews heatmap output from `DiscreteRegionData`. |
| `bsx2.viz.box_plot(...)` | Supported | Returns HoloViews distribution plot from `DiscreteRegionData`. |
| `bsx2.viz.violin_plot(...)` | Supported | Returns HoloViews distribution plot from `DiscreteRegionData`. |
| `bsx2.viz.build_gene_embedding_data(...)` | Supported | Compute-first clustering visualization API. |
| `bsx2.viz.build_gene_dendrogram_data(...)` | Supported | Compute-first clustering visualization API. |
| `bsx2.viz.build_cluster_metagene_data(...)` | Supported | Cluster metagene summary data builder. |
| `bsx2.viz.cluster_metagene_plot(...)` | Supported | One-step cluster metagene plot wrapper. |
| `bsx2.viz.chromosome_methylation_map(...)` | Supported | One-step chromosome methylation map wrapper. |
## Intentionally Out Of Scope

These surfaces are not part of the supported plotting architecture:

- built-in HTML helper generation
- legacy `bsx2.plots`
- old Plotly-specific wrapper modules from the prototype layer
- composer `.to_html()` methods on `bsx2.viz` plotting classes
- HTML-first APIs as the primary visualization contract

The supported rendering contract is HoloViews-first. HTML export is possible, but it
is a rendering concern handled outside the main package abstraction.

## Annotation Layout Guidance

Preferred annotation-driven metagene usage:

1. Build an `AnnotProfileLayout`.
2. Compute `DiscreteRegionData` with `compute_from_annot(..., layout=...)`.
3. Render with `line_plot`, `heatmap`, `box_plot`, or `violin_plot`.

Supported `AnnotProfilePart.source` values:

- `"feature"`
- `"gene"`
- `"flank5"`
- `"flank3"`

The layout system is gene-centric:

- parts are collected per gene
- orphan feature parents without a gene anchor are ignored
- multipart profiles are concatenated in declared order
- the API is not restricted to `promoter/body/terminator`
