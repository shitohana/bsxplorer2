# Code Inventory Summary

This summary maps the main Python modules that now define the supported AW25
visualization and clustering surface.

## Public Visualization Surface

### `bsx2.viz`

Primary user-facing plotting API.

Key areas:

- metagene compute wrappers
- HoloViews renderers
- clustering plot builders
- chromosome methylation map

Important modules:

- `src/bsx2/viz/__init__.py`
- `src/bsx2/viz/metagene.py`
- `src/bsx2/viz/clustering.py`
- `src/bsx2/viz/chrmap.py`

### `bsx2.viz.compute`

Compute-only data builders and intermediate data structures.

Important modules:

- `src/bsx2/viz/compute/metagene.py`
- `src/bsx2/viz/compute/data.py`
- `src/bsx2/viz/compute/clustering.py`
- `src/bsx2/viz/compute/chrmap.py`

### `bsx2.viz.render`

Render-only HoloViews composition layer.

Important modules:

- `src/bsx2/viz/render/line.py`
- `src/bsx2/viz/render/heatmap.py`
- `src/bsx2/viz/render/box.py`
- `src/bsx2/viz/render/violin.py`
- `src/bsx2/viz/render/chrmap.py`

## Public Clustering Surface

### `bsx2.clustering`

Main clustering entrypoint and exported configuration/model types.

Important modules:

- `src/bsx2/clustering/__init__.py`
- `src/bsx2/clustering/config.py`
- `src/bsx2/clustering/models.py`
- `src/bsx2/clustering/cli.py`

### Compute and backend path

- `src/bsx2/clustering/gene_profile.py`
- `src/bsx2/clustering/backend.py`
- `src/bsx2/clustering/hierarchical.py`
- `src/bsx2/clustering/io.py`

## Examples

The repo ships focused Python examples for each supported surface:

- `examples/metagene_from_contigs.py`
- `examples/metagene_from_annot.py`
- `examples/clustering_plots.py`
- `examples/chromosome_map.py`

## Tests

The supported Python surface is primarily covered by:

- `tests/viz`
- `tests/clustering`

These suites now cover:

- metagene rendering and API surface
- annotation layout edge cases
- clustering outputs and metrics
- hierarchical skip policy
- chromosome map compute/render behavior
