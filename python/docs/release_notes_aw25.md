# AW25 Release Notes

This note summarizes the state of the `bsx2` Python layer after the AW25 visualization
and clustering work.

## Delivered

- `bsx2.viz` as the main public visualization surface
- Compute/render split for metagene plots
- Annotation-driven metagene layouts via `AnnotProfileLayout`
- HoloViews renderers for line, heatmap, box, and violin plots
- PCA embedding, dendrogram, and cluster metagene plotting helpers
- Chromosome methylation map support

## Architectural Outcomes

- Metagene construction is no longer tied to fixed `up/body/down` assembly
- Annotation-driven multipart profiles are gene-centric
- Plotting now uses reusable compute data objects instead of HTML-first wrappers
- Clustering outputs include explicit performance and scalability guidance

## Release-Facing Notes

- Main plotting API: `bsx2.viz`
- Removed legacy surface: `bsx2.plots`
- Preferred annotation path: `compute_from_annot(..., layout=...)`
- Hierarchical clustering now has an explicit scalability policy with configurable limits

## Verification Snapshot

At the current release pass:

- `pytest python/tests` -> `55 passed`
- `ruff check python/src python/tests python/examples` -> `All checks passed`

See also:

- [Upgrade notes](upgrade_notes.md)
- [Performance notes](performance_notes.md)
- [Plotting support matrix](plotting_support_matrix.md)
