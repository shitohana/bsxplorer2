# Upgrade Notes

This note summarizes the practical migration path from the old prototype plotting
layer to the current `bsx2.viz` API.

## Key Changes

- Plotting is now HoloViews-first.
- Compute and render responsibilities are separated.
- Annotation-driven metagenes use `layout=` instead of fixed `up/body/down` assembly.
- Legacy `combine_parts`-style annotation assembly is removed.
- The old `bsx2.plots` module is removed.
- Built-in HTML helper generation is removed from the package surface.

## Migration Map

| Old idea | Current approach |
| --- | --- |
| HTML-first plot helpers | Compute with `bsx2.viz`, then render/save with the HoloViews backend you choose |
| Fixed `up/body/down` combined annotation profile | `AnnotProfileLayout((AnnotProfilePart(...), ...))` |
| Plot module owns both aggregation and rendering | `compute_*` builds data, `line_plot` / `heatmap` / `box_plot` / `violin_plot` render it |
| `bsx2.plots` public surface | `bsx2.viz` |
| Plot-specific internal data coupling | Shared `DiscreteRegionData` compute contract |

## Recommended Migration Pattern

### Arbitrary regions

1. Use `compute_discrete_regions(...)`.
2. Pass the returned `DiscreteRegionData` to one or more renderers.

### Annotation-driven metagenes

1. Define an `AnnotProfileLayout`.
2. Use `compute_from_annot(..., layout=...)`.
3. Render the result with the desired plot builder.

### HTML export

Render HoloViews objects and save/export with the backend you choose.

## Non-Goals

The current architecture intentionally does not restore the old plotting layer one to
one. The goal is equivalent or improved functionality on top of a cleaner compute/render
split, not strict preservation of the old module structure.
