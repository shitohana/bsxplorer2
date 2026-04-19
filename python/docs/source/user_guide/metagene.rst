Metagene analysis
=================

BSXplorer2 supports two layout-based ways to assemble metagene-compatible
``DiscreteRegionData`` objects:

Annotation-driven layout
------------------------

This high-level path uses ``RegionReader + HcAnnotStore`` and automatically
builds up/body/down parts, flanks, labels, and orientation handling.

Manual composed layout
----------------------

This low-level path lets you prepare region parts manually, compute each part
with ``compute_discrete_regions(...)``, and compose them into one shared
metagene layout.

Shared render layer
-------------------

Both modes feed the same rendering API in ``bsx2.viz``:

- line plot
- heatmap
- box plot
- violin plot
