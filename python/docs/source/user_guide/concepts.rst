Concepts
========

This section defines the core plotting concepts used throughout BSXplorer2.

DiscreteRegionData
------------------

Most plotting entrypoints in ``bsx2.viz`` consume a ``DiscreteRegionData``
object. It is the common intermediate representation shared by the compute and
render layers.

Conceptually, ``DiscreteRegionData`` stores:

- one normalized x-axis per region or gene
- one methylation density vector aligned to that axis
- one optional label per profile

Because the render layer only needs this normalized representation, the same
plotting functions can work with multiple compute pipelines.

Scaled region profile vs metagene
---------------------------------

BSXplorer2 intentionally distinguishes two different ideas:

Scaled region profile
    Each input region is scaled independently into the interval ``[0, 1]`` and
    then aggregated on that shared relative coordinate. This is the ordinary
    arbitrary-region mode driven by ``compute_discrete_regions(...)``.

Metagene
    A composed profile in which biologically distinct parts such as upstream,
    body, and downstream are first assembled explicitly and only then combined
    into a shared normalized axis. This mode is driven by
    ``compute_from_annot(...)`` or by a manual multi-part composition.

The distinction matters biologically. A scaled region profile does not imply
that upstream, body, or downstream were encoded in the data unless those parts
were assembled explicitly before rendering.

Annotation-driven vs manual metagene assembly
---------------------------------------------

BSXplorer2 supports two metagene assembly entrypoints that produce the same
render-compatible output.

Annotation-driven layout
^^^^^^^^^^^^^^^^^^^^^^^^

This high-level mode takes ``RegionReader + HcAnnotStore`` and a declared
layout, then resolves parts gene-by-gene. It is the preferred path for
classical AW25-style metagene plots.

Manual composed layout
^^^^^^^^^^^^^^^^^^^^^^

This low-level mode starts from already prepared part-specific contigs and
labels. Each part is normalized independently and then concatenated into a
shared metagene profile.

Both modes converge to the same ``DiscreteRegionData`` contract.

Shared render layer
-------------------

The render layer in ``bsx2.viz.render`` does not care whether the profiles came
from annotation-driven assembly or manual composition. The common flow is:

1. Build ``DiscreteRegionData`` in ``bsx2.viz.compute``.
2. Render the data with ``line_plot(...)``, ``heatmap(...)``, ``box_plot(...)``,
   ``violin_plot(...)``, or clustering/chromosome composers.
3. Export with the HoloViews/Plotly backend that fits the current workflow.

Grouped distributions
---------------------

For box and violin plots, the biologically meaningful default metagene
interpretation is group-level aggregation across many genes or regions.

``per_region=False``
    The default metagene-friendly mode. Values are grouped by metagene window
    or segment across the full set of genes.

``per_region=True``
    A region-summary mode. It is still useful, but it should not be treated as
    the primary metagene box/violin interpretation.
