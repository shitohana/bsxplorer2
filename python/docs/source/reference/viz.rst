Visualization
=============

``bsx2.viz`` is the public plotting surface.

Recommended flow:

1. Prepare profile or clustering data in ``bsx2.viz.compute``.
2. Render HoloViews objects through ``bsx2.viz`` wrappers or composers.
3. Export through the backend that matches the current workflow.

Main public entrypoints include:

- ``compute_discrete_regions(...)``
- ``compute_from_annot(...)``
- ``build_manual_metagene(...)``
- ``build_annotation_metagene(...)``
- ``line_plot(...)``
- ``heatmap(...)``
- ``box_plot(...)``
- ``violin_plot(...)``
- chromosome and clustering plot builders

.. automodule:: bsx2.viz
   :members:
   :undoc-members:
   :show-inheritance:
