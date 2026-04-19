How plots are built
===================

This page documents the compute pipeline used to construct each plotting
surface in BSXplorer2.

Metagene line plot
------------------

Inputs:

- ``report.bsx``
- ``annot.gff`` or manually prepared contigs
- methylation context such as ``CG``

Typical pipeline:

1. Build a ``RegionReader`` and apply context filters.
2. Build a metagene-compatible ``DiscreteRegionData`` object either by:

   - ``compute_from_annot(...)`` for annotation-driven assembly, or
   - ``build_manual_metagene(...)`` for explicit part composition.

3. Render with ``line_plot(...)``.
4. Rebin the profiles into output windows.
5. Aggregate values per window using the selected aggregation strategy.
6. Optionally smooth mean or median profiles with Savitzky-Golay smoothing.

Metagene heatmap
----------------

The heatmap reuses the same metagene input data but renders it as a matrix.

1. Each region or gene is rebinned to the requested number of windows.
2. The result becomes a matrix of shape ``n_profiles x n_windows``.
3. Rows are ranked by a scoring strategy such as ``mean`` or ``body_mean``.
4. If requested, rows are compressed to ``rank_rows`` for visualization.
5. The matrix is rendered as a heatmap with optional segment guides.

Metagene box and violin plots
-----------------------------

These plots use grouped distributions derived from metagene windows or
segments.

``distribution_mode="windows"``
    One distribution per metagene window.

``distribution_mode="segments"``
    One distribution per named segment such as ``up``, ``body``, or ``down``.

The metagene-friendly default is ``per_region=False``, which aggregates values
across the group of genes rather than summarizing one value per gene.

Chromosome methylation map
--------------------------

Pipeline:

1. Read chromosome-scale methylation values from ``report.bsx``.
2. Aggregate values into chromosome windows using a user-defined bin size.
3. Render the result with ``build_chromosome_methylation_map(...)``.

This plot is independent of the gene limit used by metagene or clustering
steps.

Gene PCA
--------

Pipeline:

1. Build a gene profile matrix from methylation data and annotation.
2. Apply preprocessing and feature filtering.
3. Run PCA in the clustering backend.
4. Convert the result with ``build_gene_embedding_data(...)``.
5. Render with ``GeneEmbeddingPlotComposer``.

Gene dendrogram
---------------

Pipeline:

1. Build the same gene profile matrix used by PCA.
2. Run hierarchical clustering when the effective gene count is within policy.
3. Convert the result with ``build_gene_dendrogram_data(...)``.
4. Render with ``GeneDendrogramPlotComposer``.

This plot is more computationally expensive than PCA and is therefore more
sensitive to limit settings and scalability policy.

Cluster metagene
----------------

Pipeline:

1. Build or load a ``GeneClusterResult``.
2. Select clusters or genes of interest.
3. Aggregate cluster-level mean profiles with
   ``build_cluster_metagene_data(...)``.
4. Render with ``build_cluster_metagene_plot(...)``.

Why the architecture stays transparent
--------------------------------------

The key design principle is that compute and render are separated:

- compute code prepares normalized or clustered data
- render code only visualizes that prepared data

This keeps the plotting API transparent, debuggable, and reusable across
annotation-driven and manual workflows.
