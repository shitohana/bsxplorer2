Parameter glossary
==================

This glossary explains the main parameters exposed by plotting and clustering
entrypoints.

``segments``
    Ordered list of named profile segments such as ``up``, ``body``, and
    ``down``. Used for bin allocation and segment-aware rendering.

``n_windows``
    Number of output bins used when rebinding a profile for display. If omitted,
    the total number of bins is derived from the segment layout.

``smooth``
    Optional Savitzky-Golay smoothing configuration for line plots. Intended for
    visual smoothing of aggregated curves, not for changing the underlying data.

``rank_rows``
    Number of rows kept in heatmap output after rank compression.

``per_region``
    Distribution mode switch for box and violin plots. ``False`` is the default
    metagene-friendly grouped mode. ``True`` summarizes per individual region.

``distribution_mode``
    Grouping strategy for box and violin plots. ``"windows"`` groups by
    metagene windows, while ``"segments"`` groups by named segments.

``flank_bp``
    Size of synthetic flanks in base pairs for annotation-driven or manual
    metagene layouts.

``limit_regions``
    Maximum number of input regions included in region-profile or metagene
    preparation.

``limit_genes``
    Maximum number of genes included in clustering workflows.

``min_coverage``
    Minimum read coverage threshold applied before clustering profile assembly.

``query_block_merge_gap_bp``
    Distance threshold used to merge nearby gene spans into fewer BSX region
    queries for clustering.

``hierarchical_max_genes``
    Upper bound for exact hierarchical clustering before policy-based fallback
    or skipping is applied.
