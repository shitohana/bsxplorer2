"""Compute-layer exports for BSX2 visualization helpers.

The full compute API imports extension-backed and optional plotting-adjacent
modules. Lightweight cache utilities remain importable from a source checkout
when those optional dependencies are unavailable.
"""

try:
    from .chrline import ChrLineData, ChrLineTrack, compute_chr_line_data
    from .chrmap import (
        ChromosomeMethylationMapData,
        ChromosomeMethylationTrack,
        compute_chromosome_methylation_map_data,
    )
    from .clustering import (
        ClusterMetageneData,
        ClusterMetageneGroup,
        GeneDendrogramData,
        GeneEmbeddingData,
        build_cluster_metagene_data,
        build_gene_dendrogram_data,
        build_gene_embedding_data,
        cluster_profile_segments,
    )
    from .data import DiscreteRegionData, SegmentData, SortBy
    from .discrete_region_cache import (
        discrete_region_data_fingerprint,
        load_discrete_region_data,
        save_discrete_region_data,
        write_discrete_region_cache_manifest,
    )
    from .distribution import (
        DistributionPlotData,
        build_box_distribution_data,
        build_violin_distribution_data,
    )
    from .input_compat import (
        SeqnameCompatibilityReport,
        build_seqname_compatibility_report,
        normalize_seqname,
        read_bsx_seqnames,
        read_gff_seqnames,
    )
    from .metagene import (
        AnnotProfileLayout,
        AnnotProfilePart,
        MetageneProfileSegment,
        collect_contigs_from_hcannot,
        collect_layout_parts_from_hcannot,
        collect_parts_from_hcannot,
        compose_layout_drd,
        compute_discrete_regions,
        compute_from_annot,
        segment_boundaries,
        segments_total_bins,
    )
    from .metagene_layouts import build_annotation_metagene, build_manual_metagene
    from .studio_workspaces import (
        ManualRegionSpec,
        PreparedClusterDendrogramFamily,
        PreparedClusterFamily,
        PreparedClusterMatrixWorkspace,
        PreparedMetageneFamily,
        build_studio_cluster_config,
        cluster_dendrogram_family_key,
        cluster_matrix_family_key,
        cluster_plot_family_key,
        metagene_family_key,
        prepare_cluster_dendrogram_family,
        prepare_cluster_family,
        prepare_cluster_matrix_workspace,
        prepare_metagene_family,
    )
    from .windowing import _bin_points_windows_fast, _rank_compress
except ModuleNotFoundError as exc:
    if exc.name not in {"bsx2._bsx2", "beartype", "holoviews", "panel"}:
        raise
    from .discrete_region_cache import (
        DiscreteRegionData,
        discrete_region_data_fingerprint,
        load_discrete_region_data,
        save_discrete_region_data,
        write_discrete_region_cache_manifest,
    )

    __all__ = [
        "DiscreteRegionData",
        "discrete_region_data_fingerprint",
        "load_discrete_region_data",
        "save_discrete_region_data",
        "write_discrete_region_cache_manifest",
    ]
else:
    __all__ = [
        "ChrLineData",
        "ChrLineTrack",
        "compute_chr_line_data",
        "ChromosomeMethylationMapData",
        "ChromosomeMethylationTrack",
        "compute_chromosome_methylation_map_data",
        "ClusterMetageneData",
        "ClusterMetageneGroup",
        "GeneDendrogramData",
        "GeneEmbeddingData",
        "build_cluster_metagene_data",
        "build_gene_dendrogram_data",
        "build_gene_embedding_data",
        "cluster_profile_segments",
        "DiscreteRegionData",
        "DistributionPlotData",
        "SegmentData",
        "SortBy",
        "build_box_distribution_data",
        "build_violin_distribution_data",
        "SeqnameCompatibilityReport",
        "normalize_seqname",
        "read_bsx_seqnames",
        "read_gff_seqnames",
        "build_seqname_compatibility_report",
        "AnnotProfilePart",
        "AnnotProfileLayout",
        "MetageneProfileSegment",
        "collect_contigs_from_hcannot",
        "collect_layout_parts_from_hcannot",
        "collect_parts_from_hcannot",
        "compose_layout_drd",
        "compute_discrete_regions",
        "compute_from_annot",
        "build_annotation_metagene",
        "build_manual_metagene",
        "PreparedMetageneFamily",
        "ManualRegionSpec",
        "PreparedClusterMatrixWorkspace",
        "PreparedClusterFamily",
        "PreparedClusterDendrogramFamily",
        "build_studio_cluster_config",
        "metagene_family_key",
        "cluster_matrix_family_key",
        "cluster_plot_family_key",
        "cluster_dendrogram_family_key",
        "prepare_metagene_family",
        "prepare_cluster_matrix_workspace",
        "prepare_cluster_family",
        "prepare_cluster_dendrogram_family",
        "segment_boundaries",
        "segments_total_bins",
        "_bin_points_windows_fast",
        "_rank_compress",
        "discrete_region_data_fingerprint",
        "load_discrete_region_data",
        "save_discrete_region_data",
        "write_discrete_region_cache_manifest",
    ]
