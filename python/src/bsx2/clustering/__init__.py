from .backend import cluster_gene_profiles, run_pca_kmeans
from .config import (
    AnnotationFormat,
    BackendConfig,
    BlockCacheConfig,
    ClusterConfig,
    ClusterSource,
    GeneProfileConfig,
    HierarchicalConfig,
    HierarchicalDistance,
    HierarchicalLinkage,
    NormalizationMode,
    OutputConfig,
    ReadConfig,
    TableFormat,
)
from .gene_annotation import infer_annotation_format, load_annotation_store, load_gene_annotations
from .gene_profile import (
    build_feature_bins,
    build_gene_profile_matrix,
    build_metagene_bins,
    clear_region_reader_cache,
)
from .hierarchical import run_hierarchical
from .metagene import cluster_metagene_summary
from .models import (
    ClusterArtifacts,
    FeatureBin,
    GeneAnnotation,
    GeneClusterResult,
    GeneProfileMatrix,
    ProfileBin,
)


def build_cluster_artifacts(*args, **kwargs):
    from .io import build_cluster_artifacts as _impl

    return _impl(*args, **kwargs)


def build_cluster_metrics(*args, **kwargs):
    from .io import build_cluster_metrics as _impl

    return _impl(*args, **kwargs)


def build_cluster_plot_data(*args, **kwargs):
    from .io import build_cluster_plot_data as _impl

    return _impl(*args, **kwargs)


def build_cluster_plots(*args, **kwargs):
    from .io import build_cluster_plots as _impl

    return _impl(*args, **kwargs)


def build_cluster_tables(*args, **kwargs):
    from .io import build_cluster_tables as _impl

    return _impl(*args, **kwargs)


def build_dendrogram_plot(*args, **kwargs):
    from .io import build_dendrogram_plot as _impl

    return _impl(*args, **kwargs)


def build_pca_plot(*args, **kwargs):
    from .io import build_pca_plot as _impl

    return _impl(*args, **kwargs)


def main(*args, **kwargs):
    from .cli import main as _impl

    return _impl(*args, **kwargs)


def write_cluster_outputs(*args, **kwargs):
    from .io import write_cluster_outputs as _impl

    return _impl(*args, **kwargs)


__all__ = [
    "AnnotationFormat",
    "BackendConfig",
    "BlockCacheConfig",
    "ClusterConfig",
    "ClusterSource",
    "ClusterArtifacts",
    "FeatureBin",
    "GeneAnnotation",
    "GeneClusterResult",
    "GeneProfileConfig",
    "GeneProfileMatrix",
    "HierarchicalConfig",
    "HierarchicalDistance",
    "HierarchicalLinkage",
    "NormalizationMode",
    "OutputConfig",
    "ProfileBin",
    "ReadConfig",
    "TableFormat",
    "build_feature_bins",
    "build_cluster_artifacts",
    "build_cluster_metrics",
    "build_cluster_plot_data",
    "build_cluster_plots",
    "build_cluster_tables",
    "build_dendrogram_plot",
    "build_gene_profile_matrix",
    "build_metagene_bins",
    "clear_region_reader_cache",
    "build_pca_plot",
    "cluster_gene_profiles",
    "cluster_metagene_summary",
    "infer_annotation_format",
    "load_annotation_store",
    "load_gene_annotations",
    "main",
    "run_hierarchical",
    "run_pca_kmeans",
    "write_cluster_outputs",
]
