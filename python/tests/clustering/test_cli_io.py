from __future__ import annotations

import json
import shutil
import uuid
from pathlib import Path

import numpy as np

from bsx2.clustering.cli import parse_args
from bsx2.clustering.config import (
    AnnotationFormat,
    BlockCacheMode,
    ClusterConfig,
    ClusterSource,
    GeneProfileConfig,
    OutputConfig,
    TableFormat,
)
from bsx2.clustering.io import (
    build_cluster_artifacts,
    build_cluster_plot_data,
    write_cluster_outputs,
)
from bsx2.clustering.models import (
    ClusterArtifacts,
    FeatureBin,
    GeneAnnotation,
    GeneClusterResult,
    GeneProfileMatrix,
)
from bsx2.viz.clustering import ClusterMetageneData, GeneDendrogramData, GeneEmbeddingData


def _is_holoviews_object(obj) -> bool:
    return obj.__class__.__module__.startswith("holoviews")


def _result() -> GeneClusterResult:
    genes = [
        GeneAnnotation("gene_a", "chr1", 0, 100, "+", "GeneA"),
        GeneAnnotation("gene_b", "chr1", 100, 200, "-", "GeneB"),
    ]
    feature_bins = [
        FeatureBin("up_1", "up", 0, 0),
        FeatureBin("body_1", "body", 0, 1),
        FeatureBin("down_1", "down", 0, 2),
    ]
    values = np.array([[0.1, 0.2, 0.3], [0.8, 0.7, 0.6]], dtype=float)
    matrix = GeneProfileMatrix(
        genes=genes,
        feature_bins=feature_bins,
        values=values,
        gene_missing_rate=np.zeros(2, dtype=float),
        gene_variance=np.var(values, axis=1),
        feature_missing_rate=np.zeros(3, dtype=float),
        feature_variance=np.var(values, axis=0),
        metadata={"timings_s": {"total": 0.01}},
    )
    return GeneClusterResult(
        feature_matrix=matrix,
        labels=np.array([0, 1], dtype=np.int64),
        cluster_source="kmeans",
        kmeans_labels=np.array([0, 1], dtype=np.int64),
        hierarchical_labels=np.array([0, 1], dtype=np.int64),
        embedding=np.array([[1.0, 0.0], [-1.0, 0.0]], dtype=float),
        components=np.array([[0.5, 0.5, 0.5], [0.1, 0.2, 0.3]], dtype=float),
        centroids=np.array([[1.0, 0.0], [-1.0, 0.0]], dtype=float),
        explained_variance_ratio=np.array([0.9, 0.1], dtype=float),
        inertia=0.5,
        silhouette_score=None,
        linkage_matrix=np.array([[0.0, 1.0, 0.7, 2.0]], dtype=float),
        leaf_order=np.array([1, 0], dtype=np.int64),
        metadata={"timings_s": {"total": 0.02}, "pca_kmeans": {}, "hierarchical": {}},
    )


def test_parse_args_supports_gene_profile_cli() -> None:
    config = parse_args(
        [
            "--bsx",
            "sample.bsx",
            "--genes",
            "genes.gff",
            "--annotation-format",
            "gff",
            "--normalization",
            "row_zscore",
            "--cluster-source",
            "hierarchical",
            "--pca-solver",
            "truncated",
            "--pca-exact-max-matrix-size",
            "250000",
            "--query-block-cache-mode",
            "uncompressed",
            "--query-block-merge-gap-bp",
            "250",
            "--table-format",
            "csv",
            "--table-format",
            "json",
            "--silhouette-mode",
            "sampled",
            "--silhouette-max-samples",
            "1200",
            "--hierarchical-mode",
            "subsample",
            "--hierarchical-subsample-genes",
            "900",
            "--hierarchical-max-genes",
            "2500",
            "--output",
            "outdir",
        ]
    )

    assert config.bsx_path == Path("sample.bsx")
    assert config.annotation_path == Path("genes.gff")
    assert config.annotation_format is AnnotationFormat.GFF
    assert config.cluster_source is ClusterSource.HIERARCHICAL
    assert config.block_cache.mode is BlockCacheMode.UNCOMPRESSED
    assert config.read.query_block_merge_gap_bp == 250
    assert config.backend.pca_solver.value == "truncated"
    assert config.backend.pca_exact_max_matrix_size == 250000
    assert config.silhouette.mode.value == "sampled"
    assert config.silhouette.max_samples == 1200
    assert config.hierarchical.mode.value == "subsample"
    assert config.hierarchical.subsample_genes == 900
    assert config.hierarchical.max_genes == 2500
    assert config.output.table_formats == (TableFormat.CSV, TableFormat.JSON)


def test_write_cluster_outputs_emits_gene_level_files() -> None:
    outdir = Path("python/tests/clustering") / f"_tmp_io_{uuid.uuid4().hex}"
    outdir.mkdir(parents=True, exist_ok=False)
    try:
        config = ClusterConfig(
            bsx_path=Path("sample.bsx"),
            annotation_path=Path("genes.gff"),
            gene_profile=GeneProfileConfig(),
            output=OutputConfig(
                output_dir=outdir,
                table_formats=(TableFormat.TSV, TableFormat.JSON),
            ),
        )
        artifacts = write_cluster_outputs(_result(), config)

        expected = {
            "gene_profiles.tsv",
            "gene_profiles.json",
            "gene_clusters.tsv",
            "gene_clusters.json",
            "gene_embedding.tsv",
            "gene_embedding.json",
            "gene_bins.tsv",
            "gene_bins.json",
            "cluster_metagenes.tsv",
            "cluster_metagenes.json",
            "dendrogram_linkage.tsv",
            "dendrogram_linkage.json",
            "dendrogram_leaves.tsv",
            "dendrogram_leaves.json",
            "pca_loadings.tsv",
            "pca_loadings.json",
            "metrics.json",
        }
        assert isinstance(artifacts, ClusterArtifacts)
        assert {path.name for path in artifacts.paths.values()} == expected
        assert isinstance(artifacts.plot_data["pca"], GeneEmbeddingData)
        assert isinstance(artifacts.plot_data["dendrogram"], GeneDendrogramData)
        assert isinstance(artifacts.plot_data["cluster_metagene"], ClusterMetageneData)
        assert _is_holoviews_object(artifacts.plots["pca"])
        assert _is_holoviews_object(artifacts.plots["dendrogram"])
        assert _is_holoviews_object(artifacts.plots["cluster_metagene"])
        assert "profiles" in artifacts.tables
        assert "performance_analysis" in artifacts.metrics
        assert "scalability_guidance" in artifacts.metrics
        assert "scalability_limits" in artifacts.metrics
        assert "optimization_hotspots" in artifacts.metrics
        assert "quadratic_steps_detected" in artifacts.metrics
        assert "recommended_large_run_policy" in artifacts.metrics
        assert "perf_non_regression_contract" in artifacts.metrics

        metrics = json.loads((outdir / "metrics.json").read_text(encoding="utf-8"))
        assert metrics["n_genes"] == 2
        assert metrics["n_profile_bins"] == 3
        assert "io_timings_s" in metrics
        assert "scalability_guidance" in metrics
        assert "scalability_limits" in metrics
        assert "quadratic_steps_detected" in metrics
        assert "recommended_large_run_policy" in metrics
        assert "perf_non_regression_contract" in metrics
        assert "query_block_cache_mode" in metrics["scalability_guidance"]
        assert "query_block_merge_gap_bp" in metrics["scalability_guidance"]
        assert metrics["performance_analysis"]["version"] == 3
        assert "perf_harness" in metrics["performance_analysis"]
        assert "reference_scenarios" in metrics["scalability_limits"]
        assert metrics["plots_are_serialized_by_library"] is False
    finally:
        shutil.rmtree(outdir, ignore_errors=True)


def test_build_cluster_artifacts_returns_plot_objects_without_writing() -> None:
    artifacts = build_cluster_artifacts(
        _result(),
        ClusterConfig(
            bsx_path=Path("sample.bsx"),
            annotation_path=Path("genes.gff"),
            gene_profile=GeneProfileConfig(),
            output=OutputConfig(
                output_dir=Path("outdir"),
                write_table_files=False,
                write_metrics_file=False,
            ),
        ),
    )

    assert isinstance(artifacts.plot_data["pca"], GeneEmbeddingData)
    assert _is_holoviews_object(artifacts.plots["pca"])
    assert _is_holoviews_object(artifacts.plots["dendrogram"])
    assert "cluster_metagene" in artifacts.plots
    assert "profiles" in artifacts.tables
    assert artifacts.paths == {}


def test_build_cluster_artifacts_supports_missing_hierarchical_outputs() -> None:
    result = _result()
    result.hierarchical_labels = None
    result.linkage_matrix = None
    result.leaf_order = None
    result.metadata["hierarchical"] = {}

    artifacts = build_cluster_artifacts(
        result,
        ClusterConfig(
            bsx_path=Path("sample.bsx"),
            annotation_path=Path("genes.gff"),
            gene_profile=GeneProfileConfig(),
            output=OutputConfig(
                output_dir=Path("outdir"),
                write_table_files=False,
                write_metrics_file=False,
            ),
        ),
    )

    assert "linkage" not in artifacts.tables
    assert "leaves" not in artifacts.tables
    assert artifacts.tables["clusters"]["hierarchical_cluster"].null_count() == 2
    assert "dendrogram" not in artifacts.plot_data
    assert "dendrogram" not in artifacts.plots


def test_build_cluster_plot_data_returns_dataclass_layer() -> None:
    plot_data = build_cluster_plot_data(_result())

    assert isinstance(plot_data["pca"], GeneEmbeddingData)
    assert isinstance(plot_data["dendrogram"], GeneDendrogramData)
    assert isinstance(plot_data["cluster_metagene"], ClusterMetageneData)
