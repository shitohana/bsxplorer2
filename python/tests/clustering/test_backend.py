from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from bsx2.clustering.backend import cluster_gene_profiles, run_pca_kmeans
from bsx2.clustering.config import (
    BackendConfig,
    ClusterConfig,
    ClusterSource,
    GeneProfileConfig,
    HierarchicalConfig,
    HierarchicalMode,
    NormalizationMode,
    OutputConfig,
    PcaSolver,
    SilhouetteConfig,
    SilhouetteMode,
)
from bsx2.clustering.models import FeatureBin, GeneAnnotation, GeneProfileMatrix


def _matrix() -> GeneProfileMatrix:
    genes = [
        GeneAnnotation("gene_a", "chr1", 0, 100, "+"),
        GeneAnnotation("gene_b", "chr1", 100, 200, "+"),
        GeneAnnotation("gene_c", "chr1", 200, 300, "-"),
        GeneAnnotation("gene_d", "chr1", 300, 400, "-"),
    ]
    feature_bins = [
        FeatureBin("up_1", "up", 0, 0),
        FeatureBin("body_1", "body", 0, 1),
        FeatureBin("body_2", "body", 1, 2),
        FeatureBin("down_1", "down", 0, 3),
    ]
    values = np.array(
        [
            [0.10, 0.12, 0.14, 0.16],
            [0.11, 0.13, 0.15, 0.17],
            [0.82, 0.78, 0.76, 0.72],
            [0.80, 0.79, 0.75, 0.71],
        ],
        dtype=float,
    )
    return GeneProfileMatrix(
        genes=genes,
        feature_bins=feature_bins,
        values=values,
        gene_missing_rate=np.zeros(4, dtype=float),
        gene_variance=np.var(values, axis=1),
        feature_missing_rate=np.zeros(4, dtype=float),
        feature_variance=np.var(values, axis=0),
    )


def _scaled_matrix(n_genes: int = 12) -> GeneProfileMatrix:
    genes = [
        GeneAnnotation(f"gene_{idx}", "chr1", idx * 100, (idx + 1) * 100, "+" if idx % 2 == 0 else "-")
        for idx in range(n_genes)
    ]
    feature_bins = [
        FeatureBin("up_1", "up", 0, 0),
        FeatureBin("body_1", "body", 0, 1),
        FeatureBin("body_2", "body", 1, 2),
        FeatureBin("down_1", "down", 0, 3),
    ]
    values = np.zeros((n_genes, 4), dtype=float)
    for idx in range(n_genes):
        cluster = idx % 3
        base = np.array(
            [
                [0.10, 0.12, 0.14, 0.16],
                [0.45, 0.50, 0.52, 0.56],
                [0.82, 0.78, 0.76, 0.72],
            ][cluster],
            dtype=float,
        )
        values[idx] = base + (0.005 * (idx // 3))
    return GeneProfileMatrix(
        genes=genes,
        feature_bins=feature_bins,
        values=values,
        gene_missing_rate=np.zeros(n_genes, dtype=float),
        gene_variance=np.var(values, axis=1),
        feature_missing_rate=np.zeros(4, dtype=float),
        feature_variance=np.var(values, axis=0),
    )


def test_run_pca_kmeans_is_deterministic() -> None:
    matrix = _matrix()
    config = BackendConfig(n_components=2, n_clusters=2, seed=7, n_init=4)
    left = run_pca_kmeans(matrix, config)
    right = run_pca_kmeans(matrix, config)

    assert np.array_equal(left["labels"], right["labels"])
    assert np.allclose(left["embedding"], right["embedding"])
    assert left["inertia"] == right["inertia"]


def test_run_pca_kmeans_reduces_effective_components_for_two_genes() -> None:
    genes = [
        GeneAnnotation("gene_a", "chr1", 0, 100, "+"),
        GeneAnnotation("gene_b", "chr1", 100, 200, "+"),
    ]
    feature_bins = [
        FeatureBin("body_1", "body", 0, 0),
        FeatureBin("body_2", "body", 1, 1),
        FeatureBin("body_3", "body", 2, 2),
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
    )

    result = run_pca_kmeans(matrix, BackendConfig(n_components=2, n_clusters=1, seed=3))

    assert result["embedding"].shape == (2, 1)
    assert result["components"].shape == (1, 3)
    assert result["metadata"]["effective_n_components"] == 1


def test_run_pca_kmeans_uses_truncated_pca_when_auto_threshold_is_exceeded() -> None:
    result = run_pca_kmeans(
        _scaled_matrix(200),
        BackendConfig(
            n_components=2,
            n_clusters=3,
            seed=7,
            n_init=4,
            pca_solver=PcaSolver.AUTO,
            pca_exact_max_matrix_size=500,
        ),
        silhouette_config=SilhouetteConfig(enabled=False),
    )

    assert result["embedding"].shape == (200, 2)
    assert result["components"].shape == (2, 4)
    assert result["metadata"]["pca_solver"] == "truncated"
    assert result["metadata"]["pca_solver_reason"] == "matrix_size_exceeds_exact_limit:500"


def test_run_pca_kmeans_uses_sampled_silhouette_in_auto_mode_above_threshold() -> None:
    result = run_pca_kmeans(
        _scaled_matrix(12),
        BackendConfig(n_components=2, n_clusters=3, seed=7, n_init=4),
        silhouette_config=SilhouetteConfig(
            mode=SilhouetteMode.AUTO,
            exact_max_genes=6,
            max_samples=9,
            seed=7,
        ),
    )

    assert result["metadata"]["silhouette_mode"] == "sampled"
    assert result["metadata"]["silhouette_n_samples"] == 9
    assert result["silhouette_score"] is not None


def test_run_pca_kmeans_can_disable_silhouette() -> None:
    result = run_pca_kmeans(
        _scaled_matrix(12),
        BackendConfig(n_components=2, n_clusters=3, seed=7, n_init=4),
        silhouette_config=SilhouetteConfig(enabled=False),
    )

    assert result["silhouette_score"] is None
    assert result["metadata"]["silhouette_mode"] == "disabled"
    assert result["metadata"]["silhouette_reason"] == "disabled_by_config"


def test_cluster_gene_profiles_returns_kmeans_and_hierarchical_labels() -> None:
    result = cluster_gene_profiles(
        _matrix(),
        ClusterConfig(
            bsx_path=Path("sample.bsx"),
            annotation_path=Path("genes.gff"),
            gene_profile=GeneProfileConfig(normalization=NormalizationMode.COLUMN_ZSCORE),
            backend=BackendConfig(n_components=2, n_clusters=2, seed=5, n_init=4),
            hierarchical=HierarchicalConfig(),
            cluster_source=ClusterSource.HIERARCHICAL,
            output=OutputConfig(output_dir=Path(".")),
        ),
    )

    assert result.embedding.shape == (4, 2)
    assert result.linkage_matrix is not None
    assert result.leaf_order is not None
    assert result.hierarchical_labels is not None
    assert result.cluster_source == "hierarchical"
    assert np.array_equal(result.labels, result.hierarchical_labels)


def test_cluster_gene_profiles_skips_hierarchical_above_threshold() -> None:
    result = cluster_gene_profiles(
        _matrix(),
        ClusterConfig(
            bsx_path=Path("sample.bsx"),
            annotation_path=Path("genes.gff"),
            backend=BackendConfig(n_components=2, n_clusters=2, seed=5, n_init=4),
            hierarchical=HierarchicalConfig(max_genes=3, subsample_genes=None),
            cluster_source=ClusterSource.KMEANS,
            output=OutputConfig(output_dir=Path(".")),
        ),
    )

    assert result.hierarchical_labels is None
    assert result.linkage_matrix is None
    assert result.leaf_order is None
    assert result.metadata["hierarchical"]["enabled"] is False
    assert result.metadata["hierarchical"]["skipped_reason"] == "n_genes_exceeds_max_genes:3"


def test_cluster_gene_profiles_runs_subsampled_hierarchical_diagnostics() -> None:
    result = cluster_gene_profiles(
        _scaled_matrix(12),
        ClusterConfig(
            bsx_path=Path("sample.bsx"),
            annotation_path=Path("genes.gff"),
            backend=BackendConfig(n_components=2, n_clusters=3, seed=5, n_init=4),
            hierarchical=HierarchicalConfig(
                mode=HierarchicalMode.AUTO,
                max_genes=6,
                subsample_genes=5,
            ),
            cluster_source=ClusterSource.KMEANS,
            output=OutputConfig(output_dir=Path(".")),
        ),
    )

    assert result.hierarchical_labels is None
    assert result.linkage_matrix is not None
    assert result.leaf_order is not None
    assert len(result.leaf_order) == 5
    assert result.metadata["hierarchical"]["enabled"] is True
    assert result.metadata["hierarchical"]["mode"] == "subsample"
    assert result.metadata["hierarchical"]["diagnostic_only"] is True
    assert result.metadata["hierarchical"]["effective_n_genes"] == 5
    assert result.metadata["hierarchical"]["subsampled"] is True
    assert len(result.metadata["hierarchical"]["leaf_gene_ids"]) == 5


def test_cluster_gene_profiles_rejects_hierarchical_source_when_policy_is_not_exact() -> None:
    with pytest.raises(ValueError, match="Hierarchical cluster_source requested"):
        cluster_gene_profiles(
            _scaled_matrix(12),
            ClusterConfig(
                bsx_path=Path("sample.bsx"),
                annotation_path=Path("genes.gff"),
                backend=BackendConfig(n_components=2, n_clusters=3, seed=5, n_init=4),
                hierarchical=HierarchicalConfig(
                    mode=HierarchicalMode.AUTO,
                    max_genes=6,
                    subsample_genes=5,
                ),
                cluster_source=ClusterSource.HIERARCHICAL,
                output=OutputConfig(output_dir=Path(".")),
            ),
        )
