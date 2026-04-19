from __future__ import annotations

import numpy as np
import pytest

bsx2 = pytest.importorskip("bsx2")
if not hasattr(bsx2, "Context"):
    pytest.skip("bsx2 extension is unavailable", allow_module_level=True)

viz = pytest.importorskip("bsx2.viz")
clustering_models = pytest.importorskip("bsx2.clustering.models")

FeatureBin = clustering_models.FeatureBin
GeneAnnotation = clustering_models.GeneAnnotation
GeneClusterResult = clustering_models.GeneClusterResult
GeneProfileMatrix = clustering_models.GeneProfileMatrix

build_cluster_metagene_data = viz.build_cluster_metagene_data
cluster_metagene_plot = viz.cluster_metagene_plot


def _is_holoviews_object(obj) -> bool:
    return obj.__class__.__module__.startswith("holoviews")


def _result() -> GeneClusterResult:
    genes = [
        GeneAnnotation("gene_a", "chr1", 0, 100, "+", "GeneA"),
        GeneAnnotation("gene_b", "chr1", 100, 200, "+", "GeneB"),
        GeneAnnotation("gene_c", "chr1", 200, 300, "-", "GeneC"),
        GeneAnnotation("gene_d", "chr1", 300, 400, "-", "GeneD"),
    ]
    feature_bins = [
        FeatureBin("up_1", "up", 0, 0),
        FeatureBin("body_1", "body", 0, 1),
        FeatureBin("body_2", "body", 1, 2),
        FeatureBin("down_1", "down", 0, 3),
    ]
    values = np.array(
        [
            [0.10, 0.20, 0.30, 0.40],
            [0.20, 0.30, 0.40, 0.50],
            [0.70, 0.60, 0.50, 0.40],
            [0.90, 0.80, 0.70, 0.60],
        ],
        dtype=float,
    )
    matrix = GeneProfileMatrix(
        genes=genes,
        feature_bins=feature_bins,
        values=values,
        gene_missing_rate=np.zeros(4, dtype=float),
        gene_variance=np.var(values, axis=1),
        feature_missing_rate=np.zeros(4, dtype=float),
        feature_variance=np.var(values, axis=0),
        metadata={"timings_s": {"total": 0.01}},
    )
    return GeneClusterResult(
        feature_matrix=matrix,
        labels=np.array([0, 0, 1, 1], dtype=np.int64),
        cluster_source="kmeans",
        kmeans_labels=np.array([0, 0, 1, 1], dtype=np.int64),
        hierarchical_labels=np.array([0, 0, 1, 1], dtype=np.int64),
        embedding=np.array(
            [[1.0, 0.1], [0.9, 0.2], [-0.9, -0.2], [-1.0, -0.1]],
            dtype=float,
        ),
        components=np.array([[0.4, 0.3, 0.2, 0.1], [0.1, 0.2, 0.3, 0.4]], dtype=float),
        centroids=np.array([[0.95, 0.15], [-0.95, -0.15]], dtype=float),
        explained_variance_ratio=np.array([0.8, 0.2], dtype=float),
        inertia=0.4,
        silhouette_score=0.7,
        linkage_matrix=np.array(
            [[0.0, 1.0, 0.1, 2.0], [2.0, 3.0, 0.2, 2.0], [4.0, 5.0, 0.9, 4.0]],
            dtype=float,
        ),
        leaf_order=np.array([0, 1, 2, 3], dtype=np.int64),
        metadata={"timings_s": {"total": 0.02}, "pca_kmeans": {}, "hierarchical": {}},
    )


def test_build_cluster_metagene_data_returns_group_metadata() -> None:
    data = build_cluster_metagene_data(_result())

    assert len(data.groups) == 2
    assert [group.cluster_id for group in data.groups] == [0, 1]
    assert [group.gene_ids for group in data.groups] == [
        ["gene_a", "gene_b"],
        ["gene_c", "gene_d"],
    ]
    assert data.profiles.labels == ["kmeans cluster 0", "kmeans cluster 1"]
    np.testing.assert_allclose(data.profiles.densities[0], np.array([0.15, 0.25, 0.35, 0.45]))
    np.testing.assert_allclose(data.profiles.densities[1], np.array([0.80, 0.70, 0.60, 0.50]))


def test_build_cluster_metagene_data_filters_clusters_and_genes() -> None:
    data = build_cluster_metagene_data(
        _result(),
        cluster_ids=[1],
        gene_ids=["gene_c"],
    )

    assert len(data.groups) == 1
    assert data.groups[0].cluster_id == 1
    assert data.groups[0].gene_ids == ["gene_c"]
    np.testing.assert_allclose(data.profiles.densities[0], np.array([0.70, 0.60, 0.50, 0.40]))


def test_build_cluster_metagene_data_can_collapse_cross_cluster_selection() -> None:
    data = build_cluster_metagene_data(
        _result(),
        gene_ids=["gene_a", "gene_c"],
        collapse=True,
        label="focus set",
    )

    assert len(data.groups) == 1
    assert data.groups[0].cluster_id is None
    assert data.groups[0].gene_ids == ["gene_a", "gene_c"]
    assert data.profiles.labels == ["focus set"]
    np.testing.assert_allclose(data.profiles.densities[0], np.array([0.40, 0.40, 0.40, 0.40]))


def test_cluster_metagene_plot_returns_holoviews_object() -> None:
    plot = cluster_metagene_plot(
        _result(),
        cluster_ids=[0],
    )

    assert _is_holoviews_object(plot)


def test_build_cluster_metagene_data_rejects_unknown_selection() -> None:
    with pytest.raises(ValueError, match="Unknown cluster_ids"):
        build_cluster_metagene_data(_result(), cluster_ids=[99])

    with pytest.raises(ValueError, match="Unknown gene_ids"):
        build_cluster_metagene_data(_result(), gene_ids=["missing_gene"])

    with pytest.raises(ValueError, match="produced no genes"):
        build_cluster_metagene_data(
            _result(),
            cluster_ids=[0],
            gene_ids=["gene_c"],
        )
