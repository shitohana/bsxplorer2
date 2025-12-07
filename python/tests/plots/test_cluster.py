import numpy as np
import pytest

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.cluster import (
    prepare_matrix,
    run_pca,
    run_kmeans,
    run_linkage,
    reorder_matrix,
    cluster_subset,
    metagene_for_cluster,
)


def _make_drd(n_regions: int = 5, n_bins: int = 8) -> DiscreteRegionData:
    drd = DiscreteRegionData()
    rng = np.random.default_rng(0)
    for i in range(n_regions):
        x = np.linspace(0, 1, n_bins + i % 3)  # разные длины
        y = rng.random(len(x))
        drd.insert(x, y, f"r{i}")
    return drd


def test_prepare_matrix_norms():
    drd = _make_drd()
    mat_none = prepare_matrix(drd, norm="none")
    assert mat_none.matrix.shape[0] == len(drd)
    mat_z = prepare_matrix(drd, norm="zscore")
    row = mat_z.matrix[0]
    assert abs(row.mean()) < 1e-9
    mat_mm = prepare_matrix(drd, norm="minmax")
    assert mat_mm.matrix.min() >= 0 and mat_mm.matrix.max() <= 1


def test_run_pca_shapes():
    drd = _make_drd()
    mat = prepare_matrix(drd, norm="zscore")
    res = run_pca(mat, n_components=3)
    assert res.scores.shape == (len(drd), 3)
    assert res.loadings.shape[1] == 3
    assert res.explained_variance.shape[0] == 3


def test_run_kmeans_deterministic():
    drd = _make_drd()
    mat = prepare_matrix(drd)
    res1 = run_kmeans(mat, n_clusters=2, random_state=42)
    res2 = run_kmeans(mat, n_clusters=2, random_state=42)
    assert np.array_equal(res1.labels, res2.labels)
    assert res1.centroids.shape[0] == 2


def test_run_linkage_optional():
    drd = _make_drd()
    mat = prepare_matrix(drd)
    try:
        res = run_linkage(mat)
    except ImportError:
        pytest.skip("scipy not installed")
    assert res.linkage.shape[0] == len(drd) - 1
    assert res.order.shape[0] == len(drd)


def test_reorder_and_subset():
    drd = _make_drd()
    mat = prepare_matrix(drd)
    order = list(reversed(range(len(drd))))
    reordered = reorder_matrix(mat, order)
    assert reordered.region_ids[0] == mat.region_ids[-1]
    labels = np.array([0, 1, 0, 1, 0])
    subset = cluster_subset(mat, labels, cluster_id=1)
    assert subset.matrix.shape[0] == np.sum(labels == 1)


def test_metagene_for_cluster_smoke():
    # synthetic contigs as simple tuples; RegionReader mocked via duck-typing
    class DummyBatch:
        def __init__(self, xs, ys):
            self._xs = xs
            self._ys = ys

        def discretise(self, n, agg):
            return self._xs, self._ys

    class DummyReader:
        def iter_contigs(self, contigs):
            for _ in contigs:
                yield DummyBatch([0.0, 0.5, 1.0], [0.1, 0.2, 0.3])

    from bsx2.plots.metagene import Segment

    contigs = ["c1", "c2", "c3"]
    labels = np.array([0, 1, 0])
    drd = metagene_for_cluster(
        DummyReader(),
        contigs,
        labels,
        cluster_id=0,
        segments=[Segment("region", 3)],
        agg_method=None,
    )
    assert len(drd.positions) == 2
