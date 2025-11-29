import numpy as np
import pytest

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.cluster import prepare_matrix, run_pca, run_kmeans, run_linkage
from bsx2.plots.cluster_vis import (
    pca_scatter,
    kmeans_centroids_heatmap,
    heatmap_ordered,
    dendrogram_plot,
)


def _drd(n_regions: int = 6, n_bins: int = 5) -> DiscreteRegionData:
    drd = DiscreteRegionData()
    rng = np.random.default_rng(0)
    for i in range(n_regions):
        x = np.linspace(0, 1, n_bins + i % 2)
        y = rng.random(len(x))
        drd.insert(x, y, f"r{i}")
    return drd


def test_pca_scatter_html():
    drd = _drd()
    mat = prepare_matrix(drd, norm="zscore")
    pca = run_pca(mat, n_components=2)
    html = pca_scatter(pca, labels=mat.region_ids, clusters=[0] * len(mat.region_ids), full_html=False)
    assert "<html" not in html.lower()
    assert "plotly" in html.lower()


def test_kmeans_heatmap_html():
    drd = _drd()
    mat = prepare_matrix(drd)
    km = run_kmeans(mat, n_clusters=2, random_state=0)
    html = kmeans_centroids_heatmap(km, bins=mat.bins, full_html=False)
    assert "heatmap" in html.lower()


def test_heatmap_ordered_html():
    drd = _drd()
    mat = prepare_matrix(drd)
    html = heatmap_ordered(mat, order=list(reversed(range(len(mat.region_ids)))), full_html=False)
    assert "heatmap" in html.lower()


def test_dendrogram_html():
    drd = _drd()
    mat = prepare_matrix(drd)
    try:
        link = run_linkage(mat)
    except ImportError:
        pytest.skip("scipy not installed")
    html = dendrogram_plot(link, labels=mat.region_ids, full_html=False)
    assert "plotly" in html.lower()
