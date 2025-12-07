from __future__ import annotations

from typing import Any, Iterable, Sequence

import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from bsx2.plots.cluster import MethylationMatrix, PCAResult, KMeansResult, LinkageResult, reorder_matrix


def pca_scatter(
    pca: PCAResult,
    *,
    labels: Sequence[str] | None = None,
    clusters: Sequence[int] | None = None,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    """Plotly scatter for PCA scores (first two components)."""
    scores = pca.scores
    if scores.shape[1] < 2:
        raise ValueError("PCA scores must have at least 2 components for scatter")
    xs, ys = scores[:, 0], scores[:, 1]
    df = {
        "PC1": xs,
        "PC2": ys,
        "label": labels if labels is not None else [f"r{i}" for i in range(len(xs))],
    }
    if clusters is not None:
        df["cluster"] = list(clusters)
    fig = px.scatter(df, x="PC1", y="PC2", color="cluster" if clusters is not None else None, hover_name="label")
    if pca.explained_variance is not None and len(pca.explained_variance) >= 2:
        total_var = np.sum(pca.explained_variance)
        if total_var > 0:
            pct = pca.explained_variance / total_var * 100.0
            fig.update_layout(xaxis_title=f"PC1 ({pct[0]:.1f}%)", yaxis_title=f"PC2 ({pct[1]:.1f}%)")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def kmeans_centroids_heatmap(
    km: KMeansResult,
    *,
    bins: Iterable[Any] | None = None,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    """Plotly heatmap of k-means centroids (clusters x bins)."""
    z = km.centroids
    x = list(bins) if bins is not None else list(range(z.shape[1]))
    y = [f"cluster_{i}" for i in range(km.n_clusters)]
    fig = go.Figure(data=go.Heatmap(z=z, x=x, y=y, colorscale="Viridis", colorbar=dict(title="centroid")))
    fig.update_layout(xaxis_title="bin", yaxis_title="cluster", yaxis_autorange="reversed")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def heatmap_ordered(
    mat: MethylationMatrix,
    *,
    order: Sequence[int] | None = None,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    """Plotly heatmap for methylation matrix; optionally reorder rows."""
    mat_use = reorder_matrix(mat, order) if order is not None else mat
    z = mat_use.matrix
    x = mat_use.bins.tolist()
    y = mat_use.region_ids
    fig = go.Figure(data=go.Heatmap(z=z, x=x, y=y, colorscale="Viridis", colorbar=dict(title="density")))
    fig.update_layout(xaxis_title="bin", yaxis_title="region", yaxis_autorange="reversed")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def dendrogram_plot(
    link: LinkageResult,
    *,
    labels: Sequence[str] | None = None,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    """Plotly dendrogram from a linkage result (requires scipy)."""
    try:
        import scipy.cluster.hierarchy as sch  # type: ignore
    except ImportError as e:  # pragma: no cover - optional dependency
        raise ImportError("scipy is required for dendrogram_plot") from e

    dendro = sch.dendrogram(link.linkage, labels=list(labels) if labels is not None else None, no_plot=True)
    icoord = np.array(dendro["icoord"])
    dcoord = np.array(dendro["dcoord"])
    xlabs = dendro.get("ivl", [])

    data = []
    for xs, ys in zip(icoord, dcoord):
        data.append(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color="black"), hoverinfo="none"))

    fig = go.Figure(data=data)
    fig.update_layout(
        xaxis=dict(ticktext=xlabs, tickvals=list(range(5, 10 * len(xlabs), 10)), showgrid=False, zeroline=False),
        yaxis=dict(title="distance"),
        showlegend=False,
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)
