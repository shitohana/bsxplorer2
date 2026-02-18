from __future__ import annotations

from typing import Optional

import numpy as np
import plotly.graph_objects as go

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import Segment, segment_boundaries, segments_total_bins
from ._html_common import (
    _bin_points_windows,
    _hv_init,
    _rank_compress,
    _segment_decor_bin,
)


def _heatmap_matrix(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    nan_fill: Optional[float] = None,
):
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40

    rows = []
    labels = []
    for pos, dens, lbl in zip(drd.positions, drd.densities, drd.labels):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        if nan_fill is not None:
            y = np.where(np.isnan(y), nan_fill, y)
        row = _bin_points_windows(
            x,
            y,
            n_windows=n_windows,
            agg="mean",
            nan_policy="keep",
        )
        rows.append(row)
        labels.append(lbl if lbl is not None else f"region_{len(labels)+1}")

    if not rows:
        return np.empty((0, 0)), [], []

    mat = np.vstack(rows)
    bins = (np.arange(n_windows, dtype=float) + 0.5) / float(n_windows)
    return mat, labels, bins


def heatmap_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: int | None = None,
    rank_rows: int = 100,
    rank_score: str = "mean",
    sort_order: str = "desc",
    colorscale: str = "Viridis",
    empty_bin_fill: float = 0.0,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
    vmax_q: float | None = 0.995,
) -> str:
    _hv_init()
    if rank_score not in {"mean", "body_mean"}:
        raise ValueError("rank_score must be 'mean' or 'body_mean'")
    if sort_order not in {"asc", "desc"}:
        raise ValueError("sort_order must be 'asc' or 'desc'")
    if segments is None:
        segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    if n_windows is None:
        n_windows = segments_total_bins(segments)

    z, _, bins = _heatmap_matrix(
        drd,
        segments=segments,
        n_windows=n_windows,
        nan_fill=None,
    )
    if z.size == 0:
        return go.Figure().to_html(full_html=full_html, include_plotlyjs=include_js)

    if rank_score == "body_mean":
        bins_arr = np.asarray(bins, dtype=float)
        bounds = segment_boundaries(segments)
        b0, b1 = bounds[0], bounds[1]
        mask = (bins_arr >= b0) & (bins_arr < b1)
    else:
        mask = np.ones(z.shape[1], dtype=bool)

    denom = float(np.sum(mask) + 1.0)
    scores = np.nansum(z[:, mask], axis=1) / denom
    order_idx = np.argsort(scores)
    if sort_order == "desc":
        order_idx = order_idx[::-1]
    z = z[order_idx]

    z = _rank_compress(z, rank_rows, fill=empty_bin_fill)
    regions = [str(i) for i in range(z.shape[0])]

    z_vis = np.asarray(z, dtype=float)
    zmin = 0.0
    zmax = 1.0
    if vmax_q is not None:
        vals = z_vis[np.isfinite(z_vis)]
        if vals.size:
            zmax = float(np.nanquantile(vals, vmax_q))
            if not np.isfinite(zmax) or zmax <= zmin:
                zmax = 1.0

    n_bins = int(z_vis.shape[1])
    fig_width = width if width is not None else max(900, min(1400, 2 * n_bins))
    fig_height = height if height is not None else max(550, min(900, 5 * int(z_vis.shape[0])))

    x_plot = np.arange(n_bins, dtype=int)
    fig = go.Figure(
        data=go.Heatmap(
            z=z_vis,
            x=x_plot,
            y=regions,
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            zsmooth=False,
            xgap=0,
            ygap=0,
            hoverongaps=False,
            colorbar=dict(title="Methylation density"),
        )
    )
    fig.update_yaxes(autorange="reversed", showticklabels=False)
    _segment_decor_bin(fig, segments, annotate_tss_tes=True)
    fig.update_layout(
        title=title or "Metagene profile - Heatmap (BSX1 ranked)",
        width=fig_width,
        height=fig_height,
        margin=dict(l=70, r=30, t=60, b=70),
        xaxis_title="Position (bin)",
        yaxis_title="Rank",
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


