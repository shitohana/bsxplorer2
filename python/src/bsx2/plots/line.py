from __future__ import annotations

from typing import Optional
import heapq

import holoviews as hv
import numpy as np

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    Segment,
    _apply_savgol_smoothing,
    _clip_profile,
    _coerce_smooth_config,
    segments_total_bins,
)
from ._html_common import (
    _bin_points_windows,
    _ensure_plotly,
    _hv_init,
    _segment_decor_rel,
)


def _line_profile(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    nan_fill: Optional[float] = None,
):
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    x_out = (np.arange(n_windows, dtype=float) + 0.5) / float(n_windows)

    streams = []
    for pos, dens in zip(drd.positions, drd.densities):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        if nan_fill is not None:
            y = np.where(np.isnan(y), nan_fill, y)
        pts = list(zip(x.tolist(), y.tolist()))
        pts.sort(key=lambda t: t[0])
        streams.append(pts)

    if not streams:
        return np.array([]), np.array([])

    merged = heapq.merge(*streams, key=lambda t: t[0])
    xs = []
    ys = []
    for x, y in merged:
        xs.append(x)
        ys.append(y)
    y_out = _bin_points_windows(
        np.asarray(xs),
        np.asarray(ys),
        n_windows=n_windows,
        agg="mean",
        nan_policy="keep",
    )
    return x_out, y_out


def line_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    smooth: dict | int | None = 50,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
) -> str:
    _hv_init()
    if segments is None:
        segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    if n_windows is None:
        n_windows = segments_total_bins(segments)

    x, y = _line_profile(
        drd,
        segments=segments,
        n_windows=n_windows,
        nan_fill=None,
    )

    if smooth is not None and y.size > 0:
        total_bins = segments_total_bins(segments)
        smooth_cfg = _coerce_smooth_config(smooth, total_bins=total_bins)
        if smooth_cfg is not None:
            y_scaled = y.astype(float, copy=True)
            y_scaled = _apply_savgol_smoothing(y_scaled, smooth_cfg, segments=segments)
            y_scaled = _clip_profile(y_scaled)
            y = y_scaled

    if x.size == 0 or y.size == 0:
        return _ensure_plotly(hv.render(hv.Curve([]), backend="plotly")).to_html(
            full_html=full_html,
            include_plotlyjs=include_js,
        )

    curve = hv.Curve((x, y), kdims="relative position", vdims="density").opts(
        xlabel="Metagene position (relative)",
        ylabel="Mean methylation density",
        show_legend=False,
    )
    fig = _ensure_plotly(hv.render(curve, backend="plotly"))
    _segment_decor_rel(fig, segments, annotate_tss_tes=True)
    fig.update_layout(
        height=600 if height is None else int(height),
        width=1000 if width is None else int(width),
        margin=dict(l=70, r=30, t=60, b=70),
        title=title or "Metagene profile - Line",
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


