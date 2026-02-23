from __future__ import annotations

from typing import Callable

import holoviews as hv
import numpy as np
import plotly.graph_objects as go

from bsx2.guards import require_equal_length
from bsx2.plots.metagene import Segment, segments_total_bins
from bsx2.validation import (
    validate_matrix_shape,
    validate_n_windows,
    validate_nan_policy,
    validate_segments,
    validate_window_agg,
    NanPolicy,
)


def _hv_init() -> None:
    hv.extension("plotly")


def _ensure_plotly(fig):
    if isinstance(fig, go.Figure):
        return fig
    return go.Figure(fig)


def _bin_points_windows(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    *,
    n_windows: int,
    agg: str,
    nan_policy: NanPolicy,
) -> np.ndarray:
    x_vals = validate_matrix_shape(x_vals, 1, name="x_vals")
    y_vals = validate_matrix_shape(y_vals, 1, name="y_vals")
    require_equal_length(x_vals, y_vals, left_name="x_vals", right_name="y_vals")
    n_windows = validate_n_windows(n_windows)
    validate_nan_policy(nan_policy)
    agg = validate_window_agg(agg)

    bins_values = [[] for _ in range(n_windows)]
    for x, y in zip(x_vals, y_vals):
        if not np.isfinite(x):
            continue
        if x < 0.0 or x > 1.0:
            continue
        if nan_policy is NanPolicy.ZERO and not np.isfinite(y):
            y = 0.0
        if nan_policy is NanPolicy.DROP and not np.isfinite(y):
            continue
        idx = int(x * n_windows)
        if idx == n_windows:
            idx = n_windows - 1
        if 0 <= idx < n_windows:
            bins_values[idx].append(y)

    out = []
    for vals in bins_values:
        if not vals:
            out.append(np.nan)
            continue
        arr = np.asarray(vals, dtype=float)
        if agg == "mean":
            if nan_policy == "drop":
                out.append(float(np.nanmean(arr)))
            elif nan_policy == "zero":
                out.append(float(np.mean(arr)))
            else:
                finite = np.isfinite(arr)
                out.append(float(np.mean(arr[finite])) if finite.any() else np.nan)
        elif agg == "median":
            if nan_policy == "drop":
                out.append(float(np.nanmedian(arr)))
            elif nan_policy == "zero":
                out.append(float(np.median(arr)))
            else:
                finite = np.isfinite(arr)
                out.append(float(np.median(arr[finite])) if finite.any() else np.nan)
        elif agg == "max":
            if nan_policy == "drop":
                out.append(float(np.nanmax(arr)))
            elif nan_policy == "zero":
                out.append(float(np.max(arr)))
            else:
                finite = np.isfinite(arr)
                out.append(float(np.max(arr[finite])) if finite.any() else np.nan)
        elif agg == "min":
            if nan_policy == "drop":
                out.append(float(np.nanmin(arr)))
            elif nan_policy == "zero":
                out.append(float(np.min(arr)))
            else:
                finite = np.isfinite(arr)
                out.append(float(np.min(arr[finite])) if finite.any() else np.nan)
    return np.asarray(out, dtype=float)


def _rank_compress(z_sorted: np.ndarray, rank_rows: int, *, fill: float | None = 0.0) -> np.ndarray:
    if z_sorted.ndim != 2:
        return z_sorted
    n_rows, n_bins = z_sorted.shape
    if n_rows == 0:
        return np.empty((0, n_bins), dtype=float)
    rows = max(int(rank_rows), 1)
    sums = np.zeros((rows, n_bins), dtype=float)
    cnts = np.zeros((rows, n_bins), dtype=np.int32)
    for i in range(n_rows):
        ridx = int(i * rows / n_rows)
        row = z_sorted[i]
        finite = np.isfinite(row)
        if np.any(finite):
            sums[ridx, finite] += row[finite]
            cnts[ridx, finite] += 1
    if fill is None:
        out = np.full((rows, n_bins), np.nan, dtype=float)
    else:
        out = np.full((rows, n_bins), float(fill), dtype=float)
    np.divide(sums, cnts, out=out, where=(cnts > 0))
    return out


def _segment_decor(
    fig,
    segments: list[Segment] | None,
    *,
    annotate_tss_tes: bool = False,
    scale: Callable[[float], float],
) -> None:
    if not segments:
        return
    validate_segments(segments)
    boundaries = []
    centers = []
    labels = []
    cum = 0
    for seg in segments:
        start = cum
        end = cum + seg.n_bins
        mid = (start + end) / 2
        boundaries.append(end)
        centers.append(mid)
        labels.append(seg.name)
        cum = end
    shapes = []
    for b in boundaries[:-1]:
        xb = scale(b)
        shapes.append(
            dict(type="line", x0=xb, x1=xb, y0=0, y1=1, yref="paper", line=dict(dash="dash", width=1, color="gray"))
        )
    tickvals = [scale(c) for c in centers]
    fig.update_xaxes(tickmode="array", tickvals=tickvals, ticktext=labels)
    fig.update_layout(shapes=shapes)
    if annotate_tss_tes and len(boundaries) >= 2:
        first = scale(boundaries[0])
        last = scale(boundaries[-2])
        fig.add_annotation(x=first, y=1.02, xref="x", yref="paper", text="TSS", showarrow=False, font=dict(size=10))
        fig.add_annotation(x=last, y=1.02, xref="x", yref="paper", text="TES", showarrow=False, font=dict(size=10))


def _segment_decor_rel(fig, segments: list[Segment] | None, *, annotate_tss_tes: bool = False) -> None:
    if not segments:
        return
    total = float(segments_total_bins(segments))
    _segment_decor(fig, segments, annotate_tss_tes=annotate_tss_tes, scale=lambda v: v / total)


def _segment_decor_bin(fig, segments: list[Segment] | None, *, annotate_tss_tes: bool = False) -> None:
    _segment_decor(fig, segments, annotate_tss_tes=annotate_tss_tes, scale=float)
