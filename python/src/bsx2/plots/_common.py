from __future__ import annotations

from beartype.typing import Callable

import holoviews as hv
import numpy as np
import plotly.graph_objects as go

from bsx2 import AggMethod

from bsx2.guards import require_equal_length
from bsx2.plots.metagene import MetageneProfileSegment, segments_total_bins
from bsx2.validation import (
    validate_matrix_shape,
    validate_n_windows,
    validate_nan_policy,
    validate_segments,
    NanPolicy,
)


def _hv_init() -> None:
    hv.extension("plotly")


def _ensure_plotly(fig):
    if isinstance(fig, go.Figure):
        return fig
    return go.Figure(fig)


def _bin_points_windows_fast(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    *,
    n_windows: int,
    agg: AggMethod,
    nan_policy: NanPolicy,
) -> np.ndarray:
    x_vals = validate_matrix_shape(x_vals, 1, name="x_vals")
    y_vals = validate_matrix_shape(y_vals, 1, name="y_vals")
    require_equal_length(x_vals, y_vals, left_name="x_vals", right_name="y_vals")
    n_windows = validate_n_windows(n_windows)
    nan_policy = validate_nan_policy(nan_policy)

    x = np.asarray(x_vals, dtype=np.float64)
    y = np.asarray(y_vals, dtype=np.float64)

    # NaN/Inf handling once
    if nan_policy is NanPolicy.ZERO:
        y = np.where(np.isfinite(y), y, 0.0)
    else:  # KEEP / DROP -> same effective behavior as before
        m = np.isfinite(y)
        x = x[m]
        y = y[m]

    out = np.full(n_windows, np.nan, dtype=np.float64)
    if y.size == 0:
        return out

    # x is assumed in [0, 1]; x == 1 maps to last bin
    idx = (x * n_windows).astype(np.int64)
    idx[idx == n_windows] = n_windows - 1

    if agg is AggMethod.Mean:
        counts = np.bincount(idx, minlength=n_windows)
        sums = np.bincount(idx, weights=y, minlength=n_windows)
        nonempty = counts > 0
        out[nonempty] = sums[nonempty] / counts[nonempty]
        return out

    if agg is AggMethod.Min:
        tmp = np.full(n_windows, np.inf, dtype=np.float64)
        np.minimum.at(tmp, idx, y)
        tmp[tmp == np.inf] = np.nan   # empty bins stayed untouched
        return tmp

    if agg is AggMethod.Max:
        tmp = np.full(n_windows, -np.inf, dtype=np.float64)
        np.maximum.at(tmp, idx, y)
        tmp[tmp == -np.inf] = np.nan  # empty bins stayed untouched
        return tmp

    if agg is AggMethod.Median:
        # Median path relies on grouped equal-bin runs (caller sorts x for median when needed).
        cuts = np.flatnonzero(np.diff(idx)) + 1
        y_groups = np.split(y, cuts)
        bin_ids = idx[np.r_[0, cuts]]

        for b, g in zip(bin_ids, y_groups):
            out[int(b)] = float(np.median(g))
        return out

    raise ValueError(f"unsupported agg: {agg}")


def _rank_compress(z_sorted: np.ndarray, rank_rows: int, *, fill: float | None = 0.0) -> np.ndarray:
    if z_sorted.ndim != 2:
        return z_sorted

    n_rows, n_bins = z_sorted.shape
    if n_rows == 0:
        return np.empty((0, n_bins), dtype=float)

    rows = max(int(rank_rows), 1)

    sums = np.zeros((rows, n_bins), dtype=float)
    cnts = np.zeros((rows, n_bins), dtype=np.int32)

    k = np.arange(rows + 1, dtype=np.int64)
    bounds = (k * n_rows + rows - 1) // rows  # ceil(k * n_rows / rows)

    for b in range(rows):
        s = int(bounds[b])
        e = int(bounds[b + 1])
        if s >= e:  
            continue

        block = z_sorted[s:e] 
        finite = np.isfinite(block)
        if not finite.any():
            continue

        sums[b] = np.where(finite, block, 0.0).sum(axis=0)
        cnts[b] = finite.sum(axis=0, dtype=np.int32)

    if fill is None:
        out = np.full((rows, n_bins), np.nan, dtype=float)
    else:
        out = np.full((rows, n_bins), float(fill), dtype=float)

    np.divide(sums, cnts, out=out, where=(cnts > 0))
    return out


def _segment_decor(
    fig,
    segments: list[MetageneProfileSegment] | None,
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


def _segment_decor_rel(fig, segments: list[MetageneProfileSegment] | None, *, annotate_tss_tes: bool = False) -> None:
    if not segments:
        return
    total = float(segments_total_bins(segments))
    _segment_decor(fig, segments, annotate_tss_tes=annotate_tss_tes, scale=lambda v: v / total)


def _segment_decor_bin(fig, segments: list[MetageneProfileSegment] | None, *, annotate_tss_tes: bool = False) -> None:
    _segment_decor(fig, segments, annotate_tss_tes=annotate_tss_tes, scale=float)
