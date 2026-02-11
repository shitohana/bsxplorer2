from __future__ import annotations
from typing import Tuple, Sequence, List, Optional
import heapq

import numpy as np
import holoviews as hv
import plotly.graph_objects as go

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    Segment,
    segment_boundaries,
    segments_total_bins,
    compute_discrete_regions,
    collect_contigs_from_hcannot,
    collect_parts_from_hcannot,
    combine_parts_drd,
    _coerce_smooth_config,
    _apply_savgol_smoothing,
    _clip_profile,
)


def _hv_init():
    # Гарантируем наличие plotly backend, даже если ранее загружался другой
    hv.extension("plotly")


def _ensure_plotly(fig):
    """hv.render может вернуть dict; приводим к plotly Figure."""
    if isinstance(fig, go.Figure):
        return fig
    return go.Figure(fig)


def discrete_to_long(drd: DiscreteRegionData, *, as_percent: bool = False, nan_fill: Optional[float] = None) -> List[tuple]:
    """Возвращает список записей (region, bin, x, density) без pandas/polars."""
    rows: List[tuple] = []
    for ridx, (pos, dens, lbl) in enumerate(zip(drd.positions, drd.densities, drd.labels)):
        region = lbl if lbl else f"region_{ridx+1}"
        vals = np.asarray(dens, dtype=float)
        if as_percent:
            vals = vals * 100.0
        vals[~np.isfinite(vals)] = np.nan
        if nan_fill is not None:
            vals = np.where(np.isnan(vals), nan_fill, vals)
        bins = np.arange(len(vals), dtype=int)
        for b, x, v in zip(bins, pos, vals):
            rows.append((region, int(b), float(x), float(v)))
    return rows

# Backwards compatibility: old name used in tests
discrete_to_long_pl = discrete_to_long_pd = lambda drd, as_percent=False, nan_fill=None: np.array(discrete_to_long(drd, as_percent=as_percent, nan_fill=nan_fill), dtype=object)


def _bin_points_windows(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    *,
    weights: Optional[np.ndarray] = None,
    n_windows: int,
    agg: str,
    nan_policy: str,
    value_mode: str = "density",
) -> np.ndarray:
    if n_windows <= 0:
        raise ValueError("n_windows must be > 0")
    if nan_policy not in {"drop", "zero", "keep"}:
        raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
    if value_mode not in {"density", "weighted"}:
        raise ValueError("value_mode must be 'density' or 'weighted'")
    if value_mode == "weighted" and weights is None:
        raise ValueError("value_mode='weighted' requires weights")

    bins_values = [[] for _ in range(n_windows)]
    bins_weighted_sum = np.zeros(n_windows, dtype=float)
    bins_weighted_total = np.zeros(n_windows, dtype=float)
    for i, (x, y) in enumerate(zip(x_vals, y_vals)):
        if not np.isfinite(x):
            continue
        if x < 0.0 or x > 1.0:
            continue
        w = None
        if weights is not None:
            if i >= len(weights):
                continue
            w = float(weights[i])
            if not np.isfinite(w) or w <= 0:
                continue
        if nan_policy == "zero" and not np.isfinite(y):
            y = 0.0
        if nan_policy == "drop" and not np.isfinite(y):
            continue
        idx = int(x * n_windows)
        if idx == n_windows:
            idx = n_windows - 1
        if 0 <= idx < n_windows:
            if value_mode == "weighted":
                if not np.isfinite(y):
                    if nan_policy == "zero":
                        y = 0.0
                    else:
                        continue
                w_use = 1.0 if w is None else w
                bins_weighted_sum[idx] += float(y) * w_use
                bins_weighted_total[idx] += w_use
            else:
                bins_values[idx].append(y)

    if value_mode == "weighted":
        out = np.full(n_windows, np.nan, dtype=float)
        mask = bins_weighted_total > 0
        out[mask] = bins_weighted_sum[mask] / bins_weighted_total[mask]
        return out

    out = []
    for vals in bins_values:
        if not vals:
            out.append(np.nan)
            continue
        arr = np.asarray(vals, dtype=float)
        if agg == "mean":
            fn = np.mean if nan_policy in {"keep", "zero"} else np.nanmean
            out.append(float(fn(arr)))
        elif agg == "median":
            fn = np.median if nan_policy in {"keep", "zero"} else np.nanmedian
            out.append(float(fn(arr)))
        elif agg == "max":
            out.append(float(np.nanmax(arr)) if nan_policy == "drop" else float(np.max(arr)))
        elif agg == "min":
            out.append(float(np.nanmin(arr)) if nan_policy == "drop" else float(np.min(arr)))
        else:
            raise ValueError(f"Unsupported agg: {agg}")
    return np.asarray(out, dtype=float)


def line_df(
    drd: DiscreteRegionData,
    agg: str = "median",
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    as_percent: bool = False,
    order: Sequence[str] | None = None,
    nan_fill: Optional[float] = None,
    agg_scope: str = "points",
    nan_policy: str = "drop",
    x_mode: str = "relative",
    value_mode: str = "density",
    max_nan_frac: Optional[float] = None,
):
    """Return x, y arrays for line plot (nan-aware)."""
    if agg_scope not in {"points", "genes"}:
        raise ValueError("agg_scope must be 'points' or 'genes'")
    if nan_policy not in {"drop", "zero", "keep"}:
        raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
    if x_mode not in {"relative", "absolute"}:
        raise ValueError("x_mode must be 'relative' or 'absolute'")

    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    x_out = (np.arange(n_windows, dtype=float) + 0.5) / float(n_windows)

    def _to_relative(x: np.ndarray) -> tuple[np.ndarray, float | None]:
        if x_mode == "relative":
            return x, None
        if x.size == 0:
            return x, None
        x_min = float(np.nanmin(x))
        x_max = float(np.nanmax(x))
        length = x_max - x_min
        if not np.isfinite(length) or length <= 0:
            return np.array([], dtype=float), None
        return (x - x_min) / length, length

    if agg_scope == "points":
        streams = []
        lengths = []
        for pos, dens, w in zip(drd.positions, drd.densities, drd.weights or [None] * len(drd.positions)):
            x = np.asarray(pos, dtype=float)
            y = np.asarray(dens, dtype=float)
            weights = None
            if w is not None:
                weights = np.asarray(w, dtype=float)
            if as_percent:
                y = y * 100.0
            if nan_fill is not None:
                y = np.where(np.isnan(y), nan_fill, y)
            x_rel, length = _to_relative(x)
            if length is not None:
                lengths.append(length)
            x = x_rel
            if weights is None:
                pts = list(zip(x.tolist(), y.tolist(), [None] * len(x)))
            else:
                pts = list(zip(x.tolist(), y.tolist(), weights.tolist()))
            pts.sort(key=lambda t: t[0])
            streams.append(pts)

        if not streams:
            return np.array([]), np.array([])

        merged = heapq.merge(*streams, key=lambda t: t[0])
        xs = []
        ys = []
        ws = []
        for x, y, w in merged:
            xs.append(x)
            ys.append(y)
            ws.append(w)
        w_arr = None
        if any(w is not None for w in ws):
            w_arr = np.asarray([0.0 if w is None else w for w in ws], dtype=float)
        y_out = _bin_points_windows(
            np.asarray(xs),
            np.asarray(ys),
            weights=w_arr,
            n_windows=n_windows,
            agg=agg,
            nan_policy=nan_policy,
            value_mode=value_mode,
        )
        if x_mode == "absolute":
            if lengths:
                scale = float(np.median(np.asarray(lengths)))
                x_out = x_out * scale
        return x_out, y_out

    rows = []
    lengths = []
    for pos, dens, w in zip(drd.positions, drd.densities, drd.weights or [None] * len(drd.positions)):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        weights = None
        if w is not None:
            weights = np.asarray(w, dtype=float)
        if as_percent:
            y = y * 100.0
        if nan_fill is not None:
            y = np.where(np.isnan(y), nan_fill, y)
        x_rel, length = _to_relative(x)
        if length is not None:
            lengths.append(length)
        rows.append(
            _bin_points_windows(
                x_rel,
                y,
                weights=weights,
                n_windows=n_windows,
                agg=agg,
                nan_policy=nan_policy,
                value_mode=value_mode,
            )
        )
    if not rows:
        return np.array([]), np.array([])
    mat = np.vstack(rows)
    if max_nan_frac is not None and mat.size > 0:
        keep = np.mean(~np.isfinite(mat), axis=1) <= max_nan_frac
        mat = mat[keep]
    if agg == "mean":
        y_out = np.nanmean(mat, axis=0)
    elif agg == "median":
        y_out = np.nanmedian(mat, axis=0)
    elif agg == "max":
        y_out = np.nanmax(mat, axis=0)
    elif agg == "min":
        y_out = np.nanmin(mat, axis=0)
    else:
        raise ValueError("agg must be one of: mean, median, max, min")
    if x_mode == "absolute":
        if lengths:
            scale = float(np.median(np.asarray(lengths)))
            x_out = x_out * scale
    return x_out, y_out



def heatmap_df(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    as_percent: bool = False,
    order: Sequence[str] | None = None,
    nan_fill: Optional[float] = None,
    agg: str = "median",
    nan_policy: str = "drop",
    x_mode: str = "relative",
    value_mode: str = "density",
    max_nan_frac: Optional[float] = None,
):
    """Return heatmap matrix (regions x bins)."""
    if nan_policy not in {"drop", "zero", "keep"}:
        raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
    if x_mode not in {"relative", "absolute"}:
        raise ValueError("x_mode must be 'relative' or 'absolute'")

    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40

    def _to_relative(x: np.ndarray) -> tuple[np.ndarray, float | None]:
        if x_mode == "relative":
            return x, None
        if x.size == 0:
            return x, None
        x_min = float(np.nanmin(x))
        x_max = float(np.nanmax(x))
        length = x_max - x_min
        if not np.isfinite(length) or length <= 0:
            return np.array([], dtype=float), None
        return (x - x_min) / length, length

    rows = []
    labels = []
    lengths = []
    for pos, dens, w, lbl in zip(
        drd.positions,
        drd.densities,
        drd.weights or [None] * len(drd.positions),
        drd.labels,
    ):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        weights = None
        if w is not None:
            weights = np.asarray(w, dtype=float)
        if as_percent:
            y = y * 100.0
        if nan_fill is not None:
            y = np.where(np.isnan(y), nan_fill, y)
        x_rel, length = _to_relative(x)
        if length is not None:
            lengths.append(length)
        row = _bin_points_windows(
            x_rel,
            y,
            weights=weights,
            n_windows=n_windows,
            agg=agg,
            nan_policy=nan_policy,
            value_mode=value_mode,
        )
        rows.append(row)
        labels.append(lbl if lbl is not None else f"region_{len(labels)+1}")

    if not rows:
        return np.empty((0, 0)), [], []

    mat = np.vstack(rows)
    if max_nan_frac is not None and mat.size > 0:
        keep = np.mean(~np.isfinite(mat), axis=1) <= max_nan_frac
        mat = mat[keep]
        labels = [l for l, k in zip(labels, keep) if k]

    if x_mode == "absolute" and lengths:
        scale = float(np.median(np.asarray(lengths)))
        bins = (np.arange(n_windows, dtype=float) + 0.5) / float(n_windows)
        bins = bins * scale
    else:
        bins = list(range(n_windows))
    return mat, labels, bins



def dist_df(
    drd: DiscreteRegionData,
    *,
    as_percent: bool = False,
    nan_fill: Optional[float] = None,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    nan_policy: str = "drop",
):
    """Collect samples for box/violin.
    Each point = one region in one window.
    """
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    rows = []
    for pos, dens, w, lbl in zip(
        drd.positions,
        drd.densities,
        drd.weights or [None] * len(drd.positions),
        drd.labels,
    ):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        weights = None
        if w is not None:
            weights = np.asarray(w, dtype=float)
        if as_percent:
            y = y * 100.0
        if nan_fill is not None:
            y = np.where(np.isnan(y), nan_fill, y)
        binned = _bin_points_windows(
            x,
            y,
            weights=weights,
            n_windows=n_windows,
            agg="mean",
            nan_policy=nan_policy,
            value_mode="density",
        )
        for b, v in enumerate(binned):
            if np.isfinite(v):
                rows.append((b, float(v), lbl if lbl is not None else f"region_{len(rows)+1}"))
    return rows
def _segment_decor(
    fig,
    segments: list[Segment] | None,
    *,
    annotate_tss_tes: bool = False,
    x_mode: str = "rel",
    x_values: Optional[Sequence[float]] = None,
) -> None:
    """Применить границы сегментов/тики к plotly Figure."""
    if not segments:
        return
    total = segments_total_bins(segments)
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
    if x_mode == "rel":
        scale = lambda v: v / total
    elif x_mode == "bin":
        scale = float
    elif x_mode == "abs":
        vals = np.asarray(x_values if x_values is not None else [], dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size >= 2:
            vals.sort()
            diffs = np.diff(vals)
            step = float(np.nanmedian(diffs)) if diffs.size else 0.0
            if np.isfinite(step) and step > 0 and total > 0:
                total_span = (vals[-1] - vals[0]) + step
                factor = total_span / float(total)
                scale = lambda v: v * factor
            else:
                scale = lambda v: v / total
        else:
            scale = lambda v: v / total
    else:
        raise ValueError("x_mode must be 'rel', 'bin', or 'abs'")
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


def line_html(
    drd: DiscreteRegionData,
    agg: str = "median",
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    agg_scope: str = "points",
    nan_policy: str = "drop",
    x_mode: str = "relative",
    value_mode: str = "density",
    max_nan_frac: Optional[float] = None,
    drop_nan_rows: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
    smooth: dict | int | None = None,
    connectgaps: bool = False,
) -> str:
    _hv_init()
    x, y = line_df(
        drd,
        agg=agg,
        segments=segments,
        n_windows=n_windows,
        as_percent=as_percent,
        order=order,
        nan_fill=nan_fill,
        agg_scope=agg_scope,
        nan_policy=nan_policy,
        x_mode=x_mode,
        value_mode=value_mode,
        max_nan_frac=max_nan_frac,
    )
    if smooth is not None and y.size > 0:
        total_bins = segments_total_bins(segments) if segments else int(y.size)
        smooth_cfg = _coerce_smooth_config(smooth, total_bins=total_bins)
        if smooth_cfg is not None:
            y_scaled = y.astype(float, copy=True)
            if as_percent:
                y_scaled = y_scaled / 100.0
            y_scaled = _apply_savgol_smoothing(y_scaled, smooth_cfg, segments=segments)
            y_scaled = _clip_profile(y_scaled)
            y = y_scaled * (100.0 if as_percent else 1.0)
    if x.size == 0 or y.size == 0:
        return _ensure_plotly(hv.render(hv.Curve([]), backend="plotly")).to_html(
            full_html=full_html,
            include_plotlyjs=include_js,
        )
    x_label = "Metagene position (relative)" if x_mode == "relative" else "Metagene position (bp)"
    curve = hv.Curve((x, y), kdims="relative position", vdims="density").opts(
        xlabel=x_label,
        ylabel=f"{agg} density" + (" (%)" if as_percent else ""),
        show_legend=False,
    )
    fig = _ensure_plotly(hv.render(curve, backend="plotly"))
    if connectgaps:
        fig.update_traces(connectgaps=True)
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    decor_mode = "rel" if x_mode == "relative" else "abs"
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode=decor_mode, x_values=x)
    fig.update_layout(
        height=600,
        width=1000,
        margin=dict(l=70, r=30, t=60, b=70),
        title=title or "Metagene profile - Line",
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def heatmap_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    agg: str = "median",
    nan_policy: str = "drop",
    x_mode: str = "relative",
    heatmap_mode: str = "genes",
    value_mode: str = "density",
    max_nan_frac: Optional[float] = None,
    drop_nan_rows: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
) -> str:
    _hv_init()
    if heatmap_mode not in {"genes", "aggregate"}:
        raise ValueError("heatmap_mode must be 'genes' or 'aggregate'")
    if heatmap_mode == "aggregate":
        x, y = line_df(
            drd,
            agg=agg,
            segments=segments,
            n_windows=n_windows,
            as_percent=as_percent,
            order=order,
            nan_fill=nan_fill,
            agg_scope="genes",
            nan_policy=nan_policy,
            x_mode=x_mode,
            value_mode=value_mode,
            max_nan_frac=max_nan_frac,
        )
        if y.size == 0:
            return _ensure_plotly(hv.render(hv.Curve([]), backend="plotly")).to_html(
                full_html=full_html,
                include_plotlyjs=include_js,
            )
        z = np.asarray(y, dtype=float)[None, :]
        regions = [f"{agg} profile"]
        bins = list(range(len(y)))
    else:
        z, regions, bins = heatmap_df(
            drd,
            segments=segments,
            n_windows=n_windows,
            as_percent=as_percent,
            order=order,
            nan_fill=nan_fill,
            agg=agg,
            nan_policy=nan_policy,
            x_mode=x_mode,
            value_mode=value_mode,
            max_nan_frac=max_nan_frac,
        )
    if z.size == 0:
        return _ensure_plotly(hv.render(hv.Curve([]), backend="plotly")).to_html(
            full_html=full_html,
            include_plotlyjs=include_js,
        )
    max_regions = 800
    if len(regions) > max_regions:
        regions = regions[:max_regions]
        z = z[:max_regions, :]
    data = [(b, r, z[i, j]) for i, r in enumerate(regions) for j, b in enumerate(bins) if np.isfinite(z[i, j]) or np.isnan(z[i, j])]
    hm = hv.HeatMap(data, kdims=["bin", "region"], vdims=["density"]).opts(
        colorbar=True,
        colorbar_opts={"title": "density (%)" if as_percent else "density"},
        invert_yaxis=True,
    )
    fig = _ensure_plotly(hv.render(hm, backend="plotly"))
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    decor_mode = "rel" if x_mode == "relative" else "abs"
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode=decor_mode, x_values=bins)
    fig.update_layout(
        xaxis_title="Metagene position (relative)" if x_mode == "relative" else "Metagene position (bp)",
        yaxis_title="Feature",
        height=700,
        width=1000,
        margin=dict(l=90, r=30, t=60, b=90),
        title=title or "Metagene profile - Heatmap",
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def box_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    nan_policy: str = "drop",
    per_region: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
) -> str:
    _hv_init()
    if per_region:
        if nan_policy not in {"drop", "zero", "keep"}:
            raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
        data = []
        for dens, lbl in zip(drd.densities, drd.labels):
            y = np.asarray(dens, dtype=float)
            if as_percent:
                y = y * 100.0
            if nan_fill is not None:
                y = np.where(np.isnan(y), nan_fill, y)
            if nan_policy == "zero":
                y = np.where(np.isnan(y), 0.0, y)
            if nan_policy == "drop":
                y = y[np.isfinite(y)]
            if y.size == 0:
                continue
            label = lbl if lbl is not None else f"region_{len(data)+1}"
            data.append((label, float(np.nanmean(y))))
        kdims = ["region"]
    else:
        dist = dist_df(
            drd,
            as_percent=as_percent,
            nan_fill=nan_fill,
            segments=segments,
            n_windows=n_windows,
            nan_policy=nan_policy,
        )
        data = [(str(b), v) for b, v, _ in dist]
        kdims = [hv.Dimension("zone", type=str)]
    y_label = "Mean methylation per feature (%)" if as_percent else "Mean methylation per feature"
    box = hv.BoxWhisker(data, kdims=kdims, vdims=["density"]).opts(
        ylabel=y_label,
        show_legend=False,
    )
    fig = _ensure_plotly(hv.render(box, backend="plotly"))
    fig.update_yaxes(range=[0, 100] if as_percent else None)
    fig.update_layout(
        margin=dict(l=90, r=30, t=60, b=80),
        title=title or "Metagene profile - Boxplot",
        showlegend=False,
        xaxis_title=("Region" if per_region else "Metagene position (relative)"),
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def violin_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    nan_policy: str = "drop",
    per_region: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
) -> str:
    _hv_init()
    if per_region:
        if nan_policy not in {"drop", "zero", "keep"}:
            raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
        data = []
        for dens, lbl in zip(drd.densities, drd.labels):
            y = np.asarray(dens, dtype=float)
            if as_percent:
                y = y * 100.0
            if nan_fill is not None:
                y = np.where(np.isnan(y), nan_fill, y)
            if nan_policy == "zero":
                y = np.where(np.isnan(y), 0.0, y)
            if nan_policy == "drop":
                y = y[np.isfinite(y)]
            if y.size == 0:
                continue
            label = lbl if lbl is not None else f"region_{len(data)+1}"
            for val in y:
                data.append((label, float(val)))
        kdims = ["region"]
    else:
        dist = dist_df(
            drd,
            as_percent=as_percent,
            nan_fill=nan_fill,
            segments=segments,
            n_windows=n_windows,
            nan_policy=nan_policy,
        )
        data = [(str(b), v) for b, v, _ in dist]
        kdims = [hv.Dimension("zone", type=str)]
    y_label = "Mean methylation per feature (%)" if as_percent else "Mean methylation per feature"
    viol = hv.Violin(data, kdims=kdims, vdims=["density"]).opts(
        ylabel=y_label,
        show_legend=False,
        box=True,
    )
    fig = _ensure_plotly(hv.render(viol, backend="plotly"))
    fig.update_yaxes(range=[0, 100] if as_percent else None)
    fig.update_layout(
        height=600,
        width=900,
        margin=dict(l=90, r=30, t=60, b=80),
        title=title or "Metagene profile - Violin",
        showlegend=False,
        xaxis_title=("Region" if per_region else "Metagene position (relative)"),
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


# -------- Wrappers: RegionReader + HcAnnotStore → HTML --------

def _drd_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None,
    agg_method=None,
    feature_type: str | None,
    reverse_negative: bool,
    labels: list[str] | None,
    limit: int | None,
    mode: str = "discretise",
    x_mode: str = "relative",
    add_flanks: bool = False,
    flank_bp: int = 2000,
    combine_parts: bool = False,
    parts: Sequence[str] | None = None,
) -> DiscreteRegionData:
    if segments is None:
        if combine_parts:
            segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
        else:
            segments = [Segment("region", 100)]
    if add_flanks:
        try:
            ft_map = annot.get_feature_types()
        except Exception:
            ft_map = {}
        gene_ids = ft_map.get("gene", []) if isinstance(ft_map, dict) else []
        if gene_ids:
            flank = int(abs(flank_bp))
            if "upstream_gene" not in ft_map:
                annot.add_flanks(gene_ids, -flank, "upstream_")
            if "downstream_gene" not in ft_map:
                annot.add_flanks(gene_ids, flank, "downstream_")
    if combine_parts:
        if mode != "raw":
            raise ValueError("combine_parts=True requires mode='raw'")
        if x_mode != "relative":
            raise ValueError("combine_parts=True requires x_mode='relative'")
        parts_order = list(parts) if parts is not None else ["upstream_gene", "gene", "downstream_gene"]
        parts_data = collect_parts_from_hcannot(annot, parts=parts_order, limit=limit)
        drd_map: dict[str, DiscreteRegionData] = {}
        for part in parts_order:
            contigs, auto_labels = parts_data.get(part, ([], []))
            if not contigs:
                continue
            drd_part = compute_discrete_regions(
                reader,
                contigs,
                segments=segments,
                agg_method=agg_method,
                reverse_negative=reverse_negative,
                labels=auto_labels,
                mode=mode,
                x_mode=x_mode,
            )
            drd_map[part] = drd_part
        if not drd_map:
            return DiscreteRegionData()
        return combine_parts_drd(drd_map, segments=segments, parts_order=parts_order)

    contigs, auto_labels = collect_contigs_from_hcannot(annot, feature_type=feature_type, limit=limit)
    use_labels = labels if labels is not None else auto_labels
    return compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        agg_method=agg_method,
        reverse_negative=reverse_negative,
        labels=use_labels,
        mode=mode,
        x_mode=x_mode,
    )


def line_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    agg: str | None = None,
    agg_method=None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    combine_parts: bool = False,
    parts: Sequence[str] | None = None,
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    agg_scope: str = "points",
    nan_policy: str | None = None,
    mode: str = "raw",
    x_mode: str = "relative",
    value_mode: str | None = None,
    max_nan_frac: Optional[float] = None,
    full_html: bool = False,
    include_js: str = "cdn",
    smooth: dict | int | None = None,
    connectgaps: bool = False,
    bsx1_compat: bool = False,
) -> str:
    if bsx1_compat:
        combine_parts = True
        if segments is None:
            segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
        elif len(segments) != 3:
            raise ValueError("bsx1_compat requires exactly 3 segments (up/body/down)")
        if value_mode is None:
            value_mode = "weighted"
        elif value_mode != "weighted":
            raise ValueError("bsx1_compat requires value_mode='weighted'")
        if agg is None:
            agg = "mean"
        elif agg != "mean":
            raise ValueError("bsx1_compat requires agg='mean'")
        if nan_policy is None:
            nan_policy = "keep"
        elif nan_policy != "keep":
            raise ValueError("bsx1_compat requires nan_policy='keep'")
        if smooth is None:
            smooth = 50
        connectgaps = True
    else:
        if value_mode is None:
            value_mode = "density"
        if agg is None:
            agg = "median"
        if nan_policy is None:
            nan_policy = "drop"

    drd = _drd_from_annot(
        reader,
        annot,
        segments=segments,
        agg_method=agg_method,
        feature_type=feature_type,
        reverse_negative=reverse_negative,
        labels=labels,
        limit=limit,
        mode=mode,
        x_mode=x_mode,
        add_flanks=add_flanks,
        flank_bp=flank_bp,
        combine_parts=combine_parts,
        parts=parts,
    )
    return line_html(
        drd,
        agg=agg,
        segments=segments,
        n_windows=n_windows,
        order=order,
        as_percent=as_percent,
        agg_scope=agg_scope,
        nan_policy=nan_policy,
        x_mode=x_mode,
        value_mode=value_mode,
        max_nan_frac=max_nan_frac,
        full_html=full_html,
        include_js=include_js,
        smooth=smooth,
        connectgaps=connectgaps,
    )


def heatmap_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    agg: str = "median",
    agg_method=None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    combine_parts: bool = False,
    parts: Sequence[str] | None = None,
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    nan_policy: str = "drop",
    mode: str = "raw",
    x_mode: str = "relative",
    heatmap_mode: str = "genes",
    value_mode: str = "density",
    max_nan_frac: Optional[float] = None,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    drd = _drd_from_annot(
        reader,
        annot,
        segments=segments,
        agg_method=agg_method,
        feature_type=feature_type,
        reverse_negative=reverse_negative,
        labels=labels,
        limit=limit,
        mode=mode,
        x_mode=x_mode,
        add_flanks=add_flanks,
        flank_bp=flank_bp,
        combine_parts=combine_parts,
        parts=parts,
    )
    return heatmap_html(
        drd,
        segments=segments,
        n_windows=n_windows,
        order=order,
        as_percent=as_percent,
        agg=agg,
        nan_policy=nan_policy,
        x_mode=x_mode,
        heatmap_mode=heatmap_mode,
        value_mode=value_mode,
        max_nan_frac=max_nan_frac,
        full_html=full_html,
        include_js=include_js,
    )


def box_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    nan_fill: Optional[float] = None,
    nan_policy: str = "drop",
    agg_method=None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    per_region: bool = False,
    as_percent: bool = True,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    drd = _drd_from_annot(
        reader,
        annot,
        segments=segments,
        agg_method=agg_method,
        feature_type=feature_type,
        reverse_negative=reverse_negative,
        labels=labels,
        limit=limit,
        add_flanks=add_flanks,
        flank_bp=flank_bp,
    )
    return box_html(drd, segments=segments, n_windows=n_windows, as_percent=as_percent, nan_fill=nan_fill, nan_policy=nan_policy, per_region=per_region, full_html=full_html, include_js=include_js)


def violin_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    nan_fill: Optional[float] = None,
    nan_policy: str = "drop",
    agg_method=None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    per_region: bool = False,
    as_percent: bool = True,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    drd = _drd_from_annot(
        reader,
        annot,
        segments=segments,
        agg_method=agg_method,
        feature_type=feature_type,
        reverse_negative=reverse_negative,
        labels=labels,
        limit=limit,
        add_flanks=add_flanks,
        flank_bp=flank_bp,
    )
    return violin_html(drd, segments=segments, n_windows=n_windows, as_percent=as_percent, nan_fill=nan_fill, nan_policy=nan_policy, per_region=per_region, full_html=full_html, include_js=include_js)
