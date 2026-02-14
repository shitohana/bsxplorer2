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
    if value_mode != "density":
        raise ValueError("value_mode must be 'density'")

    bins_values = [[] for _ in range(n_windows)]
    for x, y in zip(x_vals, y_vals):
        if not np.isfinite(x):
            continue
        if x < 0.0 or x > 1.0:
            continue
        if nan_policy == "zero" and not np.isfinite(y):
            y = 0.0
        if nan_policy == "drop" and not np.isfinite(y):
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
        for pos, dens in zip(drd.positions, drd.densities):
            x = np.asarray(pos, dtype=float)
            y = np.asarray(dens, dtype=float)
            if as_percent:
                y = y * 100.0
            if nan_fill is not None:
                y = np.where(np.isnan(y), nan_fill, y)
            x_rel, length = _to_relative(x)
            if length is not None:
                lengths.append(length)
            x = x_rel
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
    for pos, dens in zip(drd.positions, drd.densities):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
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
    for pos, dens, lbl in zip(
        drd.positions,
        drd.densities,
        drd.labels,
    ):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
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
        bins = (np.arange(n_windows, dtype=float) + 0.5) / float(n_windows)
    return mat, labels, bins


def _rank_compress(z_sorted: np.ndarray, rank_rows: int, *, fill: float | None = 0.0) -> np.ndarray:
    """BSX1-style rank compression: mean over finite values per bucket; fill only empty bins."""
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
    for pos, dens, lbl in zip(
        drd.positions,
        drd.densities,
        drd.labels,
    ):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        if as_percent:
            y = y * 100.0
        if nan_fill is not None:
            y = np.where(np.isnan(y), nan_fill, y)
        binned = _bin_points_windows(
            x,
            y,
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
    """BSX1-style line: mean density, relative axis, up/body/down."""
    _hv_init()
    if segments is None:
        segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    if n_windows is None:
        n_windows = segments_total_bins(segments)

    x, y = line_df(
        drd,
        agg="mean",
        segments=segments,
        n_windows=n_windows,
        as_percent=False,
        order=None,
        nan_fill=None,
        agg_scope="points",
        nan_policy="keep",
        x_mode="relative",
        value_mode="density",
        max_nan_frac=None,
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
    annotate = True
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="rel", x_values=x)
    fig.update_layout(
        height=600 if height is None else int(height),
        width=1000 if width is None else int(width),
        margin=dict(l=70, r=30, t=60, b=70),
        title=title or "Metagene profile - Line",
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


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
    """BSX1-style ranked heatmap (no modes)."""
    _hv_init()
    if rank_score not in {"mean", "body_mean"}:
        raise ValueError("rank_score must be 'mean' or 'body_mean'")
    if sort_order not in {"asc", "desc"}:
        raise ValueError("sort_order must be 'asc' or 'desc'")
    if segments is None:
        segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    if n_windows is None:
        n_windows = segments_total_bins(segments)

    z, _, bins = heatmap_df(
        drd,
        segments=segments,
        n_windows=n_windows,
        as_percent=False,
        order=None,
        nan_fill=None,
        agg="mean",
        nan_policy="keep",
        x_mode="relative",
        value_mode="density",
        max_nan_frac=None,
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
    fig_height = height if height is not None else max(
        550, min(900, 5 * int(z_vis.shape[0]))
    )

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
    _segment_decor(fig, segments, annotate_tss_tes=True, x_mode="bin", x_values=bins)
    fig.update_layout(
        title=title or "Metagene profile - Heatmap (BSX1 ranked)",
        width=fig_width,
        height=fig_height,
        margin=dict(l=70, r=30, t=60, b=70),
        xaxis_title="Position (bin)",
        yaxis_title="Rank",
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
                r = annot.add_flanks(gene_ids, -flank, "upstream_")
                if r is not None:
                    annot = r
            if "downstream_gene" not in ft_map:
                r = annot.add_flanks(gene_ids, flank, "downstream_")
                if r is not None:
                    annot = r
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
    add_flanks: bool = False,
    flank_bp: int = 2000,
    smooth: dict | int | None = 50,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
) -> str:
    """BSX1-style line from annotation: mean density, relative axis, up/body/down."""
    if segments is None:
        segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    drd = _drd_from_annot(
        reader,
        annot,
        segments=segments,
        agg_method=None,
        feature_type=None,
        reverse_negative=True,
        labels=None,
        limit=None,
        mode="raw",
        x_mode="relative",
        add_flanks=add_flanks,
        flank_bp=flank_bp,
        combine_parts=True,
        parts=None,
    )
    return line_html(
        drd,
        segments=segments,
        n_windows=segments_total_bins(segments),
        smooth=smooth,
        full_html=full_html,
        include_js=include_js,
        title=title,
        width=width,
        height=height,
    )


def heatmap_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    rank_rows: int = 100,
    rank_score: str = "mean",
    sort_order: str = "desc",
    colorscale: str = "Viridis",
    vmax_q: float | None = 0.995,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
) -> str:
    """BSX1-style ranked heatmap from annotation (no modes)."""
    if segments is None:
        segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    drd = _drd_from_annot(
        reader,
        annot,
        segments=segments,
        agg_method=None,
        feature_type=None,
        reverse_negative=True,
        labels=None,
        limit=None,
        mode="raw",
        x_mode="relative",
        add_flanks=add_flanks,
        flank_bp=flank_bp,
        combine_parts=True,
        parts=None,
    )
    return heatmap_html(
        drd,
        segments=segments,
        n_windows=segments_total_bins(segments),
        rank_rows=rank_rows,
        rank_score=rank_score,
        sort_order=sort_order,
        colorscale=colorscale,
        empty_bin_fill=0.0,
        vmax_q=vmax_q,
        full_html=full_html,
        include_js=include_js,
        title=title,
        width=width,
        height=height,
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
