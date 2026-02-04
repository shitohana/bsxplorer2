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


def _bin_points(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    *,
    weights: Optional[np.ndarray] = None,
    segments: list[Segment],
    agg: str,
    nan_policy: str,
    value_mode: str = "density",
) -> np.ndarray:
    bounds = segment_boundaries(segments)
    total_bins = segments_total_bins(segments)
    starts = np.array([0.0] + bounds[:-1], dtype=float)
    ends = np.array(bounds, dtype=float)
    seg_nbins = np.array([s.n_bins for s in segments], dtype=int)
    seg_offsets = np.concatenate(([0], np.cumsum(seg_nbins)[:-1]))

    if value_mode not in {"density", "weighted"}:
        raise ValueError("value_mode must be 'density' or 'weighted'")

    bins_values = [[] for _ in range(total_bins)]
    bins_weighted_sum = np.zeros(total_bins, dtype=float)
    bins_weighted_total = np.zeros(total_bins, dtype=float)
    for i, (x, y) in enumerate(zip(x_vals, y_vals)):
        if x < 0.0 or x > 1.0:
            continue
        if not np.isfinite(x):
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
        seg_idx = int(np.searchsorted(ends, x, side="right"))
        seg_idx = min(seg_idx, len(ends) - 1)
        width = ends[seg_idx] - starts[seg_idx]
        if width <= 0:
            continue
        x_seg = (x - starts[seg_idx]) / width
        local = int(x_seg * seg_nbins[seg_idx])
        local = min(local, seg_nbins[seg_idx] - 1)
        global_bin = int(seg_offsets[seg_idx] + local)
        if 0 <= global_bin < total_bins:
            if value_mode == "weighted":
                if not np.isfinite(y):
                    if nan_policy == "zero":
                        y = 0.0
                    else:
                        continue
                w_use = 1.0 if w is None else w
                bins_weighted_sum[global_bin] += float(y) * w_use
                bins_weighted_total[global_bin] += w_use
            else:
                bins_values[global_bin].append(y)

    if value_mode == "weighted":
        out = np.full(total_bins, np.nan, dtype=float)
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
    agg: str = "mean",
    *,
    segments: list[Segment] | None = None,
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
    if segments:
        bounds = segment_boundaries(segments)
        total_bins = segments_total_bins(segments)
        starts = np.array([0.0] + bounds[:-1], dtype=float)
        ends = np.array(bounds, dtype=float)
        seg_nbins = np.array([s.n_bins for s in segments], dtype=int)
        seg_offsets = np.concatenate(([0], np.cumsum(seg_nbins)[:-1]))

        centers = []
        for s, start, end in zip(segments, starts, ends):
            width = end - start
            if width <= 0:
                continue
            idx = np.arange(s.n_bins, dtype=float)
            centers.append(start + (idx + 0.5) / s.n_bins * width)
        x_out = np.concatenate(centers) if centers else np.array([], dtype=float)

        if agg_scope not in {"points", "genes"}:
            raise ValueError("agg_scope must be 'points' or 'genes'")
        if nan_policy not in {"drop", "zero", "keep"}:
            raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
        if x_mode not in {"relative", "absolute"}:
            raise ValueError("x_mode must be 'relative' or 'absolute'")

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
            y_out = _bin_points(
                np.asarray(xs),
                np.asarray(ys),
                weights=w_arr,
                segments=segments,
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
                _bin_points(
                    x_rel,
                    y,
                    weights=weights,
                    segments=segments,
                    agg=agg,
                    nan_policy=nan_policy,
                    value_mode=value_mode,
                )
            )
        if not rows:
            return np.array([]), np.array([])
        mat = np.vstack(rows)
        if max_nan_frac is not None and mat.size > 0:
            max_nan_frac = float(max_nan_frac)
            nan_frac = np.mean(~np.isfinite(mat), axis=1)
            keep = nan_frac <= max_nan_frac
            mat = mat[keep]
        if nan_policy == "keep":
            y_out = np.mean(mat, axis=0)
        else:
            y_out = np.nanmean(mat, axis=0)
        if x_mode == "absolute":
            if lengths:
                scale = float(np.median(np.asarray(lengths)))
                x_out = x_out * scale
        return x_out, y_out

    rows = discrete_to_long(drd, as_percent=as_percent, nan_fill=nan_fill)
    if not rows:
        return np.array([]), np.array([])
    data = np.array(rows, dtype=object)
    bins = data[:, 1].astype(int)
    x_vals = data[:, 2].astype(float)
    y_vals = data[:, 3].astype(float)
    max_bin = bins.max() + 1
    y_out = []
    x_out = []
    global_mean = np.nanmean(y_vals) if np.isfinite(y_vals).any() else 0.0
    for b in range(max_bin):
        mask = bins == b
        if not mask.any():
            continue
        y_bin = y_vals[mask]
        x_bin = x_vals[mask][0]
        if agg == "mean":
            val = np.nanmean(y_bin)
        elif agg == "median":
            val = np.nanmedian(y_bin)
        elif agg == "max":
            val = np.nanmax(y_bin)
        elif agg == "min":
            val = np.nanmin(y_bin)
        else:
            raise ValueError(f"Unsupported agg: {agg}")
        if np.isnan(val):
            val = global_mean
        x_out.append(x_bin)
        y_out.append(val)
    return np.array(x_out, dtype=float), np.array(y_out, dtype=float)


def heatmap_df(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    as_percent: bool = False,
    order: Sequence[str] | None = None,
    nan_fill: Optional[float] = None,
    agg: str = "mean",
    nan_policy: str = "drop",
    x_mode: str = "relative",
    value_mode: str = "density",
    max_nan_frac: Optional[float] = None,
):
    """Sobrat (z, regions, bins), gde z - np.ndarray shape (n_regions, n_bins)."""
    if segments is not None:
        if nan_policy not in {"drop", "zero", "keep"}:
            raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
        if x_mode not in {"relative", "absolute"}:
            raise ValueError("x_mode must be 'relative' or 'absolute'")
        total_bins = segments_total_bins(segments)
        rows = []
        labels = []
        for ridx, (pos, dens, w, lbl) in enumerate(zip(drd.positions, drd.densities, drd.weights or [None] * len(drd.positions), drd.labels)):
            x = np.asarray(pos, dtype=float)
            y = np.asarray(dens, dtype=float)
            weights = None
            if w is not None:
                weights = np.asarray(w, dtype=float)
            if as_percent:
                y = y * 100.0
            if nan_fill is not None:
                y = np.where(np.isnan(y), nan_fill, y)
            if x_mode == "absolute":
                if x.size == 0:
                    continue
                x_min = float(np.nanmin(x))
                x_max = float(np.nanmax(x))
                length = x_max - x_min
                if not np.isfinite(length) or length <= 0:
                    continue
                x = (x - x_min) / length
            binned = _bin_points(
                x,
                y,
                weights=weights,
                segments=segments,
                agg=agg,
                nan_policy=nan_policy,
                value_mode=value_mode,
            )
            if nan_policy == "zero":
                binned = np.where(np.isnan(binned), 0.0, binned)
            rows.append(binned)
            labels.append(lbl if lbl else f"region_{ridx+1}")

        if not rows:
            return np.empty((0, 0)), [], []
        if order:
            idx_map = {lab: i for i, lab in enumerate(labels)}
            order_idx = [idx_map[r] for r in order if r in idx_map]
            rows = [rows[i] for i in order_idx]
            labels = [labels[i] for i in order_idx]
        z = np.vstack(rows)
        if max_nan_frac is not None and z.size > 0:
            max_nan_frac = float(max_nan_frac)
            nan_frac = np.mean(~np.isfinite(z), axis=1)
            keep = nan_frac <= max_nan_frac
            z = z[keep]
            labels = [l for l, k in zip(labels, keep) if k]
        bins = list(range(total_bins))
        return z, labels, bins

    rows = discrete_to_long(drd, as_percent=as_percent, nan_fill=nan_fill)
    if not rows:
        return np.empty((0, 0)), [], []
    data = np.array(rows, dtype=object)
    regions = list(dict.fromkeys(data[:, 0]))  # preserve insertion
    bins = sorted(set(data[:, 1].astype(int).tolist()))
    if order:
        regions = [r for r in order if r in regions]
    z = np.full((len(regions), len(bins)), np.nan, dtype=float)
    region_index = {r: i for i, r in enumerate(regions)}
    bin_index = {b: i for i, b in enumerate(bins)}
    for r, b, _, v in data:
        if r in region_index:
            z[region_index[r], bin_index[int(b)]] = float(v)
    return z, regions, bins


def dist_df(
    drd: DiscreteRegionData,
    *,
    as_percent: bool = False,
    nan_fill: Optional[float] = None,
    segments: list[Segment] | None = None,
):
    """
    Sbor vyboriki dlya box/violin.
    Esli peredany segments - agregiruem po zonam (segment.name), inache po binam.
    Kazhdaya tochka = odin region v odnoj zone.
    """
    if segments:
        mat, labels, _ = drd._stack_to_common_grid()
        if mat.size == 0:
            return []
        vals = mat.astype(float)
        if as_percent:
            vals = vals * 100.0
        vals[~np.isfinite(vals)] = np.nan
        if nan_fill is not None:
            vals = np.where(np.isnan(vals), nan_fill, vals)
        out: list[tuple] = []
        start = 0
        for seg in segments:
            end = min(start + seg.n_bins, vals.shape[1])
            if end <= start:
                start = end
                continue
            zone_slice = vals[:, start:end]
            for rid, row in zip(labels, zone_slice):
                finite = np.isfinite(row)
                if not finite.any():
                    continue
                val = np.nanmean(row[finite])
                if np.isfinite(val):
                    out.append((seg.name, float(val), rid))
            start = end
        return out

    rows = discrete_to_long(drd, as_percent=as_percent, nan_fill=nan_fill)
    out = []
    for r, b, _, v in rows:
        if np.isfinite(v):
            out.append((b, v, r))
    return out

def _segment_decor(fig, segments: list[Segment] | None, *, annotate_tss_tes: bool = False, x_mode: str = "rel") -> None:
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
    scale = (lambda v: v / total) if x_mode == "rel" else float
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
    agg: str = "mean",
    *,
    segments: list[Segment] | None = None,
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
) -> str:
    _hv_init()
    x, y = line_df(
        drd,
        agg=agg,
        segments=segments,
        as_percent=as_percent,
        order=order,
        nan_fill=nan_fill,
        agg_scope=agg_scope,
        nan_policy=nan_policy,
        x_mode=x_mode,
        value_mode=value_mode,
        max_nan_frac=max_nan_frac,
    )
    x_label = "Metagene position (relative)" if x_mode == "relative" else "Metagene position (bp)"
    curve = hv.Curve((x, y), kdims="relative position", vdims="density").opts(
        xlabel=x_label,
        ylabel=f"{agg} density" + (" (%)" if as_percent else ""),
        show_legend=False,
    )
    fig = _ensure_plotly(hv.render(curve, backend="plotly"))
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="rel")
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
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    agg: str = "mean",
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
            return hv.render(hv.Curve([]), backend="plotly").to_html(full_html=full_html, include_plotlyjs=include_js)
        z = np.asarray(y, dtype=float)[None, :]
        regions = [f"{agg} profile"]
        bins = list(range(len(y)))
    else:
        z, regions, bins = heatmap_df(
            drd,
            segments=segments,
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
        return hv.render(hv.Curve([]), backend="plotly").to_html(full_html=full_html, include_plotlyjs=include_js)
    max_regions = 800
    if len(regions) > max_regions:
        regions = regions[:max_regions]
        z = z[:max_regions, :]
    data = [ (b, r, z[i,j]) for i,r in enumerate(regions) for j,b in enumerate(bins) if np.isfinite(z[i,j]) or np.isnan(z[i,j]) ]
    hm = hv.HeatMap(data, kdims=["bin", "region"], vdims=["density"]).opts(
        colorbar=True,
        colorbar_opts={"title": "density (%)" if as_percent else "density"},
        invert_yaxis=True,
    )
    fig = _ensure_plotly(hv.render(hm, backend="plotly"))
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="bin")
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
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    per_region: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
) -> str:
    _hv_init()
    dist = dist_df(drd, as_percent=as_percent, nan_fill=nan_fill, segments=segments)
    if per_region:
        # агрегат средний по бинам на регион
        by_region = {}
        for b, v, r in dist:
            by_region.setdefault(r, []).append(v)
        data = [(r, np.nanmean(vals)) for r, vals in by_region.items()]
        kdims = ["region"]
    else:
        # Для совместимости с hv plotly backend ключи делаем строками
        data = [(str(b), v) for b, v, _ in dist]
        kdims = [hv.Dimension("zone", type=str)]
    y_label = "Mean methylation per feature (%)" if as_percent else "Mean methylation per feature"
    box = hv.BoxWhisker(data, kdims=kdims, vdims=["density"]).opts(
        ylabel=y_label,
        show_legend=False,
    )
    fig = _ensure_plotly(hv.render(box, backend="plotly"))
    fig.update_yaxes(range=[0, 100] if as_percent else None)
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    # zone-level axis, without segment dividers
    fig.update_layout(
        margin=dict(l=90, r=30, t=60, b=80),
        title=title or "Metagene profile - Boxplot",
        showlegend=False,
        xaxis_title=('Region' if per_region else 'Metagene position (relative)'),
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def violin_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    per_region: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
    title: Optional[str] = None,
) -> str:
    _hv_init()
    dist = dist_df(drd, as_percent=as_percent, nan_fill=nan_fill, segments=segments)
    if per_region:
        by_region = {}
        for b, v, r in dist:
            by_region.setdefault(r, []).append(v)
        data = [(r, val) for r, vals in by_region.items() for val in vals]
        kdims = ["region"]
    else:
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
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    # zone-level axis, without segment dividers
    fig.update_layout(
        height=600,
        width=900,
        margin=dict(l=90, r=30, t=60, b=80),
        title=title or "Metagene profile - Violin",
        showlegend=False,
        xaxis_title=('Region' if per_region else 'Metagene position (relative)'),
    )
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


# ---------------- Wrappers from annot/BSX (used in end-to-end tests) ----------------

def _drd_from_annot(
    rr,
    annot,
    *,
    segments: list[Segment],
    agg_method,
    feature_type: str | None,
    limit: int | None,
):
    contigs, labels = collect_contigs_from_hcannot(annot, feature_type=feature_type, limit=limit)
    if not contigs:
        return DiscreteRegionData()
    try:
        contigs = rr.index().sort(list(contigs))
    except Exception:
        pass
    return compute_discrete_regions(rr, contigs, segments=segments, agg_method=agg_method, labels=labels)


def line_html_from_annot(
    rr,
    annot,
    *,
    segments: list[Segment],
    agg: str = "mean",
    agg_method=None,
    feature_type: str | None = None,
    limit: int | None = None,
    full_html: bool = False,
) -> str:
    drd = _drd_from_annot(rr, annot, segments=segments, agg_method=agg_method, feature_type=feature_type, limit=limit)
    return line_html(drd, agg=agg, segments=segments, full_html=full_html)


def heatmap_html_from_annot(
    rr,
    annot,
    *,
    segments: list[Segment],
    agg_method=None,
    feature_type: str | None = None,
    limit: int | None = None,
    full_html: bool = False,
) -> str:
    drd = _drd_from_annot(rr, annot, segments=segments, agg_method=agg_method, feature_type=feature_type, limit=limit)
    return heatmap_html(drd, segments=segments, full_html=full_html)


def box_html_from_annot(
    rr,
    annot,
    *,
    segments: list[Segment],
    agg_method=None,
    feature_type: str | None = None,
    limit: int | None = None,
    full_html: bool = False,
) -> str:
    drd = _drd_from_annot(rr, annot, segments=segments, agg_method=agg_method, feature_type=feature_type, limit=limit)
    return box_html(drd, segments=segments, full_html=full_html)


def violin_html_from_annot(
    rr,
    annot,
    *,
    segments: list[Segment],
    agg_method=None,
    feature_type: str | None = None,
    limit: int | None = None,
    full_html: bool = False,
) -> str:
    drd = _drd_from_annot(rr, annot, segments=segments, agg_method=agg_method, feature_type=feature_type, limit=limit)
    return violin_html(drd, segments=segments, full_html=full_html)


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
) -> DiscreteRegionData:
    if segments is None:
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
    agg: str = "mean",
    agg_method=None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    agg_scope: str = "points",
    nan_policy: str = "drop",
    mode: str = "raw",
    x_mode: str = "relative",
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
    )
    return line_html(
        drd,
        agg=agg,
        segments=segments,
        order=order,
        as_percent=as_percent,
        agg_scope=agg_scope,
        nan_policy=nan_policy,
        x_mode=x_mode,
        value_mode=value_mode,
        max_nan_frac=max_nan_frac,
        full_html=full_html,
        include_js=include_js,
    )


def heatmap_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    agg: str = "mean",
    agg_method=None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
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
    )
    return heatmap_html(
        drd,
        segments=segments,
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
    return box_html(drd, segments=segments, as_percent=as_percent, per_region=per_region, full_html=full_html, include_js=include_js)


def violin_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
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
    return violin_html(drd, segments=segments, as_percent=as_percent, per_region=per_region, full_html=full_html, include_js=include_js)
