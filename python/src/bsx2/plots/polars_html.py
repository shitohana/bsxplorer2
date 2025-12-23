from __future__ import annotations
from typing import Tuple, Sequence, List, Optional

import numpy as np
import holoviews as hv
import plotly.graph_objects as go

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    Segment,
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


def line_df(drd: DiscreteRegionData, agg: str = "mean", *, as_percent: bool = False, order: Sequence[str] | None = None, nan_fill: Optional[float] = None):
    """Возвращает два массива x, y (nan-aware)."""
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


def heatmap_df(drd: DiscreteRegionData, *, as_percent: bool = False, order: Sequence[str] | None = None, nan_fill: Optional[float] = None):
    """Возвращает (z, regions, bins) где z — np.ndarray shape (n_regions, n_bins)."""
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


def dist_df(drd: DiscreteRegionData, *, as_percent: bool = False, nan_fill: Optional[float] = None):
    """Возвращает список (group, density) с фильтром isfinite."""
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
    drop_nan_rows: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    _hv_init()
    x, y = line_df(drd, agg=agg, as_percent=as_percent, order=order, nan_fill=nan_fill)
    curve = hv.Curve((x, y), kdims="relative position", vdims="density")
    curve = curve.opts(
        xlabel="relative position",
        ylabel=f"{agg} density" + (" (%)" if as_percent else ""),
        show_legend=False,
    )
    fig = _ensure_plotly(hv.render(curve, backend="plotly"))
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="rel")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def heatmap_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    drop_nan_rows: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    _hv_init()
    z, regions, bins = heatmap_df(drd, as_percent=as_percent, order=order, nan_fill=nan_fill)
    if z.size == 0:
        return hv.render(hv.Curve([]), backend="plotly").to_html(full_html=full_html, include_plotlyjs=include_js)
    data = [ (b, r, z[i,j]) for i,r in enumerate(regions) for j,b in enumerate(bins) if np.isfinite(z[i,j]) or np.isnan(z[i,j]) ]
    hm = hv.HeatMap(data, kdims=["bin","region"], vdims=["density"]).opts(
        colorbar=True,
        colorbar_opts={"title": "density (%)" if as_percent else "density"},
        invert_yaxis=True,
    )
    fig = _ensure_plotly(hv.render(hm, backend="plotly"))
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="bin")
    fig.update_layout(xaxis_title="bin", yaxis_title="region")
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
) -> str:
    _hv_init()
    dist = dist_df(drd, as_percent=as_percent, nan_fill=nan_fill)
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
        kdims = [hv.Dimension("bin", type=str)]
    y_label = "density (%)" if as_percent else "density"
    box = hv.BoxWhisker(data, kdims=kdims, vdims=["density"]).opts(
        ylabel=y_label,
        show_legend=False,
    )
    fig = _ensure_plotly(hv.render(box, backend="plotly"))
    fig.update_yaxes(range=[0, 100] if as_percent else None)
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    if not per_region:
        _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="bin")
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
) -> str:
    _hv_init()
    dist = dist_df(drd, as_percent=as_percent, nan_fill=nan_fill)
    if per_region:
        by_region = {}
        for b, v, r in dist:
            by_region.setdefault(r, []).append(v)
        data = [(r, val) for r, vals in by_region.items() for val in vals]
        kdims = ["region"]
    else:
        data = [(str(b), v) for b, v, _ in dist]
        kdims = [hv.Dimension("bin", type=str)]
    y_label = "density (%)" if as_percent else "density"
    viol = hv.Violin(data, kdims=kdims, vdims=["density"]).opts(
        ylabel=y_label,
        show_legend=False,
        box=True,
    )
    fig = _ensure_plotly(hv.render(viol, backend="plotly"))
    fig.update_yaxes(range=[0, 100] if as_percent else None)
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    if not per_region:
        _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="bin")
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
    return compute_discrete_regions(rr, contigs, labels, segments)


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
) -> DiscreteRegionData:
    if segments is None:
        segments = [Segment("region", 100)]
    contigs, auto_labels = collect_contigs_from_hcannot(annot, feature_type=feature_type, limit=limit)
    use_labels = labels if labels is not None else auto_labels
    return compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        agg_method=agg_method,
        reverse_negative=reverse_negative,
        labels=use_labels,
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
    order: Sequence[str] | None = None,
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
    )
    return line_html(drd, agg=agg, segments=segments, order=order, as_percent=as_percent, full_html=full_html, include_js=include_js)


def heatmap_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    agg_method=None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    order: Sequence[str] | None = None,
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
    )
    return heatmap_html(drd, segments=segments, order=order, as_percent=as_percent, full_html=full_html, include_js=include_js)


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
    )
    return violin_html(drd, segments=segments, as_percent=as_percent, per_region=per_region, full_html=full_html, include_js=include_js)
