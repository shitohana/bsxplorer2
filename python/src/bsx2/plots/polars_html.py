from __future__ import annotations
from typing import Tuple, Sequence

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    Segment,
    collect_contigs_from_hcannot,
    compute_discrete_regions,
    segments_total_bins,
)


def discrete_to_long_pd(drd: DiscreteRegionData, *, as_percent: bool = False) -> pd.DataFrame:
    rows = []
    for ridx, (pos, dens, lbl) in enumerate(zip(drd.positions, drd.densities, drd.labels)):
        region = lbl if lbl else f"region_{ridx+1}"
        vals = dens * 100.0 if as_percent else dens
        vals = np.asarray(vals, dtype=float)
        vals[~np.isfinite(vals)] = np.nan
        bins = np.arange(len(dens), dtype=int)
        rows.append(pd.DataFrame({"region": region, "bin": bins, "x": pos, "density": vals}))
    if not rows:
        return pd.DataFrame({"region": [], "bin": [], "x": [], "density": []})
    return pd.concat(rows, ignore_index=True)

# Backwards compatibility: old name used in tests
discrete_to_long_pl = discrete_to_long_pd


def line_df(drd: DiscreteRegionData, agg: str = "mean", *, as_percent: bool = False, order: Sequence[str] | None = None) -> pd.DataFrame:
    df = discrete_to_long_pd(drd, as_percent=as_percent)
    if df.empty:
        return pd.DataFrame({"x": [], "y": []})
    agg_fn = {"mean": "mean", "median": "median", "max": "max", "min": "min"}[agg]
    grouped = df.groupby("bin", as_index=False).agg({"x": "first", "density": agg_fn})
    if order is not None and len(order) == len(grouped):
        grouped["bin"] = list(order)
    grouped = grouped.rename(columns={"density": "y"}).sort_values("bin")
    return grouped[["x", "y"]]


def heatmap_df(drd: DiscreteRegionData, *, as_percent: bool = False, order: Sequence[str] | None = None) -> Tuple[pd.DataFrame, Sequence[str], Sequence[int]]:
    df = discrete_to_long_pd(drd, as_percent=as_percent)
    if df.empty:
        return df, [], []
    regions = df["region"].unique().tolist()
    if order:
        regions = [r for r in order if r in regions]
    else:
        regions = sorted(regions)
    bins = sorted(df["bin"].unique().tolist())
    return df[["region", "bin", "density"]], regions, bins


def dist_df(drd: DiscreteRegionData, *, as_percent: bool = False) -> pd.DataFrame:
    df = discrete_to_long_pd(drd, as_percent=as_percent)
    if df.empty:
        return df
    out = df[["region", "bin", "density"]].copy()
    out = out[np.isfinite(out["density"])]
    return out


def _segment_decor(fig, segments: list[Segment] | None, *, annotate_tss_tes: bool = False, x_mode: str) -> None:
    """Apply segment boundaries and ticks to a plotly Figure.

    x_mode: "rel" for 0..1 positions, "bin" for bin indices.
    """
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
    # Convert to x-axis scale
    if x_mode == "rel":
        scale = lambda v: v / total
    else:
        scale = float
    shapes = []
    for b in boundaries[:-1]:
        xb = scale(b)
        shapes.append(
            dict(type="line", x0=xb, x1=xb, y0=0, y1=1, yref="paper", line=dict(dash="dash", width=1, color="gray"))
        )
    tickvals = [scale(c) for c in centers]
    fig.update_xaxes(tickmode="array", tickvals=tickvals, ticktext=labels)
    fig.update_layout(shapes=shapes)
    # Optional TSS/TES annotations for classic layout (>=3 segments)
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
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    pdf = line_df(drd, agg=agg, as_percent=as_percent, order=order)
    y_label = f"{agg} density" + (" (%)" if as_percent else "")
    fig = px.line(pdf, x="x", y="y", labels={"x": "relative position", "y": y_label})
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="rel")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def heatmap_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    order: Sequence[str] | None = None,
    as_percent: bool = True,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    long_df, regions, bins = heatmap_df(drd, as_percent=as_percent, order=order)
    if long_df.empty:
        return go.Figure().to_html(full_html=full_html, include_plotlyjs=include_js)
    pivot = long_df.pivot(index="region", columns="bin", values="density").reindex(index=regions)
    z = pivot.values
    y = pivot.index.tolist()
    x = sorted(bins)
    bar_title = "density (%)" if as_percent else "density"
    fig = go.Figure(data=go.Heatmap(z=z, x=x, y=y, colorscale="Viridis", colorbar=dict(title=bar_title)))
    fig.update_layout(xaxis_title="bin", yaxis_title="region", yaxis_autorange="reversed")
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="bin")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def box_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    as_percent: bool = True,
    per_region: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    pdf = dist_df(drd, as_percent=as_percent)
    if per_region:
        # aggregate per region (mean over bins)
        pdf = pdf.groupby("region", as_index=False)["density"].mean()
        x_col = "region"
    else:
        x_col = "bin"
    y_label = "density (%)" if as_percent else "density"
    fig = px.box(pdf, x=x_col, y="density", labels={x_col: x_col, "density": y_label}, points=False)
    fig.update_yaxes(range=[0, 100] if as_percent else None)
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    if not per_region:
        fig.update_xaxes(type="linear")
        _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="bin")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def violin_html(
    drd: DiscreteRegionData,
    *,
    segments: list[Segment] | None = None,
    as_percent: bool = True,
    per_region: bool = False,
    full_html: bool = False,
    include_js: str = "cdn",
) -> str:
    pdf = dist_df(drd, as_percent=as_percent)
    if per_region:
        pdf = pdf.groupby("region", as_index=False)["density"].mean()
        x_col = "region"
    else:
        x_col = "bin"
    y_label = "density (%)" if as_percent else "density"
    fig = px.violin(pdf, x=x_col, y="density", box=True, points=False, labels={x_col: x_col, "density": y_label})
    fig.update_yaxes(range=[0, 100] if as_percent else None)
    annotate = bool(segments) and any(s.name.lower() == "body" for s in segments)
    if not per_region:
        fig.update_xaxes(type="linear")
        _segment_decor(fig, segments, annotate_tss_tes=annotate, x_mode="bin")
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
