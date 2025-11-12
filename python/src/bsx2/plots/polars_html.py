from __future__ import annotations
from typing import Tuple, Sequence

import polars as pl
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    Segment,
    compute_discrete_regions,
    collect_contigs_from_hcannot,
)


def discrete_to_long_pl(drd: DiscreteRegionData) -> pl.DataFrame:
    rows = []
    for ridx, (pos, dens, lbl) in enumerate(zip(drd.positions, drd.densities, drd.labels)):
        region = lbl if lbl else f"region_{ridx+1}"
        n = len(dens)
        rows.append(pl.DataFrame({"region": [region] * n, "bin": np.arange(n), "x": pos, "density": dens}))
    return pl.concat(rows, how="vertical_relaxed") if rows else pl.DataFrame({"region": [], "bin": [], "x": [], "density": []})


def line_df(drd: DiscreteRegionData, agg: str = "mean") -> pl.DataFrame:
    df = discrete_to_long_pl(drd)
    if df.is_empty():
        return pl.DataFrame({"x": [], "y": []})
    agg_expr = {
        "mean": pl.col("density").mean(),
        "median": pl.col("density").median(),
        "max": pl.col("density").max(),
        "min": pl.col("density").min(),
    }[agg]
    return df.group_by("bin").agg(x=pl.col("x").first(), y=agg_expr).sort("bin").select(["x", "y"])


def heatmap_df(drd: DiscreteRegionData) -> Tuple[pl.DataFrame, Sequence[str], Sequence[int]]:
    df = discrete_to_long_pl(drd)
    if df.is_empty():
        return df, [], []
    regions = df.select(pl.col("region")).unique().sort("region").to_series().to_list()
    bins = df.select(pl.col("bin")).unique().sort("bin").to_series().to_list()
    return df.select(["region", "bin", "density"]), regions, bins


def dist_df(drd: DiscreteRegionData) -> pl.DataFrame:
    df = discrete_to_long_pl(drd)
    return df.select(["bin", "density"]) if not df.is_empty() else df


def line_html(drd: DiscreteRegionData, agg: str = "mean", *, full_html: bool = False, include_js: str = "cdn") -> str:
    pdf = line_df(drd, agg=agg).to_pandas()
    fig = px.line(pdf, x="x", y="y", labels={"x": "relative position", "y": f"{agg} density"})
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def heatmap_html(drd: DiscreteRegionData, *, full_html: bool = False, include_js: str = "cdn") -> str:
    long_df, regions, bins = heatmap_df(drd)
    if long_df.is_empty():
        return go.Figure().to_html(full_html=full_html, include_plotlyjs=include_js)
    pivot = long_df.pivot(values="density", index="region", columns="bin", aggregate_function="first").sort("region")
    z = pivot.select(pl.all().exclude("region")).to_numpy()
    y = pivot.select("region").to_series().to_list()
    x = sorted(bins)
    fig = go.Figure(data=go.Heatmap(z=z, x=x, y=y, colorscale="Viridis", colorbar=dict(title="density")))
    fig.update_layout(xaxis_title="bin", yaxis_title="region", yaxis_autorange="reversed")
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def box_html(drd: DiscreteRegionData, *, full_html: bool = False, include_js: str = "cdn") -> str:
    pdf = dist_df(drd).to_pandas()
    fig = px.box(pdf, x="bin", y="density", labels={"bin": "bin", "density": "density"})
    return fig.to_html(full_html=full_html, include_plotlyjs=include_js)


def violin_html(drd: DiscreteRegionData, *, full_html: bool = False, include_js: str = "cdn") -> str:
    pdf = dist_df(drd).to_pandas()
    fig = px.violin(pdf, x="bin", y="density", box=True, points=False, labels={"bin": "bin", "density": "density"})
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
    return line_html(drd, agg=agg, full_html=full_html, include_js=include_js)


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
    return heatmap_html(drd, full_html=full_html, include_js=include_js)


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
    return box_html(drd, full_html=full_html, include_js=include_js)


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
    return violin_html(drd, full_html=full_html, include_js=include_js)
