from __future__ import annotations

from typing import Optional

import holoviews as hv
import numpy as np

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import Segment, segments_total_bins
from ._html_common import _bin_points_windows, _drd_from_annot, _ensure_plotly, _hv_init


def _dist_rows(
    drd: DiscreteRegionData,
    *,
    as_percent: bool = False,
    nan_fill: Optional[float] = None,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    nan_policy: str = "drop",
):
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    rows = []
    for pos, dens, lbl in zip(drd.positions, drd.densities, drd.labels):
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
        )
        for b, v in enumerate(binned):
            if np.isfinite(v):
                rows.append((b, float(v), lbl if lbl is not None else f"region_{len(rows)+1}"))
    return rows


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
        dist = _dist_rows(
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


def box_html_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None = None,
    n_windows: Optional[int] = None,
    nan_fill: Optional[float] = None,
    nan_policy: str = "drop",
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
        feature_type=feature_type,
        reverse_negative=reverse_negative,
        labels=labels,
        limit=limit,
        add_flanks=add_flanks,
        flank_bp=flank_bp,
    )
    return box_html(
        drd,
        segments=segments,
        n_windows=n_windows,
        as_percent=as_percent,
        nan_fill=nan_fill,
        nan_policy=nan_policy,
        per_region=per_region,
        full_html=full_html,
        include_js=include_js,
    )
