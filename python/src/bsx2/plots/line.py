from __future__ import annotations

from typing import Optional

import holoviews as hv

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    Segment,
    _apply_savgol_smoothing,
    _clip_profile,
    _coerce_smooth_config,
    segments_total_bins,
)
from ._html_common import (
    _drd_from_annot,
    _ensure_plotly,
    _hv_init,
    _segment_decor_rel,
    line_df,
)


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

    x, y = line_df(
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
    if segments is None:
        segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
    drd = _drd_from_annot(
        reader,
        annot,
        segments=segments,
        feature_type=None,
        reverse_negative=True,
        labels=None,
        limit=None,
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
