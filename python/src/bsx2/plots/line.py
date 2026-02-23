from __future__ import annotations

from typing import Optional
import holoviews as hv
import numpy as np
import warnings
from bsx2.plots.data import DiscreteRegionData
from bsx2.validation import validate_n_windows, validate_window_agg
from bsx2.plots.metagene import (
    MetageneProfileSegment,
    _apply_savgol_smoothing,
    _clip_profile,
    _coerce_smooth_config,
    segments_total_bins,
)
from ._common import _bin_points_windows_fast


def _line_profile(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    nan_fill: Optional[float] = None,
    agg: str = "mean",
):
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    else:
        n_windows = validate_n_windows(n_windows)

    agg = validate_window_agg(agg)

    x_out = (np.arange(n_windows, dtype=float) + 0.5) / float(n_windows)

    xs_parts: list[np.ndarray] = []
    ys_parts: list[np.ndarray] = []

    for pos, dens in zip(drd.positions, drd.densities):
        x = np.asarray(pos, dtype=np.float64)
        y = np.asarray(dens, dtype=np.float64)

        if x.size == 0:
            continue

        if nan_fill is not None:
           
            y = np.where(np.isnan(y), nan_fill, y)

        xs_parts.append(x)
        ys_parts.append(y)

    if not xs_parts:
        empty = np.array([], dtype=np.float64)
        return empty, empty


    x_all = np.concatenate(xs_parts)
    y_all = np.concatenate(ys_parts)

    
    if agg == "median":
        order = np.argsort(x_all, kind="mergesort")  
        x_all = x_all[order]
        y_all = y_all[order]

    y_out = _bin_points_windows_fast(
        x_all,
        y_all,
        n_windows=n_windows,
        agg=agg,
        nan_policy="keep",
    )
    return x_out, y_out

def line_plot(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    agg: str = "mean",
    smooth: dict | int | None = 50,
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
) -> hv.Curve:
    agg = validate_window_agg(agg)
    if segments is None:
        segments = [MetageneProfileSegment("up", 100), MetageneProfileSegment("body", 200), MetageneProfileSegment("down", 100)]
    if n_windows is None:
        n_windows = segments_total_bins(segments)

    x, y = _line_profile(
        drd,
        segments=segments,
        n_windows=n_windows,
        nan_fill=None,
        agg=agg,
    )
    if smooth is not None and agg in {"min", "max"}:
        warnings.warn(
            f"smooth={smooth!r} ignored for agg={agg!r}; smoothing is only applied to mean/median",
            RuntimeWarning,
            stacklevel=2,
        )
    # Smoothing is applied only to mean/median profiles.
    if smooth is not None and y.size > 0 and agg in {"mean", "median"}:
        total_bins = segments_total_bins(segments)
        smooth_cfg = _coerce_smooth_config(smooth, total_bins=total_bins)
        if smooth_cfg is not None:
            y_scaled = y.astype(float, copy=True)
            y_scaled = _apply_savgol_smoothing(y_scaled, smooth_cfg, segments=segments)
            y_scaled = _clip_profile(y_scaled)
            y = y_scaled

    if x.size == 0 or y.size == 0:
        curve = hv.Curve([])
    else:
        curve = hv.Curve((x, y), kdims="relative position", vdims="density")

    opts_kwargs = dict(
        xlabel="Metagene position (relative)",
        ylabel=f"{agg.capitalize()} methylation density",
        show_legend=False,
        title=title or "Metagene profile - Line",
    )
    if width is not None:
        opts_kwargs["width"] = int(width)
    if height is not None:
        opts_kwargs["height"] = int(height)

    return curve.opts(**opts_kwargs)

