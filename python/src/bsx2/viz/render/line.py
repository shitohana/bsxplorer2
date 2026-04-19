from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import holoviews as hv
import numpy as np
from beartype.typing import Optional

from bsx2 import AggMethod
from bsx2.validation import (
    validate_n_windows,
    validate_nan_policy,
    validate_positive_int,
    validate_window_agg,
)

from ..compute.data import DiscreteRegionData
from ..compute.metagene import (
    MetageneProfileSegment,
    _apply_savgol_smoothing,
    _clip_profile,
    _coerce_smooth_config,
    segments_total_bins,
)
from ..compute.windowing import _bin_points_windows_fast
from ._common import NanPolicy, _hv_init


def _agg_label(agg: AggMethod) -> str:
    name = getattr(agg, "name", None)
    if isinstance(name, str) and name:
        token = name
    else:
        token = str(agg).rsplit(".", 1)[-1]
    return token.lower().capitalize()


def _profile_title() -> str:
    return "Scaled region profile - Line"


def _position_axis_label() -> str:
    return "Relative feature position"


def _line_profile(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    nan_policy: NanPolicy = NanPolicy.KEEP,
    agg: AggMethod = AggMethod.Mean,
) -> tuple[np.ndarray, np.ndarray]:
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    else:
        n_windows = validate_n_windows(n_windows)

    x_out = (np.arange(n_windows, dtype=float) + 0.5) / float(n_windows)

    if not drd.positions:
        empty = np.array([], dtype=np.float64)
        return empty, empty

    if agg is AggMethod.Mean:
        sums = np.zeros(n_windows, dtype=np.float64)
        counts = np.zeros(n_windows, dtype=np.int64)
        for pos, dens in zip(drd.positions, drd.densities, strict=True):
            x = np.asarray(pos, dtype=np.float64)
            y = np.asarray(dens, dtype=np.float64)
            if x.size == 0:
                continue
            if nan_policy is NanPolicy.ZERO:
                y = np.where(np.isfinite(y), y, 0.0)
            else:
                finite = np.isfinite(y)
                if not np.any(finite):
                    continue
                x = x[finite]
                y = y[finite]

            idx = (x * n_windows).astype(np.int64)
            idx[idx == n_windows] = n_windows - 1
            counts += np.bincount(idx, minlength=n_windows)
            sums += np.bincount(idx, weights=y, minlength=n_windows)

        y_out = np.full(n_windows, np.nan, dtype=np.float64)
        nonempty = counts > 0
        y_out[nonempty] = sums[nonempty] / counts[nonempty]
        return x_out, y_out

    if agg in (AggMethod.Min, AggMethod.Max):
        y_out = (
            np.full(n_windows, np.inf, dtype=np.float64)
            if agg is AggMethod.Min
            else np.full(n_windows, -np.inf, dtype=np.float64)
        )
        touched = np.zeros(n_windows, dtype=bool)
        for pos, dens in zip(drd.positions, drd.densities, strict=True):
            x = np.asarray(pos, dtype=np.float64)
            y = np.asarray(dens, dtype=np.float64)
            if x.size == 0:
                continue
            if nan_policy is NanPolicy.ZERO:
                y = np.where(np.isfinite(y), y, 0.0)
            else:
                finite = np.isfinite(y)
                if not np.any(finite):
                    continue
                x = x[finite]
                y = y[finite]

            idx = (x * n_windows).astype(np.int64)
            idx[idx == n_windows] = n_windows - 1
            touched[idx] = True
            if agg is AggMethod.Min:
                np.minimum.at(y_out, idx, y)
            else:
                np.maximum.at(y_out, idx, y)

        y_out[~touched] = np.nan
        return x_out, y_out

    total_points = 0
    for dens in drd.densities:
        y = np.asarray(dens, dtype=np.float64)
        if nan_policy is NanPolicy.ZERO:
            total_points += int(y.size)
        else:
            total_points += int(np.isfinite(y).sum())

    if total_points == 0:
        return x_out, np.full(n_windows, np.nan, dtype=np.float64)

    x_all = np.empty(total_points, dtype=np.float64)
    y_all = np.empty(total_points, dtype=np.float64)
    offset = 0
    for pos, dens in zip(drd.positions, drd.densities, strict=True):
        x = np.asarray(pos, dtype=np.float64)
        y = np.asarray(dens, dtype=np.float64)
        if x.size == 0:
            continue
        if nan_policy is NanPolicy.ZERO:
            y = np.where(np.isfinite(y), y, 0.0)
        else:
            finite = np.isfinite(y)
            if not np.any(finite):
                continue
            x = x[finite]
            y = y[finite]
        size = int(y.size)
        x_all[offset: offset + size] = x
        y_all[offset: offset + size] = y
        offset += size

    x_all = x_all[:offset]
    y_all = y_all[:offset]
    order = np.argsort(x_all, kind="mergesort")
    y_out = _bin_points_windows_fast(
        x_all[order],
        y_all[order],
        n_windows=n_windows,
        agg=agg,
        nan_policy=NanPolicy.KEEP,
    )
    return x_out, y_out


@dataclass
class LinePlotComposer:
    segments: list[MetageneProfileSegment] | None = None
    n_windows: Optional[int] = None
    agg: AggMethod = field(default_factory=lambda: AggMethod.Mean)
    nan_policy: NanPolicy = NanPolicy.KEEP
    smooth: dict | int | None = 50
    title: Optional[str] = None
    width: int | None = None
    height: int | None = None
    _total_bins: int = field(init=False, repr=False)
    _smooth_cfg: object | None = field(init=False, repr=False, default=None)
    _n_windows_auto: bool = field(init=False, repr=False, default=False)
    x: list[np.ndarray] = field(default_factory=list)
    y: list[np.ndarray] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)

 
    borders: list[float] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.set_agg(self.agg)
        self.set_nan_policy(self.nan_policy)
        self._n_windows_auto = self.n_windows is None
        self.set_segments(self.segments)
        self.set_n_windows(None if self._n_windows_auto else self.n_windows)
        self.set_smooth(self.smooth)
        self.set_width(self.width)
        self.set_height(self.height)

    @staticmethod
    def _default_segments() -> list[MetageneProfileSegment]:
        return [MetageneProfileSegment("region", 100)]

    def _refresh_borders(self) -> None:
        total = float(self._total_bins)
        cum = 0
        self.borders = []
        for seg in self.segments:
            cum += seg.n_bins
            self.borders.append(cum / total)

    def _refresh_smooth_cfg(self) -> None:
        self._smooth_cfg = (
            _coerce_smooth_config(self.smooth, total_bins=self._total_bins)
            if self.smooth is not None
            else None
        )

    def set_segments(
        self,
        segments: list[MetageneProfileSegment] | None,
    ) -> LinePlotComposer:
        self.segments = list(self._default_segments() if segments is None else segments)
        self._total_bins = segments_total_bins(self.segments)
        if self._n_windows_auto:
            self.n_windows = self._total_bins
        self._refresh_borders()
        self._refresh_smooth_cfg()
        return self

    def set_n_windows(self, n_windows: Optional[int]) -> LinePlotComposer:
        self._n_windows_auto = n_windows is None
        self.n_windows = self._total_bins if n_windows is None else validate_n_windows(n_windows)
        return self

    def set_agg(self, agg: AggMethod) -> LinePlotComposer:
        self.agg = validate_window_agg(agg)
        return self

    def set_nan_policy(self, nan_policy: NanPolicy) -> LinePlotComposer:
        self.nan_policy = validate_nan_policy(nan_policy)
        return self

    def set_smooth(self, smooth: dict | int | None) -> LinePlotComposer:
        self.smooth = smooth
        self._refresh_smooth_cfg()
        return self

    def set_width(self, width: int | None) -> LinePlotComposer:
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    def set_height(self, height: int | None) -> LinePlotComposer:
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    def add_data(
        self,
        drd: DiscreteRegionData,
        *,
        name: str,
    ) -> LinePlotComposer:
        x_vals, y_vals = _line_profile(
            drd,
            segments=self.segments,
            n_windows=self.n_windows,
            nan_policy=self.nan_policy,
            agg=self.agg,
        )


        if self.smooth is not None and y_vals.size > 0 and self.agg in (AggMethod.Min, AggMethod.Max):
            warnings.warn(
                f"smooth={self.smooth!r} ignored for agg={self.agg!r}; "
                "smoothing is only applied to mean/median",
                RuntimeWarning,
                stacklevel=2,
            )

        if self._smooth_cfg is not None and y_vals.size > 0 and self.agg in (AggMethod.Mean, AggMethod.Median):
            y_scaled = y_vals.astype(float, copy=True)
            y_scaled = _apply_savgol_smoothing(y_scaled, self._smooth_cfg, segments=self.segments)
            y_scaled = _clip_profile(y_scaled)
            y_vals = y_scaled

        self.x.append(x_vals)
        self.y.append(y_vals)
        self.labels.append(name)
        return self

    def set_label(self, curve_index: int, label: str) -> LinePlotComposer:
        self.labels[curve_index] = label
        return self

    def finish(self):
        _hv_init()

        if not self.x:
            curve = hv.Curve([])
            return curve.opts(
                xlabel=_position_axis_label(),
                ylabel=f"{_agg_label(self.agg)} methylation density",
                show_legend=False,
                title=self.title or _profile_title(),
                **({} if self.width is None else {"width": int(self.width)}),
                **({} if self.height is None else {"height": int(self.height)}),
            )

        curves = []
        for x_vals, y_vals, label in zip(self.x, self.y, self.labels):
            if x_vals.size == 0 or y_vals.size == 0:
                c = hv.Curve([]).relabel(label)
            else:
                c = hv.Curve((x_vals, y_vals), kdims="relative position", vdims="density").relabel(label)
            curves.append(c)

        plot = curves[0]
        for c in curves[1:]:
            plot *= c  # Overlay

        opts_kwargs = dict(
            xlabel=_position_axis_label(),
            ylabel=f"{_agg_label(self.agg)} methylation density",
            show_legend=len(curves) > 1,
            title=self.title or _profile_title(),
        )
        if self.width is not None:
            opts_kwargs["width"] = int(self.width)
        if self.height is not None:
            opts_kwargs["height"] = int(self.height)

        return plot.opts(**opts_kwargs)


def line_plot(
    drd: DiscreteRegionData,
    *,
    name: str = "sample",
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    agg: AggMethod = AggMethod.Mean,
    nan_policy: NanPolicy = NanPolicy.KEEP,
    smooth: dict | int | None = 50,
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
):
    """
    Build a HoloViews normalized-profile line plot from precomputed discrete regions.

    Parameters
    ----------
    drd
        Discrete normalized profiles to render.
    name
        Dataset label used in the legend.
    segments
        Optional normalized-profile segment layout used for binning.
    n_windows
        Number of output windows. Defaults to the total segment bin count.
    agg
        Aggregation method used while rebinding values into display windows.
    nan_policy
        Policy controlling how NaN values are handled during windowing.
    smooth
        Optional Savitzky-Golay smoothing configuration.
    title
        Optional plot title.
    width, height
        Optional plot size in pixels.
    """
    composer = LinePlotComposer(
        segments=segments,
        n_windows=n_windows,
        agg=agg,
        nan_policy=nan_policy,
        smooth=smooth,
        title=title,
        width=width,
        height=height,
    )
    return composer.add_data(drd, name=name).finish()
