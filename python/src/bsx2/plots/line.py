from __future__ import annotations

from dataclasses import dataclass, field
from beartype.typing import Optional
import warnings
from bsx2 import AggMethod
import holoviews as hv
import numpy as np

from bsx2.plots.data import DiscreteRegionData
from bsx2.validation import (
    validate_n_windows,
    validate_nan_policy,
    validate_positive_int,
    validate_window_agg,
)
from bsx2.plots.metagene import (
    MetageneProfileSegment,
    _apply_savgol_smoothing,
    _clip_profile,
    _coerce_smooth_config,
    segments_total_bins,
)
from ._common import _bin_points_windows_fast, _hv_init, NanPolicy


def _agg_label(agg: AggMethod) -> str:
    name = getattr(agg, "name", None)
    if isinstance(name, str) and name:
        token = name
    else:
        token = str(agg).rsplit(".", 1)[-1]
    return token.lower().capitalize()


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

    xs_parts: list[np.ndarray] = []
    ys_parts: list[np.ndarray] = []

    for pos, dens in zip(drd.positions, drd.densities):
        x = np.asarray(pos, dtype=np.float64)
        y = np.asarray(dens, dtype=np.float64)

        if x.size == 0:
            continue

        xs_parts.append(x)
        ys_parts.append(y)

    if not xs_parts:
        empty = np.array([], dtype=np.float64)
        return empty, empty

    x_all = np.concatenate(xs_parts)
    y_all = np.concatenate(ys_parts)

    if agg is AggMethod.Median:
        order = np.argsort(x_all, kind="mergesort")
        x_all = x_all[order]
        y_all = y_all[order]

    y_out = _bin_points_windows_fast(
        x_all,
        y_all,
        n_windows=n_windows,
        agg=agg,
        nan_policy=nan_policy,
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
        return [
            MetageneProfileSegment("up", 100),
            MetageneProfileSegment("body", 200),
            MetageneProfileSegment("down", 100),
        ]

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
    ) -> "LinePlotComposer":
        self.segments = list(self._default_segments() if segments is None else segments)
        self._total_bins = segments_total_bins(self.segments)
        if self._n_windows_auto:
            self.n_windows = self._total_bins
        self._refresh_borders()
        self._refresh_smooth_cfg()
        return self

    def set_n_windows(self, n_windows: Optional[int]) -> "LinePlotComposer":
        self._n_windows_auto = n_windows is None
        self.n_windows = self._total_bins if n_windows is None else validate_n_windows(n_windows)
        return self

    def set_agg(self, agg: AggMethod) -> "LinePlotComposer":
        self.agg = validate_window_agg(agg)
        return self

    def set_nan_policy(self, nan_policy: NanPolicy) -> "LinePlotComposer":
        self.nan_policy = validate_nan_policy(nan_policy)
        return self

    def set_smooth(self, smooth: dict | int | None) -> "LinePlotComposer":
        self.smooth = smooth
        self._refresh_smooth_cfg()
        return self

    def set_width(self, width: int | None) -> "LinePlotComposer":
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    def set_height(self, height: int | None) -> "LinePlotComposer":
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    def add_data(
        self,
        drd: DiscreteRegionData,
        *,
        name: str,
    ) -> "LinePlotComposer":
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

    def set_label(self, curve_index: int, label: str) -> "LinePlotComposer":
        self.labels[curve_index] = label
        return self

    def finish(self):
        _hv_init()

        if not self.x:
            curve = hv.Curve([])
            return curve.opts(
                xlabel="Metagene position (relative)",
                ylabel=f"{_agg_label(self.agg)} methylation density",
                show_legend=False,
                title=self.title or "Metagene profile - Line",
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
            xlabel="Metagene position (relative)",
            ylabel=f"{_agg_label(self.agg)} methylation density",
            show_legend=len(curves) > 1,
            title=self.title or "Metagene profile - Line",
        )
        if self.width is not None:
            opts_kwargs["width"] = int(self.width)
        if self.height is not None:
            opts_kwargs["height"] = int(self.height)

        return plot.opts(**opts_kwargs)
