from __future__ import annotations

from dataclasses import dataclass, field

import holoviews as hv
import numpy as np
from beartype.typing import Optional

from bsx2 import AggMethod
from bsx2.validation import (
    validate_n_windows,
    validate_positive_int,
    validate_rank_score,
    validate_sort_order,
)

from ..compute.data import DiscreteRegionData
from ..compute.metagene import (
    MetageneProfileSegment,
    segment_boundaries,
    segments_total_bins,
)
from ..compute.windowing import _bin_points_windows_fast, _rank_compress
from ._common import NanPolicy, _hv_init


def _profile_title() -> str:
    return "Scaled region profile - Heatmap"


def _position_axis_label() -> str:
    return "Relative feature position"


def _default_segments() -> list[MetageneProfileSegment]:
    return [
        MetageneProfileSegment("up", 100),
        MetageneProfileSegment("body", 200),
        MetageneProfileSegment("down", 100),
    ]


def _show_segment_guides(segments: list[MetageneProfileSegment]) -> bool:
    return len(segments) > 1


def _heatmap_matrix(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    nan_policy: NanPolicy = NanPolicy.KEEP,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    else:
        n_windows = validate_n_windows(n_windows)

    n_regions = len(drd.positions)
    if n_regions == 0:
        return np.empty((0, 0), dtype=np.float64), [], np.empty((0,), dtype=np.float64)

    mat = np.empty((n_regions, n_windows), dtype=np.float64)
    labels = [f"region_{idx + 1}" for idx in range(n_regions)]

    for row_idx, (pos, dens, lbl) in enumerate(zip(drd.positions, drd.densities, drd.labels, strict=True)):
        x = np.asarray(pos, dtype=np.float64)
        y = np.asarray(dens, dtype=np.float64)

        mat[row_idx] = _bin_points_windows_fast(
            x,
            y,
            n_windows=n_windows,
            agg=AggMethod.Mean,
            nan_policy=nan_policy,
        )
        if lbl is not None:
            labels[row_idx] = lbl

    bins_rel = (np.arange(n_windows, dtype=np.float64) + 0.5) / float(n_windows)
    return mat, labels, bins_rel


@dataclass
class HeatmapPlotComposer:
    segments: list[MetageneProfileSegment] | None = None
    n_windows: int | None = None
    rank_rows: int = 100
    rank_score: str = "mean"      # "mean" | "body_mean"
    sort_order: str = "desc"      # "asc" | "desc"
    colorscale: str = "Viridis"
    empty_bin_fill: float = 0.0
    nan_policy: NanPolicy = NanPolicy.KEEP
    title: Optional[str] = None
    width: int | None = None
    height: int | None = None
    vmax_q: float | None = 0.995

    # accumulated raw rows
    z_parts: list[np.ndarray] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)

    # cached internals
    _total_bins: int = field(init=False, repr=False)
    _bins_rel: np.ndarray = field(init=False, repr=False)
    _rank_mask: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.rank_score = validate_rank_score(self.rank_score)
        self.sort_order = validate_sort_order(self.sort_order)
        self.rank_rows = validate_positive_int(self.rank_rows, name="rank_rows")

        if self.segments is None:
            self.segments = _default_segments()

        self._total_bins = segments_total_bins(self.segments)

        if self.n_windows is None:
            self.n_windows = self._total_bins
        else:
            self.n_windows = validate_n_windows(self.n_windows)

        self._bins_rel = (np.arange(self.n_windows, dtype=np.float64) + 0.5) / float(self.n_windows)

        if self.rank_score == "body_mean" and len(self.segments) == 3:
            bounds = segment_boundaries(self.segments)
            b0, b1 = bounds[0], bounds[1]
            self._rank_mask = (self._bins_rel >= b0) & (self._bins_rel < b1)
        else:
            self._rank_mask = np.ones(self.n_windows, dtype=bool)

    def add_data(
        self,
        drd: DiscreteRegionData,
        *,
        label_prefix: str | None = None,
    ) -> HeatmapPlotComposer:
        z, row_labels, _ = _heatmap_matrix(
            drd,
            segments=self.segments,
            n_windows=self.n_windows,
            nan_policy=self.nan_policy,
        )

        if z.size == 0:
            return self

        if label_prefix is not None:
            row_labels = [f"{label_prefix}:{lbl}" for lbl in row_labels]

        self.z_parts.append(z)
        self.labels.extend(row_labels)
        return self

    def finish(self):
        _hv_init()
        if not self.z_parts:
            return hv.Image(
                np.full((1, 1), np.nan, dtype=np.float64),
                bounds=(0.0, -0.5, 1.0, 0.5),
                kdims=["Relative feature position", "Rank"],
                vdims=["Methylation density"],
            ).opts(
                cmap=self.colorscale,
                colorbar=True,
                clim=(0.0, 1.0),
                invert_yaxis=True,
                xlabel=_position_axis_label(),
                ylabel="Rank",
                title=self.title or _profile_title(),
                **({} if self.width is None else {"width": int(self.width)}),
                **({} if self.height is None else {"height": int(self.height)}),
            )

        z = self.z_parts[0] if len(self.z_parts) == 1 else np.vstack(self.z_parts)

        # Ranking score (preserve existing behavior)
        mask = self._rank_mask
        denom = float(np.sum(mask) + 1.0)
        scores = np.nansum(z[:, mask], axis=1) / denom

        order_idx = np.argsort(scores)
        if self.sort_order == "desc":
            order_idx = order_idx[::-1]
        z = z[order_idx]

        z = _rank_compress(z, self.rank_rows, fill=self.empty_bin_fill)
        z_vis = np.asarray(z, dtype=np.float64)

        zmin = 0.0
        zmax = 1.0
        if self.vmax_q is not None:
            vals = z_vis[np.isfinite(z_vis)]
            if vals.size:
                zmax = float(np.nanquantile(vals, self.vmax_q))
                if not np.isfinite(zmax) or zmax <= zmin:
                    zmax = 1.0

        n_bins = int(z_vis.shape[1])
        fig_width = self.width if self.width is not None else max(900, min(1400, 2 * n_bins))
        fig_height = self.height if self.height is not None else max(550, min(900, 5 * int(z_vis.shape[0])))

        heatmap = hv.Image(
            z_vis,
            bounds=(0.0, -0.5, 1.0, z_vis.shape[0] - 0.5),
            kdims=["Relative feature position", "Rank"],
            vdims=["Methylation density"],
        ).opts(
            cmap=self.colorscale,
            colorbar=True,
            clim=(zmin, zmax),
            invert_yaxis=True,
            xlabel=_position_axis_label(),
            ylabel="Rank",
            title=self.title or _profile_title(),
            width=fig_width,
            height=fig_height,
        )

        if not _show_segment_guides(self.segments):
            return heatmap

        xticks = []
        start = 0.0
        for end, seg in zip(segment_boundaries(self.segments), self.segments, strict=True):
            xticks.append(((start + end) / 2.0, seg.name))
            start = end
        plot = heatmap.opts(xticks=xticks)
        for boundary in segment_boundaries(self.segments)[:-1]:
            plot *= hv.VLine(float(boundary)).opts(
                line_dash="dash",
                line_color="gray",
                line_width=1,
            )
        return plot

def heatmap(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: int | None = None,
    rank_rows: int = 100,
    rank_score: str = "mean",
    sort_order: str = "desc",
    colorscale: str = "Viridis",
    empty_bin_fill: float = 0.0,
    nan_policy: NanPolicy = NanPolicy.KEEP,
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
    vmax_q: float | None = 0.995,
    label_prefix: str | None = None,
):
    """
    Build a ranked HoloViews heatmap from precomputed discrete regions.

    Parameters
    ----------
    drd
        Discrete normalized profiles to render.
    segments
        Optional normalized-profile segment layout used for binning.
    n_windows
        Number of output windows. Defaults to the total segment bin count.
    rank_rows
        Number of output rows after rank compression.
    rank_score
        Row ranking strategy. Supported values are `"mean"` and `"body_mean"`.
        `"body_mean"` is only meaningful for explicit three-part layouts.
    sort_order
        Sort direction for ranked rows.
    colorscale
        HoloViews colormap name.
    empty_bin_fill
        Fill value used after rank compression for empty bins.
    nan_policy
        Policy controlling how NaN values are handled during windowing.
    title
        Optional plot title.
    width, height
        Optional plot size in pixels.
    vmax_q
        Optional upper quantile used to derive the heatmap color limit.
    label_prefix
        Optional prefix applied to row labels before combining datasets.

    Returns
    -------
    object
        HoloViews heatmap-like object that can be rendered with the Plotly
        backend.

    Notes
    -----
    The heatmap is built by rebinding each discrete profile to a shared window
    grid, ranking rows, and optionally compressing them to ``rank_rows``.
    """
    composer = HeatmapPlotComposer(
        segments=segments,
        n_windows=n_windows,
        rank_rows=rank_rows,
        rank_score=rank_score,
        sort_order=sort_order,
        colorscale=colorscale,
        empty_bin_fill=empty_bin_fill,
        nan_policy=nan_policy,
        title=title,
        width=width,
        height=height,
        vmax_q=vmax_q,
    )
    return composer.add_data(drd, label_prefix=label_prefix).finish()
