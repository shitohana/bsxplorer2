from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import plotly.graph_objects as go

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    MetageneProfileSegment,
    segment_boundaries,
    segments_total_bins,
)
from bsx2.validation import (
    validate_n_windows,
    validate_positive_int,
    validate_rank_score,
    validate_sort_order,
)
from ._common import (
    _bin_points_windows_fast,
    _rank_compress,
    _segment_decor_bin,
    NanPolicy,
)


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

    rows: list[np.ndarray] = []
    labels: list[str] = []

    for pos, dens, lbl in zip(drd.positions, drd.densities, drd.labels):
        x = np.asarray(pos, dtype=np.float64)
        y = np.asarray(dens, dtype=np.float64)

        row = _bin_points_windows_fast(
            x,
            y,
            n_windows=n_windows,
            agg="mean",
            nan_policy=nan_policy,
        )
        rows.append(row)
        labels.append(lbl if lbl is not None else f"region_{len(labels) + 1}")

    if not rows:
        return np.empty((0, 0), dtype=np.float64), [], np.empty((0,), dtype=np.float64)

    mat = np.vstack(rows)
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
            self.segments = [
                MetageneProfileSegment("up", 100),
                MetageneProfileSegment("body", 200),
                MetageneProfileSegment("down", 100),
            ]

        self._total_bins = segments_total_bins(self.segments)

        if self.n_windows is None:
            self.n_windows = self._total_bins
        else:
            self.n_windows = validate_n_windows(self.n_windows)

        self._bins_rel = (np.arange(self.n_windows, dtype=np.float64) + 0.5) / float(self.n_windows)

        if self.rank_score == "body_mean":
            bounds = segment_boundaries(self.segments)  # relative boundaries
            b0, b1 = bounds[0], bounds[1]
            self._rank_mask = (self._bins_rel >= b0) & (self._bins_rel < b1)
        else:
            self._rank_mask = np.ones(self.n_windows, dtype=bool)

    def add_data(
        self,
        drd: DiscreteRegionData,
        *,
        label_prefix: str | None = None,
    ) -> "HeatmapPlotComposer":
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

    def finish(self) -> go.Figure:
        if not self.z_parts:
            return go.Figure()

        z = np.vstack(self.z_parts)

        # Ranking score (preserve existing behavior)
        mask = self._rank_mask
        denom = float(np.sum(mask) + 1.0)
        scores = np.nansum(z[:, mask], axis=1) / denom

        order_idx = np.argsort(scores)
        if self.sort_order == "desc":
            order_idx = order_idx[::-1]
        z = z[order_idx]

        z = _rank_compress(z, self.rank_rows, fill=self.empty_bin_fill)
        regions = [str(i) for i in range(z.shape[0])]

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

        x_plot = np.arange(n_bins, dtype=int)

        fig = go.Figure(
            data=go.Heatmap(
                z=z_vis,
                x=x_plot,
                y=regions,
                colorscale=self.colorscale,
                zmin=zmin,
                zmax=zmax,
                zsmooth=False,
                xgap=0,
                ygap=0,
                hoverongaps=False,
                colorbar=dict(title="Methylation density"),
            )
        )
        fig.update_yaxes(autorange="reversed", showticklabels=False)
        _segment_decor_bin(fig, self.segments, annotate_tss_tes=True)
        fig.update_layout(
            title=self.title or "Metagene profile - Heatmap (BSX1 ranked)",
            width=fig_width,
            height=fig_height,
            margin=dict(l=70, r=30, t=60, b=70),
            xaxis_title="Position (bin)",
            yaxis_title="Rank",
        )
        return fig

    def to_html(self, *, full_html: bool = False, include_js: str = "cdn") -> str:
        return self.finish().to_html(full_html=full_html, include_plotlyjs=include_js)