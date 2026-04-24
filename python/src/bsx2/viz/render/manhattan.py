from __future__ import annotations

from dataclasses import dataclass, field

import holoviews as hv
import numpy as np
from beartype.typing import Optional

from bsx2 import Context
from bsx2.validation import ChrEmptyPolicy, ChrLineStat, validate_positive_int

from ..compute.chrline import ChrLineData
from ..compute.chrmap import compute_chromosome_methylation_map_data
from ._common import _hv_init


def _default_palette() -> tuple[str, ...]:
    return (
        "#2563eb",
        "#0f766e",
        "#7c3aed",
        "#059669",
        "#dc2626",
        "#0891b2",
    )


@dataclass
class ChromosomeManhattanPlotComposer:
    title: Optional[str] = None
    width: int | None = None
    height: int | None = None
    show_chr_boundaries: bool = True
    show_chr_labels: bool = True
    x_label: str = "Genomic coordinate (bp)"
    y_label: Optional[str] = None
    point_size: int = 5
    palette: tuple[str, ...] = field(default_factory=_default_palette)
    datasets: list[ChrLineData] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)

    def set_width(self, width: int | None) -> ChromosomeManhattanPlotComposer:
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    def set_height(self, height: int | None) -> ChromosomeManhattanPlotComposer:
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    def set_point_size(self, point_size: int) -> ChromosomeManhattanPlotComposer:
        self.point_size = validate_positive_int(point_size, name="point_size")
        return self

    def add_data(
        self,
        data: ChrLineData,
        *,
        name: str = "sample",
    ) -> ChromosomeManhattanPlotComposer:
        self.datasets.append(data)
        self.labels.append(name)
        return self

    def finish(self):
        _hv_init()
        if not self.datasets:
            return hv.Scatter([]).opts(
                title=self.title or "Chromosome methylation Manhattan plot",
                xlabel=self.x_label,
                ylabel=self.y_label or "Methylation",
                show_legend=False,
                **({} if self.width is None else {"width": int(self.width)}),
                **({} if self.height is None else {"height": int(self.height)}),
            )

        overlay = None
        first = self.datasets[0]
        multi_dataset = len(self.datasets) > 1

        for data_index, (data, label) in enumerate(zip(self.datasets, self.labels, strict=True)):
            for track_index, track in enumerate(data.tracks):
                mask = np.isfinite(track.x_global_bp) & np.isfinite(track.value)
                if not np.any(mask):
                    continue

                color = self.palette[(data_index + track_index) % len(self.palette)]
                track_label = f"{label} · {track.seqname}" if multi_dataset else track.seqname
                points = hv.Scatter(
                    (track.x_global_bp[mask], track.value[mask]),
                    kdims=["genomic bp"],
                    vdims=["methylation"],
                ).relabel(track_label).opts(
                    color=color,
                    size=int(self.point_size),
                    marker="circle",
                    alpha=0.82,
                    show_legend=False,
                )
                overlay = points if overlay is None else overlay * points

        if overlay is None:
            overlay = hv.Scatter([])

        if self.show_chr_boundaries:
            for boundary in first.chr_boundaries_bp():
                overlay *= hv.VLine(float(boundary)).opts(
                    line_dash="dash",
                    line_color="gray",
                    line_width=1,
                )

        opts_kwargs = dict(
            title=self.title or f"Chromosome methylation Manhattan plot ({first.context.name})",
            xlabel=self.x_label,
            ylabel=self.y_label
            or ("Mean methylation" if first.stat is ChrLineStat.MEAN else "Weighted methylation"),
            show_legend=False,
        )
        if self.width is not None:
            opts_kwargs["width"] = int(self.width)
        if self.height is not None:
            opts_kwargs["height"] = int(self.height)
        if self.show_chr_labels and first.tracks:
            opts_kwargs["xticks"] = list(zip(first.chr_centers_bp(), first.chromosomes(), strict=True))
        return overlay.opts(**opts_kwargs)


@dataclass
class ChromosomeManhattanTrack:
    seqname: str
    x_global_bp: np.ndarray
    value: np.ndarray


@dataclass
class ChromosomeManhattanData:
    data: ChrLineData


@dataclass
class ManhattanMethylationPlotComposer(ChromosomeManhattanPlotComposer):
    pass


def build_chromosome_manhattan_plot(
    data: ChrLineData,
    *,
    name: str = "sample",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
    point_size: int = 5,
    show_chr_boundaries: bool = True,
    show_chr_labels: bool = True,
    x_label: str = "Genomic coordinate (bp)",
    y_label: Optional[str] = None,
):
    """
    Build a Manhattan-style chromosome methylation plot from precomputed data.
    """
    composer = ChromosomeManhattanPlotComposer(
        title=title,
        width=width,
        height=height,
        point_size=point_size,
        show_chr_boundaries=show_chr_boundaries,
        show_chr_labels=show_chr_labels,
        x_label=x_label,
        y_label=y_label,
    )
    return composer.add_data(data, name=name).finish()


def chromosome_manhattan_plot(
    reader: object,
    *,
    context: Context,
    bin_size_bp: int = 50_000,
    chr_lengths: dict[str, int] | None = None,
    min_coverage: int | float = 0,
    stat: ChrLineStat | str = ChrLineStat.WEIGHTED_MEAN,
    smooth: dict | int | None = None,
    empty_policy: ChrEmptyPolicy | str = ChrEmptyPolicy.NAN,
    name: str = "sample",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
    point_size: int = 5,
    show_chr_boundaries: bool = True,
    show_chr_labels: bool = True,
    x_label: str = "Genomic coordinate (bp)",
    y_label: Optional[str] = None,
):
    """
    Compute and build a Manhattan-style chromosome methylation plot in one step.
    """
    data = compute_chromosome_methylation_map_data(
        reader,
        context=context,
        bin_size_bp=bin_size_bp,
        chr_lengths=chr_lengths,
        min_coverage=min_coverage,
        stat=stat,
        smooth=smooth,
        empty_policy=empty_policy,
    )
    return build_chromosome_manhattan_plot(
        data,
        name=name,
        title=title,
        width=width,
        height=height,
        point_size=point_size,
        show_chr_boundaries=show_chr_boundaries,
        show_chr_labels=show_chr_labels,
        x_label=x_label,
        y_label=y_label,
    )


__all__ = [
    "ChromosomeManhattanPlotComposer",
    "ManhattanMethylationPlotComposer",
    "build_chromosome_manhattan_plot",
    "chromosome_manhattan_plot",
]
