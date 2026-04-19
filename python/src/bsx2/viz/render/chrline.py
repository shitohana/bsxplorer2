from __future__ import annotations

from dataclasses import dataclass, field

import holoviews as hv
from beartype.typing import Optional

from bsx2.validation import ChrLineStat, validate_positive_int

from ..compute.chrline import ChrLineData
from ._common import _hv_init


@dataclass
class ChrLinePlotComposer:
    title: Optional[str] = None
    width: int | None = None
    height: int | None = None
    show_chr_boundaries: bool = True
    show_chr_labels: bool = True
    x_label: str = "Genomic coordinate (bp)"
    y_label: Optional[str] = None
    datasets: list[ChrLineData] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)

    def set_width(self, width: int | None) -> ChrLinePlotComposer:
        self.width = None if width is None else validate_positive_int(width, name="width")
        return self

    def set_height(self, height: int | None) -> ChrLinePlotComposer:
        self.height = None if height is None else validate_positive_int(height, name="height")
        return self

    def add_data(
        self,
        data: ChrLineData,
        *,
        name: str = "sample",
    ) -> ChrLinePlotComposer:
        self.datasets.append(data)
        self.labels.append(name)
        return self

    def finish(self):
        _hv_init()
        if not self.datasets:
            return hv.Curve([]).opts(
                title=self.title or "Chromosome-wide methylation profile",
                xlabel=self.x_label,
                ylabel=self.y_label or "Methylation",
                show_legend=False,
                **({} if self.width is None else {"width": int(self.width)}),
                **({} if self.height is None else {"height": int(self.height)}),
            )

        curves = []
        for data, label in zip(self.datasets, self.labels, strict=True):
            x_vals, y_vals = data.global_curve()
            curve = hv.Curve(
                (x_vals, y_vals),
                kdims=["genomic bp"],
                vdims=["methylation"],
            ).relabel(label)
            curves.append(curve)

        plot = curves[0]
        for curve in curves[1:]:
            plot *= curve

        first = self.datasets[0]
        if self.show_chr_boundaries:
            for boundary in first.chr_boundaries_bp():
                plot *= hv.VLine(float(boundary)).opts(
                    line_dash="dash",
                    line_color="gray",
                    line_width=1,
                )

        opts_kwargs = dict(
            title=self.title or f"Chromosome-wide methylation profile ({first.context.name})",
            xlabel=self.x_label,
            ylabel=self.y_label
            or ("Mean methylation" if first.stat is ChrLineStat.MEAN else "Weighted methylation"),
            show_legend=len(curves) > 1,
        )
        if self.width is not None:
            opts_kwargs["width"] = int(self.width)
        if self.height is not None:
            opts_kwargs["height"] = int(self.height)
        if self.show_chr_labels and first.tracks:
            opts_kwargs["xticks"] = list(zip(first.chr_centers_bp(), first.chromosomes(), strict=True))
        return plot.opts(**opts_kwargs)


__all__ = ["ChrLinePlotComposer"]
