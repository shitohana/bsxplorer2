from __future__ import annotations

from typing import Iterable

import holoviews as hv
import numpy as np

from bsx2.plots.chrmap import ChrLineData, ChrBoxData


def chr_line_hv(line: ChrLineData, *, label: str | None = None) -> hv.Overlay:
    """Holoviews line + CI + chromosome borders. Returns an Overlay."""
    curves = [hv.Curve((line.x, line.y), label=label or "density").opts(ylabel="Methylation density (%)")]
    if line.lower is not None and line.upper is not None:
        area = hv.Area((line.x, line.lower, line.upper)).opts(alpha=0.2, color="lightgray")
        curves.append(area)
    if len(line.borders) > 0:
        borders = [hv.VLine(int(b)) for b in line.borders]
        curves.extend(borders)
    overlay = hv.Overlay(curves)
    if line.x_ticks and line.x_labels:
        overlay = overlay.opts(
            xaxis="bottom",
            xticks=list(zip(line.x_ticks, line.x_labels)),
            xlabel="Chromosome",
        )
    return overlay


def chr_box_hv(box: ChrBoxData, *, kind: str = "box") -> hv.Element:
    """Holoviews BoxWhisker/Violin for per-chromosome densities."""
    data = {"chr": np.repeat(box.labels, [len(v) for v in box.values]), "density": np.concatenate(box.values)}
    if kind == "violin":
        return hv.Violin(data, kdims=["chr"], vdims=["density"]).opts(yticks=5, xlabel="Chromosome", ylabel="Methylation density (%)")
    return hv.BoxWhisker(data, kdims=["chr"], vdims=["density"]).opts(yticks=5, xlabel="Chromosome", ylabel="Methylation density (%)")
