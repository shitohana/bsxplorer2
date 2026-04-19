from __future__ import annotations

from beartype import beartype
from beartype.typing import Optional

from bsx2 import Context
from bsx2.validation import ChrEmptyPolicy, ChrLineStat

from ..compute.chrline import ChrLineData
from ..compute.chrmap import compute_chromosome_methylation_map_data
from .chrline import ChrLinePlotComposer

ChromosomeMethylationMapComposer = ChrLinePlotComposer


@beartype
def build_chromosome_methylation_map(
    data: ChrLineData,
    *,
    name: str = "sample",
    title: Optional[str] = None,
    width: int | None = None,
    height: int | None = None,
    show_chr_boundaries: bool = True,
    show_chr_labels: bool = True,
    x_label: str = "Genomic coordinate (bp)",
    y_label: Optional[str] = None,
):
    """
    Build a HoloViews chromosome methylation map from precomputed data.
    """
    composer = ChrLinePlotComposer(
        title=title,
        width=width,
        height=height,
        show_chr_boundaries=show_chr_boundaries,
        show_chr_labels=show_chr_labels,
        x_label=x_label,
        y_label=y_label,
    )
    return composer.add_data(data, name=name).finish()


@beartype
def chromosome_methylation_map(
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
    show_chr_boundaries: bool = True,
    show_chr_labels: bool = True,
    x_label: str = "Genomic coordinate (bp)",
    y_label: Optional[str] = None,
):
    """
    Compute and build a chromosome methylation map in one step.
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
    return build_chromosome_methylation_map(
        data,
        name=name,
        title=title,
        width=width,
        height=height,
        show_chr_boundaries=show_chr_boundaries,
        show_chr_labels=show_chr_labels,
        x_label=x_label,
        y_label=y_label,
    )


__all__ = [
    "ChrLinePlotComposer",
    "ChromosomeMethylationMapComposer",
    "build_chromosome_methylation_map",
    "chromosome_methylation_map",
]
