from __future__ import annotations

"""Compute-only chromosome methylation map API built on top of `chrline`."""

from beartype import beartype
from beartype.typing import Mapping

from bsx2 import Context
from bsx2.validation import ChrEmptyPolicy, ChrLineStat

from .chrline import ChrLineData, ChrLineTrack, compute_chr_line_data

ChromosomeMethylationTrack = ChrLineTrack
ChromosomeMethylationMapData = ChrLineData


@beartype
def compute_chromosome_methylation_map_data(
    reader: object,
    *,
    context: Context,
    bin_size_bp: int = 50_000,
    chr_lengths: Mapping[str, int] | None = None,
    min_coverage: int | float = 0,
    stat: ChrLineStat | str = ChrLineStat.WEIGHTED_MEAN,
    smooth: dict | int | None = None,
    empty_policy: ChrEmptyPolicy | str = ChrEmptyPolicy.NAN,
) -> ChrLineData:
    """
    Compute chromosome-wide methylation map data.

    Parameters
    ----------
    reader
        `BsxFileReader`-like iterable or `RegionReader`-like object.
    context
        Single methylation context filter.
    bin_size_bp
        Window size in base pairs.
    chr_lengths
        Optional real chromosome lengths keyed by sequence name.
    min_coverage
        Minimum per-site coverage threshold.
    stat
        Aggregation statistic for each genomic window.
    smooth
        Optional smoothing configuration passed through to `chrline`.
    empty_policy
        Policy for empty genomic windows.
    """
    return compute_chr_line_data(
        reader,
        context=context,
        bin_size_bp=bin_size_bp,
        chr_lengths=None if chr_lengths is None else dict(chr_lengths),
        min_coverage=min_coverage,
        stat=stat,
        smooth=smooth,
        empty_policy=empty_policy,
    )


__all__ = [
    "ChrLineData",
    "ChrLineTrack",
    "ChromosomeMethylationMapData",
    "ChromosomeMethylationTrack",
    "compute_chromosome_methylation_map_data",
]
