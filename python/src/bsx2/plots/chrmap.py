from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional

import numpy as np
import polars as pl
from beartype import beartype
from scipy.signal import savgol_filter


def _smooth_series(values: np.ndarray, window: int) -> np.ndarray:
    """Savitzky-Golay smoothing; if window <= 1 or too short, returns input."""
    if window is None or window <= 1:
        return values
    n = values.shape[0]
    if n < 3:
        return values
    window = int(window)
    if window % 2 == 0:
        window += 1  # savgol_filter needs an odd window length
    # cap to available points while keeping window odd
    window = min(window, n if n % 2 == 1 else n - 1)
    if window < 3:
        return values
    polyorder = min(3, window - 1)
    return savgol_filter(values, window_length=window, polyorder=polyorder, mode="interp")


@dataclass(frozen=True)
class ChrLineData:
    """Prepared line data for chromosome methylation map."""

    x: np.ndarray
    y: np.ndarray
    x_ticks: List[int]
    x_labels: List[str]
    borders: np.ndarray
    lower: Optional[np.ndarray] = None
    upper: Optional[np.ndarray] = None


@dataclass(frozen=True)
class ChrBoxData:
    """Prepared box/violin data per chromosome."""

    labels: List[str]
    values: List[np.ndarray]  # per-chromosome densities


@beartype
def prepare_chr_line_data(
    df: pl.DataFrame,
    *,
    smooth: int = 0,
) -> ChrLineData:
    """Convert window-level report into line-plot-ready data (percent density)."""
    required = {"chr", "window", "sum", "count"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    has_ci = "upper" in df.columns and "lower" in df.columns

    agg_exprs = [
        pl.col("sum").sum().alias("sum"),
        pl.col("count").sum().alias("count"),
    ]
    if has_ci:
        agg_exprs.extend(
            [
                pl.col("upper").mean().alias("upper"),
                pl.col("lower").mean().alias("lower"),
            ]
        )

    agg_df = (
        df.sort(["chr", "window"])
        .group_by(["chr", "window"], maintain_order=True)
        .agg(agg_exprs)
        .with_columns(
            (pl.col("sum") / pl.col("count")).alias("density")
        )
        .filter(pl.col("count") > 0)
    )

    # x-axis: sequential index over windows
    y = agg_df["density"].to_numpy() * 100.0
    x = np.arange(len(y), dtype=int)

    lower = upper = None
    if has_ci:
        lower = agg_df["lower"].to_numpy() * 100.0
        upper = agg_df["upper"].to_numpy() * 100.0

    if smooth and smooth > 1:
        y = _smooth_series(y, smooth)
        if has_ci:
            lower = _smooth_series(lower, smooth)
            upper = _smooth_series(upper, smooth)

    # ticks and borders per chromosome
    # Build mapping chr -> first index
    chr_order = agg_df["chr"].to_list()
    chr_unique: List[str] = []
    first_idx: List[int] = []
    seen = set()
    for i, c in enumerate(chr_order):
        if c not in seen:
            seen.add(c)
            chr_unique.append(c)
            first_idx.append(i)

    borders = np.array(first_idx + [len(x)], dtype=int)
    x_ticks = [int((borders[i] + borders[i + 1]) // 2) for i in range(len(borders) - 1)]

    return ChrLineData(
        x=x,
        y=y,
        x_ticks=x_ticks,
        x_labels=chr_unique,
        borders=borders,
        lower=lower,
        upper=upper,
    )


@beartype
def prepare_chr_box_data(df: pl.DataFrame) -> ChrBoxData:
    """Prepare per-chromosome density distributions for box/violin plots."""
    required = {"chr", "sum", "count"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    dens_df = (
        df.sort("chr")
        .group_by("chr", maintain_order=True)
        .agg((pl.col("sum") / pl.col("count")).alias("density"))
    )
    labels = dens_df["chr"].to_list()
    values = [np.array(d, dtype=float) * 100.0 for d in dens_df["density"].to_list()]
    return ChrBoxData(labels=labels, values=values)
