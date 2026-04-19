from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from beartype.typing import Optional

from bsx2 import AggMethod
from bsx2.validation import NanPolicy, validate_n_windows, validate_nan_policy

from .data import DiscreteRegionData
from .metagene import MetageneProfileSegment, segments_total_bins
from .windowing import _bin_points_windows_fast


@dataclass(frozen=True)
class DistributionPlotData:
    rows: list[tuple[str, float]]
    x_label: str
    y_label: str


def _region_label(
    label: str | None,
    *,
    index: int,
) -> str:
    return label if label is not None else f"region_{index}"


def _prepare_region_values(
    values: np.ndarray,
    *,
    as_percent: bool,
    nan_fill: Optional[float],
    nan_policy: NanPolicy,
) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if as_percent:
        out = out * 100.0
    if nan_fill is not None:
        out = np.where(np.isnan(out), nan_fill, out)
    if nan_policy is NanPolicy.ZERO:
        out = np.where(np.isnan(out), 0.0, out)
    if nan_policy is NanPolicy.DROP:
        out = out[np.isfinite(out)]
    return out


def _distribution_rows(
    drd: DiscreteRegionData,
    *,
    as_percent: bool,
    nan_fill: Optional[float],
    segments: list[MetageneProfileSegment] | None,
    n_windows: Optional[int],
    nan_policy: NanPolicy,
) -> list[tuple[str, float]]:
    if n_windows is None:
        n_windows = segments_total_bins(segments) if segments else 40
    else:
        n_windows = validate_n_windows(n_windows)

    rows: list[tuple[str, float]] = []
    for pos, dens in zip(drd.positions, drd.densities):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        if as_percent:
            y = y * 100.0
        if nan_fill is not None:
            y = np.where(np.isnan(y), nan_fill, y)
        if nan_policy is NanPolicy.ZERO:
            y = np.where(np.isnan(y), 0.0, y)
        elif nan_policy is NanPolicy.DROP:
            finite = np.isfinite(y)
            x = x[finite]
            y = y[finite]
        if x.size == 0:
            continue
        binned = _bin_points_windows_fast(
            x,
            y,
            n_windows=n_windows,
            agg=AggMethod.Mean,
            nan_policy=nan_policy,
        )
        for bin_index, value in enumerate(binned):
            if np.isfinite(value):
                rows.append((str(bin_index), float(value)))
    return rows


def _position_axis_label() -> str:
    return "Relative feature position"


def build_box_distribution_data(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    nan_policy: NanPolicy = NanPolicy.DROP,
    per_region: bool = False,
) -> DistributionPlotData:
    nan_policy = validate_nan_policy(nan_policy)

    if per_region:
        rows: list[tuple[str, float]] = []
        for index, (dens, label) in enumerate(zip(drd.densities, drd.labels), start=1):
            values = _prepare_region_values(
                dens,
                as_percent=as_percent,
                nan_fill=nan_fill,
                nan_policy=nan_policy,
            )
            if values.size == 0:
                continue
            rows.append((_region_label(label, index=index), float(np.nanmean(values))))
        x_label = "Region"
    else:
        rows = _distribution_rows(
            drd,
            as_percent=as_percent,
            nan_fill=nan_fill,
            segments=segments,
            n_windows=n_windows,
            nan_policy=nan_policy,
        )
        x_label = _position_axis_label()

    y_label = (
        "Mean methylation per feature (%)"
        if as_percent
        else "Mean methylation per feature"
    )
    return DistributionPlotData(rows=rows, x_label=x_label, y_label=y_label)


def build_violin_distribution_data(
    drd: DiscreteRegionData,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    n_windows: Optional[int] = None,
    as_percent: bool = True,
    nan_fill: Optional[float] = None,
    nan_policy: NanPolicy = NanPolicy.DROP,
    per_region: bool = False,
) -> DistributionPlotData:
    nan_policy = validate_nan_policy(nan_policy)

    if per_region:
        rows: list[tuple[str, float]] = []
        for index, (dens, label) in enumerate(zip(drd.densities, drd.labels), start=1):
            values = _prepare_region_values(
                dens,
                as_percent=as_percent,
                nan_fill=nan_fill,
                nan_policy=nan_policy,
            )
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue
            region_label = _region_label(label, index=index)
            rows.extend((region_label, float(value)) for value in values)
        x_label = "Region"
    else:
        rows = _distribution_rows(
            drd,
            as_percent=as_percent,
            nan_fill=nan_fill,
            segments=segments,
            n_windows=n_windows,
            nan_policy=nan_policy,
        )
        x_label = _position_axis_label()

    y_label = (
        "Mean methylation per feature (%)"
        if as_percent
        else "Mean methylation per feature"
    )
    return DistributionPlotData(rows=rows, x_label=x_label, y_label=y_label)
