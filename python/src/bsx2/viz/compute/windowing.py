from __future__ import annotations

import numpy as np

from bsx2 import AggMethod
from bsx2.guards import require_equal_length
from bsx2.validation import (
    NanPolicy,
    validate_matrix_shape,
    validate_n_windows,
    validate_nan_policy,
)


def _bin_points_windows_fast(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    *,
    n_windows: int,
    agg: AggMethod,
    nan_policy: NanPolicy,
) -> np.ndarray:
    x_vals = validate_matrix_shape(x_vals, 1, name="x_vals")
    y_vals = validate_matrix_shape(y_vals, 1, name="y_vals")
    require_equal_length(x_vals, y_vals, left_name="x_vals", right_name="y_vals")
    n_windows = validate_n_windows(n_windows)
    nan_policy = validate_nan_policy(nan_policy)

    x = np.asarray(x_vals, dtype=np.float64)
    y = np.asarray(y_vals, dtype=np.float64)

    if nan_policy is NanPolicy.ZERO:
        y = np.where(np.isfinite(y), y, 0.0)
    else:
        finite = np.isfinite(y)
        x = x[finite]
        y = y[finite]

    out = np.full(n_windows, np.nan, dtype=np.float64)
    if y.size == 0:
        return out

    idx = (x * n_windows).astype(np.int64)
    idx[idx == n_windows] = n_windows - 1

    if agg is AggMethod.Mean:
        counts = np.bincount(idx, minlength=n_windows)
        sums = np.bincount(idx, weights=y, minlength=n_windows)
        nonempty = counts > 0
        out[nonempty] = sums[nonempty] / counts[nonempty]
        return out

    if agg is AggMethod.Min:
        tmp = np.full(n_windows, np.inf, dtype=np.float64)
        np.minimum.at(tmp, idx, y)
        tmp[tmp == np.inf] = np.nan
        return tmp

    if agg is AggMethod.Max:
        tmp = np.full(n_windows, -np.inf, dtype=np.float64)
        np.maximum.at(tmp, idx, y)
        tmp[tmp == -np.inf] = np.nan
        return tmp

    if agg is AggMethod.Median:
        cuts = np.flatnonzero(np.diff(idx)) + 1
        y_groups = np.split(y, cuts)
        bin_ids = idx[np.r_[0, cuts]]
        for bin_id, group in zip(bin_ids, y_groups, strict=True):
            out[int(bin_id)] = float(np.median(group))
        return out

    raise ValueError(f"unsupported agg: {agg}")


def _rank_compress(
    z_sorted: np.ndarray,
    rank_rows: int,
    *,
    fill: float | None = 0.0,
) -> np.ndarray:
    if z_sorted.ndim != 2:
        return z_sorted

    n_rows, n_bins = z_sorted.shape
    if n_rows == 0:
        return np.empty((0, n_bins), dtype=float)

    rows = max(int(rank_rows), 1)
    sums = np.zeros((rows, n_bins), dtype=float)
    counts = np.zeros((rows, n_bins), dtype=np.int32)

    k = np.arange(rows + 1, dtype=np.int64)
    bounds = (k * n_rows + rows - 1) // rows

    for bucket in range(rows):
        start = int(bounds[bucket])
        end = int(bounds[bucket + 1])
        if start >= end:
            continue
        block = z_sorted[start:end]
        finite = np.isfinite(block)
        if not finite.any():
            continue
        sums[bucket] = np.where(finite, block, 0.0).sum(axis=0)
        counts[bucket] = finite.sum(axis=0, dtype=np.int32)

    out = (
        np.full((rows, n_bins), np.nan, dtype=float)
        if fill is None
        else np.full((rows, n_bins), float(fill), dtype=float)
    )
    np.divide(sums, counts, out=out, where=counts > 0)
    return out
