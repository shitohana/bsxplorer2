"""Small statistical helpers shared by validation checks."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd


def safe_ratio(numerator: float, denominator: float, default: float = float("nan")) -> float:
    if denominator == 0 or np.isnan(denominator):
        return default
    return float(numerator) / float(denominator)


def sign(value: float, *, zero: int = 0) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return zero


def direction_from_delta(delta: float, threshold: float = 0.0) -> str:
    if delta > threshold:
        return "hyper"
    if delta < -threshold:
        return "hypo"
    return "near_zero"


def bootstrap_ci(values: np.ndarray, alpha: float = 0.05) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan"), float("nan")
    return (
        float(np.quantile(values, alpha / 2.0)),
        float(np.quantile(values, 1.0 - alpha / 2.0)),
    )

def bh_qvalues(p_values: Iterable[float]) -> list[float]:
    values = np.asarray([np.nan if pd.isna(p) else float(p) for p in p_values], dtype=float)
    q = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return q.tolist()
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    valid_indices = np.where(valid)[0]
    q[valid_indices[order]] = adjusted
    return q.tolist()
