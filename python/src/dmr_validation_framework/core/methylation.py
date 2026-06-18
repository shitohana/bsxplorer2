"""Methylation effect-size helpers."""

from __future__ import annotations

import numpy as np


def methylation_ratio(methylated: float, total: float) -> float:
    if total <= 0 or np.isnan(total):
        return float("nan")
    return float(methylated) / float(total)


def delta_from_group_means(group_a: float, group_b: float) -> float:
    return float(group_b) - float(group_a)

