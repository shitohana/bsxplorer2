"""Interval matching helpers for the validation framework.

Scalar interval arithmetic lives in :mod:`bsx2.analysis.intervals` and is the
single source of truth (BED-style 0-based half-open). This module only adds the
``matched_pairs`` sweep on top of it, plus thin pandas-row wrappers kept for
backward compatibility.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from bsx2.analysis.intervals import (
    interval_length,
    overlap_length,
    reciprocal_overlap as reciprocal_overlap_scalar,
)


def interval_len(start: int | float, end: int | float) -> float:
    return interval_length(start, end)


def overlap_len(start_a: int | float, end_a: int | float, start_b: int | float, end_b: int | float) -> float:
    return overlap_length(start_a, end_a, start_b, end_b)


def interval_overlap(start_a: int | float, end_a: int | float, start_b: int | float, end_b: int | float) -> float:
    return overlap_length(start_a, end_a, start_b, end_b)


def min_length_overlap(row_a: pd.Series, row_b: pd.Series) -> float:
    """Legacy diagnostic overlap fraction relative to the shorter interval."""
    overlap = overlap_length(row_a.start, row_a.end, row_b.start, row_b.end)
    if overlap <= 0:
        return 0.0
    len_a = interval_length(row_a.start, row_a.end)
    len_b = interval_length(row_b.start, row_b.end)
    return overlap / max(1.0, min(len_a, len_b))


def reciprocal_overlap(row_a: pd.Series, row_b: pd.Series) -> float:
    """Strict reciprocal overlap for two pandas rows with start/end attributes."""
    return reciprocal_overlap_scalar(row_a.start, row_a.end, row_b.start, row_b.end)


def matched_pairs(a: pd.DataFrame, b: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows: list[dict] = []
    if a.empty or b.empty:
        return pd.DataFrame(rows)
    required = {"chrom", "context", "start", "end"}
    if not required.issubset(a.columns) or not required.issubset(b.columns):
        return pd.DataFrame(rows)

    b_groups = {
        key: sub.sort_values("start").reset_index(drop=True)
        for key, sub in b.groupby(["chrom", "context"], dropna=False, sort=False)
    }
    for key, sub_a in a.groupby(["chrom", "context"], dropna=False, sort=False):
        sub_b = b_groups.get(key)
        if sub_b is None or sub_b.empty:
            continue
        starts = sub_b["start"].to_numpy(dtype=np.int64, copy=False)
        ends = sub_b["end"].to_numpy(dtype=np.int64, copy=False)
        max_len_b = int((ends - starts).max()) if len(starts) else 0
        records_b = list(sub_b.itertuples(index=False))
        for idx_a, row_a in enumerate(sub_a.itertuples(index=False)):
            a_map = row_a._asdict()
            a_start = int(a_map["start"])
            a_end = int(a_map["end"])
            candidate_start = int(np.searchsorted(starts, a_start - max_len_b - 1, side="left"))
            candidate_end = int(np.searchsorted(starts, a_end, side="right"))
            if candidate_start >= candidate_end:
                continue
            for rel_idx_b in range(candidate_start, candidate_end):
                if int(ends[rel_idx_b]) <= a_start:
                    continue
                row_b = records_b[rel_idx_b]
                b_map = row_b._asdict()
                overlap = overlap_length(a_start, a_end, b_map["start"], b_map["end"])
                if overlap <= 0:
                    continue
                len_a = interval_length(a_start, a_end)
                len_b = interval_length(b_map["start"], b_map["end"])
                if len_a <= 0 or len_b <= 0:
                    continue
                ro = min(overlap / len_a, overlap / len_b)
                if ro >= threshold:
                    rows.append(
                        {
                            "region_id_a": a_map.get("region_id", idx_a),
                            "region_id_b": b_map.get("region_id", rel_idx_b),
                            "chrom": key[0],
                            "context": key[1],
                            "start_a": a_start,
                            "end_a": a_end,
                            "start_b": b_map["start"],
                            "end_b": b_map["end"],
                            "overlap_bp": overlap,
                            "length_a": len_a,
                            "length_b": len_b,
                            "reciprocal_overlap": ro,
                            "min_length_overlap": overlap / max(1.0, min(len_a, len_b)),
                            "reciprocal_overlap_definition": "strict_min_overlap_over_each_interval",
                            "delta_a": a_map.get("delta", np.nan),
                            "delta_b": b_map.get("delta", np.nan),
                        }
                    )
    return pd.DataFrame(rows)
