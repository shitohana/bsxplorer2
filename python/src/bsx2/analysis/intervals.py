"""Canonical interval arithmetic for DMR coordinates.

Single source of truth for interval length and reciprocal overlap. Every DMR
coordinate in the BSX2 / validation stack is normalized to BED-style
**0-based half-open** intervals at import time (see
:mod:`bsx2.analysis.external_dmr_callers`), so all helpers here assume that
convention. Do not reimplement these with ``end - start + 1`` style math
elsewhere: mixing closed/1-based and half-open/0-based definitions silently
shifts reciprocal-overlap fractions for short regions.
"""

from __future__ import annotations

from typing import Any


def interval_length(start: Any, end: Any) -> float:
    """Length of a 0-based half-open interval ``[start, end)``."""
    return max(float(end) - float(start), 0.0)


def overlap_length(a_start: Any, a_end: Any, b_start: Any, b_end: Any) -> float:
    """Overlap length (bp) of two 0-based half-open intervals."""
    return max(0.0, min(float(a_end), float(b_end)) - max(float(a_start), float(b_start)))


def reciprocal_overlap(a_start: Any, a_end: Any, b_start: Any, b_end: Any) -> float:
    """Strict reciprocal overlap: ``min(overlap / len_a, overlap / len_b)``.

    Returns ``0.0`` for empty or zero-length intervals.
    """
    overlap = overlap_length(a_start, a_end, b_start, b_end)
    len_a = interval_length(a_start, a_end)
    len_b = interval_length(b_start, b_end)
    if overlap <= 0.0 or len_a <= 0.0 or len_b <= 0.0:
        return 0.0
    return min(overlap / len_a, overlap / len_b)


def min_length_overlap(a_start: Any, a_end: Any, b_start: Any, b_end: Any) -> float:
    """Diagnostic overlap fraction relative to the shorter interval."""
    overlap = overlap_length(a_start, a_end, b_start, b_end)
    if overlap <= 0.0:
        return 0.0
    len_a = interval_length(a_start, a_end)
    len_b = interval_length(b_start, b_end)
    return overlap / max(1.0, min(len_a, len_b))
