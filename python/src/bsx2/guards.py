from __future__ import annotations

import numpy as np
from beartype.typing import Callable, Sequence


def require_equal_length(
    left: Sequence[object],
    right: Sequence[object],
    *,
    left_name: str,
    right_name: str,
    message: str | None = None,
) -> None:
    if len(left) != len(right):
        if message is None:
            message = f"{left_name} and {right_name} must have the same length"
        raise ValueError(message)


def require_nan_policy_compatible(y: np.ndarray, *, nan_policy: str) -> None:
    if np is None:
        raise ModuleNotFoundError("numpy is required for NaN-policy checks")
    if nan_policy == "raise" and np.isnan(y).any():
        raise ValueError("NaN values present; use nan_policy='interp' or 'mask'")


def require_per_segment_profile(
    segments: Sequence[object] | None,
    *,
    profile_len: int,
) -> list[int]:
    if not segments:
        raise ValueError("per_segment=True requires segments")
    seg_nbins: list[int] = []
    for seg in segments:
        n_bins = getattr(seg, "n_bins", None)
        if isinstance(n_bins, bool) or not isinstance(n_bins, int):
            raise ValueError("per_segment=True requires integer segment n_bins")
        seg_nbins.append(n_bins)
    if sum(seg_nbins) != profile_len:
        raise ValueError("per_segment=True requires profile length equal to total_bins")
    return seg_nbins


def resolve_reader_accessors(
    reader: object,
) -> tuple[
    Callable[[object], object] | None,
    Callable[[list[object]], object] | None,
]:
    query_fn = getattr(reader, "query", None)
    iter_fn = getattr(reader, "iter_contigs", None)
    query = query_fn if callable(query_fn) else None
    iter_contigs = iter_fn if callable(iter_fn) else None
    if query is None and iter_contigs is None:
        raise AttributeError("reader must implement query(contig) or iter_contigs(contigs)")
    return query, iter_contigs


def require_single_choice(
    value: str | Sequence[str],
    *,
    name: str,
    allowed_hint: str,
) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, list | tuple | set):
        if len(value) != 1:
            raise ValueError(f"{name} must be a single value: {allowed_hint}")
        if isinstance(value, set):
            return next(iter(value))
        return value[0]
    return str(value)


def require_savgol_filter():
    try:
        from scipy.signal import savgol_filter  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("scipy is required for Savitzky-Golay smoothing") from e
    return savgol_filter
