from __future__ import annotations

import math
from typing import Sequence

try:
    import numpy as np  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - optional for non-plot usage
    np = None  # type: ignore

_NAN_POLICIES = {"drop", "zero", "keep"}
_SMOOTH_MODES = {"interp", "nearest", "mirror", "constant", "wrap"}
_SMOOTH_NAN_POLICIES = {"interp", "mask", "raise"}
_CONTEXT_VALUES = {"CG", "CHG", "CHH"}
_WINDOW_AGG_VALUES = {"mean", "median", "max", "min"}
_RANK_SCORE_VALUES = {"mean", "body_mean"}
_SORT_ORDER_VALUES = {"asc", "desc"}
_STRAND_ALIASES = {
    "+": "+",
    "-": "-",
    "both": "both",
    "all": "both",
    "forward": "+",
    "reverse": "-",
    "f": "+",
    "r": "-",
}


def _make_odd(n: int) -> int:
    return n + 1 if n % 2 == 0 else n


def validate_nan_policy(nan_policy: str) -> str:
    if nan_policy not in _NAN_POLICIES:
        raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")
    return nan_policy


def validate_positive_int(value: object, *, name: str, allow_zero: bool = False) -> int:
    if np is not None and isinstance(value, np.integer):
        value = int(value)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be >= 0" if allow_zero else f"{name} must be > 0")
    if allow_zero:
        if value < 0:
            raise ValueError(f"{name} must be >= 0")
    else:
        if value <= 0:
            raise ValueError(f"{name} must be > 0")
    return value


def validate_n_windows(n_windows: object) -> int:
    return validate_positive_int(n_windows, name="n_windows")


def validate_window_agg(agg: object) -> str:
    text = str(agg).strip().lower()
    if text not in _WINDOW_AGG_VALUES:
        raise ValueError("agg must be one of: mean, median, max, min")
    return text


def validate_rank_score(rank_score: object) -> str:
    text = str(rank_score).strip().lower()
    if text not in _RANK_SCORE_VALUES:
        raise ValueError("rank_score must be 'mean' or 'body_mean'")
    return text


def validate_sort_order(sort_order: object) -> str:
    text = str(sort_order).strip().lower()
    if text not in _SORT_ORDER_VALUES:
        raise ValueError("sort_order must be 'asc' or 'desc'")
    return text


def validate_smoothing(
    smooth: dict | int | None,
    *,
    total_bins: int | None = None,
    series_len: int | None = None,
    label: str = "series",
) -> dict | None:
    if smooth is None:
        return None
    if isinstance(smooth, bool):
        raise ValueError("smooth must be a dict, int, or None")

    if isinstance(smooth, int):
        if smooth == 0:
            return None
        if smooth < 0:
            raise ValueError("smooth int must be >= 0")
        if total_bins is None or total_bins <= 0:
            raise ValueError("total_bins must be > 0 when smooth is int")
        window_length = _make_odd(max(3, total_bins // smooth))
        cfg = {
            "method": "savgol",
            "window_length": window_length,
            "polyorder": 2,
            "apply": "post",
            "mode": "interp",
            "nan_policy": "interp",
            "per_segment": False,
            "cval": 0.0,
        }
    elif isinstance(smooth, dict):
        method = smooth.get("method")
        if method != "savgol":
            raise ValueError("smooth.method must be 'savgol'")
        cfg = {
            "method": method,
            "window_length": smooth.get("window_length"),
            "polyorder": smooth.get("polyorder"),
            "apply": smooth.get("apply", "post"),
            "mode": smooth.get("mode", "interp"),
            "nan_policy": smooth.get("nan_policy", "interp"),
            "per_segment": smooth.get("per_segment", False),
            "cval": smooth.get("cval", 0.0),
        }
    else:
        raise ValueError("smooth must be a dict, int, or None")

    window_length = cfg.get("window_length")
    polyorder = cfg.get("polyorder")
    mode = cfg.get("mode", "interp")
    nan_policy = cfg.get("nan_policy", "interp")

    if not isinstance(window_length, int) or window_length < 3:
        raise ValueError("smooth.window_length must be an odd integer >= 3")
    if window_length % 2 == 0:
        raise ValueError("smooth.window_length must be an odd integer")
    if not isinstance(polyorder, int) or polyorder < 0:
        raise ValueError("smooth.polyorder must be an integer >= 0")
    if polyorder >= window_length:
        raise ValueError("smooth.polyorder must be < smooth.window_length")

    if mode not in _SMOOTH_MODES:
        raise ValueError(
            "smooth.mode must be one of: interp, nearest, mirror, constant, wrap"
        )
    if nan_policy not in _SMOOTH_NAN_POLICIES:
        raise ValueError("smooth.nan_policy must be one of: interp, mask, raise")

    if series_len is not None and mode == "interp" and window_length > series_len:
        raise ValueError(
            f"smooth.window_length ({window_length}) must be <= {label} "
            f"length ({series_len}) when mode='interp'"
        )

    return cfg


def validate_segments(
    segments: Sequence[object],
    *,
    expected_len: int | None = None,
    flank_bp: int | float | None = None,
) -> int:
    if not segments:
        raise ValueError("segments must not be empty")

    if expected_len is not None and len(segments) != expected_len:
        if expected_len == 3:
            raise ValueError(
                "combined metagene requires exactly 3 segments (up/body/down)"
            )
        raise ValueError(f"segments must contain exactly {expected_len} items")

    total_bins = 0
    for seg in segments:
        name = getattr(seg, "name", None)
        n_bins = getattr(seg, "n_bins", None)
        if not isinstance(name, str) or not name.strip():
            raise ValueError("each segment name must be a non-empty string")
        if isinstance(n_bins, bool) or not isinstance(n_bins, int) or n_bins <= 0:
            raise ValueError("each segment n_bins must be > 0")
        total_bins += n_bins

    if flank_bp is not None:
        if isinstance(flank_bp, bool) or not isinstance(flank_bp, (int, float)):
            raise ValueError("flank_bp must be >= 0")
        if flank_bp < 0:
            raise ValueError("flank_bp must be >= 0")

    return total_bins


def validate_matrix_shape(
    arr: np.ndarray,
    expected_dims: int | Sequence[int | None],
    *,
    name: str = "array",
) -> np.ndarray:
    if np is None:
        raise ModuleNotFoundError("numpy is required for matrix validation")
    arr = np.asarray(arr)

    if isinstance(expected_dims, int):
        if arr.ndim != expected_dims:
            raise ValueError(f"{name} must be {expected_dims}D")
        return arr

    expected_dims = tuple(expected_dims)
    if arr.ndim != len(expected_dims):
        raise ValueError(f"{name} must be {len(expected_dims)}D")

    for axis, expected in enumerate(expected_dims):
        if expected is not None and arr.shape[axis] != expected:
            raise ValueError(
                f"{name} has invalid shape: expected axis {axis} size "
                f"{expected}, got {arr.shape[axis]}"
            )
    return arr


def _normalize_context_value(value: object) -> str:
    if hasattr(value, "name"):
        value = getattr(value, "name")
    text = str(value).strip().upper()
    if text not in _CONTEXT_VALUES:
        raise ValueError("context must be CG, CHG, CHH or a list of these values")
    return text


def validate_context(context: object) -> str | list[str]:
    if isinstance(context, (list, tuple, set)):
        values = [_normalize_context_value(v) for v in context]
        if not values:
            raise ValueError("context list must not be empty")
        return values
    return _normalize_context_value(context)


def validate_strand(strand: object, *, allow_both: bool = True) -> str:
    if hasattr(strand, "name"):
        strand = getattr(strand, "name")
    text = str(strand).strip().lower()
    normalized = _STRAND_ALIASES.get(text)
    if normalized is None:
        raise ValueError("strand must be '+', '-', or 'both'")
    if not allow_both and normalized == "both":
        raise ValueError("strand must be '+' or '-'")
    return normalized


def validate_min_coverage(
    min_cov: int | float,
    *,
    integer: bool = False,
) -> int | float:
    if isinstance(min_cov, bool) or not isinstance(min_cov, (int, float)):
        raise ValueError("min_coverage must be >= 0")
    value = float(min_cov)
    if not math.isfinite(value) or value < 0:
        raise ValueError("min_coverage must be >= 0")
    if integer:
        if not value.is_integer():
            raise ValueError("min_coverage must be an integer >= 0")
        return int(value)
    return min_cov
