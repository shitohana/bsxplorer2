from __future__ import annotations
from enum import StrEnum
import math
from typing import Sequence
from bsx2 import Strand, Context, AggMethod

import numpy as np

class NanPolicy(StrEnum):
    KEEP = "keep"
    ZERO = "zero"
    DROP = "drop"


class ChrEmptyPolicy(StrEnum):
    NAN = "nan"
    ZERO = "zero"
    DROP = "drop"


class ChrLineStat(StrEnum):
    WEIGHTED_MEAN = "weighted_mean"
    MEAN = "mean"


_SMOOTH_MODES = {"interp", "nearest", "mirror", "constant", "wrap"}
_SMOOTH_NAN_POLICIES = {"interp", "mask", "raise"}
_RANK_SCORE_VALUES = {"mean", "body_mean"}
_SORT_ORDER_VALUES = {"asc", "desc"}
_CHR_EMPTY_POLICY_VALUES = {p.value for p in ChrEmptyPolicy}
_CHR_LINE_STAT_VALUES = {p.value for p in ChrLineStat}


def _make_odd(n: int) -> int:
    return n + 1 if n % 2 == 0 else n


def validate_nan_policy(nan_policy: object) -> NanPolicy:
    if not isinstance(nan_policy, NanPolicy):
        raise ValueError("nan_policy must be a NanPolicy value")
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


def validate_bin_size_bp(bin_size_bp: object) -> int:
    return validate_positive_int(bin_size_bp, name="bin_size_bp")


def validate_window_agg(agg: object) -> AggMethod:
    if not isinstance(agg, AggMethod) or agg is AggMethod.GeometricMean:
        raise ValueError(
            "agg must be one of: AggMethod.Mean, AggMethod.Median, "
            "AggMethod.Max, AggMethod.Min"
        )
    return agg


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


def _validate_context_value(value: object) -> Context:
    if not isinstance(value, Context):
        raise ValueError(
            "context must be a Context value or a list of Context values"
        )
    return value


def validate_context(context: object) -> Context | list[Context]:
    if isinstance(context, (list, tuple, set)):
        values = [_validate_context_value(v) for v in context]
        if not values:
            raise ValueError("context list must not be empty")
        return values
    return _validate_context_value(context)


def validate_single_context(context: object) -> Context:
    value = validate_context(context)
    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError("context must contain exactly one value")
        return value[0]
    return value


def validate_chr_empty_policy(policy: object) -> ChrEmptyPolicy:
    if isinstance(policy, ChrEmptyPolicy):
        return policy
    text = str(policy).strip().lower()
    if text not in _CHR_EMPTY_POLICY_VALUES:
        raise ValueError("empty_policy must be one of: nan, zero, drop")
    return ChrEmptyPolicy(text)


def validate_chr_line_stat(stat: object) -> ChrLineStat:
    if isinstance(stat, ChrLineStat):
        return stat
    text = str(stat).strip().lower()
    if text in {"weighted", "wmean"}:
        text = ChrLineStat.WEIGHTED_MEAN.value
    elif text in {"unweighted", "unweighted_mean", "site_mean"}:
        text = ChrLineStat.MEAN.value
    if text not in _CHR_LINE_STAT_VALUES:
        raise ValueError(
            "stat must be one of: weighted_mean, mean "
            "(aliases: weighted, wmean, unweighted, unweighted_mean, site_mean)"
        )
    return ChrLineStat(text)


def validate_chr_lengths(
    chr_lengths: dict[str, int] | None,
) -> dict[str, int] | None:
    if chr_lengths is None:
        return None
    if not isinstance(chr_lengths, dict):
        raise ValueError("chr_lengths must be a mapping {chr_name: length_bp}")

    validated: dict[str, int] = {}
    for seqname, length in chr_lengths.items():
        if not isinstance(seqname, str) or not seqname.strip():
            raise ValueError("chr_lengths keys must be non-empty strings")
        validated[seqname] = validate_positive_int(
            length,
            name=f"chr_lengths[{seqname}]",
        )
    return validated


def validate_strand(strand: object, *, allow_both: bool = True) -> Strand:
    if not isinstance(strand, Strand):
        raise ValueError(
            "strand must be a Strand value"
            if allow_both
            else "strand must be Strand.Forward or Strand.Reverse"
        )

    if not allow_both and strand is Strand.Null:
        raise ValueError("strand must be Strand.Forward or Strand.Reverse")
    return strand


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
