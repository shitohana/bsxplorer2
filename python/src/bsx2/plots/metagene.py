from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict, deque
from typing import TYPE_CHECKING, Callable, List, Optional, Sequence, Tuple, cast

import numpy as np
from bsx2 import RegionReader
from bsx2.guards import (
    require_nan_policy_compatible,
    require_per_segment_profile,
    require_savgol_filter,
    resolve_reader_accessors,
)
from bsx2.plots.data import DiscreteRegionData
from bsx2.validation import validate_segments, validate_smoothing

if TYPE_CHECKING:
    from bsx2.types import Contig


# ============================================================
# Segment definition
# ============================================================

@dataclass(frozen=True)
class MetageneProfileSegment:
    name: str
    n_bins: int


_DEFAULT_SEGMENTS: tuple[MetageneProfileSegment, ...] = (
    MetageneProfileSegment("region", 100),
)


def segments_total_bins(segments: Sequence[MetageneProfileSegment]) -> int:
    return validate_segments(segments)


def segment_boundaries(segments: Sequence[MetageneProfileSegment]) -> List[float]:
    total = segments_total_bins(segments)
    cum = 0
    bounds: List[float] = []
    for s in segments:
        cum += s.n_bins
        bounds.append(cum / total)
    return bounds


# ============================================================
# Smoothing helpers
# ============================================================

def _coerce_smooth_config(
    smooth: dict | int | None,
    *,
    total_bins: int,
) -> dict | None:
    return validate_smoothing(smooth, total_bins=total_bins)


def _validate_savgol_config(
    cfg: dict,
    *,
    series_len: int,
    label: str = "series",
) -> None:
    validate_smoothing(cfg, series_len=series_len, label=label)


def _savgol_filter_1d(y: np.ndarray, **kwargs) -> np.ndarray:
    savgol_filter = require_savgol_filter()

    if y.size == 0:
        return y

    window_length = int(kwargs["window_length"])
    polyorder = int(kwargs["polyorder"])
    mode = kwargs.get("mode", "interp")
    cval = kwargs.get("cval", 0.0)

    nan_policy = kwargs.get("nan_policy", "interp")
    _validate_savgol_config(kwargs, series_len=len(y))
    require_nan_policy_compatible(y, nan_policy=nan_policy)

    if nan_policy == "interp":
        if np.all(np.isnan(y)):
            return y
        good = np.isfinite(y)
        if int(good.sum()) < window_length:
            return y
        idx = np.arange(y.size)
        y_filled = y.copy()
        y_filled[~good] = np.interp(idx[~good], idx[good], y[good])
        return savgol_filter(y_filled, window_length, polyorder, mode=mode, cval=cval)

    if nan_policy == "mask":
        if not np.isnan(y).any():
            if y.size < window_length:
                return y
            return savgol_filter(y, window_length, polyorder, mode=mode, cval=cval)

        out = y.copy()
        good = np.isfinite(y)
        i = 0
        while i < y.size:
            if not good[i]:
                i += 1
                continue
            j = i
            while j < y.size and good[j]:
                j += 1
            if j - i >= window_length:
                out[i:j] = savgol_filter(y[i:j], window_length, polyorder, mode=mode, cval=cval)
            i = j
        return out

    if y.size < window_length:
        return y

    return savgol_filter(y, window_length, polyorder, mode=mode, cval=cval)


def _apply_savgol_smoothing(
    y: np.ndarray,
    cfg: dict,
    *,
    segments: Sequence[MetageneProfileSegment] | None = None,
) -> np.ndarray:
    if not isinstance(y, np.ndarray):
        y = np.asarray(y, dtype=float)

    if y.size == 0:
        return y

    if cfg.get("per_segment"):
        seg_nbins = require_per_segment_profile(segments, profile_len=y.size)
        out = []
        offset = 0
        for n in seg_nbins:
            seg = y[offset: offset + n]
            out.append(_savgol_filter_1d(seg, **cfg))
            offset += n
        return np.concatenate(out) if out else y

    return _savgol_filter_1d(y, **cfg)


def _clip_profile(y: np.ndarray) -> np.ndarray:
    if y.size == 0:
        return y
    y_out = y.astype(float, copy=True)
    mask = np.isfinite(y_out)
    if np.any(mask):
        np.clip(y_out, 0.0, 1.0, out=y_out, where=mask)
    return y_out


# ============================================================
# Strand helper
# ============================================================

def _is_negative_strand(contig: Contig) -> bool:
    val = None
    if hasattr(contig, "strand_str"):
        try:
            v = contig.strand_str
            val = v() if callable(v) else v
        except Exception:
            val = None
    if val is None and hasattr(contig, "strand"):
        try:
            val = str(contig.strand)
        except Exception:
            pass
    return str(val).strip() == "-"


# ============================================================
# Fast numeric conversion
# ============================================================

def _to_np_float64(col) -> np.ndarray:
    if hasattr(col, "to_numpy"):
        return np.asarray(col.to_numpy(), dtype=np.float64)
    if hasattr(col, "to_numpy_array"):
        return np.asarray(col.to_numpy_array(), dtype=np.float64)
    if hasattr(col, "to_array"):
        return np.asarray(col.to_array(), dtype=np.float64)
    return np.asarray(col.to_list(), dtype=np.float64)


# ============================================================
# Main compute function (optimized)
# ============================================================

def compute_discrete_regions(
    reader: RegionReader,
    contigs: Sequence[Contig],
    *,
    segments: Sequence[MetageneProfileSegment] = _DEFAULT_SEGMENTS,
    reverse_negative: bool = True,
    labels: Optional[Sequence[str]] = None,
    progress: bool = False,
    progress_every: int = 1,
) -> DiscreteRegionData:

    try:
        sorted_contigs = reader.index().sort(list(contigs))
    except Exception:
        sorted_contigs = list(contigs)

    contigs_list = sorted_contigs
    data = DiscreteRegionData()

    query_fn, iter_contigs_fn = resolve_reader_accessors(reader)

    total = len(contigs_list)
    step = progress_every if progress_every > 0 else 1

    np_isfinite = np.isfinite
    np_argsort = np.argsort
    np_clip = np.clip

    if query_fn is not None:
        last_seqname = None
        for idx, contig in enumerate(contigs_list):

            seqname = getattr(contig, "seqname", None)
            try:
                seqname = seqname() if callable(seqname) else seqname
            except Exception:
                seqname = None

            if seqname is not None and seqname != last_seqname:
                reset_fn = getattr(reader, "reset", None)
                if callable(reset_fn):
                    try:
                        reset_fn()
                    except Exception:
                        pass
                last_seqname = seqname

            try:
                batch = query_fn(contig)
            except Exception:
                continue
            if batch is None:
                continue

            start = getattr(contig, "start", None)
            end = getattr(contig, "end", None)
            if start is None or end is None or end <= start:
                continue

            try:
                pos = _to_np_float64(batch.position())
                dens = _to_np_float64(batch.density())
            except Exception:
                continue

            if pos.size == 0 or dens.size == 0:
                continue

            x = (pos - float(start)) / float(end - start)
            y = dens

            if reverse_negative and _is_negative_strand(contig):
                x = 1.0 - x

            mask = np_isfinite(x) & (x >= 0.0) & (x <= 1.0)
            x = x[mask]
            y = y[mask]
            if x.size == 0:
                continue

            if x.size >= 2 and not np.all(x[:-1] <= x[1:]):
                order = np_argsort(x, kind="mergesort")
                x = x[order]
                y = y[order]

            bad = ~np_isfinite(y)
            if bad.any():
                y = y.copy()
                y[bad] = np.nan

            good = np_isfinite(y)
            if good.any():
                np_clip(y, 0.0, 1.0, out=y, where=good)

            label = labels[idx] if labels and idx < len(labels) else None
            data.insert(x, y, label)

    else:
        batch_iter = cast(Callable[[list[object]], object], iter_contigs_fn)(contigs_list)

        for idx, contig in enumerate(contigs_list):
            try:
                batch = next(batch_iter)
            except StopIteration:
                break
            except Exception:
                continue

            if batch is None:
                continue

            start = getattr(contig, "start", None)
            end = getattr(contig, "end", None)
            if start is None or end is None or end <= start:
                continue

            try:
                pos = _to_np_float64(batch.position())
                dens = _to_np_float64(batch.density())
            except Exception:
                continue

            if pos.size == 0 or dens.size == 0:
                continue

            x = (pos - float(start)) / float(end - start)
            y = dens

            if reverse_negative and _is_negative_strand(contig):
                x = 1.0 - x

            mask = np_isfinite(x) & (x >= 0.0) & (x <= 1.0)
            x = x[mask]
            y = y[mask]
            if x.size == 0:
                continue

            if x.size >= 2 and not np.all(x[:-1] <= x[1:]):
                order = np_argsort(x, kind="mergesort")
                x = x[order]
                y = y[order]

            bad = ~np_isfinite(y)
            if bad.any():
                y = y.copy()
                y[bad] = np.nan

            good = np_isfinite(y)
            if good.any():
                np_clip(y, 0.0, 1.0, out=y, where=good)

            label = labels[idx] if labels and idx < len(labels) else None
            data.insert(x, y, label)

    return data
