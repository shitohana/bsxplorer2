from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple
import heapq

import numpy as np
import pandas as pd

from bsx2 import io as _io
from bsx2.plots.data import DiscreteRegionData, LinePlotData


@dataclass(frozen=True)
class Segment:
    """Metagene segment definition (name + bin count)."""

    name: str
    n_bins: int


_DEFAULT_SEGMENTS: tuple[Segment, ...] = (Segment("region", 100),)


def segments_total_bins(segments: Sequence[Segment]) -> int:
    """Return total bin count for all segments.

    Raises
    ------
    ValueError
        If segments are empty or contain non-positive bin counts.
    """
    if not segments:
        raise ValueError("segments must not be empty")
    if any(s.n_bins <= 0 for s in segments):
        raise ValueError("each segment n_bins must be > 0")
    return sum(s.n_bins for s in segments)


def segment_boundaries(segments: Sequence[Segment]) -> List[float]:
    total = segments_total_bins(segments)
    cum = 0
    bounds: List[float] = []
    for s in segments:
        cum += s.n_bins
        bounds.append(cum / total)
    return bounds


def segment_ticks(segments: Sequence[Segment]) -> Tuple[List[float], List[str]]:
    bounds = segment_boundaries(segments)
    labels = [s.name for s in segments]
    return bounds, labels


def _line_from_points(drd: DiscreteRegionData, segments: Sequence[Segment]) -> Tuple[np.ndarray, np.ndarray]:
    bounds = segment_boundaries(segments)
    total_bins = segments_total_bins(segments)
    starts = np.array([0.0] + bounds[:-1], dtype=float)
    ends = np.array(bounds, dtype=float)
    seg_nbins = np.array([s.n_bins for s in segments], dtype=int)
    seg_offsets = np.concatenate(([0], np.cumsum(seg_nbins)[:-1]))

    centers = []
    for s, start, end in zip(segments, starts, ends):
        width = end - start
        if width <= 0:
            continue
        idx = np.arange(s.n_bins, dtype=float)
        centers.append(start + (idx + 0.5) / s.n_bins * width)
    x_out = np.concatenate(centers) if centers else np.array([], dtype=float)

    streams = []
    for pos, dens in zip(drd.positions, drd.densities):
        x = np.asarray(pos, dtype=float)
        y = np.asarray(dens, dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            continue
        pts = list(zip(x[mask].tolist(), y[mask].tolist()))
        pts.sort(key=lambda t: t[0])
        streams.append(pts)

    if not streams:
        return np.array([]), np.array([])

    bins_values = [[] for _ in range(total_bins)]
    for x, y in heapq.merge(*streams, key=lambda t: t[0]):
        if x < 0.0 or x > 1.0:
            continue
        seg_idx = int(np.searchsorted(ends, x, side="right"))
        seg_idx = min(seg_idx, len(ends) - 1)
        width = ends[seg_idx] - starts[seg_idx]
        if width <= 0:
            continue
        x_seg = (x - starts[seg_idx]) / width
        local = int(x_seg * seg_nbins[seg_idx])
        local = min(local, seg_nbins[seg_idx] - 1)
        global_bin = int(seg_offsets[seg_idx] + local)
        if 0 <= global_bin < total_bins:
            bins_values[global_bin].append(y)

    y_out = []
    for vals in bins_values:
        if not vals:
            y_out.append(np.nan)
            continue
        y_out.append(float(np.nanmean(np.asarray(vals, dtype=float))))

    return x_out, np.array(y_out, dtype=float)


def _make_odd(n: int) -> int:
    if n % 2 == 0:
        return n + 1
    return n


def _coerce_smooth_config(
    smooth: dict | int | None,
    *,
    total_bins: int,
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
        window_length = _make_odd(max(3, total_bins // smooth))
        return {
            "method": "savgol",
            "window_length": window_length,
            "polyorder": 2,
            "apply": "post",
            "mode": "interp",
            "nan_policy": "interp",
            "per_segment": False,
        }
    if not isinstance(smooth, dict):
        raise ValueError("smooth must be a dict, int, or None")

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
        "cval": smooth.get("cval"),
    }
    return cfg


def _validate_savgol_config(
    cfg: dict,
    *,
    series_len: int,
    label: str = "series",
) -> None:
    window_length = cfg.get("window_length")
    polyorder = cfg.get("polyorder")
    mode = cfg.get("mode", "interp")

    if not isinstance(window_length, int) or window_length < 3:
        raise ValueError("smooth.window_length must be an odd integer >= 3")
    if window_length % 2 == 0:
        raise ValueError("smooth.window_length must be an odd integer")
    if not isinstance(polyorder, int) or polyorder < 0:
        raise ValueError("smooth.polyorder must be an integer >= 0")
    if polyorder >= window_length:
        raise ValueError("smooth.polyorder must be < smooth.window_length")

    if mode not in {"interp", "nearest", "mirror", "constant", "wrap"}:
        raise ValueError("smooth.mode must be one of: interp, nearest, mirror, constant, wrap")
    if cfg.get("nan_policy") not in {"interp", "mask", "raise"}:
        raise ValueError("smooth.nan_policy must be one of: interp, mask, raise")

    if mode == "interp" and window_length > series_len:
        raise ValueError(
            f"smooth.window_length ({window_length}) must be <= {label} length ({series_len}) when mode='interp'"
        )


def _savgol_filter_1d(y: np.ndarray, cfg: dict) -> np.ndarray:
    try:
        from scipy.signal import savgol_filter  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("scipy is required for Savitzky–Golay smoothing") from e

    window_length = int(cfg["window_length"])
    polyorder = int(cfg["polyorder"])
    mode = cfg.get("mode", "interp")
    cval = cfg.get("cval", 0.0)

    _validate_savgol_config(cfg, series_len=len(y))

    nan_policy = cfg.get("nan_policy", "interp")
    if nan_policy == "raise" and np.isnan(y).any():
        raise ValueError("NaN values present; use nan_policy='interp' or 'mask'")

    if y.size == 0:
        return y

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
    segments: Sequence[Segment] | None = None,
) -> np.ndarray:
    if not isinstance(y, np.ndarray):
        y = np.asarray(y, dtype=float)
    if y.size == 0:
        return y

    if cfg.get("per_segment"):
        if not segments:
            raise ValueError("per_segment=True требует segments")
        seg_nbins = [s.n_bins for s in segments]
        if sum(seg_nbins) != y.size:
            raise ValueError("per_segment=True требует профиль длиной total_bins")
        out = []
        offset = 0
        for idx, n in enumerate(seg_nbins):
            seg = y[offset : offset + n]
            _validate_savgol_config(cfg, series_len=len(seg), label=f"segment {idx}")
            out.append(_savgol_filter_1d(seg, cfg))
            offset += n
        return np.concatenate(out) if out else y

    return _savgol_filter_1d(y, cfg)


def _clip_profile(y: np.ndarray) -> np.ndarray:
    if y.size == 0:
        return y
    y_out = y.astype(float, copy=True)
    mask = np.isfinite(y_out)
    if np.any(mask):
        y_out[mask] = np.clip(y_out[mask], 0.0, 1.0)
    return y_out


def _is_negative_strand(contig) -> bool:
    val = None
    if hasattr(contig, "strand_str"):
        try:
            v = contig.strand_str
            val = v() if callable(v) else v
        except (AttributeError, TypeError, ValueError):
            val = None
    if val is None and hasattr(contig, "strand"):
        try:
            val = str(contig.strand)
        except (AttributeError, TypeError, ValueError):
            pass
    return str(val).strip() == "-"


def compute_discrete_regions(
    reader: _io.RegionReader,
    contigs: Sequence,
    *,
    segments: Sequence[Segment] = _DEFAULT_SEGMENTS,
    agg_method=None,
    reverse_negative: bool = True,
    labels: Optional[Sequence[str]] = None,
    mode: str = "raw",
    x_mode: str = "relative",
    progress: bool = False,
    progress_every: int = 1,
) -> DiscreteRegionData:
    """Compute discretised metagene profiles for contigs.

    Parameters
    ----------
    reader
        RegionReader that implements ``iter_contigs``.
    contigs
        Contigs to extract.
    segments
        Segmentation scheme (default single ``Segment("region", 100)``).
    agg_method
        Aggregation method for ``BsxBatch.discretise``.
    reverse_negative
        Flip negative-strand profiles if True.
    labels
        Optional labels for resulting regions.
    """
    try:
        sorted_contigs = reader.index().sort(list(contigs))
    except Exception:
        sorted_contigs = list(contigs)
    else:
        if labels is not None:
            id_to_label = {id(c): lbl for c, lbl in zip(contigs, labels)}
            labels = [id_to_label.get(id(c)) for c in sorted_contigs]
    contigs_list = sorted_contigs
    total_bins = segments_total_bins(segments)

    if agg_method is None:
        from importlib import import_module

        agg_method = getattr(import_module("bsx2._bsx2"), "AggMethod").Mean

    data = DiscreteRegionData()
    batches_iter = getattr(reader, "iter_contigs", None)
    if not callable(batches_iter):
        raise AttributeError("reader must implement iter_contigs(contigs)")

    total = len(contigs_list)
    step = progress_every if progress_every > 0 else 1

    def _progress(idx: int) -> None:
        if not progress or total <= 0:
            return
        if idx % step != 0 and idx + 1 != total:
            return
        width = 30
        frac = (idx + 1) / total
        filled = int(round(frac * width))
        bar = "#" * filled + "-" * (width - filled)
        print(f"\r[{bar}] {idx + 1}/{total}", end="", flush=True)
        if idx + 1 == total:
            print()

    for idx, batch in enumerate(reader.iter_contigs(contigs_list)):
        contig = contigs_list[idx] if idx < len(contigs_list) else None
        if mode != "raw":
            raise ValueError("Only mode='raw' is supported (discretise mode removed)")
        if mode == "raw":
            if contig is None:
                _progress(idx)
                continue
            start = getattr(contig, "start", None)
            end = getattr(contig, "end", None)
            if start is None or end is None or end <= start:
                _progress(idx)
                continue
            try:
                pos = np.asarray(batch.position().to_list(), dtype=np.float64)
                dens = np.asarray(batch.density().to_list(), dtype=np.float64)
                weights = np.asarray(batch.count_total().to_list(), dtype=np.float64)
            except Exception:
                _progress(idx)
                continue
            if pos.size == 0 or dens.size == 0:
                _progress(idx)
                continue
            if weights.size != dens.size:
                _progress(idx)
                continue
            if x_mode == "absolute":
                x = pos - float(start)
            else:
                x = (pos - float(start)) / float(end - start)
            y = dens
            if reverse_negative and _is_negative_strand(contig):
                if x_mode == "absolute":
                    x = float(end - start) - x
                else:
                    x = 1.0 - x
            mask = np.isfinite(x) & np.isfinite(y)
            if x_mode != "absolute":
                mask = mask & (x >= 0.0) & (x <= 1.0)
            x = x[mask]
            y = y[mask]
            w = weights[mask]
            if x.size == 0:
                _progress(idx)
                continue
            order = np.argsort(x, kind="mergesort")
            x = x[order]
            y = y[order]
            w = w[order]
        else:
            # unreachable due to guard above
            _progress(idx)
            continue
        # Clean only infinities; keep NaN as "no data"
        if np.any(~np.isfinite(y)):
            y = y.astype(float, copy=True)
            y[~np.isfinite(y)] = np.nan
        mask = np.isfinite(y)
        if np.any(mask):
            y = y.astype(float, copy=False)
            y[mask] = np.clip(y[mask], 0.0, 1.0)
        label = labels[idx] if labels and idx < len(labels) else None
        data.insert(x, y, label, weights=w if mode == "raw" else None)
        _progress(idx)

    return data


def _get_contig_start(contig) -> Optional[int]:
    v = getattr(contig, "start", None)
    if v is None:
        return None
    try:
        return int(v() if callable(v) else v)
    except (TypeError, ValueError):
        return None


def collect_contigs_from_hcannot(
    annot,
    *,
    feature_type: Optional[str] = None,
    limit: Optional[int] = None,
    label_getter: Optional[Callable[[object, int], str]] = None,
) -> Tuple[List[object], List[str]]:
    """Extract contigs and labels from an HcAnnotStore-like object."""

    def _default_label(entry, idx: int) -> str:
        for attr in ("id", "get_id"):
            v = getattr(entry, attr, None)
            if v is not None:
                lab = v() if callable(v) else v
                if lab:
                    return str(lab)
        ft = getattr(entry, "feature_type", None)
        return f"{ft}_{idx}" if ft else f"entry_{idx}"

    label_getter = label_getter or _default_label

    entries = []
    it = getattr(annot, "iter", None)
    iterator = it() if callable(it) else iter(annot)

    idx = 0
    while True:
        try:
            item = next(iterator)
        except StopIteration:
            break
        entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item

        ft = getattr(entry, "feature_type", None)
        if feature_type and ft != feature_type:
            continue

        contig = getattr(entry, "contig", None)
        if contig is None:
            getter = getattr(entry, "get_contig", None)
            contig = getter() if callable(getter) else None
        if contig is None:
            continue
        start = _get_contig_start(contig)
        if start is not None and start < 1:
            continue

        label = label_getter(entry, idx)
        entries.append((contig, label))
        idx += 1
        if limit and idx >= limit:
            break

    contigs, labels = zip(*entries) if entries else ([], [])
    return list(contigs), list(labels)


def _entry_feature_type(entry) -> Optional[str]:
    v = getattr(entry, "feature_type", None)
    if v is not None:
        try:
            return str(v() if callable(v) else v)
        except (TypeError, ValueError):
            pass
    getter = getattr(entry, "get_feature_type", None)
    if callable(getter):
        try:
            return str(getter())
        except (TypeError, ValueError):
            return None
    return None


def _entry_id(entry) -> Optional[str]:
    for attr in ("id", "get_id"):
        v = getattr(entry, attr, None)
        if v is None:
            continue
        try:
            val = v() if callable(v) else v
        except (TypeError, ValueError):
            continue
        if val:
            return str(val)
    return None


def _entry_parents(entry) -> List[str]:
    attrs = getattr(entry, "attributes", None)
    if attrs is not None and callable(attrs):
        try:
            attrs = attrs()
        except Exception:
            attrs = None
    if attrs is not None:
        p = getattr(attrs, "parent", None)
        if p is not None:
            try:
                parents = p() if callable(p) else p
                if parents:
                    return [str(x) for x in parents]
            except Exception:
                pass
        getter = getattr(attrs, "get_parent", None)
        if callable(getter):
            try:
                parents = getter()
                if parents:
                    return [str(x) for x in parents]
            except Exception:
                pass
    getter = getattr(entry, "get_attributes", None)
    if callable(getter):
        try:
            attrs = getter()
        except Exception:
            attrs = None
        if attrs is not None:
            getter = getattr(attrs, "get_parent", None)
            if callable(getter):
                try:
                    parents = getter()
                    if parents:
                        return [str(x) for x in parents]
                except Exception:
                    pass
            p = getattr(attrs, "parent", None)
            if p is not None:
                try:
                    parents = p() if callable(p) else p
                    if parents:
                        return [str(x) for x in parents]
                except Exception:
                    pass
    return []


def collect_parts_from_hcannot(
    annot,
    *,
    parts: Sequence[str],
    limit: Optional[int] = None,
) -> dict[str, tuple[list[object], list[str]]]:
    parts_set = set(parts)
    out: dict[str, list[tuple[object, str]]] = {p: [] for p in parts_set}

    it = getattr(annot, "iter", None)
    iterator = it() if callable(it) else iter(annot)

    gene_ids: list[str] = []
    if "gene" in parts_set:
        idx = 0
        for item in iterator:
            entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item
            if _entry_feature_type(entry) != "gene":
                continue
            gid = _entry_id(entry)
            if not gid:
                continue
            gene_ids.append(gid)
            idx += 1
            if limit and idx >= limit:
                break
        allowed = set(gene_ids)
    else:
        allowed = None

    it = getattr(annot, "iter", None)
    iterator = it() if callable(it) else iter(annot)

    for item in iterator:
        entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item
        ft = _entry_feature_type(entry)
        if ft is None or ft not in parts_set:
            continue

        contig = getattr(entry, "contig", None)
        if contig is None:
            getter = getattr(entry, "get_contig", None)
            contig = getter() if callable(getter) else None
        if contig is None:
            continue
        start = _get_contig_start(contig)
        if start is not None and start < 1:
            continue

        if ft == "gene":
            gid = _entry_id(entry)
        else:
            parents = _entry_parents(entry)
            gid = parents[0] if parents else None
        if not gid:
            continue
        if allowed is not None and gid not in allowed:
            continue

        out[ft].append((contig, gid))

    result: dict[str, tuple[list[object], list[str]]] = {}
    for ft, items in out.items():
        if not items:
            result[ft] = ([], [])
        else:
            contigs, labels = zip(*items)
            result[ft] = (list(contigs), list(labels))
    return result


def combine_parts_drd(
    drd_map: dict[str, DiscreteRegionData],
    *,
    segments: Sequence[Segment],
    parts_order: Sequence[str],
) -> DiscreteRegionData:
    if len(segments) != 3:
        raise ValueError("combined metagene requires exactly 3 segments (up/body/down)")
    bounds = segment_boundaries(segments)
    starts = [0.0, bounds[0], bounds[1]]
    ends = [bounds[0], bounds[1], bounds[2]]

    by_part: dict[str, dict[str, tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]]] = {}
    labels_all: set[str] = set()
    for part, drd in drd_map.items():
        lookup: dict[str, tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]] = {}
        for pos, dens, w, lbl in zip(drd.positions, drd.densities, drd.weights, drd.labels):
            if lbl is None:
                continue
            lookup[lbl] = (np.asarray(pos, dtype=float), np.asarray(dens, dtype=float), w)
            labels_all.add(lbl)
        by_part[part] = lookup

    out = DiscreteRegionData()
    for lbl in sorted(labels_all):
        xs = []
        ys = []
        ws = []
        weights_ok = True
        for idx, part in enumerate(parts_order):
            lookup = by_part.get(part, {})
            if lbl not in lookup:
                continue
            x, y, w = lookup[lbl]
            width = ends[idx] - starts[idx]
            if width <= 0:
                continue
            x_mapped = starts[idx] + x * width
            xs.append(x_mapped)
            ys.append(y)
            if w is None:
                weights_ok = False
            ws.append(w)
        if not xs:
            continue
        x_all = np.concatenate(xs)
        y_all = np.concatenate(ys)
        if weights_ok and ws:
            w_all = np.concatenate([np.asarray(w, dtype=float) for w in ws if w is not None])
        else:
            w_all = None
        out.insert(x_all, y_all, lbl, weights=w_all)
    return out


def compute_from_annot(
    reader: _io.RegionReader,
    annot,
    *,
    segments: Sequence[Segment] = _DEFAULT_SEGMENTS,
    agg_method=None,
    feature_type: Optional[str] = None,
    reverse_negative: bool = True,
    labels: Optional[Sequence[str]] = None,
) -> DiscreteRegionData:
    """Build DiscreteRegionData from an annotation store."""
    contigs, auto_labels = collect_contigs_from_hcannot(annot, feature_type=feature_type)
    if labels is None:
        labels = auto_labels
    return compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        agg_method=agg_method,
        reverse_negative=reverse_negative,
        labels=labels,
    )


def line_plot(
    reader: _io.RegionReader,
    *,
    contigs: Sequence,
    segments: Sequence[Segment] | None = None,
    agg_method=None,
    mode: str = "raw",
    x_mode: str = "relative",
    smooth: dict | int | None = None,
):
    """HoloViews Curve for averaged metagene profile.

    smooth
        Optional Savitzky–Golay smoothing configuration. Accepts:
        - dict with keys:
          - method: "savgol" (required)
          - window_length: int (odd, >=3)
          - polyorder: int (>=0, < window_length)
          - apply: "post" or "pre" (default "post")
          - mode: "interp" | "nearest" | "mirror" | "constant" | "wrap" (default "interp")
          - nan_policy: "interp" | "mask" | "raise" (default "interp")
          - per_segment: bool (default False)
          - cval: float (optional, for mode="constant")
        - int (legacy): number of windows; 0 disables smoothing.
    """
    segments = segments or _DEFAULT_SEGMENTS
    bounds, names = segment_ticks(segments)
    smooth_cfg = _coerce_smooth_config(smooth, total_bins=segments_total_bins(segments))
    if smooth_cfg is not None:
        apply_mode = smooth_cfg.get("apply", "post")
        if apply_mode not in {"pre", "post"}:
            raise ValueError("smooth.apply must be 'pre' or 'post'")
    drd = compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        agg_method=agg_method,
        mode=mode,
        x_mode=x_mode,
    )
    if smooth_cfg is not None and smooth_cfg.get("apply", "post") == "pre":
        if smooth_cfg.get("per_segment") and mode == "raw":
            raise ValueError("per_segment=True не поддерживается для mode='raw' с apply='pre'")
        smoothed = DiscreteRegionData()
        for pos, dens, w, lbl in zip(drd.positions, drd.densities, drd.weights, drd.labels):
            y_sm = _apply_savgol_smoothing(np.asarray(dens, dtype=float), smooth_cfg, segments=segments)
            y_sm = _clip_profile(y_sm)
            smoothed.insert(np.asarray(pos, dtype=float), y_sm, lbl, weights=w)
        drd = smoothed
    x, y = _line_from_points(drd, segments)
    if smooth_cfg is not None and smooth_cfg.get("apply", "post") == "post":
        y = _apply_savgol_smoothing(np.asarray(y, dtype=float), smooth_cfg, segments=segments)
        y = _clip_profile(y)
    return LinePlotData(x=x, y=y, x_ticks=bounds, x_labels=names).to_curve()


def _stack_for_heatmap(drd: DiscreteRegionData) -> Tuple[pd.DataFrame, int]:
    mat, row_labels = drd.stack_matrix()
    if mat.size == 0:
        return pd.DataFrame(columns=["bin", "region", "density"]), 0
    n_regions, n_bins = mat.shape
    df = pd.DataFrame(
        {"bin": np.tile(np.arange(n_bins, dtype=int), n_regions), "region": np.repeat(row_labels, n_bins), "density": mat.reshape(-1)}
    )
    return df, n_bins


def heatmap(reader: _io.RegionReader, *, contigs: Sequence, segments: Sequence[Segment] | None = None, agg_method=None):
    """HoloViews HeatMap (regions x bins)."""
    try:
        import holoviews as hv  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("holoviews is required for heatmap; install with 'pip install holoviews'") from e
    segments = segments or _DEFAULT_SEGMENTS
    drd = compute_discrete_regions(reader, contigs, segments=segments, agg_method=agg_method)
    df, n_bins = _stack_for_heatmap(drd)
    if df.empty:
        return hv.HeatMap([])
    hm = hv.HeatMap(df, kdims=["bin", "region"], vdims=["density"]).opts(invert_yaxis=True, colorbar=True)
    bounds, names = segment_ticks(segments)
    if n_bins > 0 and bounds:
        xticks = [(int(round(b * (n_bins - 1))), label) for b, label in zip(bounds[:-1], names[:-1])]
        if xticks:
            hm = hm.opts(xticks=xticks)
    return hm


def box_plot(reader: _io.RegionReader, *, contigs: Sequence, segments: Sequence[Segment] | None = None, agg_method=None):
    """HoloViews BoxWhisker per-bin distributions."""
    try:
        import holoviews as hv  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("holoviews is required for box_plot; install with 'pip install holoviews'") from e
    segments = segments or _DEFAULT_SEGMENTS
    drd = compute_discrete_regions(reader, contigs, segments=segments, agg_method=agg_method)
    mat, _ = drd.stack_matrix()
    if mat.size == 0:
        return hv.BoxWhisker([])
    n_regions, n_bins = mat.shape
    df = pd.DataFrame({"bin": np.tile(np.arange(n_bins), n_regions), "density": mat.reshape(-1)})
    return hv.BoxWhisker(df, kdims=["bin"], vdims=["density"])


def violin_plot(reader: _io.RegionReader, *, contigs: Sequence, segments: Sequence[Segment] | None = None, agg_method=None):
    """HoloViews Violin per-bin distributions."""
    try:
        import holoviews as hv  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("holoviews is required for violin_plot; install with 'pip install holoviews'") from e
    segments = segments or _DEFAULT_SEGMENTS
    drd = compute_discrete_regions(reader, contigs, segments=segments, agg_method=agg_method)
    mat, _ = drd.stack_matrix()
    if mat.size == 0:
        return hv.Violin([])
    n_regions, n_bins = mat.shape
    df = pd.DataFrame({"bin": np.tile(np.arange(n_bins), n_regions), "density": mat.reshape(-1)})
    return hv.Violin(df, kdims=["bin"], vdims=["density"])
