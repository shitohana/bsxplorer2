from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict, deque
from typing import Callable, List, Optional, Sequence, Tuple, cast

import numpy as np

from bsx2 import io as _io
from bsx2.guards import (
    require_nan_policy_compatible,
    require_per_segment_profile,
    require_savgol_filter,
    resolve_reader_accessors,
)
from bsx2.plots.data import DiscreteRegionData
from bsx2.validation import validate_segments, validate_smoothing


@dataclass(frozen=True)
class MetageneProfileSegment:
    """Metagene segment definition (name + bin count)."""

    name: str
    n_bins: int


_DEFAULT_SEGMENTS: tuple[MetageneProfileSegment, ...] = (MetageneProfileSegment("region", 100),)


def segments_total_bins(segments: Sequence[MetageneProfileSegment]) -> int:
    """Return total bin count for all segments.

    Raises
    ------
    ValueError
        If segments are empty or contain non-positive bin counts.
    """
    return validate_segments(segments)


def segment_boundaries(segments: Sequence[MetageneProfileSegment]) -> List[float]:
    total = segments_total_bins(segments)
    cum = 0
    bounds: List[float] = []
    for s in segments:
        cum += s.n_bins
        bounds.append(cum / total)
    return bounds


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


def _savgol_filter_1d(y: np.ndarray, cfg: dict) -> np.ndarray:
    savgol_filter = require_savgol_filter()

    window_length = int(cfg["window_length"])
    polyorder = int(cfg["polyorder"])
    mode = cfg.get("mode", "interp")
    cval = cfg.get("cval", 0.0)

    _validate_savgol_config(cfg, series_len=len(y))

    nan_policy = cfg.get("nan_policy", "interp")
    require_nan_policy_compatible(y, nan_policy=nan_policy)

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
    segments: Sequence[MetageneProfileSegment] = _DEFAULT_SEGMENTS,
    reverse_negative: bool = True,
    labels: Optional[Sequence[str]] = None,
    progress: bool = False,
    progress_every: int = 1,
) -> DiscreteRegionData:
    """Compute relative metagene profiles (0..1) for contigs.

    Parameters
    ----------
    reader
        RegionReader that implements ``iter_contigs``.
    contigs
        Contigs to extract.
    segments
        Segmentation scheme (default single ``Segment("region", 100)``).
    reverse_negative
        Flip negative-strand profiles if True.
    labels
        Optional labels for resulting regions.
    """
    def _v(val):
        try:
            return val() if callable(val) else val
        except Exception:
            return None

    def _contig_key(c):
        seq = _v(getattr(c, "seqname", None))
        start = _v(getattr(c, "start", None))
        end = _v(getattr(c, "end", None))
        strand = getattr(c, "strand_str", None)
        if strand is not None:
            strand = _v(strand)
        else:
            strand = _v(getattr(c, "strand", None))
        return (
            str(seq),
            int(start) if start is not None else None,
            int(end) if end is not None else None,
            str(strand),
        )

    try:
        sorted_contigs = reader.index().sort(list(contigs))
    except Exception:
        sorted_contigs = list(contigs)
    else:
        if labels is not None:
            buckets: dict[tuple, deque] = defaultdict(deque)
            for c, lbl in zip(contigs, labels):
                buckets[_contig_key(c)].append(lbl)
            new_labels = []
            for c in sorted_contigs:
                q = buckets.get(_contig_key(c))
                new_labels.append(q.popleft() if q else None)
            labels = new_labels
    contigs_list = sorted_contigs
    segments_total_bins(segments)

    data = DiscreteRegionData()
    query_fn, iter_contigs_fn = resolve_reader_accessors(reader)

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

    if query_fn is not None:
        contig_iter = list(enumerate(contigs_list))
        batch_iter = None
        last_seqname = None
    else:
        contig_iter = list(enumerate(contigs_list))
        batch_iter = cast(Callable[[list[object]], object], iter_contigs_fn)(contigs_list)

    for idx, contig in contig_iter:
        if query_fn is not None:
            seqname = getattr(contig, "seqname", None)
            try:
                seqname = seqname() if callable(seqname) else seqname
            except Exception:
                seqname = None
            if seqname is not None and seqname != last_seqname:
                try:
                    reset_fn = getattr(reader, "reset", None)
                    if callable(reset_fn):
                        reset_fn()
                except Exception:
                    pass
                last_seqname = seqname
            try:
                batch = query_fn(contig)
            except Exception:
                _progress(idx)
                continue
            if batch is None:
                _progress(idx)
                continue
        else:
            try:
                batch = next(batch_iter)
            except StopIteration:
                break
            except Exception:
                _progress(idx)
                continue
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
        except Exception:
            _progress(idx)
            continue
        if pos.size == 0 or dens.size == 0:
            _progress(idx)
            continue
        x = (pos - float(start)) / float(end - start)
        y = dens
        if reverse_negative and _is_negative_strand(contig):
            x = 1.0 - x
        # Keep NaN in y as "no data"; drop only non-finite/out-of-range x.
        mask = np.isfinite(x) & (x >= 0.0) & (x <= 1.0)
        x = x[mask]
        y = y[mask]
        if x.size == 0:
            _progress(idx)
            continue
        order = np.argsort(x, kind="mergesort")
        x = x[order]
        y = y[order]
        # Clean only infinities; keep NaN as "no data"
        if np.any(~np.isfinite(y)):
            y = y.astype(float, copy=True)
            y[~np.isfinite(y)] = np.nan
        mask = np.isfinite(y)
        if np.any(mask):
            y = y.astype(float, copy=False)
            y[mask] = np.clip(y[mask], 0.0, 1.0)
        label = labels[idx] if labels and idx < len(labels) else None
        data.insert(x, y, label)
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
        if start is not None and start < 0:
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
        if start is not None and start < 0:
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
    segments: Sequence[MetageneProfileSegment],
    parts_order: Sequence[str],
) -> DiscreteRegionData:
    validate_segments(segments, expected_len=3)
    bounds = segment_boundaries(segments)
    starts = [0.0, bounds[0], bounds[1]]
    ends = [bounds[0], bounds[1], bounds[2]]

    by_part: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    labels_all: set[str] = set()
    for part, drd in drd_map.items():
        lookup: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for pos, dens, lbl in zip(drd.positions, drd.densities, drd.labels):
            if lbl is None:
                continue
            lookup[lbl] = (np.asarray(pos, dtype=float), np.asarray(dens, dtype=float))
            labels_all.add(lbl)
        by_part[part] = lookup

    out = DiscreteRegionData()
    for lbl in sorted(labels_all):
        xs = []
        ys = []
        for idx, part in enumerate(parts_order):
            lookup = by_part.get(part, {})
            if lbl not in lookup:
                continue
            x, y = lookup[lbl]
            width = ends[idx] - starts[idx]
            if width <= 0:
                continue
            x_mapped = starts[idx] + x * width
            xs.append(x_mapped)
            ys.append(y)
        if not xs:
            continue
        x_all = np.concatenate(xs)
        y_all = np.concatenate(ys)
        out.insert(x_all, y_all, lbl)
    return out


def compute_from_annot(
    reader: _io.RegionReader,
    annot,
    *,
    segments: list[MetageneProfileSegment] | None = None,
    feature_type: str | None = None,
    reverse_negative: bool = True,
    labels: list[str] | None = None,
    limit: int | None = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    combine_parts: bool = False,
    parts: Sequence[str] | None = None,
) -> DiscreteRegionData:
    """Build DiscreteRegionData from an annotation store.

    When ``combine_parts=True`` this will build a BSX1-like metagene by
    collecting gene parts (upstream/gene/downstream) and stitching them
    into a single profile per gene.
    """
    if segments is None:
        if combine_parts:
            segments = [MetageneProfileSegment("up", 100), MetageneProfileSegment("body", 200), MetageneProfileSegment("down", 100)]
        else:
            segments = [MetageneProfileSegment("region", 100)]
    validate_segments(
        segments,
        expected_len=3 if combine_parts else None,
        flank_bp=flank_bp,
    )

    if add_flanks:
        try:
            ft_map = annot.get_feature_types()
        except Exception:
            ft_map = {}
        gene_ids = ft_map.get("gene", []) if isinstance(ft_map, dict) else []
        if gene_ids:
            flank = int(abs(flank_bp))
            if "upstream_gene" not in ft_map:
                r = annot.add_flanks(gene_ids, -flank, "upstream_")
                if r is not None:
                    annot = r
            if "downstream_gene" not in ft_map:
                r = annot.add_flanks(gene_ids, flank, "downstream_")
                if r is not None:
                    annot = r

    if combine_parts:
        parts_order = list(parts) if parts is not None else ["upstream_gene", "gene", "downstream_gene"]
        parts_data = collect_parts_from_hcannot(annot, parts=parts_order, limit=limit)
        drd_map: dict[str, DiscreteRegionData] = {}
        for part in parts_order:
            contigs, auto_labels = parts_data.get(part, ([], []))
            if not contigs:
                continue
            drd_part = compute_discrete_regions(
                reader,
                contigs,
                segments=segments,
                reverse_negative=reverse_negative,
                labels=auto_labels,
            )
            drd_map[part] = drd_part
        if not drd_map:
            return DiscreteRegionData()
        return combine_parts_drd(drd_map, segments=segments, parts_order=parts_order)

    contigs, auto_labels = collect_contigs_from_hcannot(annot, feature_type=feature_type, limit=limit)
    use_labels = labels if labels is not None else auto_labels
    return compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        reverse_negative=reverse_negative,
        labels=use_labels,
    )
