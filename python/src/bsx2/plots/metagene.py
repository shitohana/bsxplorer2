from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict, deque
from typing import TYPE_CHECKING, Callable, List, Optional, Sequence, Tuple, cast

import numpy as np
from bsx2 import Contig, RegionReader
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


def _insert_region_data(
    data: DiscreteRegionData,
    pos: np.ndarray,
    dens: np.ndarray,
    *,
    contig: Contig,
    start: int | float,
    end: int | float,
    label: Optional[str],
    reverse_negative: bool,
) -> None:
    x = (pos - float(start)) / float(end - start)
    y = dens

    if reverse_negative and _is_negative_strand(contig):
        x = 1.0 - x[::-1]
        y = y[::-1]

    data.insert_unchecked(x, y, label)


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

            label = labels[idx] if labels and idx < len(labels) else None
            _insert_region_data(
                data,
                pos,
                dens,
                contig=contig,
                start=start,
                end=end,
                label=label,
                reverse_negative=reverse_negative,
            )

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

            label = labels[idx] if labels and idx < len(labels) else None
            _insert_region_data(
                data,
                pos,
                dens,
                contig=contig,
                start=start,
                end=end,
                label=label,
                reverse_negative=reverse_negative,
            )

    return data


def _get_contig_start(contig: object) -> Optional[int]:
    v = getattr(contig, "start", None)
    if v is None:
        return None
    try:
        return int(v() if callable(v) else v)
    except (TypeError, ValueError):
        return None


def collect_contigs_from_hcannot(
    annot: object,
    *,
    feature_type: Optional[str] = None,
    limit: Optional[int] = None,
    label_getter: Optional[Callable[[object, int], str]] = None,
) -> Tuple[List[object], List[str]]:
    def _default_label(entry: object, idx: int) -> str:
        for attr in ("id", "get_id"):
            v = getattr(entry, attr, None)
            if v is None:
                continue
            try:
                label = v() if callable(v) else v
            except (TypeError, ValueError):
                continue
            if label:
                return str(label)

        ft = _entry_feature_type(entry)
        return f"{ft}_{idx}" if ft else f"entry_{idx}"

    label_getter = label_getter or _default_label

    entries: List[Tuple[object, str]] = []
    iterator_factory = getattr(annot, "iter", None)
    iterator = iterator_factory() if callable(iterator_factory) else iter(annot)

    idx = 0
    while True:
        try:
            item = next(iterator)
        except StopIteration:
            break

        entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item

        ft = _entry_feature_type(entry)
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

        entries.append((contig, label_getter(entry, idx)))
        idx += 1
        if limit and idx >= limit:
            break

    contigs, labels = zip(*entries) if entries else ([], [])
    return list(contigs), list(labels)


def _entry_feature_type(entry: object) -> Optional[str]:
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


def _entry_id(entry: object) -> Optional[str]:
    for attr in ("id", "get_id"):
        v = getattr(entry, attr, None)
        if v is None:
            continue
        try:
            value = v() if callable(v) else v
        except (TypeError, ValueError):
            continue
        if value:
            return str(value)
    return None


def _entry_parents(entry: object) -> List[str]:
    attrs = getattr(entry, "attributes", None)
    if attrs is not None and callable(attrs):
        try:
            attrs = attrs()
        except Exception:
            attrs = None

    if attrs is not None:
        parent = getattr(attrs, "parent", None)
        if parent is not None:
            try:
                parents = parent() if callable(parent) else parent
            except Exception:
                parents = None
            if parents:
                return [str(x) for x in parents]

        getter = getattr(attrs, "get_parent", None)
        if callable(getter):
            try:
                parents = getter()
            except Exception:
                parents = None
            if parents:
                return [str(x) for x in parents]

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
                except Exception:
                    parents = None
                if parents:
                    return [str(x) for x in parents]

            parent = getattr(attrs, "parent", None)
            if parent is not None:
                try:
                    parents = parent() if callable(parent) else parent
                except Exception:
                    parents = None
                if parents:
                    return [str(x) for x in parents]

    return []


def _entry_contig(entry: object) -> Optional[Contig]:
    contig = getattr(entry, "contig", None)
    if contig is None:
        getter = getattr(entry, "get_contig", None)
        contig = getter() if callable(getter) else None
    return cast(Optional[Contig], contig)


def _make_flank_contig(
    contig: Contig,
    *,
    flank_bp: int,
    kind: str,
) -> Optional[Contig]:
    if flank_bp <= 0:
        return None

    start = int(getattr(contig, "start"))
    end = int(getattr(contig, "end"))
    strand = getattr(contig, "strand")
    is_negative = _is_negative_strand(contig)

    if kind == "upstream_gene":
        if is_negative:
            flank_start = end
            flank_end = end + flank_bp
        else:
            flank_start = max(0, start - flank_bp)
            flank_end = start
    elif kind == "downstream_gene":
        if is_negative:
            flank_start = max(0, start - flank_bp)
            flank_end = start
        else:
            flank_start = end
            flank_end = end + flank_bp
    else:
        return None

    if flank_end <= flank_start:
        return None

    return Contig(str(getattr(contig, "seqname")), flank_start, flank_end, strand)


def _collect_gene_contigs(
    annot: object,
    *,
    limit: Optional[int] = None,
) -> Tuple[List[Contig], List[str]]:
    contigs: List[Contig] = []
    labels: List[str] = []

    iterator_factory = getattr(annot, "iter", None)
    iterator = iterator_factory() if callable(iterator_factory) else iter(annot)

    for item in iterator:
        entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item
        if _entry_feature_type(entry) != "gene":
            continue

        contig = _entry_contig(entry)
        if contig is None:
            continue

        start = _get_contig_start(contig)
        if start is not None and start < 0:
            continue

        label = _entry_id(entry)
        if not label:
            label = f"{contig.seqname}:{contig.start}-{contig.end}"

        contigs.append(contig)
        labels.append(label)

        if limit and len(contigs) >= limit:
            break

    return contigs, labels


def _collect_synthetic_flanks(
    annot: object,
    *,
    flank_bp: int,
    kind: str,
    limit: Optional[int] = None,
) -> Tuple[List[Contig], List[str]]:
    gene_contigs, gene_labels = _collect_gene_contigs(annot, limit=limit)

    contigs: List[Contig] = []
    labels: List[str] = []
    for gene_contig, gene_label in zip(gene_contigs, gene_labels):
        flank_contig = _make_flank_contig(
            gene_contig,
            flank_bp=flank_bp,
            kind=kind,
        )
        if flank_contig is None:
            continue
        contigs.append(flank_contig)
        labels.append(gene_label)

    return contigs, labels


def collect_parts_from_hcannot(
    annot: object,
    *,
    parts: Sequence[str],
    limit: Optional[int] = None,
) -> dict[str, tuple[list[object], list[str]]]:
    parts_set = set(parts)
    out: dict[str, list[tuple[object, str]]] = {part: [] for part in parts_set}

    iterator_factory = getattr(annot, "iter", None)
    iterator = iterator_factory() if callable(iterator_factory) else iter(annot)

    gene_ids: List[str] = []
    if "gene" in parts_set:
        idx = 0
        for item in iterator:
            entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item
            if _entry_feature_type(entry) != "gene":
                continue
            gene_id = _entry_id(entry)
            if not gene_id:
                continue
            gene_ids.append(gene_id)
            idx += 1
            if limit and idx >= limit:
                break
        allowed = set(gene_ids)
    else:
        allowed = None

    iterator_factory = getattr(annot, "iter", None)
    iterator = iterator_factory() if callable(iterator_factory) else iter(annot)

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
            gene_id = _entry_id(entry)
        else:
            parents = _entry_parents(entry)
            gene_id = parents[0] if parents else None

        if not gene_id:
            continue
        if allowed is not None and gene_id not in allowed:
            continue

        out[ft].append((contig, gene_id))

    result: dict[str, tuple[list[object], list[str]]] = {}
    for ft, items in out.items():
        if not items:
            result[ft] = ([], [])
            continue
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
        for pos, dens, label in zip(drd.positions, drd.densities, drd.labels):
            if label is None:
                continue
            lookup[label] = (np.asarray(pos, dtype=float), np.asarray(dens, dtype=float))
            labels_all.add(label)
        by_part[part] = lookup

    out = DiscreteRegionData()
    for label in sorted(labels_all):
        xs = []
        ys = []
        for idx, part in enumerate(parts_order):
            lookup = by_part.get(part, {})
            if label not in lookup:
                continue
            x, y = lookup[label]
            width = ends[idx] - starts[idx]
            if width <= 0:
                continue
            xs.append(starts[idx] + x * width)
            ys.append(y)

        if not xs:
            continue

        x_all = np.concatenate(xs)
        y_all = np.concatenate(ys)
        out.insert(x_all, y_all, label)

    return out


def compute_from_annot(
    reader: RegionReader,
    annot: object,
    *,
    segments: Optional[Sequence[MetageneProfileSegment]] = None,
    feature_type: Optional[str] = None,
    reverse_negative: bool = True,
    labels: Optional[Sequence[str]] = None,
    limit: Optional[int] = None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    combine_parts: bool = False,
    parts: Optional[Sequence[str]] = None,
) -> DiscreteRegionData:
    if segments is None:
        if combine_parts:
            segments = (
                MetageneProfileSegment("up", 100),
                MetageneProfileSegment("body", 200),
                MetageneProfileSegment("down", 100),
            )
        else:
            segments = _DEFAULT_SEGMENTS

    validate_segments(
        segments,
        expected_len=3 if combine_parts else None,
        flank_bp=flank_bp,
    )

    if combine_parts:
        parts_order = list(parts) if parts is not None else ["upstream_gene", "gene", "downstream_gene"]

        drd_map: dict[str, DiscreteRegionData] = {}
        for part in parts_order:
            if add_flanks and part in {"upstream_gene", "downstream_gene"}:
                contigs, auto_labels = _collect_synthetic_flanks(
                    annot,
                    flank_bp=int(flank_bp),
                    kind=part,
                    limit=limit,
                )
            elif part == "gene":
                contigs, auto_labels = _collect_gene_contigs(annot, limit=limit)
            else:
                contigs, auto_labels = collect_contigs_from_hcannot(
                    annot,
                    feature_type=part,
                    limit=limit,
                )
            if not contigs:
                continue
            drd_map[part] = compute_discrete_regions(
                reader,
                contigs,
                segments=segments,
                reverse_negative=reverse_negative,
                labels=auto_labels,
            )

        if not drd_map:
            return DiscreteRegionData()

        return combine_parts_drd(drd_map, segments=segments, parts_order=parts_order)

    if add_flanks and feature_type in {"upstream_gene", "downstream_gene"}:
        contigs, auto_labels = _collect_synthetic_flanks(
            annot,
            flank_bp=int(flank_bp),
            kind=cast(str, feature_type),
            limit=limit,
        )
    elif feature_type == "gene":
        contigs, auto_labels = _collect_gene_contigs(annot, limit=limit)
    else:
        contigs, auto_labels = collect_contigs_from_hcannot(
            annot,
            feature_type=feature_type,
            limit=limit,
        )
    use_labels = list(labels) if labels is not None else auto_labels
    return compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        reverse_negative=reverse_negative,
        labels=use_labels,
    )
