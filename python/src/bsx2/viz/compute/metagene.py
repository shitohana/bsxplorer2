from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from beartype.typing import Callable, List, Optional, Sequence, Tuple, cast

from bsx2 import Contig, RegionReader
from bsx2.guards import (
    require_nan_policy_compatible,
    require_per_segment_profile,
    require_savgol_filter,
    resolve_reader_accessors,
)
from bsx2.validation import validate_segments, validate_smoothing

from .data import DiscreteRegionData

if TYPE_CHECKING:
    from bsx2.types import Contig


# ============================================================
# Segment definition
# ============================================================

@dataclass(frozen=True)
class MetageneProfileSegment:
    """One named segment of a normalized profile."""

    name: str
    n_bins: int


_ANNOT_LAYOUT_SOURCES = {"feature", "gene", "flank5", "flank3"}


@dataclass(frozen=True)
class AnnotProfilePart:
    """
    Describe one ordered section of an annotation-driven metagene layout.

    Parameters
    ----------
    name
        Public segment name used in the composed profile.
    n_bins
        Number of bins allocated to this section in the final normalized profile.
    source
        Part source kind. Supported values are ``"feature"``, ``"gene"``,
        ``"flank5"``, and ``"flank3"``.
    feature_type
        Annotation feature type to collect when ``source="feature"``. If omitted,
        ``name`` is reused.
    flank_bp
        Flank length in base pairs for synthetic 5' or 3' flanks.
    required
        Whether genes missing this part should be skipped entirely.
    """

    name: str
    n_bins: int
    source: str = "feature"
    feature_type: Optional[str] = None
    flank_bp: Optional[int] = None
    required: bool = False

    def to_segment(self) -> MetageneProfileSegment:
        return MetageneProfileSegment(self.name, self.n_bins)


@dataclass(frozen=True)
class AnnotProfileLayout:
    """
    Ordered layout definition for annotation-driven metagene composition.

    Notes
    -----
    The layout is gene-centric: parts are collected per gene, resolved in the
    declared order, and then concatenated into one normalized profile.
    """

    parts: tuple[AnnotProfilePart, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "parts", tuple(self.parts))

    @property
    def segments(self) -> tuple[MetageneProfileSegment, ...]:
        return tuple(part.to_segment() for part in self.parts)


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


def _validate_annot_layout(layout: AnnotProfileLayout) -> tuple[AnnotProfilePart, ...]:
    if not isinstance(layout, AnnotProfileLayout):
        raise ValueError("layout must be an AnnotProfileLayout instance")

    parts = tuple(layout.parts)
    if not parts:
        raise ValueError("layout.parts must not be empty")

    names: set[str] = set()
    for idx, part in enumerate(parts):
        validate_segments([part.to_segment()])
        source = str(part.source).strip().lower()
        if source not in _ANNOT_LAYOUT_SOURCES:
            raise ValueError(
                f"layout.parts[{idx}].source must be one of: {sorted(_ANNOT_LAYOUT_SOURCES)}"
            )
        if part.name in names:
            raise ValueError(f"layout.parts names must be unique; duplicate: {part.name!r}")
        names.add(part.name)

        if source == "feature":
            feature_type = part.feature_type or part.name
            if not isinstance(feature_type, str) or not feature_type.strip():
                raise ValueError(
                    f"layout.parts[{idx}].feature_type must be a non-empty string "
                    "for source='feature'"
                )

        if source in {"flank5", "flank3"}:
            if part.flank_bp is None:
                raise ValueError(
                    f"layout.parts[{idx}].flank_bp must be set for source={source!r}"
                )
            if isinstance(part.flank_bp, bool) or not isinstance(part.flank_bp, int):
                raise ValueError(
                    f"layout.parts[{idx}].flank_bp must be a positive integer"
                )
            if part.flank_bp <= 0:
                raise ValueError(
                    f"layout.parts[{idx}].flank_bp must be a positive integer"
                )

    return parts


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


def _get_contig_end(contig: object) -> Optional[int]:
    v = getattr(contig, "end", None)
    if v is None:
        return None
    try:
        return int(v() if callable(v) else v)
    except (TypeError, ValueError):
        return None


def _get_contig_seqname(contig: object) -> Optional[str]:
    v = getattr(contig, "seqname", None)
    if v is None:
        return None
    try:
        value = v() if callable(v) else v
    except (TypeError, ValueError):
        return None
    if value is None:
        return None
    return str(value)


def _get_contig_strand(contig: object):
    v = getattr(contig, "strand", None)
    if v is None:
        return None
    try:
        return v() if callable(v) else v
    except (TypeError, ValueError):
        return None


def _fallback_contig_label(contig: object) -> str:
    seqname = _get_contig_seqname(contig) or "unknown"
    start = _get_contig_start(contig)
    end = _get_contig_end(contig)
    return f"{seqname}:{start}-{end}"


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

    start = _get_contig_start(contig)
    end = _get_contig_end(contig)
    strand = _get_contig_strand(contig)
    is_negative = _is_negative_strand(contig)
    seqname = _get_contig_seqname(contig)

    if start is None or end is None or strand is None or seqname is None:
        return None

    if kind in {"upstream_gene", "flank5"}:
        if is_negative:
            flank_start = end
            flank_end = end + flank_bp
        else:
            flank_start = max(0, start - flank_bp)
            flank_end = start
    elif kind in {"downstream_gene", "flank3"}:
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

    return Contig(seqname, flank_start, flank_end, strand)


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
            label = _fallback_contig_label(contig)

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


def _remember_gene_id(gene_id: str, gene_order: list[str], gene_seen: set[str]) -> None:
    if gene_id in gene_seen:
        return
    gene_order.append(gene_id)
    gene_seen.add(gene_id)


def _merge_contigs_span(contigs: Sequence[Contig]) -> Optional[Contig]:
    if not contigs:
        return None
    if len(contigs) == 1:
        return contigs[0]

    seqname = _get_contig_seqname(contigs[0])
    strand = _get_contig_strand(contigs[0])
    starts: list[int] = []
    ends: list[int] = []

    if seqname is None or strand is None:
        return None

    for contig in contigs:
        current_seqname = _get_contig_seqname(contig)
        current_strand = _get_contig_strand(contig)
        start = _get_contig_start(contig)
        end = _get_contig_end(contig)
        if (
            current_seqname != seqname
            or current_strand != strand
            or start is None
            or end is None
            or end <= start
        ):
            return None
        starts.append(start)
        ends.append(end)

    return Contig(seqname, min(starts), max(ends), strand)


def _resolve_layout_part_contig(
    gene_id: str,
    part: AnnotProfilePart,
    *,
    gene_contigs: dict[str, Contig],
    feature_contigs: dict[str, dict[str, list[Contig]]],
) -> Optional[Contig]:
    source = str(part.source).strip().lower()

    if source == "gene":
        return gene_contigs.get(gene_id)

    if source in {"flank5", "flank3"}:
        gene_contig = gene_contigs.get(gene_id)
        if gene_contig is None or part.flank_bp is None:
            return None
        return _make_flank_contig(gene_contig, flank_bp=part.flank_bp, kind=source)

    contigs = feature_contigs.get(gene_id, {}).get(part.name, [])
    return _merge_contigs_span(contigs)


def collect_layout_parts_from_hcannot(
    annot: object,
    *,
    layout: AnnotProfileLayout,
    limit: Optional[int] = None,
) -> dict[str, tuple[list[object], list[str]]]:
    """
    Collect per-gene annotation parts for a declared metagene layout.

    Parameters
    ----------
    annot
        `HcAnnotStore`-like object exposing ``iter()`` or regular iteration.
    layout
        Ordered annotation layout describing which parts should be resolved for
        each gene.
    limit
        Optional maximum number of genes to keep in the observed annotation order.

    Returns
    -------
    dict[str, tuple[list[object], list[str]]]
        Mapping from part name to ``(contigs, gene_ids)`` pairs. Optional parts
        may be empty; required missing parts cause the whole gene to be skipped.
    """

    parts = _validate_annot_layout(layout)
    feature_parts = [part for part in parts if str(part.source).strip().lower() == "feature"]
    feature_types = {part.name: str(part.feature_type or part.name) for part in feature_parts}

    gene_contigs: dict[str, Contig] = {}
    feature_contigs: dict[str, dict[str, list[Contig]]] = defaultdict(lambda: defaultdict(list))
    gene_order: list[str] = []
    gene_seen: set[str] = set()

    iterator_factory = getattr(annot, "iter", None)
    iterator = iterator_factory() if callable(iterator_factory) else iter(annot)

    for item in iterator:
        entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item
        contig = _entry_contig(entry)
        if contig is None:
            continue

        start = _get_contig_start(contig)
        if start is not None and start < 0:
            continue

        ft = _entry_feature_type(entry)
        if ft == "gene":
            gene_id = _entry_id(entry) or _fallback_contig_label(contig)
            gene_contigs.setdefault(gene_id, contig)
            _remember_gene_id(gene_id, gene_order, gene_seen)

        for part in feature_parts:
            if ft != feature_types[part.name]:
                continue
            if ft == "gene":
                gene_id = _entry_id(entry) or _fallback_contig_label(contig)
            else:
                parents = _entry_parents(entry)
                gene_id = parents[0] if parents else None
            if not gene_id:
                continue
            feature_contigs[gene_id][part.name].append(contig)
            _remember_gene_id(gene_id, gene_order, gene_seen)

    anchored_gene_ids = [gene_id for gene_id in gene_order if gene_id in gene_contigs]
    selected_gene_ids = (
        anchored_gene_ids if limit is None else anchored_gene_ids[: int(limit)]
    )
    part_items: dict[str, list[tuple[object, str]]] = {part.name: [] for part in parts}

    for gene_id in selected_gene_ids:
        resolved_parts: dict[str, Contig] = {}
        missing_required = False
        for part in parts:
            contig = _resolve_layout_part_contig(
                gene_id,
                part,
                gene_contigs=gene_contigs,
                feature_contigs=feature_contigs,
            )
            if contig is None:
                if part.required:
                    missing_required = True
                    break
                continue
            resolved_parts[part.name] = contig
        if missing_required or not resolved_parts:
            continue
        for part_name, contig in resolved_parts.items():
            part_items[part_name].append((contig, gene_id))

    result: dict[str, tuple[list[object], list[str]]] = {}
    for part_name, items in part_items.items():
        if not items:
            result[part_name] = ([], [])
            continue
        contigs, labels = zip(*items)
        result[part_name] = (list(contigs), list(labels))

    return result


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


def compose_layout_drd(
    drd_map: dict[str, DiscreteRegionData],
    *,
    segments: Sequence[MetageneProfileSegment],
    parts_order: Sequence[str],
) -> DiscreteRegionData:
    """
    Compose one metagene profile from multiple per-part `DiscreteRegionData` objects.

    Parameters
    ----------
    drd_map
        Per-part data indexed by part name.
    segments
        Final segment definitions for the composed normalized profile.
    parts_order
        Ordered list of part names matching ``segments``.
    """

    validate_segments(segments)
    if len(parts_order) != len(segments):
        raise ValueError("segments and parts_order must contain the same number of items")

    bounds = segment_boundaries(segments)
    starts = [0.0, *bounds[:-1]]
    ends = list(bounds)

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


def _compute_from_annot_layout(
    reader: RegionReader,
    annot: object,
    *,
    layout: AnnotProfileLayout,
    segments: Sequence[MetageneProfileSegment],
    reverse_negative: bool,
    limit: Optional[int],
) -> DiscreteRegionData:
    parts = _validate_annot_layout(layout)
    if len(segments) != len(parts):
        raise ValueError("segments and layout.parts must contain the same number of items")

    part_map = collect_layout_parts_from_hcannot(
        annot,
        layout=layout,
        limit=limit,
    )

    drd_map: dict[str, DiscreteRegionData] = {}
    for part in parts:
        contigs, auto_labels = part_map.get(part.name, ([], []))
        if not contigs:
            continue
        drd_map[part.name] = compute_discrete_regions(
            reader,
            contigs,
            segments=[part.to_segment()],
            reverse_negative=reverse_negative,
            labels=auto_labels,
        )

    if not drd_map:
        return DiscreteRegionData()

    return compose_layout_drd(
        drd_map,
        segments=segments,
        parts_order=[part.name for part in parts],
    )


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
    layout: AnnotProfileLayout | None = None,
) -> DiscreteRegionData:
    """
    Compute metagene data from an annotation store.

    Parameters
    ----------
    reader
        `RegionReader`-like object used to query methylation data.
    annot
        `HcAnnotStore`-like annotation container.
    segments
        Final metagene segments. When ``layout`` is provided and ``segments`` is
        omitted, ``layout.segments`` is used automatically.
    feature_type
        Single annotation feature type for the legacy single-part mode.
    reverse_negative
        Whether negative-strand regions should be reversed into 5'->3' profile
        orientation before insertion.
    labels
        Optional explicit labels for the selected regions.
    limit
        Optional maximum number of regions or genes to aggregate.
    add_flanks
        Whether to synthesize upstream or downstream flanks in single-part mode.
    flank_bp
        Synthetic flank length in base pairs for single-part mode.
    layout
        Ordered multi-part annotation layout. This is the preferred API for
        promoter/gene/terminator and arbitrary annotation-driven profiles.

    Returns
    -------
    DiscreteRegionData
        Compute-only metagene representation suitable for any downstream renderer.
    """

    if layout is not None:
        if segments is None:
            segments = layout.segments
        validate_segments(segments)
        return _compute_from_annot_layout(
            reader,
            annot,
            layout=layout,
            segments=segments,
            reverse_negative=reverse_negative,
            limit=limit,
        )

    if segments is None:
        segments = _DEFAULT_SEGMENTS

    validate_segments(
        segments,
        flank_bp=flank_bp,
    )

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
