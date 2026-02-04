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
    mode: str = "discretise",
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
            try:
                if len(batch.data()) < 3:
                    _progress(idx)
                    continue
            except Exception:
                _progress(idx)
                continue
            try:
                xs, ys = batch.discretise(total_bins, agg_method)
            except Exception:
                _progress(idx)
                continue
            x = np.asarray(xs, dtype=np.float64)
            y = np.asarray(ys, dtype=np.float64)
            if reverse_negative and contig is not None and _is_negative_strand(contig):
                x = 1.0 - x[::-1]
                y = y[::-1]
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

    def _get_contig_start(contig) -> Optional[int]:
        v = getattr(contig, "start", None)
        if v is None:
            return None
        try:
            return int(v() if callable(v) else v)
        except (TypeError, ValueError):
            return None

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
):
    """HoloViews Curve for averaged metagene profile."""
    segments = segments or _DEFAULT_SEGMENTS
    bounds, names = segment_ticks(segments)
    drd = compute_discrete_regions(
        reader,
        contigs,
        segments=segments,
        agg_method=agg_method,
        mode=mode,
        x_mode=x_mode,
    )
    x, y = _line_from_points(drd, segments)
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
