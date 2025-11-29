from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from bsx2 import io as _io
from bsx2.plots.data import DiscreteRegionData, LinePlotData


@dataclass(frozen=True)
class Segment:
    """Сегмент метагена с именем и числом бинов."""
    name: str
    n_bins: int


def segments_total_bins(segments: Sequence[Segment]) -> int:
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
    segments: Sequence[Segment] | None = None,
    agg_method=None,
    reverse_negative: bool = True,
    labels: Optional[Sequence[str]] = None,
) -> DiscreteRegionData:
    """Строит дискретные профили по списку регионов (Contig) с использованием BsxBatch.discretise()."""
    if segments is None:
        segments = [Segment("region", 100)]
    total_bins = segments_total_bins(segments)

    if agg_method is None:
        from importlib import import_module
        agg_method = getattr(import_module("bsx2._bsx2"), "AggMethod").Mean

    data = DiscreteRegionData()
    batches_iter = getattr(reader, "iter_contigs", None)
    idx = 0

    if callable(batches_iter):
        for batch in reader.iter_contigs(list(contigs)):
            try:
                xs, ys = batch.discretise(total_bins, agg_method)
            except (AttributeError, TypeError, ValueError):
                idx += 1
                continue
            x = np.asarray(xs, dtype=np.float64)
            y = np.asarray(ys, dtype=np.float64)
            if reverse_negative and _is_negative_strand(contigs[idx]):
                x = 1.0 - x[::-1]
                y = y[::-1]
            label = labels[idx] if labels and idx < len(labels) else None
            data.insert(x, y, label)
            idx += 1
    else:
        for idx, contig in enumerate(contigs):
            try:
                batch = reader.query(contig)
                if batch is None:
                    continue
                xs, ys = batch.discretise(total_bins, agg_method)
            except (AttributeError, TypeError, ValueError):
                continue
            x = np.asarray(xs, dtype=np.float64)
            y = np.asarray(ys, dtype=np.float64)
            if reverse_negative and _is_negative_strand(contig):
                x = 1.0 - x[::-1]
                y = y[::-1]
            label = labels[idx] if labels and idx < len(labels) else None
            data.insert(x, y, label)

    return data


def collect_contigs_from_hcannot(
    annot,
    *,
    feature_type: Optional[str] = None,
    limit: Optional[int] = None,
    label_getter: Optional[Callable[[object, int], str]] = None,
) -> Tuple[List[object], List[str]]:
    """Извлекает список Contig и метки из HcAnnotStore (совместимо с dev-вариантом API)."""
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
    segments: Sequence[Segment] | None = None,
    agg_method=None,
    feature_type: Optional[str] = None,
    reverse_negative: bool = True,
    labels: Optional[Sequence[str]] = None,
) -> DiscreteRegionData:
    """Строит DiscreteRegionData из HcAnnotStore (фильтр по feature_type при необходимости)."""
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


def line_plot(reader: _io.RegionReader, *, contigs: Sequence, segments: Sequence[Segment] | None = None, agg_method=None):
    """Линейный метагенный профиль (HoloViews Curve)."""
    segments = segments or [Segment("region", 100)]
    bounds, names = segment_ticks(segments)
    drd = compute_discrete_regions(reader, contigs, segments=segments, agg_method=agg_method)
    lp = LinePlotData.from_discrete(drd)
    return LinePlotData(x=lp.x, y=lp.y, x_ticks=bounds, x_labels=names).to_curve()


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
    """Теплокарта (regions × bins)."""
    try:
        import holoviews as hv  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("holoviews is required for heatmap; install with 'pip install holoviews'") from e
    segments = segments or [Segment("region", 100)]
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
    """Коробчатая диаграмма распределений по бинам."""
    try:
        import holoviews as hv  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("holoviews is required for box_plot; install with 'pip install holoviews'") from e
    segments = segments or [Segment("region", 100)]
    drd = compute_discrete_regions(reader, contigs, segments=segments, agg_method=agg_method)
    mat, _ = drd.stack_matrix()
    if mat.size == 0:
        return hv.BoxWhisker([])
    n_regions, n_bins = mat.shape
    df = pd.DataFrame({"bin": np.tile(np.arange(n_bins), n_regions), "density": mat.reshape(-1)})
    return hv.BoxWhisker(df, kdims=["bin"], vdims=["density"])


def violin_plot(reader: _io.RegionReader, *, contigs: Sequence, segments: Sequence[Segment] | None = None, agg_method=None):
    """Виолин‑плот распределений по бинам."""
    try:
        import holoviews as hv  # type: ignore
    except ModuleNotFoundError as e:
        raise ImportError("holoviews is required for violin_plot; install with 'pip install holoviews'") from e
    segments = segments or [Segment("region", 100)]
    drd = compute_discrete_regions(reader, contigs, segments=segments, agg_method=agg_method)
    mat, _ = drd.stack_matrix()
    if mat.size == 0:
        return hv.Violin([])
    n_regions, n_bins = mat.shape
    df = pd.DataFrame({"bin": np.tile(np.arange(n_bins), n_regions), "density": mat.reshape(-1)})
    return hv.Violin(df, kdims=["bin"], vdims=["density"])
