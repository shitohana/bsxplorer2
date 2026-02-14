from __future__ import annotations

from typing import Callable, Optional, Sequence

import holoviews as hv
import numpy as np
import plotly.graph_objects as go

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.metagene import (
    Segment,
    combine_parts_drd,
    collect_contigs_from_hcannot,
    collect_parts_from_hcannot,
    compute_discrete_regions,
    segments_total_bins,
)


def _hv_init() -> None:
    hv.extension("plotly")


def _ensure_plotly(fig):
    if isinstance(fig, go.Figure):
        return fig
    return go.Figure(fig)


def _bin_points_windows(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    *,
    n_windows: int,
    agg: str,
    nan_policy: str,
) -> np.ndarray:
    if n_windows <= 0:
        raise ValueError("n_windows must be > 0")
    if nan_policy not in {"drop", "zero", "keep"}:
        raise ValueError("nan_policy must be 'drop', 'zero', or 'keep'")

    bins_values = [[] for _ in range(n_windows)]
    for x, y in zip(x_vals, y_vals):
        if not np.isfinite(x):
            continue
        if x < 0.0 or x > 1.0:
            continue
        if nan_policy == "zero" and not np.isfinite(y):
            y = 0.0
        if nan_policy == "drop" and not np.isfinite(y):
            continue
        idx = int(x * n_windows)
        if idx == n_windows:
            idx = n_windows - 1
        if 0 <= idx < n_windows:
            bins_values[idx].append(y)

    out = []
    for vals in bins_values:
        if not vals:
            out.append(np.nan)
            continue
        arr = np.asarray(vals, dtype=float)
        if agg == "mean":
            fn = np.mean if nan_policy in {"keep", "zero"} else np.nanmean
            out.append(float(fn(arr)))
        elif agg == "median":
            fn = np.median if nan_policy in {"keep", "zero"} else np.nanmedian
            out.append(float(fn(arr)))
        elif agg == "max":
            out.append(float(np.nanmax(arr)) if nan_policy == "drop" else float(np.max(arr)))
        elif agg == "min":
            out.append(float(np.nanmin(arr)) if nan_policy == "drop" else float(np.min(arr)))
        else:
            raise ValueError(f"Unsupported agg: {agg}")
    return np.asarray(out, dtype=float)


def _rank_compress(z_sorted: np.ndarray, rank_rows: int, *, fill: float | None = 0.0) -> np.ndarray:
    if z_sorted.ndim != 2:
        return z_sorted
    n_rows, n_bins = z_sorted.shape
    if n_rows == 0:
        return np.empty((0, n_bins), dtype=float)
    rows = max(int(rank_rows), 1)
    sums = np.zeros((rows, n_bins), dtype=float)
    cnts = np.zeros((rows, n_bins), dtype=np.int32)
    for i in range(n_rows):
        ridx = int(i * rows / n_rows)
        row = z_sorted[i]
        finite = np.isfinite(row)
        if np.any(finite):
            sums[ridx, finite] += row[finite]
            cnts[ridx, finite] += 1
    if fill is None:
        out = np.full((rows, n_bins), np.nan, dtype=float)
    else:
        out = np.full((rows, n_bins), float(fill), dtype=float)
    np.divide(sums, cnts, out=out, where=(cnts > 0))
    return out


def _segment_decor(
    fig,
    segments: list[Segment] | None,
    *,
    annotate_tss_tes: bool = False,
    scale: Callable[[float], float],
) -> None:
    if not segments:
        return
    boundaries = []
    centers = []
    labels = []
    cum = 0
    for seg in segments:
        start = cum
        end = cum + seg.n_bins
        mid = (start + end) / 2
        boundaries.append(end)
        centers.append(mid)
        labels.append(seg.name)
        cum = end
    shapes = []
    for b in boundaries[:-1]:
        xb = scale(b)
        shapes.append(
            dict(type="line", x0=xb, x1=xb, y0=0, y1=1, yref="paper", line=dict(dash="dash", width=1, color="gray"))
        )
    tickvals = [scale(c) for c in centers]
    fig.update_xaxes(tickmode="array", tickvals=tickvals, ticktext=labels)
    fig.update_layout(shapes=shapes)
    if annotate_tss_tes and len(boundaries) >= 2:
        first = scale(boundaries[0])
        last = scale(boundaries[-2])
        fig.add_annotation(x=first, y=1.02, xref="x", yref="paper", text="TSS", showarrow=False, font=dict(size=10))
        fig.add_annotation(x=last, y=1.02, xref="x", yref="paper", text="TES", showarrow=False, font=dict(size=10))


def _segment_decor_rel(fig, segments: list[Segment] | None, *, annotate_tss_tes: bool = False) -> None:
    if not segments:
        return
    total = float(segments_total_bins(segments))
    _segment_decor(fig, segments, annotate_tss_tes=annotate_tss_tes, scale=lambda v: v / total)


def _segment_decor_bin(fig, segments: list[Segment] | None, *, annotate_tss_tes: bool = False) -> None:
    _segment_decor(fig, segments, annotate_tss_tes=annotate_tss_tes, scale=float)


def _drd_from_annot(
    reader,
    annot,
    *,
    segments: list[Segment] | None,
    feature_type: str | None,
    reverse_negative: bool,
    labels: list[str] | None,
    limit: int | None,
    add_flanks: bool = False,
    flank_bp: int = 2000,
    combine_parts: bool = False,
    parts: Sequence[str] | None = None,
) -> DiscreteRegionData:
    if segments is None:
        if combine_parts:
            segments = [Segment("up", 100), Segment("body", 200), Segment("down", 100)]
        else:
            segments = [Segment("region", 100)]
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
