from __future__ import annotations

from dataclasses import dataclass, field
import math
import warnings
from beartype.typing import Iterable, Mapping, Optional

import holoviews as hv
import numpy as np

from bsx2 import Context
from bsx2.validation import (
    ChrEmptyPolicy,
    ChrLineStat,
    validate_bin_size_bp,
    validate_chr_empty_policy,
    validate_chr_line_stat,
    validate_chr_lengths,
    validate_min_coverage,
    validate_positive_int,
    validate_smoothing,
    validate_single_context,
)

from ._common import _ensure_plotly, _hv_init
from .metagene import _apply_savgol_smoothing, _clip_profile


_U32_MAX = (1 << 32) - 1


@dataclass(frozen=True)
class ChrLineTrack:
    seqname: str
    context: Context
    stat: ChrLineStat
    bin_size_bp: int
    chr_length_bp: int
    offset_bp: int
    window_start_bp: np.ndarray
    window_end_bp: np.ndarray
    x_bp: np.ndarray
    x_global_bp: np.ndarray
    value: np.ndarray
    sum_m: np.ndarray
    sum_total: np.ndarray
    site_count: np.ndarray
    sum_density: np.ndarray
    density_site_count: np.ndarray


@dataclass(frozen=True)
class ChrLineData:
    tracks: list[ChrLineTrack]
    context: Context
    stat: ChrLineStat
    bin_size_bp: int
    empty_policy: ChrEmptyPolicy

    def is_empty(self) -> bool:
        return len(self.tracks) == 0

    def chromosomes(self) -> list[str]:
        return [track.seqname for track in self.tracks]

    def genome_length_bp(self) -> int:
        if not self.tracks:
            return 0
        last = self.tracks[-1]
        return int(last.offset_bp + last.chr_length_bp)

    def chr_boundaries_bp(self) -> list[float]:
        if len(self.tracks) <= 1:
            return []
        return [
            float(track.offset_bp + track.chr_length_bp)
            for track in self.tracks[:-1]
        ]

    def chr_centers_bp(self) -> list[float]:
        return [
            float(track.offset_bp + 0.5 * track.chr_length_bp)
            for track in self.tracks
        ]

    def global_curve(self) -> tuple[np.ndarray, np.ndarray]:
        if not self.tracks:
            return (
                np.array([], dtype=np.float64),
                np.array([], dtype=np.float64),
            )

        x_parts: list[np.ndarray] = []
        y_parts: list[np.ndarray] = []

        for index, track in enumerate(self.tracks):
            x_parts.append(track.x_global_bp.astype(np.float64, copy=False))
            y_parts.append(track.value.astype(np.float64, copy=False))

            if index + 1 != len(self.tracks):
                boundary = float(track.offset_bp + track.chr_length_bp)
                x_parts.append(np.array([boundary], dtype=np.float64))
                y_parts.append(np.array([np.nan], dtype=np.float64))

        return np.concatenate(x_parts), np.concatenate(y_parts)


def _series_to_numpy(
    series: object,
    *,
    dtype: np.dtype | type | None = None,
) -> np.ndarray:
    if hasattr(series, "to_numpy"):
        values = np.asarray(series.to_numpy())
    elif hasattr(series, "to_numpy_array"):
        values = np.asarray(series.to_numpy_array())
    elif hasattr(series, "to_array"):
        values = np.asarray(series.to_array())
    elif hasattr(series, "to_list"):
        values = np.asarray(series.to_list(), dtype=object)
    else:
        values = np.asarray(series, dtype=object)

    if dtype is not None:
        values = values.astype(dtype, copy=False)
    return values


def _context_mask(
    context_values: np.ndarray,
    context: Context,
) -> np.ndarray:
    values = np.asarray(context_values, dtype=object)
    mask_cg = (
        np.equal(values, True)
        | np.equal(values, "CG")
        | np.equal(values, "cg")
    )
    mask_chg = (
        np.equal(values, False)
        | np.equal(values, "CHG")
        | np.equal(values, "chg")
    )
    if context is Context.CG:
        return mask_cg
    if context is Context.CHG:
        return mask_chg
    return ~(mask_cg | mask_chg)


def _iter_reader_batches(
    reader: object,
    *,
    chr_lengths: Mapping[str, int] | None,
) -> Iterable[object]:
    try:
        iterator = iter(reader)
    except TypeError:
        iterator = None

    if iterator is not None:
        for batch in iterator:
            if batch is not None:
                yield batch
        return

    query_fn = getattr(reader, "query", None)
    chr_order_fn = getattr(reader, "chr_order", None)
    if not callable(query_fn) or not callable(chr_order_fn):
        raise AttributeError(
            "reader must be iterable (BsxFileReader-like) or support "
            "chr_order()/query(contig) (RegionReader-like)"
        )

    from bsx2 import Contig, Strand

    strand_null = getattr(Strand, "Null", None)
    if strand_null is None:
        raise AttributeError("Strand.Null is required for RegionReader queries")

    reset_fn = getattr(reader, "reset", None)
    if callable(reset_fn):
        try:
            reset_fn()
        except Exception:
            pass

    for seqname in chr_order_fn():
        end = (
            int(chr_lengths[seqname])
            if chr_lengths is not None and seqname in chr_lengths
            else _U32_MAX
        )
        contig = Contig(seqname, 0, end, strand_null)
        batch = query_fn(contig)
        if batch is not None:
            yield batch


def _safe_seqname(batch: object) -> str | None:
    seqname_attr = getattr(batch, "seqname", None)
    if callable(seqname_attr):
        try:
            seqname = seqname_attr()
            if seqname is None:
                return None
            return str(seqname)
        except Exception:
            return None
    if isinstance(seqname_attr, str):
        return seqname_attr
    return None


def _aggregate_batch(
    batch: object,
    *,
    context: Context,
    bin_size_bp: int,
    min_coverage: int,
) -> tuple[str | None, dict[int, tuple[float, float, int, float, int]], int]:
    seqname = _safe_seqname(batch)
    if seqname is None:
        return None, {}, 0

    positions = _series_to_numpy(getattr(batch, "position")(), dtype=np.int64)
    if positions.size == 0:
        return seqname, {}, 0

    count_m = _series_to_numpy(getattr(batch, "count_m")(), dtype=np.float64)
    count_total = _series_to_numpy(getattr(batch, "count_total")(), dtype=np.float64)
    context_values = _series_to_numpy(getattr(batch, "context")(), dtype=object)

    mask = _context_mask(context_values, context)
    if min_coverage > 0:
        mask &= np.isfinite(count_total) & (count_total >= min_coverage)
    else:
        mask &= np.isfinite(count_total)
    mask &= np.isfinite(count_m)

    max_pos = int(np.max(positions))
    if not np.any(mask):
        return seqname, {}, max_pos

    pos_f = positions[mask]
    m_f = count_m[mask]
    t_f = count_total[mask]

    window_ids = ((pos_f - 1) // bin_size_bp).astype(np.int64)
    order = np.argsort(window_ids, kind="mergesort")
    window_ids = window_ids[order]
    m_f = m_f[order]
    t_f = t_f[order]

    cuts = np.flatnonzero(np.diff(window_ids)) + 1
    starts = np.r_[0, cuts]
    uniq_windows = window_ids[starts]
    sum_m = np.add.reduceat(m_f, starts)
    sum_t = np.add.reduceat(t_f, starts)
    site_count = np.diff(np.r_[starts, window_ids.size])

    payload = {
        int(win): (float(sm), float(st), int(sc), 0.0, 0)
        for win, sm, st, sc in zip(uniq_windows, sum_m, sum_t, site_count)
    }

    valid_density = t_f > 0
    if np.any(valid_density):
        win_d = window_ids[valid_density]
        dens = m_f[valid_density] / t_f[valid_density]

        order_d = np.argsort(win_d, kind="mergesort")
        win_d = win_d[order_d]
        dens = dens[order_d]

        cuts_d = np.flatnonzero(np.diff(win_d)) + 1
        starts_d = np.r_[0, cuts_d]
        uniq_d = win_d[starts_d]
        sum_d = np.add.reduceat(dens, starts_d)
        count_d = np.diff(np.r_[starts_d, win_d.size])

        for win, sd, cd in zip(uniq_d, sum_d, count_d):
            current = payload[int(win)]
            payload[int(win)] = (
                current[0],
                current[1],
                current[2],
                float(sd),
                int(cd),
            )

    return seqname, payload, max_pos


def _finalize_tracks(
    *,
    context: Context,
    stat: ChrLineStat,
    bin_size_bp: int,
    empty_policy: ChrEmptyPolicy,
    smooth: dict | int | None,
    chr_order: list[str],
    chr_lengths_bp: dict[str, int],
    aggregate: dict[str, dict[int, tuple[float, float, int, float, int]]],
) -> list[ChrLineTrack]:
    tracks: list[ChrLineTrack] = []
    offset = 0

    for seqname in chr_order:
        chr_length = int(chr_lengths_bp[seqname])
        if chr_length <= 0:
            continue

        n_windows = int(math.ceil(chr_length / float(bin_size_bp)))
        if n_windows <= 0:
            continue

        starts = (np.arange(n_windows, dtype=np.int64) * bin_size_bp) + 1
        ends = np.minimum(starts + bin_size_bp - 1, chr_length)
        centers = 0.5 * (starts + ends)

        sum_m = np.zeros(n_windows, dtype=np.float64)
        sum_t = np.zeros(n_windows, dtype=np.float64)
        sites = np.zeros(n_windows, dtype=np.int64)
        sum_density = np.zeros(n_windows, dtype=np.float64)
        density_sites = np.zeros(n_windows, dtype=np.int64)

        chr_payload = aggregate.get(seqname, {})
        if chr_payload:
            win_ids = np.array(list(chr_payload.keys()), dtype=np.int64)
            vals = np.array(list(chr_payload.values()), dtype=np.float64)
            sum_m[win_ids] = vals[:, 0]
            sum_t[win_ids] = vals[:, 1]
            sites[win_ids] = vals[:, 2].astype(np.int64)
            sum_density[win_ids] = vals[:, 3]
            density_sites[win_ids] = vals[:, 4].astype(np.int64)

        values = np.full(n_windows, np.nan, dtype=np.float64)
        if stat is ChrLineStat.MEAN:
            non_empty = density_sites > 0
            values[non_empty] = (
                sum_density[non_empty] / density_sites[non_empty]
            )
        else:
            non_empty = sum_t > 0
            values[non_empty] = sum_m[non_empty] / sum_t[non_empty]

        if empty_policy is ChrEmptyPolicy.ZERO:
            values[~non_empty] = 0.0
            keep_mask = np.ones(n_windows, dtype=bool)
        elif empty_policy is ChrEmptyPolicy.DROP:
            keep_mask = non_empty
        else:
            keep_mask = np.ones(n_windows, dtype=bool)

        if smooth is not None and values.size > 0:
            try:
                smooth_cfg = validate_smoothing(
                    smooth,
                    total_bins=int(values.size),
                    series_len=int(values.size),
                    label=f"{seqname} profile",
                )
                if smooth_cfg is not None:
                    if smooth_cfg.get("per_segment"):
                        smooth_cfg = dict(smooth_cfg)
                        smooth_cfg["per_segment"] = False
                    values = _apply_savgol_smoothing(values, smooth_cfg)
                    values = _clip_profile(values)
            except ValueError as exc:
                warnings.warn(
                    f"smoothing skipped for {seqname}: {exc}",
                    RuntimeWarning,
                    stacklevel=2,
                )

        if not np.any(keep_mask):
            offset += chr_length
            continue

        starts = starts[keep_mask]
        ends = ends[keep_mask]
        centers = centers[keep_mask]
        values = values[keep_mask]
        sum_m = sum_m[keep_mask]
        sum_t = sum_t[keep_mask]
        sites = sites[keep_mask]
        sum_density = sum_density[keep_mask]
        density_sites = density_sites[keep_mask]
        x_global = centers + offset

        tracks.append(
            ChrLineTrack(
                seqname=seqname,
                context=context,
                stat=stat,
                bin_size_bp=bin_size_bp,
                chr_length_bp=chr_length,
                offset_bp=offset,
                window_start_bp=starts.astype(np.int64, copy=False),
                window_end_bp=ends.astype(np.int64, copy=False),
                x_bp=centers.astype(np.float64, copy=False),
                x_global_bp=x_global.astype(np.float64, copy=False),
                value=values.astype(np.float64, copy=False),
                sum_m=sum_m.astype(np.float64, copy=False),
                sum_total=sum_t.astype(np.float64, copy=False),
                site_count=sites.astype(np.int64, copy=False),
                sum_density=sum_density.astype(np.float64, copy=False),
                density_site_count=density_sites.astype(np.int64, copy=False),
            )
        )
        offset += chr_length

    return tracks


def compute_chr_line_data(
    reader: object,
    *,
    context: Context,
    bin_size_bp: int = 50_000,
    chr_lengths: dict[str, int] | None = None,
    min_coverage: int | float = 0,
    stat: ChrLineStat | str = ChrLineStat.WEIGHTED_MEAN,
    smooth: dict | int | None = None,
    empty_policy: ChrEmptyPolicy | str = ChrEmptyPolicy.NAN,
) -> ChrLineData:
    """
    Compute chromosome-wide methylation profile in absolute genomic coordinates.

    Parameters
    ----------
    reader
        `BsxFileReader`-like iterable over batches or `RegionReader`-like object
        exposing `chr_order()` and `query(contig)`.
    context
        Required context filter: `Context.CG`, `Context.CHG`, or `Context.CHH`.
    bin_size_bp
        Fixed genomic window size in base pairs.
    chr_lengths
        Optional real chromosome lengths `{seqname: length_bp}`.
        If omitted, lengths are inferred from maximal observed positions.
    min_coverage
        Minimum per-site coverage threshold before a site is included.
    stat
        Aggregation statistic: `"weighted_mean"` (default) or `"mean"`.
    smooth
        Optional Savitzky-Golay smoothing config (`dict`), ratio (`int`),
        or `None`. Applied independently for each chromosome profile.
    empty_policy
        Policy for windows with no covered/context-matching sites:
        `"nan"` (default), `"zero"`, `"drop"`.
    """
    context = validate_single_context(context)
    bin_size_bp = validate_bin_size_bp(bin_size_bp)
    min_cov = int(validate_min_coverage(min_coverage, integer=True))
    stat = validate_chr_line_stat(stat)
    empty_policy = validate_chr_empty_policy(empty_policy)
    chr_lengths_validated = validate_chr_lengths(chr_lengths)

    aggregate: dict[str, dict[int, list[float]]] = {}
    observed_max_pos: dict[str, int] = {}
    chr_order: list[str] = []
    seen_chr: set[str] = set()

    for batch in _iter_reader_batches(reader, chr_lengths=chr_lengths_validated):
        seqname, payload, max_pos = _aggregate_batch(
            batch,
            context=context,
            bin_size_bp=bin_size_bp,
            min_coverage=min_cov,
        )
        if seqname is None:
            continue

        if seqname not in seen_chr:
            seen_chr.add(seqname)
            chr_order.append(seqname)

        prev_max = observed_max_pos.get(seqname, 0)
        if max_pos > prev_max:
            observed_max_pos[seqname] = max_pos

        if not payload:
            continue

        target = aggregate.setdefault(seqname, {})
        for window_id, (
            sum_m,
            sum_t,
            site_count,
            sum_density,
            density_site_count,
        ) in payload.items():
            entry = target.get(window_id)
            if entry is None:
                target[window_id] = [
                    sum_m,
                    sum_t,
                    float(site_count),
                    sum_density,
                    float(density_site_count),
                ]
            else:
                entry[0] += sum_m
                entry[1] += sum_t
                entry[2] += float(site_count)
                entry[3] += sum_density
                entry[4] += float(density_site_count)

    if chr_lengths_validated is not None:
        for seqname in chr_lengths_validated:
            if seqname not in seen_chr:
                chr_order.append(seqname)
        chr_lengths_bp = dict(chr_lengths_validated)
    else:
        chr_lengths_bp = dict(observed_max_pos)

    for seqname, max_pos in observed_max_pos.items():
        if seqname not in chr_lengths_bp:
            chr_lengths_bp[seqname] = max_pos
        elif max_pos > chr_lengths_bp[seqname]:
            raise ValueError(
                f"chr_lengths[{seqname}]={chr_lengths_bp[seqname]} is smaller than "
                f"observed max position {max_pos}"
            )

    aggregate_cast: dict[str, dict[int, tuple[float, float, int, float, int]]] = {}
    for seqname, windows in aggregate.items():
        aggregate_cast[seqname] = {
            w: (vals[0], vals[1], int(vals[2]), vals[3], int(vals[4]))
            for w, vals in windows.items()
        }

    tracks = _finalize_tracks(
        context=context,
        stat=stat,
        bin_size_bp=bin_size_bp,
        empty_policy=empty_policy,
        smooth=smooth,
        chr_order=chr_order,
        chr_lengths_bp=chr_lengths_bp,
        aggregate=aggregate_cast,
    )

    return ChrLineData(
        tracks=tracks,
        context=context,
        stat=stat,
        bin_size_bp=bin_size_bp,
        empty_policy=empty_policy,
    )


@dataclass
class ChrLinePlotComposer:
    title: Optional[str] = None
    width: int | None = None
    height: int | None = None
    show_chr_boundaries: bool = True
    show_chr_labels: bool = True
    x_label: str = "Genomic coordinate (bp)"
    y_label: Optional[str] = None
    datasets: list[ChrLineData] = field(default_factory=list)
    labels: list[str] = field(default_factory=list)

    def set_width(self, width: int | None) -> "ChrLinePlotComposer":
        self.width = None if width is None else validate_positive_int(
            width,
            name="width",
        )
        return self

    def set_height(self, height: int | None) -> "ChrLinePlotComposer":
        self.height = None if height is None else validate_positive_int(
            height,
            name="height",
        )
        return self

    def add_data(
        self,
        data: ChrLineData,
        *,
        name: str = "sample",
    ) -> "ChrLinePlotComposer":
        self.datasets.append(data)
        self.labels.append(name)
        return self

    def finish(self):
        _hv_init()
        if not self.datasets:
            return hv.Curve([]).opts(
                title=self.title or "Chromosome-wide methylation profile",
                xlabel=self.x_label,
                ylabel=self.y_label or "Methylation",
                show_legend=False,
                **({} if self.width is None else {"width": int(self.width)}),
                **({} if self.height is None else {"height": int(self.height)}),
            )

        curves = []
        for data, label in zip(self.datasets, self.labels):
            x_vals, y_vals = data.global_curve()
            curve = hv.Curve(
                (x_vals, y_vals),
                kdims=["genomic bp"],
                vdims=["methylation"],
            ).relabel(label)
            curves.append(curve)

        plot = curves[0]
        for curve in curves[1:]:
            plot *= curve

        first = self.datasets[0]
        if self.show_chr_boundaries:
            for boundary in first.chr_boundaries_bp():
                plot *= hv.VLine(float(boundary)).opts(
                    line_dash="dash",
                    line_color="gray",
                    line_width=1,
                )

        opts_kwargs = dict(
            title=self.title
            or f"Chromosome-wide methylation profile ({first.context.name})",
            xlabel=self.x_label,
            ylabel=self.y_label
            or (
                "Mean methylation"
                if first.stat is ChrLineStat.MEAN
                else "Weighted methylation"
            ),
            show_legend=len(curves) > 1,
        )
        if self.width is not None:
            opts_kwargs["width"] = int(self.width)
        if self.height is not None:
            opts_kwargs["height"] = int(self.height)

        if self.show_chr_labels and first.tracks:
            opts_kwargs["xticks"] = list(
                zip(first.chr_centers_bp(), first.chromosomes())
            )

        return plot.opts(**opts_kwargs)

    def to_html(
        self,
        *,
        full_html: bool = False,
        include_js: str = "cdn",
    ) -> str:
        fig = _ensure_plotly(hv.render(self.finish(), backend="plotly"))
        fig.update_layout(margin=dict(l=70, r=30, t=60, b=70))
        return fig.to_html(full_html=full_html, include_plotlyjs=include_js)
