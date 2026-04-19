from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np
from beartype.typing import Iterable, Mapping

from bsx2 import Context
from bsx2.validation import (
    ChrEmptyPolicy,
    ChrLineStat,
    validate_bin_size_bp,
    validate_chr_empty_policy,
    validate_chr_lengths,
    validate_chr_line_stat,
    validate_min_coverage,
    validate_single_context,
    validate_smoothing,
)

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
        return [float(track.offset_bp + track.chr_length_bp) for track in self.tracks[:-1]]

    def chr_centers_bp(self) -> list[float]:
        return [float(track.offset_bp + 0.5 * track.chr_length_bp) for track in self.tracks]

    def global_curve(self) -> tuple[np.ndarray, np.ndarray]:
        if not self.tracks:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

        n_boundaries = max(len(self.tracks) - 1, 0)
        total_points = sum(track.x_global_bp.size for track in self.tracks) + n_boundaries
        x = np.empty(total_points, dtype=np.float64)
        y = np.empty(total_points, dtype=np.float64)
        offset = 0

        for index, track in enumerate(self.tracks):
            track_size = int(track.x_global_bp.size)
            x[offset: offset + track_size] = track.x_global_bp.astype(np.float64, copy=False)
            y[offset: offset + track_size] = track.value.astype(np.float64, copy=False)
            offset += track_size
            if index + 1 != len(self.tracks):
                x[offset] = float(track.offset_bp + track.chr_length_bp)
                y[offset] = np.nan
                offset += 1

        return x, y


@dataclass
class _ChrWindowPayload:
    window_ids: list[np.ndarray]
    sum_m: list[np.ndarray]
    sum_total: list[np.ndarray]
    site_count: list[np.ndarray]
    sum_density: list[np.ndarray]
    density_site_count: list[np.ndarray]


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
    mask_cg = np.equal(values, True) | np.equal(values, "CG") | np.equal(values, "cg")
    mask_chg = np.equal(values, False) | np.equal(values, "CHG") | np.equal(values, "chg")
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
        end = int(chr_lengths[seqname]) if chr_lengths is not None and seqname in chr_lengths else _U32_MAX
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
) -> tuple[
    str | None,
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    int,
]:
    seqname = _safe_seqname(batch)
    if seqname is None:
        empty = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=np.float64)
        return None, (empty, empty_f, empty_f, empty, empty_f, empty), 0

    positions = _series_to_numpy(getattr(batch, "position")(), dtype=np.int64)
    if positions.size == 0:
        empty = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=np.float64)
        return seqname, (empty, empty_f, empty_f, empty, empty_f, empty), 0

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
        empty = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=np.float64)
        return seqname, (empty, empty_f, empty_f, empty, empty_f, empty), max_pos

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

    sum_density = np.zeros_like(sum_m, dtype=np.float64)
    density_site_count = np.zeros_like(site_count, dtype=np.int64)

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
        match_idx = np.searchsorted(uniq_d, uniq_windows)
        valid_match = (match_idx < uniq_d.size) & (uniq_d[match_idx] == uniq_windows)
        sum_density[valid_match] = sum_d[match_idx[valid_match]]
        density_site_count[valid_match] = count_d[match_idx[valid_match]]

    return seqname, (uniq_windows, sum_m, sum_t, site_count, sum_density, density_site_count), max_pos


def _finalize_tracks(
    *,
    context: Context,
    stat: ChrLineStat,
    bin_size_bp: int,
    empty_policy: ChrEmptyPolicy,
    smooth: dict | int | None,
    chr_order: list[str],
    chr_lengths_bp: dict[str, int],
    aggregate: dict[str, _ChrWindowPayload],
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

        chr_payload = aggregate.get(seqname)
        if chr_payload is not None and chr_payload.window_ids:
            win_ids = np.concatenate(chr_payload.window_ids)
            chr_sum_m = np.concatenate(chr_payload.sum_m)
            chr_sum_t = np.concatenate(chr_payload.sum_total)
            chr_sites = np.concatenate(chr_payload.site_count)
            chr_sum_density = np.concatenate(chr_payload.sum_density)
            chr_density_sites = np.concatenate(chr_payload.density_site_count)

            order = np.argsort(win_ids, kind="mergesort")
            win_ids = win_ids[order]
            chr_sum_m = chr_sum_m[order]
            chr_sum_t = chr_sum_t[order]
            chr_sites = chr_sites[order]
            chr_sum_density = chr_sum_density[order]
            chr_density_sites = chr_density_sites[order]

            cuts = np.flatnonzero(np.diff(win_ids)) + 1
            starts_idx = np.r_[0, cuts]
            uniq_windows = win_ids[starts_idx]
            sum_m[uniq_windows] = np.add.reduceat(chr_sum_m, starts_idx)
            sum_t[uniq_windows] = np.add.reduceat(chr_sum_t, starts_idx)
            sites[uniq_windows] = np.add.reduceat(chr_sites, starts_idx).astype(np.int64, copy=False)
            sum_density[uniq_windows] = np.add.reduceat(chr_sum_density, starts_idx)
            density_sites[uniq_windows] = np.add.reduceat(chr_density_sites, starts_idx).astype(
                np.int64,
                copy=False,
            )

        values = np.full(n_windows, np.nan, dtype=np.float64)
        if stat is ChrLineStat.MEAN:
            non_empty = density_sites > 0
            values[non_empty] = sum_density[non_empty] / density_sites[non_empty]
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

    aggregate: dict[str, _ChrWindowPayload] = {}
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

        window_ids, sum_m, sum_t, site_count, sum_density, density_site_count = payload
        if window_ids.size == 0:
            continue

        target = aggregate.setdefault(
            seqname,
            _ChrWindowPayload(
                window_ids=[],
                sum_m=[],
                sum_total=[],
                site_count=[],
                sum_density=[],
                density_site_count=[],
            ),
        )
        target.window_ids.append(window_ids)
        target.sum_m.append(sum_m)
        target.sum_total.append(sum_t)
        target.site_count.append(site_count)
        target.sum_density.append(sum_density)
        target.density_site_count.append(density_site_count)

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

    tracks = _finalize_tracks(
        context=context,
        stat=stat,
        bin_size_bp=bin_size_bp,
        empty_policy=empty_policy,
        smooth=smooth,
        chr_order=chr_order,
        chr_lengths_bp=chr_lengths_bp,
        aggregate=aggregate,
    )

    return ChrLineData(
        tracks=tracks,
        context=context,
        stat=stat,
        bin_size_bp=bin_size_bp,
        empty_policy=empty_policy,
    )


__all__ = [
    "ChrEmptyPolicy",
    "ChrLineData",
    "ChrLineStat",
    "ChrLineTrack",
    "compute_chr_line_data",
]
