from __future__ import annotations

import hashlib
import json
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from time import perf_counter
from uuid import uuid4

import numpy as np
from beartype import beartype

from bsx2 import Contig, RegionReader, Strand

from .agg import finalize_gene_profile_matrix
from .config import BlockCacheMode, ClusterConfig
from .gene_annotation import load_gene_annotations
from .models import FeatureBin, GeneAnnotation, GeneProfileMatrix, ProfileBin


@dataclass(frozen=True)
class _QueryBlock:
    chrom: str
    start: int
    end: int

    def to_contig(self) -> Contig:
        return Contig(self.chrom, self.start, self.end, Strand.Null)


@dataclass(frozen=True)
class _AssignedQueryBlock:
    block: _QueryBlock
    gene_indices: tuple[int, ...]


@dataclass(frozen=True)
class _PreparedRegionReader:
    reader: RegionReader
    chr_order: tuple[str, ...]
    cache_hit: bool


@dataclass
class _CachedRegionReader:
    reader: RegionReader
    chr_order: tuple[str, ...]


@dataclass(frozen=True)
class _MetageneIntervalTemplate:
    total_bins: int
    up_steps: np.ndarray
    body_steps: np.ndarray
    down_steps: np.ndarray
    up_slice: slice
    body_slice: slice
    down_slice: slice


# Process-local cache: reuses the Python RegionReader wrapper and its in-memory index.
# This is intentionally small and should be treated as non-thread-safe shared state.
_REGION_READER_CACHE_MAX_ENTRIES = 4
_REGION_READER_CACHE: OrderedDict[tuple[str, int | None, int | None], _CachedRegionReader] = (
    OrderedDict()
)
_CACHE_MISS = object()


def _region_reader_cache_key(path: Path) -> tuple[str, int | None, int | None]:
    resolved = path.resolve(strict=False)
    try:
        stat = resolved.stat()
    except OSError:
        return str(resolved), None, None
    return str(resolved), int(stat.st_size), int(stat.st_mtime_ns)


@beartype
def clear_region_reader_cache() -> None:
    _REGION_READER_CACHE.clear()


def _resolve_query_block_cache_dir(config: ClusterConfig) -> Path:
    return (
        config.block_cache.cache_dir
        if config.block_cache.cache_dir is not None
        else config.output.output_dir / ".bsx2_cache" / "queried_blocks"
    )


def _query_block_cache_key(
    config: ClusterConfig,
    block: _QueryBlock,
) -> str:
    payload = {
        "version": 1,
        "bsx": _region_reader_cache_key(config.bsx_path),
        "context": None if config.read.context is None else config.read.context.name,
        "min_coverage": config.read.min_coverage,
        "block": {
            "chrom": block.chrom,
            "start": block.start,
            "end": block.end,
        },
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.blake2b(encoded, digest_size=16).hexdigest()


def _query_block_cache_path(
    config: ClusterConfig,
    block: _QueryBlock,
) -> Path:
    return _resolve_query_block_cache_dir(config) / f"{_query_block_cache_key(config, block)}.npz"


def _read_query_block_cache(cache_path: Path):
    if not cache_path.exists():
        return _CACHE_MISS
    try:
        with np.load(cache_path, allow_pickle=False) as data:
            positions = np.asarray(data["positions"], dtype=np.int64)
            count_m = np.asarray(data["count_m"], dtype=np.int64)
            count_total = np.asarray(data["count_total"], dtype=np.int64)
    except Exception:
        try:
            cache_path.unlink(missing_ok=True)
        except OSError:
            pass
        return _CACHE_MISS
    if positions.size == 0:
        return None
    return positions, count_m, count_total


def _write_query_block_cache(
    cache_path: Path,
    block_arrays: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
    *,
    mode: BlockCacheMode,
) -> bool:
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        temp_path = cache_path.with_name(f"{cache_path.stem}.{uuid4().hex}.tmp.npz")
        savez = np.savez_compressed if mode is BlockCacheMode.COMPRESSED else np.savez
        if block_arrays is None:
            savez(
                temp_path,
                positions=np.empty(0, dtype=np.int64),
                count_m=np.empty(0, dtype=np.int64),
                count_total=np.empty(0, dtype=np.int64),
            )
        else:
            savez(
                temp_path,
                positions=block_arrays[0],
                count_m=block_arrays[1],
                count_total=block_arrays[2],
            )
        temp_path.replace(cache_path)
        return True
    except OSError:
        return False


def _load_or_query_block_arrays(
    reader: RegionReader,
    block: _QueryBlock,
    config: ClusterConfig,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray] | None, dict[str, float | bool]]:
    if config.block_cache.enabled:
        cache_path = _query_block_cache_path(config, block)
        t_cache_read = perf_counter()
        cached = _read_query_block_cache(cache_path)
        cache_read_s = perf_counter() - t_cache_read
        if cached is not _CACHE_MISS:
            return cached, {
                "cache_hit": True,
                "cache_write": False,
                "cache_read_s": cache_read_s,
                "cache_write_s": 0.0,
                "query_s": 0.0,
            }
    else:
        cache_path = None
        cache_read_s = 0.0

    t_query = perf_counter()
    block_arrays = _query_block_arrays(reader, block)
    query_s = perf_counter() - t_query

    cache_write_s = 0.0
    cache_write = False
    if config.block_cache.enabled and cache_path is not None:
        t_cache_write = perf_counter()
        cache_write = _write_query_block_cache(
            cache_path,
            block_arrays,
            mode=config.block_cache.mode,
        )
        cache_write_s = perf_counter() - t_cache_write

    return block_arrays, {
        "cache_hit": False,
        "cache_write": cache_write,
        "cache_read_s": cache_read_s,
        "cache_write_s": cache_write_s,
        "query_s": query_s,
    }


def _apply_reader_filters(reader: RegionReader, config: ClusterConfig) -> None:
    if config.read.min_coverage > 0:
        # Python bindings expose only `filter_coverage_gt`, so emulate >= min_coverage.
        reader.filter_coverage_gt(max(config.read.min_coverage - 1, 0))
    if config.read.context is not None:
        reader.filter_context(config.read.context)


def _make_region_reader(config: ClusterConfig) -> _PreparedRegionReader:
    cache_key = _region_reader_cache_key(config.bsx_path)
    cached = _REGION_READER_CACHE.pop(cache_key, None)
    if cached is not None:
        try:
            cached.reader.reset()
            cached.reader.clear_filters()
            _apply_reader_filters(cached.reader, config)
            _REGION_READER_CACHE[cache_key] = cached
            return _PreparedRegionReader(
                reader=cached.reader,
                chr_order=cached.chr_order,
                cache_hit=True,
            )
        except Exception:
            pass

    reader = RegionReader(str(config.bsx_path))
    _apply_reader_filters(reader, config)
    cached = _CachedRegionReader(reader=reader, chr_order=tuple(reader.chr_order()))
    _REGION_READER_CACHE[cache_key] = cached
    while len(_REGION_READER_CACHE) > _REGION_READER_CACHE_MAX_ENTRIES:
        _REGION_READER_CACHE.popitem(last=False)
    return _PreparedRegionReader(
        reader=cached.reader,
        chr_order=cached.chr_order,
        cache_hit=False,
    )


def _series_to_numpy(series, dtype) -> np.ndarray:
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
    return values.astype(dtype, copy=False)


def _extract_batch_arrays(batch) -> tuple[str, np.ndarray, np.ndarray, np.ndarray]:
    chrom = batch.seqname()
    if chrom is None:
        raise ValueError("Encountered batch without chromosome name")

    positions = _series_to_numpy(batch.position(), np.int64)
    count_m = _series_to_numpy(batch.count_m(), np.int64)
    count_total = _series_to_numpy(batch.count_total(), np.int64)
    return str(chrom), positions, count_m, count_total


def _gene_profile_span(
    gene: GeneAnnotation,
    config: ClusterConfig,
) -> tuple[int, int]:
    gp = config.gene_profile
    if gene.strand == "+":
        return max(gene.start - gp.upstream_bp, 0), gene.end + gp.downstream_bp
    return max(gene.start - gp.downstream_bp, 0), gene.end + gp.upstream_bp


def _collect_gene_profile_spans(
    genes: list[GeneAnnotation],
    config: ClusterConfig,
) -> list[tuple[int, int]]:
    return [_gene_profile_span(gene, config) for gene in genes]


def _build_query_blocks(
    genes: list[GeneAnnotation],
    config: ClusterConfig,
    *,
    gene_spans: list[tuple[int, int]] | None = None,
) -> list[_QueryBlock]:
    spans = _collect_gene_profile_spans(genes, config) if gene_spans is None else gene_spans
    merge_gap_bp = max(config.read.query_block_merge_gap_bp, 0)
    spans_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for gene, span in zip(genes, spans, strict=True):
        spans_by_chrom[gene.chrom].append(span)

    blocks: list[_QueryBlock] = []
    for chrom in sorted(spans_by_chrom):
        spans = sorted(spans_by_chrom[chrom])
        current_start, current_end = spans[0]
        for start, end in spans[1:]:
            if start <= current_end + merge_gap_bp:
                current_end = max(current_end, end)
            else:
                blocks.append(_QueryBlock(chrom=chrom, start=current_start, end=current_end))
                current_start, current_end = start, end
        blocks.append(_QueryBlock(chrom=chrom, start=current_start, end=current_end))
    return blocks


def _assign_genes_to_query_blocks(
    genes: list[GeneAnnotation],
    query_blocks: list[_QueryBlock],
    config: ClusterConfig,
    *,
    gene_spans: list[tuple[int, int]] | None = None,
) -> list[_AssignedQueryBlock]:
    spans = _collect_gene_profile_spans(genes, config) if gene_spans is None else gene_spans
    blocks_by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    grouped: dict[str, list[tuple[int, _QueryBlock]]] = defaultdict(list)
    for block_idx, block in enumerate(query_blocks):
        grouped[block.chrom].append((block_idx, block))

    for chrom, blocks in grouped.items():
        block_indices = np.asarray([block_idx for block_idx, _ in blocks], dtype=np.int64)
        block_starts = np.asarray([block.start for _, block in blocks], dtype=np.int64)
        block_ends = np.asarray([block.end for _, block in blocks], dtype=np.int64)
        blocks_by_chrom[chrom] = (block_indices, block_starts, block_ends)

    assignments: list[list[int]] = [[] for _ in query_blocks]
    for gene_idx, (gene, (span_start, span_end)) in enumerate(zip(genes, spans, strict=True)):
        block_arrays = blocks_by_chrom.get(gene.chrom)
        if block_arrays is None:
            raise ValueError(
                f"Gene {gene.gene_id} could not be assigned to any query block on {gene.chrom}"
            )
        block_indices, block_starts, block_ends = block_arrays
        candidate = int(np.searchsorted(block_starts, span_start, side="right") - 1)
        if candidate < 0 or span_end > int(block_ends[candidate]):
            raise ValueError(
                f"Gene {gene.gene_id} could not be assigned to any query block on {gene.chrom}"
            )
        assignments[int(block_indices[candidate])].append(gene_idx)

    return [
        _AssignedQueryBlock(block=block, gene_indices=tuple(assignments[idx]))
        for idx, block in enumerate(query_blocks)
        if assignments[idx]
    ]


def _query_block_arrays(
    reader: RegionReader,
    block: _QueryBlock,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    batch = reader.query(block.to_contig())
    if batch is None or batch.is_empty():
        return None

    _, positions, count_m, count_total = _extract_batch_arrays(batch)
    if positions.size == 0:
        return None
    if positions.size > 1 and np.any(positions[:-1] > positions[1:]):
        order = np.argsort(positions, kind="stable")
        positions = positions[order]
        count_m = count_m[order]
        count_total = count_total[order]
    return positions, count_m, count_total


def _canonical_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _resolve_gene_chromosomes(
    genes: list[GeneAnnotation],
    available_chroms: set[str],
) -> tuple[list[GeneAnnotation], int]:
    canonical_lookup: dict[str, list[str]] = defaultdict(list)
    for chrom in sorted(available_chroms):
        canonical_lookup[_canonical_chrom(chrom)].append(chrom)

    remapped: list[GeneAnnotation] = []
    skipped = 0
    for gene in genes:
        if gene.chrom in available_chroms:
            remapped.append(gene)
            continue
        matches = canonical_lookup.get(_canonical_chrom(gene.chrom), [])
        if len(matches) == 1:
            remapped.append(
                GeneAnnotation(
                    gene_id=gene.gene_id,
                    gene_name=gene.gene_name,
                    chrom=matches[0],
                    start=gene.start,
                    end=gene.end,
                    strand=gene.strand,
                )
            )
            continue
        skipped += 1
    return remapped, skipped


@beartype
def build_feature_bins(config: ClusterConfig) -> list[FeatureBin]:
    feature_bins: list[FeatureBin] = []
    segments = (
        ("up", config.gene_profile.upstream_bins),
        ("body", config.gene_profile.body_bins),
        ("down", config.gene_profile.downstream_bins),
    )
    global_index = 0
    for segment, n_bins in segments:
        for local_index in range(n_bins):
            feature_bins.append(
                FeatureBin(
                    feature_name=f"{segment}_{local_index + 1}",
                    segment=segment,
                    local_bin_index=local_index,
                    global_bin_index=global_index,
                )
            )
            global_index += 1
    return feature_bins


@lru_cache(maxsize=16)
def _build_metagene_interval_template(
    upstream_bins: int,
    body_bins: int,
    downstream_bins: int,
) -> _MetageneIntervalTemplate:
    up_slice = slice(0, upstream_bins)
    body_slice = slice(upstream_bins, upstream_bins + body_bins)
    down_slice = slice(upstream_bins + body_bins, upstream_bins + body_bins + downstream_bins)
    return _MetageneIntervalTemplate(
        total_bins=upstream_bins + body_bins + downstream_bins,
        up_steps=np.arange(upstream_bins + 1, dtype=np.int64),
        body_steps=np.arange(body_bins + 1, dtype=np.int64),
        down_steps=np.arange(downstream_bins + 1, dtype=np.int64),
        up_slice=up_slice,
        body_slice=body_slice,
        down_slice=down_slice,
    )


def _metagene_interval_template(config: ClusterConfig) -> _MetageneIntervalTemplate:
    gp = config.gene_profile
    return _build_metagene_interval_template(
        gp.upstream_bins,
        gp.body_bins,
        gp.downstream_bins,
    )


def _integer_edges(start: int, end: int, n_bins: int) -> np.ndarray:
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    clipped_start = max(start, 0)
    clipped_end = max(end, 0)
    length = max(clipped_end - clipped_start, 0)
    steps = np.arange(n_bins + 1, dtype=np.int64)
    return clipped_start + ((length * steps) // n_bins)


def _segment_profile_bins(
    *,
    gene: GeneAnnotation,
    chrom: str,
    segment: str,
    region_start: int,
    region_end: int,
    n_bins: int,
    global_offset: int,
) -> list[ProfileBin]:
    edges = _integer_edges(max(region_start, 0), max(region_end, 0), n_bins)
    raw_pairs = [(int(edges[idx]), int(edges[idx + 1])) for idx in range(n_bins)]
    ordered_pairs = raw_pairs if gene.strand == "+" else list(reversed(raw_pairs))

    bins: list[ProfileBin] = []
    for local_index, (start, end) in enumerate(ordered_pairs):
        feature_name = f"{segment}_{local_index + 1}"
        bins.append(
            ProfileBin(
                gene_id=gene.gene_id,
                chrom=chrom,
                start=start,
                end=end,
                strand=gene.strand,
                segment=segment,
                local_bin_index=local_index,
                global_bin_index=global_offset + local_index,
                feature_name=feature_name,
            )
    )
    return bins


def _metagene_segments(
    gene: GeneAnnotation,
    config: ClusterConfig,
) -> tuple[tuple[str, int, int, int, int], ...]:
    gp = config.gene_profile
    if gene.strand == "+":
        return (
            ("up", max(gene.start - gp.upstream_bp, 0), gene.start, gp.upstream_bins, 0),
            ("body", gene.start, gene.end, gp.body_bins, gp.upstream_bins),
            (
                "down",
                gene.end,
                gene.end + gp.downstream_bp,
                gp.downstream_bins,
                gp.upstream_bins + gp.body_bins,
            ),
        )
    return (
        ("up", gene.end, gene.end + gp.upstream_bp, gp.upstream_bins, 0),
        ("body", gene.start, gene.end, gp.body_bins, gp.upstream_bins),
        (
            "down",
            max(gene.start - gp.downstream_bp, 0),
            gene.start,
            gp.downstream_bins,
            gp.upstream_bins + gp.body_bins,
        ),
    )


@beartype
def build_metagene_bins(
    gene: GeneAnnotation,
    config: ClusterConfig,
) -> list[ProfileBin]:
    chrom = gene.chrom
    bins: list[ProfileBin] = []
    for segment, region_start, region_end, n_bins, global_offset in _metagene_segments(
        gene,
        config,
    ):
        bins.extend(
            _segment_profile_bins(
                gene=gene,
                chrom=chrom,
                segment=segment,
                region_start=region_start,
                region_end=region_end,
                n_bins=n_bins,
                global_offset=global_offset,
            )
        )
    return bins


def _fill_segment_interval_arrays(
    starts_out: np.ndarray,
    ends_out: np.ndarray,
    *,
    gene: GeneAnnotation,
    region_start: int,
    region_end: int,
    n_bins: int,
    global_offset: int,
) -> None:
    edges = _integer_edges(max(region_start, 0), max(region_end, 0), n_bins)
    starts = edges[:-1]
    ends = edges[1:]
    if gene.strand == "-":
        starts = starts[::-1]
        ends = ends[::-1]
    starts_out[global_offset: global_offset + n_bins] = starts
    ends_out[global_offset: global_offset + n_bins] = ends


def _fill_segment_interval_matrix(
    starts_out: np.ndarray,
    ends_out: np.ndarray,
    *,
    segment_slice: slice,
    region_starts: np.ndarray,
    region_ends: np.ndarray,
    steps: np.ndarray,
    reverse_mask: np.ndarray,
) -> None:
    n_bins = int(steps.size - 1)
    if n_bins <= 0 or starts_out.shape[0] == 0:
        return

    safe_starts = np.maximum(region_starts.astype(np.int64, copy=False), 0)
    safe_ends = np.maximum(region_ends.astype(np.int64, copy=False), 0)
    lengths = np.maximum(safe_ends - safe_starts, 0)
    edges = safe_starts[:, np.newaxis] + ((lengths[:, np.newaxis] * steps[np.newaxis, :]) // n_bins)
    seg_starts = edges[:, :-1]
    seg_ends = edges[:, 1:]

    if np.any(reverse_mask):
        seg_starts = seg_starts.copy()
        seg_ends = seg_ends.copy()
        seg_starts[reverse_mask] = seg_starts[reverse_mask, ::-1]
        seg_ends[reverse_mask] = seg_ends[reverse_mask, ::-1]

    starts_out[:, segment_slice] = seg_starts
    ends_out[:, segment_slice] = seg_ends


def _build_metagene_interval_matrices(
    genes: list[GeneAnnotation],
    config: ClusterConfig,
) -> tuple[np.ndarray, np.ndarray]:
    template = _metagene_interval_template(config)
    n_genes = len(genes)
    starts = np.empty((n_genes, template.total_bins), dtype=np.int64)
    ends = np.empty((n_genes, template.total_bins), dtype=np.int64)
    if n_genes == 0:
        return starts, ends

    gene_starts = np.asarray([gene.start for gene in genes], dtype=np.int64)
    gene_ends = np.asarray([gene.end for gene in genes], dtype=np.int64)
    reverse_mask = np.asarray([gene.strand == "-" for gene in genes], dtype=bool)
    gp = config.gene_profile

    up_region_starts = np.where(reverse_mask, gene_ends, np.maximum(gene_starts - gp.upstream_bp, 0))
    up_region_ends = np.where(reverse_mask, gene_ends + gp.upstream_bp, gene_starts)
    body_region_starts = gene_starts
    body_region_ends = gene_ends
    down_region_starts = np.where(reverse_mask, np.maximum(gene_starts - gp.downstream_bp, 0), gene_ends)
    down_region_ends = np.where(reverse_mask, gene_starts, gene_ends + gp.downstream_bp)

    _fill_segment_interval_matrix(
        starts,
        ends,
        segment_slice=template.up_slice,
        region_starts=up_region_starts,
        region_ends=up_region_ends,
        steps=template.up_steps,
        reverse_mask=reverse_mask,
    )
    _fill_segment_interval_matrix(
        starts,
        ends,
        segment_slice=template.body_slice,
        region_starts=body_region_starts,
        region_ends=body_region_ends,
        steps=template.body_steps,
        reverse_mask=reverse_mask,
    )
    _fill_segment_interval_matrix(
        starts,
        ends,
        segment_slice=template.down_slice,
        region_starts=down_region_starts,
        region_ends=down_region_ends,
        steps=template.down_steps,
        reverse_mask=reverse_mask,
    )
    return starts, ends


def _build_metagene_interval_arrays(
    gene: GeneAnnotation,
    config: ClusterConfig,
) -> tuple[np.ndarray, np.ndarray]:
    starts, ends = _build_metagene_interval_matrices([gene], config)
    return starts[0], ends[0]


def _prepare_block_cumsums(
    positions: np.ndarray,
    count_m: np.ndarray,
    count_total: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (
        positions,
        np.concatenate(([0], np.cumsum(count_m, dtype=np.int64))),
        np.concatenate(([0], np.cumsum(count_total, dtype=np.int64))),
    )


def _profile_gene(
    gene: GeneAnnotation,
    config: ClusterConfig,
    block_cumsums: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
) -> np.ndarray:
    profile = np.full(config.gene_profile.total_bins, np.nan, dtype=float)
    if block_cumsums is None:
        return profile

    positions, count_m_cumsum, count_total_cumsum = block_cumsums
    starts, ends = _build_metagene_interval_arrays(gene, config)
    left = np.searchsorted(positions, starts, side="left")
    right = np.searchsorted(positions, ends, side="left")
    valid = (ends > starts) & (right > left)
    if not np.any(valid):
        return profile

    count_m_sum = count_m_cumsum[right] - count_m_cumsum[left]
    count_total_sum = count_total_cumsum[right] - count_total_cumsum[left]
    valid &= count_total_sum >= config.gene_profile.min_bin_total_coverage
    if not np.any(valid):
        return profile

    profile[valid] = count_m_sum[valid].astype(float) / count_total_sum[valid].astype(float)
    return profile


def _profile_genes_in_block(
    genes: list[GeneAnnotation],
    gene_indices: tuple[int, ...],
    config: ClusterConfig,
    block_cumsums: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
) -> np.ndarray:
    profiles = np.full((len(gene_indices), config.gene_profile.total_bins), np.nan, dtype=float)
    if block_cumsums is None or not gene_indices:
        return profiles

    positions, count_m_cumsum, count_total_cumsum = block_cumsums
    block_genes = [genes[gene_idx] for gene_idx in gene_indices]
    starts, ends = _build_metagene_interval_matrices(block_genes, config)
    starts_flat = starts.reshape(-1)
    ends_flat = ends.reshape(-1)

    left = np.searchsorted(positions, starts_flat, side="left")
    right = np.searchsorted(positions, ends_flat, side="left")
    valid = (ends_flat > starts_flat) & (right > left)
    if not np.any(valid):
        return profiles

    count_m_sum = count_m_cumsum[right] - count_m_cumsum[left]
    count_total_sum = count_total_cumsum[right] - count_total_cumsum[left]
    valid &= count_total_sum >= config.gene_profile.min_bin_total_coverage
    if not np.any(valid):
        return profiles

    profiles_flat = profiles.reshape(-1)
    profiles_flat[valid] = (
        count_m_sum[valid].astype(float) / count_total_sum[valid].astype(float)
    )
    return profiles


@beartype
def build_gene_profile_matrix(config: ClusterConfig) -> GeneProfileMatrix:
    t0 = perf_counter()
    genes = load_gene_annotations(
        config.annotation_path,
        annotation_format=config.annotation_format,
        min_gene_length_bp=config.gene_profile.min_gene_length_bp,
        limit=config.gene_profile.limit_genes,
    )
    t_annot = perf_counter()
    prepared_reader = _make_region_reader(config)
    reader = prepared_reader.reader
    available_chroms = set(prepared_reader.chr_order)
    t_reader = perf_counter()
    genes, skipped_chrom_mismatch = _resolve_gene_chromosomes(genes, available_chroms)
    if not genes:
        raise ValueError(
            "No overlapping genes found between annotation and BSX chromosomes"
        )

    gene_spans = _collect_gene_profile_spans(genes, config)
    query_blocks = _build_query_blocks(genes, config, gene_spans=gene_spans)
    t_blocks = perf_counter()
    feature_bins = build_feature_bins(config)
    assigned_blocks = _assign_genes_to_query_blocks(
        genes,
        query_blocks,
        config,
        gene_spans=gene_spans,
    )
    values = np.full((len(genes), len(feature_bins)), np.nan, dtype=float)

    t_query_blocks = 0.0
    t_block_cache_reads = 0.0
    t_block_cache_writes = 0.0
    t_profile_blocks = 0.0
    last_chrom: str | None = None
    blocks_with_data = 0
    block_cache_hits = 0
    block_cache_misses = 0
    block_cache_writes = 0

    for assigned_block in assigned_blocks:
        block = assigned_block.block
        if block.chrom != last_chrom:
            try:
                reader.reset()
            except Exception:
                pass
            last_chrom = block.chrom

        block_arrays, block_io = _load_or_query_block_arrays(reader, block, config)
        t_query_blocks += float(block_io["query_s"])
        t_block_cache_reads += float(block_io["cache_read_s"])
        t_block_cache_writes += float(block_io["cache_write_s"])
        if bool(block_io["cache_hit"]):
            block_cache_hits += 1
        else:
            block_cache_misses += 1
        if bool(block_io["cache_write"]):
            block_cache_writes += 1

        block_cumsums = (
            None
            if block_arrays is None
            else _prepare_block_cumsums(*block_arrays)
        )
        if block_cumsums is not None:
            blocks_with_data += 1

        t_block_profile = perf_counter()
        block_profiles = _profile_genes_in_block(
            genes,
            assigned_block.gene_indices,
            config,
            block_cumsums,
        )
        values[np.asarray(assigned_block.gene_indices, dtype=np.int64)] = block_profiles
        t_profile_blocks += perf_counter() - t_block_profile

    t_values = perf_counter()

    if np.isnan(values).all():
        raise ValueError(f"No methylation data remained after filtering in {config.bsx_path}")

    matrix = finalize_gene_profile_matrix(
        genes=genes,
        feature_bins=feature_bins,
        values=values,
        max_gene_missing_rate=config.gene_profile.max_gene_missing_rate,
        max_feature_missing_rate=config.gene_profile.max_feature_missing_rate,
        min_gene_profile_variance=config.gene_profile.min_gene_profile_variance,
        min_feature_variance=config.gene_profile.min_feature_variance,
    )
    t_finalize = perf_counter()
    matrix.metadata.update(
        {
            "annotation_path": str(config.annotation_path),
            "bsx_path": str(config.bsx_path),
            "context": None if config.read.context is None else str(config.read.context),
            "min_coverage": config.read.min_coverage,
            "profile": {
                "upstream_bp": config.gene_profile.upstream_bp,
                "downstream_bp": config.gene_profile.downstream_bp,
                "upstream_bins": config.gene_profile.upstream_bins,
                "body_bins": config.gene_profile.body_bins,
                "downstream_bins": config.gene_profile.downstream_bins,
                "min_gene_length_bp": config.gene_profile.min_gene_length_bp,
                "min_bin_total_coverage": config.gene_profile.min_bin_total_coverage,
                "normalization": config.gene_profile.normalization.value,
            },
            "read_strategy": "region_reader_selective",
            "materialization_strategy": "batched_query_block_profiles",
            "region_reader_cache": {
                "mode": "process_local_python_cache",
                "cache_hit": prepared_reader.cache_hit,
                "cache_entries": len(_REGION_READER_CACHE),
                "max_entries": _REGION_READER_CACHE_MAX_ENTRIES,
            },
            "query_block_cache": {
                "enabled": config.block_cache.enabled,
                "mode": config.block_cache.mode.value,
                "cache_dir": (
                    None
                    if not config.block_cache.enabled
                    else str(_resolve_query_block_cache_dir(config))
                ),
                "hits": block_cache_hits,
                "misses": block_cache_misses,
                "writes": block_cache_writes,
            },
            "query_block_planning": {
                "strategy": "chromosome_span_merge_with_gap",
                "merge_gap_bp": config.read.query_block_merge_gap_bp,
            },
            "query_blocks": len(query_blocks),
            "query_blocks_with_data": blocks_with_data,
            "skipped_chrom_mismatch": skipped_chrom_mismatch,
            "timings_s": {
                "load_annotations": round(t_annot - t0, 6),
                "init_region_reader": round(t_reader - t_annot, 6),
                "build_query_blocks": round(t_blocks - t_reader, 6),
                "load_query_block_cache": round(t_block_cache_reads, 6),
                "query_region_blocks": round(t_query_blocks, 6),
                "write_query_block_cache": round(t_block_cache_writes, 6),
                "materialize_gene_profiles": round(t_profile_blocks, 6),
                "finalize_matrix": round(t_finalize - t_values, 6),
                "total": round(t_finalize - t0, 6),
            },
        }
    )
    return matrix
