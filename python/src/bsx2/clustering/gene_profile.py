from __future__ import annotations

import hashlib
import json
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from uuid import uuid4

import numpy as np
from beartype import beartype

from bsx2 import Contig, RegionReader, Strand

from .agg import coverage_weighted_ratio, finalize_gene_profile_matrix
from .config import ClusterConfig
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
) -> bool:
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        temp_path = cache_path.with_name(f"{cache_path.stem}.{uuid4().hex}.tmp.npz")
        if block_arrays is None:
            np.savez_compressed(
                temp_path,
                positions=np.empty(0, dtype=np.int64),
                count_m=np.empty(0, dtype=np.int64),
                count_total=np.empty(0, dtype=np.int64),
            )
        else:
            np.savez_compressed(
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
        cache_write = _write_query_block_cache(cache_path, block_arrays)
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


def _build_query_blocks(
    genes: list[GeneAnnotation],
    config: ClusterConfig,
) -> list[_QueryBlock]:
    spans_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for gene in genes:
        spans_by_chrom[gene.chrom].append(_gene_profile_span(gene, config))

    blocks: list[_QueryBlock] = []
    for chrom in sorted(spans_by_chrom):
        spans = sorted(spans_by_chrom[chrom])
        current_start, current_end = spans[0]
        for start, end in spans[1:]:
            if start <= current_end:
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
) -> list[_AssignedQueryBlock]:
    blocks_by_chrom: dict[str, list[tuple[int, _QueryBlock]]] = defaultdict(list)
    for block_idx, block in enumerate(query_blocks):
        blocks_by_chrom[block.chrom].append((block_idx, block))

    assignments: list[list[int]] = [[] for _ in query_blocks]
    for gene_idx, gene in enumerate(genes):
        span_start, span_end = _gene_profile_span(gene, config)
        matched = False
        for block_idx, block in blocks_by_chrom.get(gene.chrom, []):
            if span_start >= block.start and span_end <= block.end:
                assignments[block_idx].append(gene_idx)
                matched = True
                break
        if not matched:
            raise ValueError(
                f"Gene {gene.gene_id} could not be assigned to any query block on {gene.chrom}"
            )

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


def _integer_edges(start: int, end: int, n_bins: int) -> np.ndarray:
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    edges = np.floor(np.linspace(start, end, n_bins + 1)).astype(np.int64)
    edges[0] = start
    edges[-1] = end
    return np.maximum.accumulate(edges)


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


@beartype
def build_metagene_bins(
    gene: GeneAnnotation,
    config: ClusterConfig,
) -> list[ProfileBin]:
    gp = config.gene_profile
    chrom = gene.chrom
    bins: list[ProfileBin] = []

    if gene.strand == "+":
        upstream_region = (max(gene.start - gp.upstream_bp, 0), gene.start)
        body_region = (gene.start, gene.end)
        downstream_region = (gene.end, gene.end + gp.downstream_bp)
    else:
        upstream_region = (gene.end, gene.end + gp.upstream_bp)
        body_region = (gene.start, gene.end)
        downstream_region = (max(gene.start - gp.downstream_bp, 0), gene.start)

    bins.extend(
        _segment_profile_bins(
            gene=gene,
            chrom=chrom,
            segment="up",
            region_start=upstream_region[0],
            region_end=upstream_region[1],
            n_bins=gp.upstream_bins,
            global_offset=0,
        )
    )
    bins.extend(
        _segment_profile_bins(
            gene=gene,
            chrom=chrom,
            segment="body",
            region_start=body_region[0],
            region_end=body_region[1],
            n_bins=gp.body_bins,
            global_offset=gp.upstream_bins,
        )
    )
    bins.extend(
        _segment_profile_bins(
            gene=gene,
            chrom=chrom,
            segment="down",
            region_start=downstream_region[0],
            region_end=downstream_region[1],
            n_bins=gp.downstream_bins,
            global_offset=gp.upstream_bins + gp.body_bins,
        )
    )
    return bins


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


def _interval_ratio(
    positions: np.ndarray,
    count_m_cumsum: np.ndarray,
    count_total_cumsum: np.ndarray,
    *,
    start: int,
    end: int,
    min_bin_total_coverage: int,
) -> float:
    if end <= start:
        return float("nan")
    left = int(np.searchsorted(positions, start, side="left"))
    right = int(np.searchsorted(positions, end, side="left"))
    if right <= left:
        return float("nan")
    count_m_sum = int(count_m_cumsum[right] - count_m_cumsum[left])
    count_total_sum = int(count_total_cumsum[right] - count_total_cumsum[left])
    if count_total_sum < min_bin_total_coverage:
        return float("nan")
    return coverage_weighted_ratio(count_m_sum, count_total_sum)


def _profile_gene(
    gene: GeneAnnotation,
    config: ClusterConfig,
    block_cumsums: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
) -> np.ndarray:
    profile = np.full(config.gene_profile.total_bins, np.nan, dtype=float)
    if block_cumsums is None:
        return profile

    positions, count_m_cumsum, count_total_cumsum = block_cumsums
    for profile_bin in build_metagene_bins(gene, config):
        profile[profile_bin.global_bin_index] = _interval_ratio(
            positions,
            count_m_cumsum,
            count_total_cumsum,
            start=profile_bin.start,
            end=profile_bin.end,
            min_bin_total_coverage=config.gene_profile.min_bin_total_coverage,
        )
    return profile


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

    query_blocks = _build_query_blocks(genes, config)
    t_blocks = perf_counter()
    feature_bins = build_feature_bins(config)
    assigned_blocks = _assign_genes_to_query_blocks(genes, query_blocks, config)
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
        for gene_idx in assigned_block.gene_indices:
            values[gene_idx] = _profile_gene(genes[gene_idx], config, block_cumsums)
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
            "region_reader_cache": {
                "mode": "process_local_python_cache",
                "cache_hit": prepared_reader.cache_hit,
                "cache_entries": len(_REGION_READER_CACHE),
                "max_entries": _REGION_READER_CACHE_MAX_ENTRIES,
            },
            "query_block_cache": {
                "enabled": config.block_cache.enabled,
                "cache_dir": (
                    None
                    if not config.block_cache.enabled
                    else str(_resolve_query_block_cache_dir(config))
                ),
                "hits": block_cache_hits,
                "misses": block_cache_misses,
                "writes": block_cache_writes,
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
