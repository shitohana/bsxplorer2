from __future__ import annotations

from collections import defaultdict
from time import perf_counter
from typing import Iterable

import numpy as np

from bsx2 import BsxFileReader

from .agg import coverage_weighted_ratio, finalize_gene_profile_matrix
from .config import ClusterConfig
from .gene_annotation import load_gene_annotations
from .models import FeatureBin, GeneAnnotation, GeneProfileMatrix, ProfileBin


def _iter_filtered_batches(config: ClusterConfig) -> Iterable:
    reader = BsxFileReader(str(config.bsx_path))
    for batch in reader:
        filtered = batch.lazy()
        if config.read.min_coverage > 0:
            # Python bindings expose only `filter_coverage_gt`, so emulate >= min_coverage.
            filtered = filtered.filter_coverage_gt(max(config.read.min_coverage - 1, 0))
        if config.read.context is not None:
            filtered = filtered.filter_context(config.read.context)
        collected = filtered.collect()
        if not collected.is_empty():
            yield collected


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


def _load_bsx_arrays(
    config: ClusterConfig,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    positions_by_chrom: dict[str, list[np.ndarray]] = defaultdict(list)
    count_m_by_chrom: dict[str, list[np.ndarray]] = defaultdict(list)
    count_total_by_chrom: dict[str, list[np.ndarray]] = defaultdict(list)

    for batch in _iter_filtered_batches(config):
        chrom, positions, count_m, count_total = _extract_batch_arrays(batch)
        positions_by_chrom[chrom].append(positions)
        count_m_by_chrom[chrom].append(count_m)
        count_total_by_chrom[chrom].append(count_total)

    arrays_by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, pos_parts in positions_by_chrom.items():
        positions = np.concatenate(pos_parts) if len(pos_parts) > 1 else pos_parts[0]
        count_m = (
            np.concatenate(count_m_by_chrom[chrom])
            if len(count_m_by_chrom[chrom]) > 1
            else count_m_by_chrom[chrom][0]
        )
        count_total = (
            np.concatenate(count_total_by_chrom[chrom])
            if len(count_total_by_chrom[chrom]) > 1
            else count_total_by_chrom[chrom][0]
        )
        if positions.size > 1 and np.any(positions[:-1] > positions[1:]):
            order = np.argsort(positions, kind="stable")
            positions = positions[order]
            count_m = count_m[order]
            count_total = count_total[order]
        arrays_by_chrom[chrom] = (positions, count_m, count_total)
    return arrays_by_chrom


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


def _prepare_cumsums(
    arrays_by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    cumsums: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, (positions, count_m, count_total) in arrays_by_chrom.items():
        cumsums[chrom] = (
            positions,
            np.concatenate(([0], np.cumsum(count_m, dtype=np.int64))),
            np.concatenate(([0], np.cumsum(count_total, dtype=np.int64))),
        )
    return cumsums


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
    cumsums_by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> np.ndarray:
    profile = np.full(config.gene_profile.total_bins, np.nan, dtype=float)
    chrom_arrays = cumsums_by_chrom.get(gene.chrom)
    if chrom_arrays is None:
        return profile

    positions, count_m_cumsum, count_total_cumsum = chrom_arrays
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


def build_gene_profile_matrix(config: ClusterConfig) -> GeneProfileMatrix:
    t0 = perf_counter()
    genes = load_gene_annotations(
        config.annotation_path,
        annotation_format=config.annotation_format,
        min_gene_length_bp=config.gene_profile.min_gene_length_bp,
        limit=config.gene_profile.limit_genes,
    )
    t_annot = perf_counter()
    arrays_by_chrom = _load_bsx_arrays(config)
    t_bsx = perf_counter()
    if not arrays_by_chrom:
        raise ValueError(f"No methylation data remained after filtering in {config.bsx_path}")

    genes, skipped_chrom_mismatch = _resolve_gene_chromosomes(genes, set(arrays_by_chrom))
    if not genes:
        raise ValueError(
            "No overlapping genes found between annotation and BSX chromosomes"
        )

    cumsums_by_chrom = _prepare_cumsums(arrays_by_chrom)
    feature_bins = build_feature_bins(config)
    values = np.vstack([_profile_gene(gene, config, cumsums_by_chrom) for gene in genes])
    t_values = perf_counter()

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
            "skipped_chrom_mismatch": skipped_chrom_mismatch,
            "timings_s": {
                "load_annotations": round(t_annot - t0, 6),
                "load_bsx_arrays": round(t_bsx - t_annot, 6),
                "materialize_gene_profiles": round(t_values - t_bsx, 6),
                "finalize_matrix": round(t_finalize - t_values, 6),
                "total": round(t_finalize - t0, 6),
            },
        }
    )
    return matrix
