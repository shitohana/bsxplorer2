from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from beartype import beartype
from beartype.typing import Sequence

from bsx2.clustering.models import FeatureBin, GeneClusterResult

from .data import DiscreteRegionData
from .metagene import MetageneProfileSegment


@beartype
@dataclass(frozen=True)
class GeneEmbeddingData:
    gene_ids: list[str]
    gene_names: list[str | None]
    chromosomes: list[str]
    starts: np.ndarray
    ends: np.ndarray
    strands: list[str]
    labels: np.ndarray
    embedding: np.ndarray
    cluster_source: str

    def __post_init__(self) -> None:
        n_genes = len(self.gene_ids)
        if self.embedding.ndim != 2:
            raise ValueError("GeneEmbeddingData.embedding must be 2D")
        if self.embedding.shape[0] != n_genes:
            raise ValueError("GeneEmbeddingData.embedding row count must match gene_ids")
        if self.labels.shape != (n_genes,):
            raise ValueError("GeneEmbeddingData.labels shape must match gene_ids")
        if len(self.gene_names) != n_genes:
            raise ValueError("GeneEmbeddingData.gene_names length must match gene_ids")
        if len(self.chromosomes) != n_genes:
            raise ValueError("GeneEmbeddingData.chromosomes length must match gene_ids")
        if self.starts.shape != (n_genes,):
            raise ValueError("GeneEmbeddingData.starts shape must match gene_ids")
        if self.ends.shape != (n_genes,):
            raise ValueError("GeneEmbeddingData.ends shape must match gene_ids")
        if len(self.strands) != n_genes:
            raise ValueError("GeneEmbeddingData.strands length must match gene_ids")


@beartype
@dataclass(frozen=True)
class GeneDendrogramData:
    linkage_matrix: np.ndarray
    leaf_labels: list[str]
    leaf_order: np.ndarray

    def __post_init__(self) -> None:
        if self.linkage_matrix.ndim != 2 or self.linkage_matrix.shape[1] != 4:
            raise ValueError("GeneDendrogramData.linkage_matrix must have shape [n-1, 4]")
        if self.leaf_order.ndim != 1:
            raise ValueError("GeneDendrogramData.leaf_order must be 1D")
        if len(self.leaf_labels) == 0:
            raise ValueError("GeneDendrogramData.leaf_labels must not be empty")


@beartype
@dataclass(frozen=True)
class ClusterMetageneGroup:
    label: str
    gene_ids: list[str]
    gene_count: int
    cluster_id: int | None = None

    def __post_init__(self) -> None:
        if not self.label:
            raise ValueError("ClusterMetageneGroup.label must not be empty")
        if self.gene_count < 1:
            raise ValueError("ClusterMetageneGroup.gene_count must be >= 1")
        if len(self.gene_ids) != self.gene_count:
            raise ValueError("ClusterMetageneGroup.gene_ids length must match gene_count")


@beartype
@dataclass(frozen=True)
class ClusterMetageneData:
    profiles: DiscreteRegionData
    segments: list[MetageneProfileSegment]
    groups: list[ClusterMetageneGroup]

    def __post_init__(self) -> None:
        n_profiles = len(self.profiles.positions)
        if len(self.groups) != n_profiles:
            raise ValueError("ClusterMetageneData.groups length must match profile count")
        if len(self.profiles.densities) != n_profiles or len(self.profiles.labels) != n_profiles:
            raise ValueError("ClusterMetageneData.profiles storage is internally inconsistent")
        for group, label in zip(self.groups, self.profiles.labels, strict=True):
            if label != group.label:
                raise ValueError("ClusterMetageneData group labels must match profile labels")


def _cluster_label(cluster: int, *, cluster_source: str) -> str:
    return f"{cluster_source} cluster {cluster}"


@beartype
def cluster_profile_segments(feature_bins: list[FeatureBin]) -> list[MetageneProfileSegment]:
    """
    Collapse ordered feature bins into contiguous metagene segments.
    """
    segments: list[MetageneProfileSegment] = []
    current_name: str | None = None
    current_count = 0

    for feature_bin in feature_bins:
        segment = str(getattr(feature_bin, "segment"))
        if current_name is None:
            current_name = segment
            current_count = 1
            continue
        if segment == current_name:
            current_count += 1
            continue
        segments.append(MetageneProfileSegment(current_name, current_count))
        current_name = segment
        current_count = 1

    if current_name is not None:
        segments.append(MetageneProfileSegment(current_name, current_count))
    return segments


@beartype
def _normalise_cluster_ids(cluster_ids: Sequence[int] | None) -> list[int] | None:
    if cluster_ids is None:
        return None
    out = sorted({int(cluster_id) for cluster_id in cluster_ids})
    if not out:
        raise ValueError("cluster_ids must not be empty")
    return out


@beartype
def _normalise_gene_ids(gene_ids: Sequence[str] | None) -> list[str] | None:
    if gene_ids is None:
        return None
    out: list[str] = []
    seen: set[str] = set()
    for gene_id in gene_ids:
        value = str(gene_id)
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    if not out:
        raise ValueError("gene_ids must not be empty")
    return out


@beartype
def _resolve_gene_mask(
    result: GeneClusterResult,
    *,
    cluster_ids: list[int] | None,
    gene_ids: list[str] | None,
) -> np.ndarray:
    labels = np.asarray(result.labels, dtype=np.int64)
    mask = np.ones(result.feature_matrix.n_genes, dtype=bool)

    if cluster_ids is not None:
        available = sorted(np.unique(labels).tolist())
        missing_clusters = [cluster_id for cluster_id in cluster_ids if cluster_id not in available]
        if missing_clusters:
            raise ValueError(f"Unknown cluster_ids: {missing_clusters}; available clusters: {available}")
        mask &= np.isin(labels, np.asarray(cluster_ids, dtype=np.int64))

    if gene_ids is not None:
        available_gene_ids = result.feature_matrix.gene_ids
        missing_gene_ids = [gene_id for gene_id in gene_ids if gene_id not in available_gene_ids]
        if missing_gene_ids:
            raise ValueError(f"Unknown gene_ids: {missing_gene_ids}")
        mask &= np.isin(np.asarray(available_gene_ids, dtype=object), np.asarray(gene_ids, dtype=object))

    if not np.any(mask):
        raise ValueError("Cluster metagene selection produced no genes")
    return mask


def _aggregate_cluster_profile(
    values: np.ndarray,
    *,
    total_bins: int,
) -> np.ndarray:
    finite = np.isfinite(values)
    sums = np.where(finite, values, 0.0).sum(axis=0)
    counts = finite.sum(axis=0, dtype=np.int64)
    return np.divide(
        sums,
        counts,
        out=np.full(total_bins, np.nan, dtype=np.float64),
        where=counts > 0,
    ).astype(np.float64, copy=False)


@beartype
def _default_collapsed_label(
    result: GeneClusterResult,
    *,
    selected_labels: np.ndarray,
    label: str | None,
) -> str:
    if label is not None:
        return str(label)
    unique_labels = np.unique(selected_labels)
    if unique_labels.size == 1:
        return _cluster_label(int(unique_labels[0]), cluster_source=result.cluster_source)
    return "selected genes"


@beartype
def build_cluster_metagene_data(
    result: GeneClusterResult,
    *,
    cluster_ids: Sequence[int] | None = None,
    gene_ids: Sequence[str] | None = None,
    collapse: bool = False,
    label: str | None = None,
) -> ClusterMetageneData:
    """
    Build metagene profiles for selected gene clusters.
    """
    profiles = DiscreteRegionData()
    values = np.asarray(result.feature_matrix.values, dtype=float)
    labels = np.asarray(result.labels, dtype=np.int64)
    genes = result.feature_matrix.genes
    total_bins = result.feature_matrix.n_features
    x_vals = (np.arange(total_bins, dtype=np.float64) + 0.5) / float(total_bins)
    groups: list[ClusterMetageneGroup] = []

    selected_cluster_ids = _normalise_cluster_ids(cluster_ids)
    selected_gene_ids = _normalise_gene_ids(gene_ids)
    gene_mask = _resolve_gene_mask(
        result,
        cluster_ids=selected_cluster_ids,
        gene_ids=selected_gene_ids,
    )
    selected_values = values[gene_mask]
    selected_labels = labels[gene_mask]
    selected_genes = [gene for gene, keep in zip(genes, gene_mask, strict=True) if keep]

    if collapse:
        group_label = _default_collapsed_label(result, selected_labels=selected_labels, label=label)
        profiles.insert_unchecked(
            x_vals.copy(),
            _aggregate_cluster_profile(selected_values, total_bins=total_bins),
            group_label,
        )
        unique_clusters = np.unique(selected_labels)
        groups.append(
            ClusterMetageneGroup(
                label=group_label,
                gene_ids=[gene.gene_id for gene in selected_genes],
                gene_count=len(selected_genes),
                cluster_id=int(unique_clusters[0]) if unique_clusters.size == 1 else None,
            )
        )
    else:
        for cluster in sorted(np.unique(selected_labels).tolist()):
            cluster_mask = selected_labels == cluster
            cluster_values = selected_values[cluster_mask]
            cluster_genes = [
                gene.gene_id for gene, keep in zip(selected_genes, cluster_mask, strict=True) if keep
            ]
            cluster_label = _cluster_label(int(cluster), cluster_source=result.cluster_source)
            profiles.insert_unchecked(
                x_vals.copy(),
                _aggregate_cluster_profile(cluster_values, total_bins=total_bins),
                cluster_label,
            )
            groups.append(
                ClusterMetageneGroup(
                    label=cluster_label,
                    gene_ids=cluster_genes,
                    gene_count=len(cluster_genes),
                    cluster_id=int(cluster),
                )
            )

    return ClusterMetageneData(
        profiles=profiles,
        segments=cluster_profile_segments(result.feature_matrix.feature_bins),
        groups=groups,
    )


@beartype
def build_gene_embedding_data(result: GeneClusterResult) -> GeneEmbeddingData:
    """
    Convert a clustering result into plot-ready PCA embedding data.
    """
    genes = result.feature_matrix.genes
    return GeneEmbeddingData(
        gene_ids=[gene.gene_id for gene in genes],
        gene_names=[gene.gene_name for gene in genes],
        chromosomes=[gene.chrom for gene in genes],
        starts=np.asarray([gene.start for gene in genes], dtype=np.int64),
        ends=np.asarray([gene.end for gene in genes], dtype=np.int64),
        strands=[gene.strand for gene in genes],
        labels=np.asarray(result.labels, dtype=np.int64),
        embedding=np.asarray(result.embedding, dtype=np.float64),
        cluster_source=str(result.cluster_source),
    )


@beartype
def build_gene_dendrogram_data(result: GeneClusterResult) -> GeneDendrogramData | None:
    """
    Convert a clustering result into plot-ready dendrogram data.
    """
    if result.linkage_matrix is None or result.leaf_order is None:
        return None
    hierarchical_meta = result.metadata.get("hierarchical", {})
    leaf_labels = hierarchical_meta.get("leaf_gene_ids")
    if not isinstance(leaf_labels, list) or len(leaf_labels) != len(result.leaf_order):
        leaf_labels = result.feature_matrix.gene_ids
    return GeneDendrogramData(
        linkage_matrix=np.asarray(result.linkage_matrix, dtype=np.float64),
        leaf_labels=leaf_labels,
        leaf_order=np.asarray(result.leaf_order, dtype=np.int64),
    )


__all__ = [
    "ClusterMetageneData",
    "ClusterMetageneGroup",
    "GeneDendrogramData",
    "GeneEmbeddingData",
    "build_cluster_metagene_data",
    "build_gene_dendrogram_data",
    "build_gene_embedding_data",
    "cluster_profile_segments",
]
