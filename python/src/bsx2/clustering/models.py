from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from beartype.typing import Any

import numpy as np
from beartype import beartype


@beartype
@dataclass(frozen=True, order=True)
class GeneAnnotation:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    gene_name: str | None = None

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError("GeneAnnotation.start must be >= 0")
        if self.end <= self.start:
            raise ValueError("GeneAnnotation.end must be greater than start")
        if self.strand not in {"+", "-"}:
            raise ValueError("GeneAnnotation.strand must be '+' or '-'")

    @property
    def length_bp(self) -> int:
        return self.end - self.start


@beartype
@dataclass(frozen=True)
class FeatureBin:
    feature_name: str
    segment: str
    local_bin_index: int
    global_bin_index: int


@beartype
@dataclass(frozen=True)
class ProfileBin:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    segment: str
    local_bin_index: int
    global_bin_index: int
    feature_name: str

    def __post_init__(self) -> None:
        if self.end < self.start:
            raise ValueError("ProfileBin.end must be >= start")


@beartype
@dataclass
class GeneProfileMatrix:
    genes: list[GeneAnnotation]
    feature_bins: list[FeatureBin]
    values: np.ndarray
    gene_missing_rate: np.ndarray
    gene_variance: np.ndarray
    feature_missing_rate: np.ndarray
    feature_variance: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.values.ndim != 2:
            raise ValueError("GeneProfileMatrix.values must be 2D")
        if self.values.shape != (len(self.genes), len(self.feature_bins)):
            raise ValueError(
                "GeneProfileMatrix.values shape does not match genes/feature_bins lengths"
            )
        if self.gene_missing_rate.shape != (len(self.genes),):
            raise ValueError("gene_missing_rate shape does not match genes")
        if self.gene_variance.shape != (len(self.genes),):
            raise ValueError("gene_variance shape does not match genes")
        if self.feature_missing_rate.shape != (len(self.feature_bins),):
            raise ValueError("feature_missing_rate shape does not match feature_bins")
        if self.feature_variance.shape != (len(self.feature_bins),):
            raise ValueError("feature_variance shape does not match feature_bins")

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_features(self) -> int:
        return len(self.feature_bins)

    @property
    def gene_ids(self) -> list[str]:
        return [gene.gene_id for gene in self.genes]


@beartype
@dataclass
class GeneClusterResult:
    feature_matrix: GeneProfileMatrix
    labels: np.ndarray
    cluster_source: str
    kmeans_labels: np.ndarray
    hierarchical_labels: np.ndarray | None
    embedding: np.ndarray
    components: np.ndarray
    centroids: np.ndarray
    explained_variance_ratio: np.ndarray
    inertia: float
    silhouette_score: float | None
    linkage_matrix: np.ndarray | None = None
    leaf_order: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n_genes = self.feature_matrix.n_genes
        if self.labels.shape != (n_genes,):
            raise ValueError("labels shape does not match feature matrix")
        if self.kmeans_labels.shape != (n_genes,):
            raise ValueError("kmeans_labels shape does not match feature matrix")
        if self.hierarchical_labels is not None and self.hierarchical_labels.shape != (n_genes,):
            raise ValueError("hierarchical_labels shape does not match feature matrix")
        if self.embedding.shape[0] != n_genes:
            raise ValueError("embedding row count does not match feature matrix")


@beartype
@dataclass
class ClusterArtifacts:
    tables: dict[str, Any]
    plots: dict[str, Any]
    metrics: dict[str, Any]
    plot_data: dict[str, Any] = field(default_factory=dict)
    paths: dict[str, Path] = field(default_factory=dict)
