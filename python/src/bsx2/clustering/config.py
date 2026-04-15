from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from beartype import beartype
from bsx2 import Context


class AnnotationFormat(str, Enum):
    GFF = "gff"
    GTF = "gtf"
    BED = "bed"


class NormalizationMode(str, Enum):
    NONE = "none"
    COLUMN_ZSCORE = "column_zscore"
    ROW_ZSCORE = "row_zscore"


class ClusterSource(str, Enum):
    KMEANS = "kmeans"
    HIERARCHICAL = "hierarchical"


class HierarchicalDistance(str, Enum):
    EUCLIDEAN = "euclidean"
    CORRELATION = "correlation"


class HierarchicalLinkage(str, Enum):
    AVERAGE = "average"
    COMPLETE = "complete"
    WARD = "ward"


class TableFormat(str, Enum):
    TSV = "tsv"
    CSV = "csv"
    PARQUET = "parquet"
    JSON = "json"


@beartype
@dataclass(frozen=True)
class ReadConfig:
    context: Context | None = None
    min_coverage: int = 5


@beartype
@dataclass(frozen=True)
class BlockCacheConfig:
    enabled: bool = False
    cache_dir: Path | None = None


@beartype
@dataclass(frozen=True)
class GeneProfileConfig:
    upstream_bp: int = 2_000
    downstream_bp: int = 2_000
    upstream_bins: int = 20
    body_bins: int = 50
    downstream_bins: int = 20
    min_gene_length_bp: int = 0
    min_bin_total_coverage: int = 1
    max_gene_missing_rate: float = 0.4
    max_feature_missing_rate: float = 0.95
    min_gene_profile_variance: float = 0.0
    min_feature_variance: float = 0.0
    normalization: NormalizationMode = NormalizationMode.ROW_ZSCORE
    limit_genes: int | None = None

    @property
    def total_bins(self) -> int:
        return self.upstream_bins + self.body_bins + self.downstream_bins


@beartype
@dataclass(frozen=True)
class BackendConfig:
    n_components: int = 2
    n_clusters: int = 4
    seed: int = 0
    n_init: int = 8
    max_iter: int = 200
    tol: float = 1e-4


@beartype
@dataclass(frozen=True)
class HierarchicalConfig:
    distance: HierarchicalDistance = HierarchicalDistance.CORRELATION
    linkage: HierarchicalLinkage = HierarchicalLinkage.AVERAGE


@beartype
@dataclass(frozen=True)
class OutputConfig:
    output_dir: Path
    prefix: str = ""
    table_formats: tuple[TableFormat, ...] = (TableFormat.TSV,)
    write_table_files: bool = True
    write_metrics_file: bool = True


@beartype
@dataclass(frozen=True)
class ClusterConfig:
    bsx_path: Path
    annotation_path: Path
    annotation_format: AnnotationFormat | None = None
    read: ReadConfig = ReadConfig()
    block_cache: BlockCacheConfig = BlockCacheConfig()
    gene_profile: GeneProfileConfig = GeneProfileConfig()
    backend: BackendConfig = BackendConfig()
    hierarchical: HierarchicalConfig = HierarchicalConfig()
    cluster_source: ClusterSource = ClusterSource.KMEANS
    output: OutputConfig = OutputConfig(output_dir=Path("."))
