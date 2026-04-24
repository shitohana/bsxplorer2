from __future__ import annotations

import json
import tempfile
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from hashlib import sha256
from pathlib import Path
from typing import Any, Literal

from bsx2 import Context, Contig, HcAnnotStore, RegionReader, Strand
from bsx2.clustering import (
    AnnotationFormat,
    BackendConfig,
    BlockCacheConfig,
    ClusterConfig,
    ClusterSource,
    GeneProfileConfig,
    HierarchicalConfig,
    OutputConfig,
    ReadConfig,
    build_gene_profile_matrix,
    cluster_gene_profiles,
)
from bsx2.clustering.models import GeneClusterResult, GeneProfileMatrix

from .clustering import (
    ClusterMetageneData,
    GeneDendrogramData,
    GeneEmbeddingData,
    build_cluster_metagene_data,
    build_gene_dendrogram_data,
    build_gene_embedding_data,
)
from .data import DiscreteRegionData
from .distribution import (
    DistributionPlotData,
    build_box_distribution_data,
    build_violin_distribution_data,
)
from .metagene import (
    AnnotProfileLayout,
    MetageneProfileSegment,
    collect_layout_parts_from_hcannot,
)
from .metagene_layouts import build_annotation_metagene, build_manual_metagene

MetageneAssemblyMode = Literal["annotation-driven", "manual-composed"]
ManualRegionSpec = tuple[str, str, int, int, str]


@dataclass(frozen=True)
class PreparedMetageneFamily:
    family_key: str
    drd: DiscreteRegionData
    segments: tuple[MetageneProfileSegment, ...]
    n_regions: int
    box_segments: DistributionPlotData
    box_segments_percent: DistributionPlotData
    violin_segments: DistributionPlotData
    violin_segments_percent: DistributionPlotData


@dataclass(frozen=True)
class PreparedClusterMatrixWorkspace:
    family_key: str
    config: ClusterConfig
    feature_matrix: GeneProfileMatrix
    n_genes: int
    n_features: int


@dataclass(frozen=True)
class PreparedClusterFamily:
    family_key: str
    matrix_key: str
    feature_matrix: GeneProfileMatrix
    result: GeneClusterResult
    embedding_data: GeneEmbeddingData
    cluster_metagene_data: ClusterMetageneData
    n_genes: int
    n_features: int


@dataclass(frozen=True)
class PreparedClusterDendrogramFamily:
    family_key: str
    matrix_key: str
    result: GeneClusterResult
    dendrogram_data: GeneDendrogramData | None
    n_genes: int


def _parse_context(value: Context | str) -> Context:
    if isinstance(value, Context):
        return value
    try:
        return getattr(Context, str(value).upper())
    except AttributeError as exc:
        raise ValueError(f"Unsupported methylation context: {value}") from exc


def _new_reader(bsx_path: str | Path, context: Context | str) -> RegionReader:
    reader = RegionReader(str(bsx_path))
    reader.clear_filters()
    reader.filter_context(_parse_context(context))
    return reader


def _context_token(value: Context | str) -> str:
    context = _parse_context(value)
    return getattr(context, "value", getattr(context, "name", str(context)))


def _normalise_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if is_dataclass(value):
        return {field.name: _normalise_value(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Enum):
        return value.value
    if type(value).__name__ == "Context":
        return _context_token(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _normalise_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_normalise_value(item) for item in value]
    return value


def _stable_payload_digest(payload: dict[str, Any]) -> str:
    return sha256(
        json.dumps(
            payload,
            sort_keys=True,
            ensure_ascii=False,
            separators=(",", ":"),
            default=str,
        ).encode("utf-8")
    ).hexdigest()[:24]


def _layout_payload(layout: AnnotProfileLayout) -> list[dict[str, Any]]:
    payload: list[dict[str, Any]] = []
    for part in layout.parts:
        payload.append(
            {
                "name": part.name,
                "n_bins": part.n_bins,
                "source": part.source,
                "flank_bp": part.flank_bp,
            }
        )
    return payload


def _manual_specs_payload(
    manual_part_specs: dict[str, list[ManualRegionSpec]] | None,
) -> dict[str, list[list[Any]]] | None:
    if not manual_part_specs:
        return None
    payload: dict[str, list[list[Any]]] = {}
    for part_name, rows in manual_part_specs.items():
        payload[part_name] = [
            [label, seqname, int(start), int(end), strand]
            for label, seqname, start, end, strand in rows
        ]
    return payload


def _strand_from_symbol(symbol: str) -> Strand:
    cleaned = str(symbol).strip()
    if cleaned == "+":
        return Strand.Forward
    if cleaned == "-":
        return Strand.Reverse
    if cleaned in {".", "*", "?", "null", "none"}:
        return Strand.Null
    raise ValueError(f"Unsupported strand symbol: {symbol!r}")


def _manual_part_map_from_specs(
    manual_part_specs: dict[str, list[ManualRegionSpec]],
    *,
    layout: AnnotProfileLayout,
) -> dict[str, tuple[list[object], list[str]]]:
    part_map: dict[str, tuple[list[object], list[str]]] = {}
    for part in layout.parts:
        rows = manual_part_specs.get(part.name, [])
        contigs: list[object] = []
        labels: list[str] = []
        for label, seqname, start, end, strand_symbol in rows:
            contigs.append(Contig(str(seqname), int(start), int(end), _strand_from_symbol(strand_symbol)))
            labels.append(str(label))
        part_map[part.name] = (contigs, labels)
    return part_map


def metagene_family_key(
    *,
    bsx_path: str | Path,
    annot_path: str | Path,
    context: Context | str,
    assembly_mode: MetageneAssemblyMode,
    layout: AnnotProfileLayout,
    limit_regions: int | None,
    manual_part_specs: dict[str, list[ManualRegionSpec]] | None = None,
) -> str:
    return _stable_payload_digest(
        {
            "family": "metagene",
            "bsx_path": str(bsx_path),
            "annot_path": str(annot_path),
            "context": _context_token(context),
            "assembly_mode": assembly_mode,
            "layout": _layout_payload(layout),
            "limit_regions": limit_regions,
            "manual_part_specs": _manual_specs_payload(manual_part_specs),
        }
    )


def build_studio_cluster_config(
    *,
    bsx_path: str | Path,
    annot_path: str | Path,
    context: Context | str,
    limit_genes: int | None,
    min_coverage: int,
    query_block_merge_gap_bp: int,
    n_clusters: int,
    seed: int,
    hierarchical_enabled: bool,
    hierarchical_max_genes: int,
) -> ClusterConfig:
    return ClusterConfig(
        bsx_path=Path(bsx_path),
        annotation_path=Path(annot_path),
        annotation_format=AnnotationFormat.GFF,
        read=ReadConfig(
            context=_parse_context(context),
            min_coverage=min_coverage,
            query_block_merge_gap_bp=query_block_merge_gap_bp,
        ),
        block_cache=BlockCacheConfig(enabled=False),
        gene_profile=GeneProfileConfig(limit_genes=limit_genes),
        backend=BackendConfig(n_components=2, n_clusters=n_clusters, seed=seed),
        hierarchical=HierarchicalConfig(
            enabled=hierarchical_enabled,
            max_genes=hierarchical_max_genes,
        ),
        cluster_source=ClusterSource.KMEANS,
        output=OutputConfig(
            output_dir=Path(tempfile.gettempdir()),
            write_table_files=False,
            write_metrics_file=False,
        ),
    )


def cluster_matrix_family_key(config: ClusterConfig) -> str:
    return _stable_payload_digest(
        {
            "family": "cluster-matrix",
            "bsx_path": str(config.bsx_path),
            "annotation_path": str(config.annotation_path),
            "annotation_format": _normalise_value(config.annotation_format),
            "read": _normalise_value(config.read),
            "block_cache": _normalise_value(config.block_cache),
            "gene_profile": _normalise_value(config.gene_profile),
        }
    )


def cluster_plot_family_key(matrix_key: str, config: ClusterConfig) -> str:
    return _stable_payload_digest(
        {
            "family": "cluster-plot",
            "matrix_key": matrix_key,
            "backend": _normalise_value(config.backend),
            "cluster_source": _normalise_value(config.cluster_source),
        }
    )


def cluster_dendrogram_family_key(matrix_key: str, config: ClusterConfig) -> str:
    return _stable_payload_digest(
        {
            "family": "cluster-dendrogram",
            "matrix_key": matrix_key,
            "backend": _normalise_value(config.backend),
            "cluster_source": _normalise_value(config.cluster_source),
            "hierarchical": _normalise_value(config.hierarchical),
        }
    )


def prepare_metagene_family(
    *,
    bsx_path: str | Path,
    annot_path: str | Path,
    context: Context | str,
    assembly_mode: MetageneAssemblyMode,
    layout: AnnotProfileLayout,
    limit_regions: int | None,
    manual_part_specs: dict[str, list[ManualRegionSpec]] | None = None,
) -> PreparedMetageneFamily:
    reader = _new_reader(bsx_path, context)
    if manual_part_specs:
        drd = build_manual_metagene(
            reader,
            part_map=_manual_part_map_from_specs(manual_part_specs, layout=layout),
            layout=layout,
        )
    elif assembly_mode == "annotation-driven":
        annot = HcAnnotStore.from_gff(str(annot_path))
        drd = build_annotation_metagene(
            reader,
            annot,
            layout=layout,
            limit=limit_regions,
        )
    elif assembly_mode == "manual-composed":
        annot = HcAnnotStore.from_gff(str(annot_path))
        part_map = collect_layout_parts_from_hcannot(annot, layout=layout, limit=limit_regions)
        drd = build_manual_metagene(reader, part_map=part_map, layout=layout)
    else:
        raise ValueError(f"Unsupported metagene assembly_mode: {assembly_mode}")

    segments = tuple(layout.segments)
    return PreparedMetageneFamily(
        family_key=metagene_family_key(
            bsx_path=bsx_path,
            annot_path=annot_path,
            context=context,
            assembly_mode=assembly_mode,
            layout=layout,
            limit_regions=limit_regions,
            manual_part_specs=manual_part_specs,
        ),
        drd=drd,
        segments=segments,
        n_regions=len(drd.positions),
        box_segments=build_box_distribution_data(
            drd,
            segments=list(segments),
            per_region=False,
            distribution_mode="segments",
            as_percent=False,
        ),
        box_segments_percent=build_box_distribution_data(
            drd,
            segments=list(segments),
            per_region=False,
            distribution_mode="segments",
            as_percent=True,
        ),
        violin_segments=build_violin_distribution_data(
            drd,
            segments=list(segments),
            per_region=False,
            distribution_mode="segments",
            as_percent=False,
        ),
        violin_segments_percent=build_violin_distribution_data(
            drd,
            segments=list(segments),
            per_region=False,
            distribution_mode="segments",
            as_percent=True,
        ),
    )


def prepare_cluster_matrix_workspace(config: ClusterConfig) -> PreparedClusterMatrixWorkspace:
    feature_matrix = build_gene_profile_matrix(config)
    return PreparedClusterMatrixWorkspace(
        family_key=cluster_matrix_family_key(config),
        config=config,
        feature_matrix=feature_matrix,
        n_genes=feature_matrix.n_genes,
        n_features=feature_matrix.n_features,
    )


def prepare_cluster_family(
    matrix_workspace: PreparedClusterMatrixWorkspace,
    config: ClusterConfig,
) -> PreparedClusterFamily:
    result = cluster_gene_profiles(matrix_workspace.feature_matrix, config)
    return PreparedClusterFamily(
        family_key=cluster_plot_family_key(matrix_workspace.family_key, config),
        matrix_key=matrix_workspace.family_key,
        feature_matrix=matrix_workspace.feature_matrix,
        result=result,
        embedding_data=build_gene_embedding_data(result),
        cluster_metagene_data=build_cluster_metagene_data(result),
        n_genes=matrix_workspace.n_genes,
        n_features=matrix_workspace.n_features,
    )


def prepare_cluster_dendrogram_family(
    matrix_workspace: PreparedClusterMatrixWorkspace,
    config: ClusterConfig,
) -> PreparedClusterDendrogramFamily:
    result = cluster_gene_profiles(matrix_workspace.feature_matrix, config)
    return PreparedClusterDendrogramFamily(
        family_key=cluster_dendrogram_family_key(matrix_workspace.family_key, config),
        matrix_key=matrix_workspace.family_key,
        result=result,
        dendrogram_data=build_gene_dendrogram_data(result),
        n_genes=matrix_workspace.n_genes,
    )


__all__ = [
    "MetageneAssemblyMode",
    "ManualRegionSpec",
    "PreparedMetageneFamily",
    "PreparedClusterMatrixWorkspace",
    "PreparedClusterFamily",
    "PreparedClusterDendrogramFamily",
    "build_studio_cluster_config",
    "cluster_dendrogram_family_key",
    "cluster_matrix_family_key",
    "cluster_plot_family_key",
    "metagene_family_key",
    "prepare_cluster_dendrogram_family",
    "prepare_cluster_family",
    "prepare_cluster_matrix_workspace",
    "prepare_metagene_family",
]
