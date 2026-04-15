from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter
from beartype.typing import Any

import polars as pl
from beartype import beartype
from bsx2.plots.clustering import (
    GeneDendrogramPlotComposer,
    GeneEmbeddingPlotComposer,
    build_cluster_metagene_data,
    build_cluster_metagene_plot,
    build_gene_dendrogram_data,
    build_gene_embedding_data,
)

from .config import ClusterConfig, TableFormat
from .metagene import cluster_metagene_summary
from .models import ClusterArtifacts, GeneClusterResult


_TABLE_BASENAMES = {
    "profiles": "gene_profiles",
    "clusters": "gene_clusters",
    "embedding": "gene_embedding",
    "bins": "gene_bins",
    "cluster_metagenes": "cluster_metagenes",
    "linkage": "dendrogram_linkage",
    "leaves": "dendrogram_leaves",
    "pca_loadings": "pca_loadings",
}

def _output_path(
    output_dir: Path,
    prefix: str,
    filename: str,
) -> Path:
    stem = f"{prefix}_{filename}" if prefix else filename
    return output_dir / stem


@beartype
def build_pca_plot(result: GeneClusterResult):
    return GeneEmbeddingPlotComposer().add_data(build_gene_embedding_data(result)).finish()


@beartype
def build_dendrogram_plot(result: GeneClusterResult):
    dendrogram_data = build_gene_dendrogram_data(result)
    if dendrogram_data is None:
        return None
    return GeneDendrogramPlotComposer().add_data(dendrogram_data).finish()


@beartype
def build_cluster_plot_data(result: GeneClusterResult) -> dict[str, Any]:
    plot_data: dict[str, Any] = {
        "pca": build_gene_embedding_data(result),
        "cluster_metagene": build_cluster_metagene_data(result),
    }
    dendrogram_data = build_gene_dendrogram_data(result)
    if dendrogram_data is not None:
        plot_data["dendrogram"] = dendrogram_data
    return plot_data


@beartype
def build_cluster_tables(result: GeneClusterResult) -> dict[str, pl.DataFrame]:
    profile_data: dict[str, list[str] | list[int] | list[float | None]] = {
        "gene_id": [gene.gene_id for gene in result.feature_matrix.genes],
        "gene_name": [gene.gene_name for gene in result.feature_matrix.genes],
        "chr": [gene.chrom for gene in result.feature_matrix.genes],
        "start": [gene.start for gene in result.feature_matrix.genes],
        "end": [gene.end for gene in result.feature_matrix.genes],
        "strand": [gene.strand for gene in result.feature_matrix.genes],
    }
    for idx, feature_bin in enumerate(result.feature_matrix.feature_bins):
        profile_data[feature_bin.feature_name] = result.feature_matrix.values[:, idx].tolist()

    tables: dict[str, pl.DataFrame] = {
        "profiles": pl.DataFrame(profile_data),
        "clusters": pl.DataFrame(
            {
                "gene_id": [gene.gene_id for gene in result.feature_matrix.genes],
                "gene_name": [gene.gene_name for gene in result.feature_matrix.genes],
                "chr": [gene.chrom for gene in result.feature_matrix.genes],
                "start": [gene.start for gene in result.feature_matrix.genes],
                "end": [gene.end for gene in result.feature_matrix.genes],
                "strand": [gene.strand for gene in result.feature_matrix.genes],
                "cluster": result.labels.tolist(),
                "cluster_source": [result.cluster_source] * result.feature_matrix.n_genes,
                "kmeans_cluster": result.kmeans_labels.tolist(),
                "hierarchical_cluster": (
                    result.hierarchical_labels.tolist()
                    if result.hierarchical_labels is not None
                    else [None] * result.feature_matrix.n_genes
                ),
            }
        ),
        "embedding": pl.DataFrame(
            {
                "gene_id": [gene.gene_id for gene in result.feature_matrix.genes],
                "cluster": result.labels.tolist(),
                "kmeans_cluster": result.kmeans_labels.tolist(),
                "hierarchical_cluster": (
                    result.hierarchical_labels.tolist()
                    if result.hierarchical_labels is not None
                    else [None] * result.feature_matrix.n_genes
                ),
                **{
                    f"PC{idx + 1}": result.embedding[:, idx].tolist()
                    for idx in range(result.embedding.shape[1])
                },
            }
        ),
        "bins": pl.DataFrame(
            {
                "feature_name": [feature_bin.feature_name for feature_bin in result.feature_matrix.feature_bins],
                "segment": [feature_bin.segment for feature_bin in result.feature_matrix.feature_bins],
                "local_bin_index": [feature_bin.local_bin_index for feature_bin in result.feature_matrix.feature_bins],
                "global_bin_index": [feature_bin.global_bin_index for feature_bin in result.feature_matrix.feature_bins],
                "missing_rate": result.feature_matrix.feature_missing_rate.tolist(),
                "variance": result.feature_matrix.feature_variance.tolist(),
            }
        ),
        "cluster_metagenes": cluster_metagene_summary(result),
        "pca_loadings": pl.DataFrame(
            {
                "feature_name": [feature_bin.feature_name for feature_bin in result.feature_matrix.feature_bins],
                "segment": [feature_bin.segment for feature_bin in result.feature_matrix.feature_bins],
                "global_bin_index": [feature_bin.global_bin_index for feature_bin in result.feature_matrix.feature_bins],
                **{
                    f"PC{idx + 1}": result.components[idx].tolist()
                    for idx in range(result.components.shape[0])
                },
            }
        ),
    }

    if result.linkage_matrix is not None and result.leaf_order is not None:
        tables["linkage"] = pl.DataFrame(
            {
                "left_node": result.linkage_matrix[:, 0].astype(int).tolist(),
                "right_node": result.linkage_matrix[:, 1].astype(int).tolist(),
                "distance": result.linkage_matrix[:, 2].tolist(),
                "cluster_size": result.linkage_matrix[:, 3].astype(int).tolist(),
            }
        )
        tables["leaves"] = pl.DataFrame(
            {
                "gene_id": [result.feature_matrix.genes[idx].gene_id for idx in result.leaf_order.tolist()],
                "dendrogram_order": list(range(result.feature_matrix.n_genes)),
                "leaf_index": result.leaf_order.tolist(),
            }
        )

    return tables


@beartype
def build_cluster_plots(
    result: GeneClusterResult,
    plot_data: dict[str, Any] | None = None,
) -> dict[str, Any]:
    data = build_cluster_plot_data(result) if plot_data is None else plot_data
    plots: dict[str, Any] = {
        "pca": GeneEmbeddingPlotComposer().add_data(data["pca"]).finish(),
        "cluster_metagene": build_cluster_metagene_plot(data["cluster_metagene"]),
    }
    dendrogram_data = data.get("dendrogram")
    if dendrogram_data is not None:
        plots["dendrogram"] = GeneDendrogramPlotComposer().add_data(dendrogram_data).finish()
    return plots


def _build_efficiency_report(result: GeneClusterResult) -> dict[str, Any]:
    feature_timings = result.feature_matrix.metadata.get("timings_s", {})
    reader_cache = result.feature_matrix.metadata.get("region_reader_cache", {})
    block_cache = result.feature_matrix.metadata.get("query_block_cache", {})
    pipeline_timings = result.metadata.get("pipeline_timings_s", {})
    backend_timings = result.metadata.get("timings_s", {})
    notes: list[str] = []

    feature_total = float(feature_timings.get("total", 0.0) or 0.0)
    init_reader = float(feature_timings.get("init_region_reader", 0.0) or 0.0)
    load_bsx = float(
        feature_timings.get("query_region_blocks", feature_timings.get("load_bsx_arrays", 0.0)) or 0.0
    )
    read_cost = init_reader + load_bsx
    materialize = float(feature_timings.get("materialize_gene_profiles", 0.0) or 0.0)
    cache_hit = bool(reader_cache.get("cache_hit", False))
    block_cache_hits = int(block_cache.get("hits", 0) or 0)

    if feature_total > 0 and read_cost / feature_total >= 0.7:
        notes.append(
            "The dominant cost is RegionReader initialization plus selective BSX reads; "
            "gene-profile assembly is comparatively cheap."
        )
    if cache_hit:
        notes.append(
            "The process-local RegionReader cache was hit; repeated runs avoided rebuilding "
            "the BSX reader/index in Python."
        )
    if block_cache_hits > 0:
        notes.append(
            "Persistent queried-block cache hits avoided repeated RegionReader.query calls "
            "for previously materialized blocks."
        )
    if materialize > 0 and read_cost > 0 and materialize / read_cost < 0.2:
        notes.append(
            "Per-gene metagene aggregation is not the main bottleneck; selective BSX reads "
            "or persistent caching would yield larger gains than optimizing PCA/KMeans."
        )
    if pipeline_timings and float(pipeline_timings.get("build_gene_profile_matrix", 0.0) or 0.0) > float(
        pipeline_timings.get("cluster_gene_profiles", 0.0) or 0.0
    ):
        notes.append(
            "The pipeline is I/O-bound before clustering; optimization should focus on "
            "BSX loading, chromosome filtering, or cache reuse."
        )
    if backend_timings:
        notes.append(
            "PCA, KMeans, hierarchical clustering, and plot composition are already cheap "
            "relative to BSX loading on current matrix sizes."
        )

    return {
        "feature_matrix_timings_s": feature_timings,
        "pipeline_timings_s": pipeline_timings,
        "backend_timings_s": backend_timings,
        "recommendations": notes,
    }


def _build_optimization_hotspots(result: GeneClusterResult) -> list[dict[str, Any]]:
    timings = result.feature_matrix.metadata.get("timings_s", {})
    reader_cache = result.feature_matrix.metadata.get("region_reader_cache", {})
    block_cache = result.feature_matrix.metadata.get("query_block_cache", {})
    total_feature_time = float(timings.get("total", 0.0) or 0.0)
    init_reader = float(timings.get("init_region_reader", 0.0) or 0.0)
    load_bsx = float(
        timings.get("query_region_blocks", timings.get("load_bsx_arrays", 0.0)) or 0.0
    )
    read_cost = init_reader + load_bsx
    materialize = float(timings.get("materialize_gene_profiles", 0.0) or 0.0)
    finalize = float(timings.get("finalize_matrix", 0.0) or 0.0)
    n_genes = int(result.feature_matrix.n_genes)
    n_bins = int(result.feature_matrix.n_features)

    load_score = 10 if total_feature_time <= 0 or read_cost / max(total_feature_time, 1e-9) >= 0.5 else 8
    materialize_score = 8 if n_genes >= 5_000 else 7 if n_genes >= 1_000 else 5
    hierarchical_score = 8 if n_genes >= 5_000 else 7 if n_genes >= 1_000 else 4
    pca_score = 6 if (n_genes * n_bins) >= 1_000_000 else 4
    finalize_score = 6 if finalize > 0.2 * max(total_feature_time, 1e-9) else 4

    hotspots = [
        {
            "area": "Selective BSX region reads and block materialization",
            "score_0_to_10": load_score,
            "location": "bsx2.clustering.gene_profile._make_region_reader/_query_block_arrays",
            "reason": (
                "Even after switching to RegionReader, reading and materializing requested "
                "blocks still dominates runtime compared with clustering itself."
            ),
            "evidence": {
                "init_region_reader_s": init_reader,
                "query_region_blocks_s": load_bsx,
                "feature_matrix_total_s": total_feature_time,
                "reader_cache_hit": bool(reader_cache.get("cache_hit", False)),
                "block_cache_hits": int(block_cache.get("hits", 0) or 0),
                "block_cache_writes": int(block_cache.get("writes", 0) or 0),
            },
            "suggested_action": (
                "Reuse a reader/index across runs, persist chromosome/block-level caches, or "
                "merge nearby query blocks more aggressively."
            ),
        },
        {
            "area": "Per-gene Python loop for metagene bin aggregation",
            "score_0_to_10": materialize_score,
            "location": "bsx2.clustering.gene_profile._profile_gene/build_metagene_bins",
            "reason": (
                "Every gene rebuilds its bins and performs Python-level searchsorted-based "
                "aggregation bin by bin. This scales linearly with genes and bins."
            ),
            "evidence": {
                "materialize_gene_profiles_s": materialize,
                "n_genes": n_genes,
                "n_profile_bins": n_bins,
            },
            "suggested_action": (
                "Cache per-gene bin boundaries, vectorize interval lookup by chromosome, or move "
                "the inner loop into Rust/Python extension code."
            ),
        },
        {
            "area": "Hierarchical clustering O(n^2) distance matrix",
            "score_0_to_10": hierarchical_score,
            "location": "bsx2.clustering.hierarchical.run_hierarchical",
            "reason": (
                "Pairwise distances and linkage become quadratic in the number of genes and can "
                "overtake the rest of the backend on larger cohorts."
            ),
            "evidence": {
                "n_genes": n_genes,
                "distance_metric": result.metadata.get("hierarchical", {}).get("distance"),
            },
            "suggested_action": (
                "Gate hierarchical clustering for large n, subsample, or support approximate / "
                "user-optional hierarchical mode."
            ),
        },
        {
            "area": "Dense SVD for PCA",
            "score_0_to_10": pca_score,
            "location": "bsx2.clustering.backend._fit_pca",
            "reason": (
                "Full dense SVD is simple and correct, but it does unnecessary work when only a "
                "small number of principal components is requested."
            ),
            "evidence": {
                "n_genes": n_genes,
                "n_profile_bins": n_bins,
            },
            "suggested_action": (
                "Switch to truncated/randomized PCA once matrix sizes grow enough to justify it."
            ),
        },
        {
            "area": "Matrix finalization scans and copies",
            "score_0_to_10": finalize_score,
            "location": "bsx2.clustering.agg.finalize_gene_profile_matrix",
            "reason": (
                "The matrix is scanned multiple times for NaN rates, variances, filtering, "
                "imputation, and copying. This is moderate today but grows with matrix size."
            ),
            "evidence": {
                "finalize_matrix_s": finalize,
                "n_genes": n_genes,
                "n_profile_bins": n_bins,
            },
            "suggested_action": (
                "Fuse missing-rate/variance passes where possible and reduce intermediate copies."
            ),
        },
    ]
    return hotspots


@beartype
def build_cluster_metrics(
    result: GeneClusterResult,
    config: ClusterConfig,
) -> dict[str, Any]:
    return {
        "n_genes": result.feature_matrix.n_genes,
        "n_profile_bins": result.feature_matrix.n_features,
        "cluster_source": result.cluster_source,
        "normalization_mode": config.gene_profile.normalization.value,
        "n_components": int(result.embedding.shape[1]),
        "n_clusters": int(len(set(result.labels.tolist()))),
        "inertia": result.inertia,
        "silhouette_score": result.silhouette_score,
        "explained_variance_ratio": result.explained_variance_ratio.tolist(),
        "feature_matrix": result.feature_matrix.metadata,
        "backend": result.metadata,
        "performance_analysis": _build_efficiency_report(result),
        "optimization_hotspots": _build_optimization_hotspots(result),
        "plots_are_serialized_by_library": False,
    }


@beartype
def build_cluster_artifacts(
    result: GeneClusterResult,
    config: ClusterConfig,
) -> ClusterArtifacts:
    plot_data = build_cluster_plot_data(result)
    return ClusterArtifacts(
        tables=build_cluster_tables(result),
        plots=build_cluster_plots(result, plot_data=plot_data),
        metrics=build_cluster_metrics(result, config),
        plot_data=plot_data,
    )


def _write_dataframe(
    df: pl.DataFrame,
    path: Path,
    table_format: TableFormat,
) -> None:
    if table_format is TableFormat.TSV:
        df.write_csv(path, separator="\t")
        return
    if table_format is TableFormat.CSV:
        df.write_csv(path)
        return
    if table_format is TableFormat.PARQUET:
        df.write_parquet(path)
        return
    if table_format is TableFormat.JSON:
        path.write_text(json.dumps(df.to_dicts(), indent=2), encoding="utf-8")
        return
    raise NotImplementedError(f"Unsupported table format: {table_format.value}")

@beartype
def write_cluster_outputs(
    result: GeneClusterResult,
    config: ClusterConfig,
) -> ClusterArtifacts:
    t0 = perf_counter()
    output_dir = config.output.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = config.output.prefix

    artifacts = build_cluster_artifacts(result, config)
    io_timings: dict[str, float] = {}
    paths: dict[str, Path] = {}

    if config.output.write_table_files:
        for table_name, df in artifacts.tables.items():
            basename = _TABLE_BASENAMES[table_name]
            for table_format in config.output.table_formats:
                path = _output_path(output_dir, prefix, f"{basename}.{table_format.value}")
                t_write = perf_counter()
                _write_dataframe(df, path, table_format)
                io_timings[f"{basename}_{table_format.value}"] = round(perf_counter() - t_write, 6)
                paths[f"{table_name}.{table_format.value}"] = path

    if config.output.write_metrics_file:
        metrics_path = _output_path(output_dir, prefix, "metrics.json")
        t_metrics_start = perf_counter()
        artifacts.metrics["io_timings_s"] = io_timings
        metrics_path.write_text(json.dumps(artifacts.metrics, indent=2), encoding="utf-8")
        io_timings["metrics_json"] = round(perf_counter() - t_metrics_start, 6)
        io_timings["total"] = round(perf_counter() - t0, 6)
        artifacts.metrics["io_timings_s"] = io_timings
        metrics_path.write_text(json.dumps(artifacts.metrics, indent=2), encoding="utf-8")
        paths["metrics.json"] = metrics_path
    else:
        io_timings["total"] = round(perf_counter() - t0, 6)
        artifacts.metrics["io_timings_s"] = io_timings

    artifacts.paths = paths
    return artifacts
