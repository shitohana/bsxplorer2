from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter

import numpy as np
import polars as pl
from beartype import beartype
from beartype.typing import Any

from bsx2.viz.clustering import (
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

_REFERENCE_PERF_SCENARIOS = {
    "small": {"n_genes": 1_000, "n_bins": 100},
    "medium": {"n_genes": 5_000, "n_bins": 100},
    "large": {"n_genes": 10_000, "n_bins": 100},
}


def _output_path(
    output_dir: Path,
    prefix: str,
    filename: str,
) -> Path:
    stem = f"{prefix}_{filename}" if prefix else filename
    return output_dir / stem


def _gene_table_base_columns(result: GeneClusterResult) -> dict[str, Any]:
    genes = result.feature_matrix.genes
    return {
        "gene_id": [gene.gene_id for gene in genes],
        "gene_name": [gene.gene_name for gene in genes],
        "chr": [gene.chrom for gene in genes],
        "start": np.asarray([gene.start for gene in genes], dtype=np.int64),
        "end": np.asarray([gene.end for gene in genes], dtype=np.int64),
        "strand": [gene.strand for gene in genes],
    }


def _feature_bin_base_columns(result: GeneClusterResult) -> dict[str, Any]:
    feature_bins = result.feature_matrix.feature_bins
    return {
        "feature_name": [feature_bin.feature_name for feature_bin in feature_bins],
        "segment": [feature_bin.segment for feature_bin in feature_bins],
        "local_bin_index": np.asarray(
            [feature_bin.local_bin_index for feature_bin in feature_bins],
            dtype=np.int64,
        ),
        "global_bin_index": np.asarray(
            [feature_bin.global_bin_index for feature_bin in feature_bins],
            dtype=np.int64,
        ),
    }


def _hierarchical_cluster_column(result: GeneClusterResult) -> np.ndarray | list[None]:
    if result.hierarchical_labels is None:
        return [None] * result.feature_matrix.n_genes
    return result.hierarchical_labels.astype(np.int64, copy=False)


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
    gene_columns = _gene_table_base_columns(result)
    feature_columns = _feature_bin_base_columns(result)
    hierarchical_cluster = _hierarchical_cluster_column(result)
    profile_data: dict[str, Any] = dict(gene_columns)
    for idx, feature_bin in enumerate(result.feature_matrix.feature_bins):
        profile_data[feature_bin.feature_name] = result.feature_matrix.values[:, idx]

    tables: dict[str, pl.DataFrame] = {
        "profiles": pl.DataFrame(profile_data),
        "clusters": pl.DataFrame(
            {
                **gene_columns,
                "cluster": result.labels.astype(np.int64, copy=False),
                "cluster_source": [result.cluster_source] * result.feature_matrix.n_genes,
                "kmeans_cluster": result.kmeans_labels.astype(np.int64, copy=False),
                "hierarchical_cluster": hierarchical_cluster,
            }
        ),
        "embedding": pl.DataFrame(
            {
                "gene_id": gene_columns["gene_id"],
                "cluster": result.labels.astype(np.int64, copy=False),
                "kmeans_cluster": result.kmeans_labels.astype(np.int64, copy=False),
                "hierarchical_cluster": hierarchical_cluster,
                **{
                    f"PC{idx + 1}": result.embedding[:, idx]
                    for idx in range(result.embedding.shape[1])
                },
            }
        ),
        "bins": pl.DataFrame(
            {
                **feature_columns,
                "missing_rate": result.feature_matrix.feature_missing_rate,
                "variance": result.feature_matrix.feature_variance,
            }
        ),
        "cluster_metagenes": cluster_metagene_summary(result),
        "pca_loadings": pl.DataFrame(
            {
                "feature_name": feature_columns["feature_name"],
                "segment": feature_columns["segment"],
                "global_bin_index": feature_columns["global_bin_index"],
                **{
                    f"PC{idx + 1}": result.components[idx]
                    for idx in range(result.components.shape[0])
                },
            }
        ),
    }

    if result.linkage_matrix is not None and result.leaf_order is not None:
        hierarchical_meta = result.metadata.get("hierarchical", {})
        leaf_gene_ids = hierarchical_meta.get("leaf_gene_ids")
        if not isinstance(leaf_gene_ids, list) or len(leaf_gene_ids) != len(result.leaf_order):
            leaf_gene_ids = [
                result.feature_matrix.genes[idx].gene_id for idx in result.leaf_order.tolist()
            ]
        tables["linkage"] = pl.DataFrame(
            {
                "left_node": result.linkage_matrix[:, 0].astype(np.int64, copy=False),
                "right_node": result.linkage_matrix[:, 1].astype(np.int64, copy=False),
                "distance": result.linkage_matrix[:, 2],
                "cluster_size": result.linkage_matrix[:, 3].astype(np.int64, copy=False),
            }
        )
        tables["leaves"] = pl.DataFrame(
            {
                "gene_id": leaf_gene_ids,
                "dendrogram_order": np.arange(len(result.leaf_order), dtype=np.int64),
                "leaf_index": result.leaf_order.astype(np.int64, copy=False),
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
    block_planning = result.feature_matrix.metadata.get("query_block_planning", {})
    pipeline_timings = result.metadata.get("pipeline_timings_s", {})
    backend_timings = result.metadata.get("timings_s", {})
    hierarchical_meta = result.metadata.get("hierarchical", {})
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
    if block_cache.get("enabled") and block_cache.get("mode") == "uncompressed":
        notes.append(
            "Persistent queried-block cache used uncompressed npz payloads to trade disk "
            "space for lower cache write/read CPU overhead."
        )
    if int(block_planning.get("merge_gap_bp", 0) or 0) > 0:
        notes.append(
            "Query-block planning merged nearby gene spans with a positive merge gap to "
            "reduce RegionReader query fragmentation."
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
    if result.metadata.get("pca_kmeans", {}).get("silhouette_mode") == "sampled":
        notes.append(
            "Silhouette scoring used stratified sampling to avoid a full quadratic distance "
            "matrix on the retained embedding."
        )
    if hierarchical_meta.get("skipped_reason"):
        notes.append(
            "Hierarchical clustering was skipped by policy to avoid quadratic distance work "
            "on a matrix that exceeded the configured threshold."
        )
    elif hierarchical_meta.get("subsampled"):
        notes.append(
            "Hierarchical clustering ran on a stratified diagnostic subsample rather than the "
            "full matrix to preserve a dendrogram without full quadratic cost."
        )

    return {
        "version": 3,
        "feature_matrix_timings_s": feature_timings,
        "pipeline_timings_s": pipeline_timings,
        "backend_timings_s": backend_timings,
        "recommendations": notes,
        "perf_harness": {
            "reference_scenarios": _REFERENCE_PERF_SCENARIOS,
            "extended_perf_env_var": "BSX2_RUN_PERF",
            "extended_perf_env_value": "1",
        },
    }


def _build_scalability_limits(config: ClusterConfig) -> dict[str, Any]:
    return {
        "reference_scenarios": _REFERENCE_PERF_SCENARIOS,
        "silhouette": {
            "enabled": config.silhouette.enabled,
            "mode": config.silhouette.mode.value,
            "exact_max_genes": config.silhouette.exact_max_genes,
            "max_samples": config.silhouette.max_samples,
        },
        "hierarchical": {
            "enabled": config.hierarchical.enabled,
            "mode": config.hierarchical.mode.value,
            "max_genes": config.hierarchical.max_genes,
            "subsample_genes": config.hierarchical.subsample_genes,
        },
        "pca": {
            "solver": config.backend.pca_solver.value,
            "exact_max_matrix_size": config.backend.pca_exact_max_matrix_size,
        },
        "query_block_cache": {
            "enabled": config.block_cache.enabled,
            "mode": config.block_cache.mode.value,
            "merge_gap_bp": config.read.query_block_merge_gap_bp,
        },
    }


def _build_scalability_guidance(result: GeneClusterResult) -> dict[str, Any]:
    n_genes = int(result.feature_matrix.n_genes)
    n_bins = int(result.feature_matrix.n_features)
    hierarchical_meta = result.metadata.get("hierarchical", {})
    pca_meta = result.metadata.get("pca_kmeans", {})

    guidance: list[str] = []
    if n_genes >= 5_000:
        guidance.append(
            "Hierarchical clustering is quadratic in gene count; keep it disabled or capped "
            "for matrices above roughly 5,000 genes unless a dendrogram is required."
        )
    elif n_genes >= 2_000:
        guidance.append(
            "Hierarchical clustering and exact silhouette remain feasible here, but they are "
            "the first backend steps likely to dominate as gene count grows."
        )
    else:
        guidance.append(
            "On small and medium matrices, selective BSX reads usually dominate runtime more "
            "than PCA or KMeans."
        )

    if n_genes * n_bins >= 1_000_000:
        guidance.append(
            "PCA can now switch to truncated SVD in auto mode; for very large matrices, "
            "tuning pca_solver or reducing retained genes/bins still helps."
        )

    if not bool(result.feature_matrix.metadata.get("query_block_cache", {}).get("enabled", False)):
        guidance.append(
            "Repeated runs benefit from `--query-block-cache` because selective RegionReader "
            "queries are still the main cost center."
        )
    elif result.feature_matrix.metadata.get("query_block_cache", {}).get("mode") == "compressed":
        guidance.append(
            "If block-cache CPU time starts to matter, `--query-block-cache-mode uncompressed` "
            "can reduce cache serialization overhead at the cost of larger cache files."
        )

    return {
        "mode": "safe_auto_defaults",
        "n_genes": n_genes,
        "n_profile_bins": n_bins,
        "hierarchical_enabled": bool(hierarchical_meta.get("enabled", True)),
        "hierarchical_mode": hierarchical_meta.get("mode"),
        "hierarchical_max_genes": hierarchical_meta.get("max_genes"),
        "hierarchical_skipped_reason": hierarchical_meta.get("skipped_reason"),
        "pca_solver": pca_meta.get("pca_solver"),
        "silhouette_mode": pca_meta.get("silhouette_mode"),
        "query_block_cache_mode": result.feature_matrix.metadata.get("query_block_cache", {}).get("mode"),
        "query_block_merge_gap_bp": result.feature_matrix.metadata.get("query_block_planning", {}).get(
            "merge_gap_bp"
        ),
        "guidance": guidance,
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
    silhouette_score = 8 if n_genes >= 5_000 else 7 if n_genes >= 2_000 else 4

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
                "block_cache_mode": block_cache.get("mode"),
                "query_block_merge_gap_bp": result.feature_matrix.metadata.get(
                    "query_block_planning", {}
                ).get("merge_gap_bp"),
            },
            "suggested_action": (
                "Reuse a reader/index across runs, persist block-level caches, consider "
                "uncompressed cache mode, or merge nearby query blocks more aggressively."
            ),
        },
        {
            "area": "Per-gene Python loop for metagene bin aggregation",
            "score_0_to_10": materialize_score,
            "location": "bsx2.clustering.gene_profile._profile_genes_in_block/_build_metagene_interval_matrices",
            "reason": (
                "Gene profiles are now aggregated in query-block batches, but interval matrices "
                "and searchsorted lookups still scale with retained genes and bins."
            ),
            "evidence": {
                "materialize_gene_profiles_s": materialize,
                "n_genes": n_genes,
                "n_profile_bins": n_bins,
                "materialization_strategy": result.feature_matrix.metadata.get(
                    "materialization_strategy"
                ),
            },
            "suggested_action": (
                "Further gains now likely require lower-level interval kernels, reduced BSX "
                "block reads, or moving more of the block aggregation path into Rust."
            ),
        },
        {
            "area": "Silhouette score pairwise distances",
            "score_0_to_10": silhouette_score,
            "location": "bsx2.clustering.backend._silhouette_score_exact",
            "reason": (
                "Exact silhouette builds a full pairwise distance matrix on the embedding, which "
                "is quadratic in the number of retained genes."
            ),
            "evidence": {
                "n_genes": n_genes,
                "silhouette_mode": result.metadata.get("pca_kmeans", {}).get("silhouette_mode"),
                "silhouette_n_samples": result.metadata.get("pca_kmeans", {}).get(
                    "silhouette_n_samples"
                ),
            },
            "suggested_action": (
                "Prefer sampled or auto silhouette on larger runs and reserve exact mode for "
                "small matrices or validation."
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
                "Keep hierarchical clustering disabled or capped for large n, subsample, or "
                "support approximate linkage for exploratory runs."
            ),
        },
        {
            "area": "PCA solver selection and matrix factorization",
            "score_0_to_10": pca_score,
            "location": "bsx2.clustering.backend._fit_pca",
            "reason": (
                "PCA still dominates backend time once the retained matrix gets large, even with "
                "automatic fallback between exact and truncated solvers."
            ),
            "evidence": {
                "n_genes": n_genes,
                "n_profile_bins": n_bins,
                "pca_solver": result.metadata.get("pca_kmeans", {}).get("pca_solver"),
                "pca_solver_reason": result.metadata.get("pca_kmeans", {}).get("pca_solver_reason"),
            },
            "suggested_action": (
                "Tune pca_solver or pca_exact_max_matrix_size if PCA starts to dominate backend "
                "time on retained matrices."
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
        "scalability_guidance": _build_scalability_guidance(result),
        "scalability_limits": _build_scalability_limits(config),
        "optimization_hotspots": _build_optimization_hotspots(result),
        "quadratic_steps_detected": {
            "silhouette_exact": result.metadata.get("pca_kmeans", {}).get("silhouette_mode") == "exact",
            "hierarchical_exact": result.metadata.get("hierarchical", {}).get("mode") == "exact",
        },
        "recommended_large_run_policy": {
            "silhouette_mode": "auto_or_sampled",
            "hierarchical_mode": "auto_or_skip",
            "notes": (
                "Use sampled silhouette and avoid exact full-matrix hierarchical clustering "
                "unless a full dendrogram is explicitly required."
            ),
        },
        "perf_non_regression_contract": {
            "small_runs_execute_in_default_test_pass": True,
            "medium_large_runs_require_env": {
                "name": "BSX2_RUN_PERF",
                "value": "1",
            },
            "asserts_policy_modes_not_wall_clock": True,
        },
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
