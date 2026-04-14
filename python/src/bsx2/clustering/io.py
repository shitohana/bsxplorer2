from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter

import polars as pl
import plotly.graph_objects as go
from scipy.cluster.hierarchy import dendrogram

from .config import ClusterConfig
from .metagene import cluster_metagene_summary
from .models import GeneClusterResult


def _output_path(
    output_dir: Path,
    prefix: str,
    filename: str,
) -> Path:
    stem = f"{prefix}_{filename}" if prefix else filename
    return output_dir / stem


def _write_pca_plot(result: GeneClusterResult, path: Path) -> None:
    fig = go.Figure()
    labels = result.labels
    genes = result.feature_matrix.genes
    embedding = result.embedding
    for cluster in sorted(set(labels.tolist())):
        indices = [idx for idx, label in enumerate(labels.tolist()) if label == cluster]
        fig.add_trace(
            go.Scatter(
                x=embedding[indices, 0],
                y=embedding[indices, 1] if embedding.shape[1] > 1 else [0.0] * len(indices),
                mode="markers",
                name=f"cluster {cluster}",
                text=[genes[idx].gene_id for idx in indices],
                customdata=[
                    [
                        genes[idx].gene_name or "",
                        genes[idx].chrom,
                        genes[idx].start,
                        genes[idx].end,
                        genes[idx].strand,
                    ]
                    for idx in indices
                ],
                hovertemplate=(
                    "gene=%{text}<br>"
                    "name=%{customdata[0]}<br>"
                    "chr=%{customdata[1]}:%{customdata[2]}-%{customdata[3]}<br>"
                    "strand=%{customdata[4]}<br>"
                    "PC1=%{x:.4f}<br>"
                    "PC2=%{y:.4f}<extra></extra>"
                ),
            )
        )
    fig.update_layout(
        title="Gene-level PCA",
        xaxis_title="PC1",
        yaxis_title="PC2" if embedding.shape[1] > 1 else "PC2 (not available)",
    )
    fig.write_html(path, include_plotlyjs="cdn")


def _write_dendrogram_plot(result: GeneClusterResult, path: Path) -> None:
    if result.linkage_matrix is None:
        return
    dendro = dendrogram(
        result.linkage_matrix,
        labels=result.feature_matrix.gene_ids,
        no_plot=True,
    )
    fig = go.Figure()
    for xs, ys in zip(dendro["icoord"], dendro["dcoord"], strict=True):
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(color="#1f2937", width=1),
                hoverinfo="skip",
                showlegend=False,
            )
        )
    fig.update_layout(
        title="Gene-level dendrogram",
        xaxis=dict(
            tickmode="array",
            tickvals=[5 + 10 * idx for idx in range(len(dendro["ivl"]))],
            ticktext=dendro["ivl"],
            tickangle=90,
        ),
        yaxis_title="Distance",
    )
    fig.write_html(path, include_plotlyjs="cdn")


def write_cluster_outputs(
    result: GeneClusterResult,
    config: ClusterConfig,
) -> dict[str, Path]:
    t0 = perf_counter()
    output_dir = config.output.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = config.output.prefix

    profiles_path = _output_path(output_dir, prefix, "gene_profiles.tsv")
    clusters_path = _output_path(output_dir, prefix, "gene_clusters.tsv")
    embedding_path = _output_path(output_dir, prefix, "gene_embedding.tsv")
    bins_path = _output_path(output_dir, prefix, "gene_bins.tsv")
    cluster_metagenes_path = _output_path(output_dir, prefix, "cluster_metagenes.tsv")
    linkage_path = _output_path(output_dir, prefix, "dendrogram_linkage.tsv")
    leaves_path = _output_path(output_dir, prefix, "dendrogram_leaves.tsv")
    loadings_path = _output_path(output_dir, prefix, "pca_loadings.tsv")
    pca_plot_path = _output_path(output_dir, prefix, "pca_plot.html")
    dendrogram_plot_path = _output_path(output_dir, prefix, "dendrogram.html")
    metrics_path = _output_path(output_dir, prefix, "metrics.json")

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
    pl.DataFrame(profile_data).write_csv(profiles_path, separator="\t")
    t_profiles = perf_counter()

    pl.DataFrame(
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
    ).write_csv(clusters_path, separator="\t")
    t_clusters = perf_counter()

    embedding_data: dict[str, list[float] | list[str] | list[int | None]] = {
        "gene_id": [gene.gene_id for gene in result.feature_matrix.genes],
        "cluster": result.labels.tolist(),
        "kmeans_cluster": result.kmeans_labels.tolist(),
        "hierarchical_cluster": (
            result.hierarchical_labels.tolist()
            if result.hierarchical_labels is not None
            else [None] * result.feature_matrix.n_genes
        ),
    }
    for idx in range(result.embedding.shape[1]):
        embedding_data[f"PC{idx + 1}"] = result.embedding[:, idx].tolist()
    pl.DataFrame(embedding_data).write_csv(embedding_path, separator="\t")
    t_embedding = perf_counter()

    pl.DataFrame(
        {
            "feature_name": [feature_bin.feature_name for feature_bin in result.feature_matrix.feature_bins],
            "segment": [feature_bin.segment for feature_bin in result.feature_matrix.feature_bins],
            "local_bin_index": [feature_bin.local_bin_index for feature_bin in result.feature_matrix.feature_bins],
            "global_bin_index": [feature_bin.global_bin_index for feature_bin in result.feature_matrix.feature_bins],
            "missing_rate": result.feature_matrix.feature_missing_rate.tolist(),
            "variance": result.feature_matrix.feature_variance.tolist(),
        }
    ).write_csv(bins_path, separator="\t")
    t_bins = perf_counter()

    cluster_metagene_summary(result).write_csv(cluster_metagenes_path, separator="\t")
    t_cluster_metagenes = perf_counter()

    if result.linkage_matrix is not None:
        pl.DataFrame(
            {
                "left_node": result.linkage_matrix[:, 0].astype(int).tolist(),
                "right_node": result.linkage_matrix[:, 1].astype(int).tolist(),
                "distance": result.linkage_matrix[:, 2].tolist(),
                "cluster_size": result.linkage_matrix[:, 3].astype(int).tolist(),
            }
        ).write_csv(linkage_path, separator="\t")
        pl.DataFrame(
            {
                "gene_id": [result.feature_matrix.genes[idx].gene_id for idx in result.leaf_order.tolist()],
                "dendrogram_order": list(range(result.feature_matrix.n_genes)),
                "leaf_index": result.leaf_order.tolist(),
            }
        ).write_csv(leaves_path, separator="\t")
    t_linkage = perf_counter()

    loading_data: dict[str, list[str] | list[int] | list[float]] = {
        "feature_name": [feature_bin.feature_name for feature_bin in result.feature_matrix.feature_bins],
        "segment": [feature_bin.segment for feature_bin in result.feature_matrix.feature_bins],
        "global_bin_index": [feature_bin.global_bin_index for feature_bin in result.feature_matrix.feature_bins],
    }
    for idx in range(result.components.shape[0]):
        loading_data[f"PC{idx + 1}"] = result.components[idx].tolist()
    pl.DataFrame(loading_data).write_csv(loadings_path, separator="\t")
    t_loadings = perf_counter()

    _write_pca_plot(result, pca_plot_path)
    t_pca_plot = perf_counter()
    _write_dendrogram_plot(result, dendrogram_plot_path)
    t_dendrogram_plot = perf_counter()

    metrics = {
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
        "io_timings_s": {
            "gene_profiles_tsv": round(t_profiles - t0, 6),
            "gene_clusters_tsv": round(t_clusters - t_profiles, 6),
            "gene_embedding_tsv": round(t_embedding - t_clusters, 6),
            "gene_bins_tsv": round(t_bins - t_embedding, 6),
            "cluster_metagenes_tsv": round(t_cluster_metagenes - t_bins, 6),
            "dendrogram_tables": round(t_linkage - t_cluster_metagenes, 6),
            "pca_loadings_tsv": round(t_loadings - t_linkage, 6),
            "pca_plot_html": round(t_pca_plot - t_loadings, 6),
            "dendrogram_html": round(t_dendrogram_plot - t_pca_plot, 6),
        },
    }
    t_metrics_start = perf_counter()
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    t_metrics_end = perf_counter()
    metrics["io_timings_s"]["metrics_json"] = round(t_metrics_end - t_metrics_start, 6)
    metrics["io_timings_s"]["total"] = round(t_metrics_end - t0, 6)
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    return {
        "profiles": profiles_path,
        "clusters": clusters_path,
        "embedding": embedding_path,
        "bins": bins_path,
        "cluster_metagenes": cluster_metagenes_path,
        "linkage": linkage_path,
        "leaves": leaves_path,
        "pca_loadings": loadings_path,
        "pca_plot": pca_plot_path,
        "dendrogram_plot": dendrogram_plot_path,
        "metrics": metrics_path,
    }
