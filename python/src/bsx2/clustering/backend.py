from __future__ import annotations

import math
from time import perf_counter

import numpy as np

from .agg import normalize_matrix
from .config import BackendConfig, ClusterConfig, ClusterSource
from .hierarchical import run_hierarchical
from .models import GeneClusterResult, GeneProfileMatrix


def _resolve_effective_n_components(
    requested_n_components: int,
    n_objects: int,
    n_features: int,
) -> int:
    if n_objects < 2:
        raise ValueError("PCA/KMeans requires at least 2 genes")
    if n_features < 1:
        raise ValueError("PCA input must contain at least 1 feature")

    max_components = min(n_objects - 1, n_features)
    if max_components < 1:
        raise ValueError(
            "PCA requires at least 1 effective component; check retained gene/bin counts"
        )
    return max(1, min(requested_n_components, max_components))


def _fit_pca(
    values: np.ndarray,
    n_components: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if values.ndim != 2:
        raise ValueError("PCA input must be 2D")
    n_objects, n_features = values.shape
    n_components = _resolve_effective_n_components(
        requested_n_components=n_components,
        n_objects=n_objects,
        n_features=n_features,
    )
    centered = values - values.mean(axis=0, keepdims=True)
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)

    components = vt[:n_components]
    embedding = centered @ components.T

    explained_variance = (singular_values**2) / (n_objects - 1)
    total_variance = explained_variance.sum()
    if total_variance == 0:
        ratio = np.zeros(n_components, dtype=float)
    else:
        ratio = explained_variance[:n_components] / total_variance

    return embedding, components, ratio


def _init_centroids_plus_plus(
    data: np.ndarray,
    n_clusters: int,
    rng: np.random.Generator,
) -> np.ndarray:
    n_objects = data.shape[0]
    first_idx = int(rng.integers(0, n_objects))
    centroids = [data[first_idx]]

    while len(centroids) < n_clusters:
        distances = np.min(
            np.stack([np.sum((data - centroid) ** 2, axis=1) for centroid in centroids]),
            axis=0,
        )
        total = distances.sum()
        if total <= 0:
            next_idx = int(rng.integers(0, n_objects))
        else:
            probs = distances / total
            next_idx = int(rng.choice(n_objects, p=probs))
        centroids.append(data[next_idx])

    return np.asarray(centroids, dtype=float)


def _assign_clusters(
    data: np.ndarray,
    centroids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    distances = np.sum((data[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
    labels = np.argmin(distances, axis=1)
    min_distances = distances[np.arange(data.shape[0]), labels]
    return labels, min_distances


def _run_single_kmeans(
    data: np.ndarray,
    config: BackendConfig,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    rng = np.random.default_rng(seed)
    centroids = _init_centroids_plus_plus(data, config.n_clusters, rng)

    for _ in range(config.max_iter):
        labels, _ = _assign_clusters(data, centroids)
        new_centroids = centroids.copy()

        for cluster_idx in range(config.n_clusters):
            members = data[labels == cluster_idx]
            if len(members) == 0:
                new_centroids[cluster_idx] = data[int(rng.integers(0, data.shape[0]))]
            else:
                new_centroids[cluster_idx] = members.mean(axis=0)

        shift = np.linalg.norm(new_centroids - centroids)
        centroids = new_centroids
        if shift <= config.tol:
            break

    labels, min_distances = _assign_clusters(data, centroids)
    inertia = float(min_distances.sum())
    return labels, centroids, inertia


def _fit_kmeans(
    data: np.ndarray,
    config: BackendConfig,
) -> tuple[np.ndarray, np.ndarray, float]:
    if config.n_clusters < 1:
        raise ValueError("n_clusters must be >= 1")
    if config.n_clusters >= data.shape[0]:
        raise ValueError(
            "Degenerate clustering configuration: n_clusters must be smaller than "
            "the number of genes"
        )

    best_labels: np.ndarray | None = None
    best_centroids: np.ndarray | None = None
    best_inertia = math.inf

    for init_idx in range(config.n_init):
        labels, centroids, inertia = _run_single_kmeans(
            data,
            config,
            seed=config.seed + init_idx,
        )
        if inertia < best_inertia:
            best_labels = labels
            best_centroids = centroids
            best_inertia = inertia

    assert best_labels is not None
    assert best_centroids is not None
    return best_labels, best_centroids, best_inertia


def _silhouette_score(
    data: np.ndarray,
    labels: np.ndarray,
) -> tuple[float | None, str | None]:
    unique_labels = np.unique(labels)
    if len(unique_labels) < 2:
        return None, "silhouette requires at least 2 clusters"

    _, cluster_sizes = np.unique(labels, return_counts=True)
    if np.any(cluster_sizes < 2):
        return None, "silhouette requires at least 2 genes in every cluster"

    distances = np.sqrt(
        np.sum((data[:, None, :] - data[None, :, :]) ** 2, axis=2, dtype=float)
    )
    silhouettes = np.zeros(data.shape[0], dtype=float)

    for idx in range(data.shape[0]):
        own_cluster = labels[idx]
        own_mask = labels == own_cluster
        own_mask[idx] = False

        if own_mask.any():
            a_i = distances[idx, own_mask].mean()
        else:
            silhouettes[idx] = 0.0
            continue

        b_i = math.inf
        for cluster in unique_labels:
            if cluster == own_cluster:
                continue
            cluster_mask = labels == cluster
            if cluster_mask.any():
                b_i = min(b_i, float(distances[idx, cluster_mask].mean()))

        silhouettes[idx] = (b_i - a_i) / max(a_i, b_i)

    return float(silhouettes.mean()), None


def run_pca_kmeans(
    matrix: GeneProfileMatrix,
    config: BackendConfig,
    *,
    analysis_values: np.ndarray | None = None,
) -> dict[str, object]:
    n_objects = matrix.n_genes
    n_features = matrix.n_features

    effective_n_components = _resolve_effective_n_components(
        requested_n_components=config.n_components,
        n_objects=n_objects,
        n_features=n_features,
    )

    t0 = perf_counter()
    values = matrix.values if analysis_values is None else analysis_values
    embedding, components, ratio = _fit_pca(values, effective_n_components)
    t_pca = perf_counter()
    labels, centroids, inertia = _fit_kmeans(embedding, config)
    t_kmeans = perf_counter()
    silhouette, silhouette_reason = _silhouette_score(embedding, labels)
    t_silhouette = perf_counter()
    unique_labels, cluster_sizes = np.unique(labels, return_counts=True)

    return {
        "labels": labels,
        "embedding": embedding,
        "components": components,
        "centroids": centroids,
        "explained_variance_ratio": ratio,
        "inertia": inertia,
        "silhouette_score": silhouette,
        "metadata": {
            "n_components": int(embedding.shape[1]),
            "requested_n_components": config.n_components,
            "effective_n_components": int(embedding.shape[1]),
            "n_clusters": config.n_clusters,
            "seed": config.seed,
            "n_init": config.n_init,
            "cluster_sizes": {
                int(label): int(size)
                for label, size in zip(unique_labels.tolist(), cluster_sizes.tolist(), strict=True)
            },
            "silhouette_reason": silhouette_reason,
            "timings_s": {
                "pca": round(t_pca - t0, 6),
                "kmeans": round(t_kmeans - t_pca, 6),
                "silhouette": round(t_silhouette - t_kmeans, 6),
                "total": round(t_silhouette - t0, 6),
            },
        },
    }


def cluster_gene_profiles(
    matrix: GeneProfileMatrix,
    config: ClusterConfig,
) -> GeneClusterResult:
    t0 = perf_counter()
    analysis_values = normalize_matrix(
        matrix.values,
        config.gene_profile.normalization,
    )
    t_norm = perf_counter()
    pca_kmeans = run_pca_kmeans(matrix, config.backend, analysis_values=analysis_values)
    t_kmeans = perf_counter()
    hierarchical = run_hierarchical(
        analysis_values,
        config.hierarchical,
        n_clusters=config.backend.n_clusters,
    )
    t_hier = perf_counter()

    if config.cluster_source is ClusterSource.KMEANS:
        selected_labels = pca_kmeans["labels"]
    else:
        selected_labels = hierarchical["labels"]

    return GeneClusterResult(
        feature_matrix=matrix,
        labels=np.asarray(selected_labels, dtype=np.int64),
        cluster_source=config.cluster_source.value,
        kmeans_labels=np.asarray(pca_kmeans["labels"], dtype=np.int64),
        hierarchical_labels=np.asarray(hierarchical["labels"], dtype=np.int64),
        embedding=np.asarray(pca_kmeans["embedding"], dtype=float),
        components=np.asarray(pca_kmeans["components"], dtype=float),
        centroids=np.asarray(pca_kmeans["centroids"], dtype=float),
        explained_variance_ratio=np.asarray(pca_kmeans["explained_variance_ratio"], dtype=float),
        inertia=float(pca_kmeans["inertia"]),
        silhouette_score=pca_kmeans["silhouette_score"],
        linkage_matrix=np.asarray(hierarchical["linkage_matrix"], dtype=float),
        leaf_order=np.asarray(hierarchical["leaf_order"], dtype=np.int64),
        metadata={
            "normalization_mode": config.gene_profile.normalization.value,
            "cluster_source": config.cluster_source.value,
            "pca_kmeans": pca_kmeans["metadata"],
            "hierarchical": {
                "distance": hierarchical["distance"],
                "linkage": hierarchical["linkage"],
                "timings_s": hierarchical["timings_s"],
            },
            "timings_s": {
                "normalize_matrix": round(t_norm - t0, 6),
                "pca_kmeans": round(t_kmeans - t_norm, 6),
                "hierarchical": round(t_hier - t_kmeans, 6),
                "total": round(t_hier - t0, 6),
            },
        },
    )
