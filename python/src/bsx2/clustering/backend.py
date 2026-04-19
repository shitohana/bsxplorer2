from __future__ import annotations

import math
from dataclasses import dataclass
from time import perf_counter

import numpy as np
from beartype import beartype
from scipy.sparse.linalg import svds

from .agg import normalize_matrix
from .config import (
    BackendConfig,
    ClusterConfig,
    ClusterSource,
    HierarchicalMode,
    PcaSolver,
    SilhouetteConfig,
    SilhouetteMode,
)
from .hierarchical import run_hierarchical, run_hierarchical_subsample
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


@dataclass(frozen=True)
class _SilhouettePolicy:
    enabled: bool
    mode: str
    reason: str | None
    n_samples: int | None


@dataclass(frozen=True)
class _HierarchicalPolicy:
    mode: HierarchicalMode
    reason: str | None
    effective_n_genes: int
    diagnostic_only: bool


@dataclass(frozen=True)
class _PcaPolicy:
    solver: str
    reason: str | None


def _resolve_pca_policy(
    n_objects: int,
    n_features: int,
    n_components: int,
    config: BackendConfig,
) -> _PcaPolicy:
    matrix_size = n_objects * n_features
    truncated_supported = n_components < min(n_objects, n_features)

    if config.pca_solver is PcaSolver.EXACT:
        return _PcaPolicy(solver=PcaSolver.EXACT.value, reason="forced_exact_solver")

    if config.pca_solver is PcaSolver.TRUNCATED:
        if truncated_supported:
            return _PcaPolicy(
                solver=PcaSolver.TRUNCATED.value,
                reason="forced_truncated_solver",
            )
        return _PcaPolicy(
            solver=PcaSolver.EXACT.value,
            reason="truncated_solver_fallback_due_to_component_limit",
        )

    if config.pca_exact_max_matrix_size is None:
        return _PcaPolicy(
            solver=PcaSolver.EXACT.value,
            reason="auto_exact_without_matrix_size_cap",
        )

    if matrix_size <= config.pca_exact_max_matrix_size or not truncated_supported:
        return _PcaPolicy(
            solver=PcaSolver.EXACT.value,
            reason=(
                f"matrix_size_within_exact_limit:{config.pca_exact_max_matrix_size}"
                if truncated_supported
                else "auto_exact_due_to_component_limit"
            ),
        )

    return _PcaPolicy(
        solver=PcaSolver.TRUNCATED.value,
        reason=f"matrix_size_exceeds_exact_limit:{config.pca_exact_max_matrix_size}",
    )


def _total_variance(centered: np.ndarray, n_objects: int) -> float:
    if n_objects <= 1:
        return 0.0
    return float(np.sum(centered * centered, dtype=float) / (n_objects - 1))


def _fit_pca_exact(
    centered: np.ndarray,
    *,
    n_components: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_objects = centered.shape[0]
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    components = np.asarray(vt[:n_components], dtype=float)
    embedding = np.asarray(centered @ components.T, dtype=float)
    explained_variance = (singular_values[:n_components] ** 2) / (n_objects - 1)
    total_variance = _total_variance(centered, n_objects)
    if total_variance <= 0:
        ratio = np.zeros(n_components, dtype=float)
    else:
        ratio = explained_variance / total_variance
    return embedding, components, ratio


def _fit_pca_truncated(
    centered: np.ndarray,
    *,
    n_components: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_objects = centered.shape[0]
    u, singular_values, vt = svds(
        centered,
        k=n_components,
        return_singular_vectors=True,
    )
    order = np.argsort(singular_values)[::-1]
    singular_values = singular_values[order]
    u = u[:, order]
    vt = vt[order]

    components = np.asarray(vt, dtype=float)
    embedding = np.asarray(u * singular_values[np.newaxis, :], dtype=float)
    explained_variance = (singular_values**2) / (n_objects - 1)
    total_variance = _total_variance(centered, n_objects)
    if total_variance <= 0:
        ratio = np.zeros(n_components, dtype=float)
    else:
        ratio = explained_variance / total_variance
    return embedding, components, ratio


def _fit_pca(
    values: np.ndarray,
    config: BackendConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str, str | None]:
    if values.ndim != 2:
        raise ValueError("PCA input must be 2D")
    n_objects, n_features = values.shape
    n_components = _resolve_effective_n_components(
        requested_n_components=config.n_components,
        n_objects=n_objects,
        n_features=n_features,
    )
    centered = values - values.mean(axis=0, keepdims=True)
    policy = _resolve_pca_policy(
        n_objects=n_objects,
        n_features=n_features,
        n_components=n_components,
        config=config,
    )

    if policy.solver == PcaSolver.EXACT.value:
        embedding, components, ratio = _fit_pca_exact(centered, n_components=n_components)
        return embedding, components, ratio, policy.solver, policy.reason

    try:
        embedding, components, ratio = _fit_pca_truncated(centered, n_components=n_components)
        return embedding, components, ratio, policy.solver, policy.reason
    except Exception:
        embedding, components, ratio = _fit_pca_exact(centered, n_components=n_components)
        return (
            embedding,
            components,
            ratio,
            PcaSolver.EXACT.value,
            "truncated_solver_failed_fallback_to_exact",
        )


def _resolve_silhouette_policy(
    n_genes: int,
    config: SilhouetteConfig,
) -> _SilhouettePolicy:
    if not config.enabled:
        return _SilhouettePolicy(
            enabled=False,
            mode="disabled",
            reason="disabled_by_config",
            n_samples=None,
        )

    if config.mode is SilhouetteMode.EXACT:
        return _SilhouettePolicy(
            enabled=True,
            mode=SilhouetteMode.EXACT.value,
            reason=None,
            n_samples=None,
        )

    if config.mode is SilhouetteMode.SAMPLED:
        if config.max_samples is None:
            return _SilhouettePolicy(
                enabled=False,
                mode="disabled",
                reason="sampled_mode_requires_max_samples",
                n_samples=None,
            )
        if config.max_samples >= n_genes:
            return _SilhouettePolicy(
                enabled=True,
                mode=SilhouetteMode.EXACT.value,
                reason="sample_size_covers_all_genes",
                n_samples=n_genes,
            )
        return _SilhouettePolicy(
            enabled=True,
            mode=SilhouetteMode.SAMPLED.value,
            reason=None,
            n_samples=config.max_samples,
        )

    if config.exact_max_genes is None or n_genes <= config.exact_max_genes:
        return _SilhouettePolicy(
            enabled=True,
            mode=SilhouetteMode.EXACT.value,
            reason=(
                None
                if config.exact_max_genes is None
                else f"n_genes_within_exact_max_genes:{config.exact_max_genes}"
            ),
            n_samples=None,
        )

    if config.max_samples is None:
        return _SilhouettePolicy(
            enabled=False,
            mode="disabled",
            reason="auto_sampled_mode_requires_max_samples",
            n_samples=None,
        )

    return _SilhouettePolicy(
        enabled=True,
        mode=SilhouetteMode.SAMPLED.value,
        reason=f"n_genes_exceeds_exact_max_genes:{config.exact_max_genes}",
        n_samples=min(config.max_samples, n_genes),
    )


def _resolve_hierarchical_policy(
    n_genes: int,
    config: ClusterConfig,
) -> _HierarchicalPolicy:
    hierarchical = config.hierarchical
    if not hierarchical.enabled:
        return _HierarchicalPolicy(
            mode=HierarchicalMode.SKIP,
            reason="disabled_by_config",
            effective_n_genes=0,
            diagnostic_only=False,
        )

    if hierarchical.mode is HierarchicalMode.SKIP:
        return _HierarchicalPolicy(
            mode=HierarchicalMode.SKIP,
            reason="disabled_by_mode",
            effective_n_genes=0,
            diagnostic_only=False,
        )

    if hierarchical.mode is HierarchicalMode.EXACT:
        return _HierarchicalPolicy(
            mode=HierarchicalMode.EXACT,
            reason="forced_exact_mode",
            effective_n_genes=n_genes,
            diagnostic_only=False,
        )

    if hierarchical.mode is HierarchicalMode.SUBSAMPLE:
        if hierarchical.subsample_genes is None or hierarchical.subsample_genes >= n_genes:
            return _HierarchicalPolicy(
                mode=HierarchicalMode.EXACT,
                reason="subsample_not_needed",
                effective_n_genes=n_genes,
                diagnostic_only=False,
            )
        return _HierarchicalPolicy(
            mode=HierarchicalMode.SUBSAMPLE,
            reason=f"forced_subsample_mode:{hierarchical.subsample_genes}",
            effective_n_genes=hierarchical.subsample_genes,
            diagnostic_only=True,
        )

    if hierarchical.max_genes is None or n_genes <= hierarchical.max_genes:
        return _HierarchicalPolicy(
            mode=HierarchicalMode.EXACT,
            reason=(
                None
                if hierarchical.max_genes is None
                else f"n_genes_within_max_genes:{hierarchical.max_genes}"
            ),
            effective_n_genes=n_genes,
            diagnostic_only=False,
        )

    if hierarchical.subsample_genes is not None and hierarchical.subsample_genes >= 2:
        return _HierarchicalPolicy(
            mode=HierarchicalMode.SUBSAMPLE,
            reason=f"n_genes_exceeds_max_genes:{hierarchical.max_genes}",
            effective_n_genes=min(hierarchical.subsample_genes, n_genes),
            diagnostic_only=True,
        )

    return _HierarchicalPolicy(
        mode=HierarchicalMode.SKIP,
        reason=f"n_genes_exceeds_max_genes:{hierarchical.max_genes}",
        effective_n_genes=0,
        diagnostic_only=False,
    )


def _init_centroids_plus_plus(
    data: np.ndarray,
    n_clusters: int,
    rng: np.random.Generator,
) -> np.ndarray:
    n_objects = data.shape[0]
    first_idx = int(rng.integers(0, n_objects))
    centroids = np.empty((n_clusters, data.shape[1]), dtype=float)
    centroids[0] = data[first_idx]
    closest_sq_distances = np.sum((data - centroids[0]) ** 2, axis=1)
    closest_sq_distances[first_idx] = 0.0

    for centroid_idx in range(1, n_clusters):
        total = float(closest_sq_distances.sum())
        if total <= 0:
            next_idx = int(rng.integers(0, n_objects))
        else:
            probs = closest_sq_distances / total
            next_idx = int(rng.choice(n_objects, p=probs))
        centroids[centroid_idx] = data[next_idx]
        new_sq_distances = np.sum((data - centroids[centroid_idx]) ** 2, axis=1)
        closest_sq_distances = np.minimum(closest_sq_distances, new_sq_distances)
        closest_sq_distances[next_idx] = 0.0

    return centroids


def _assign_clusters(
    data: np.ndarray,
    centroids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    data_sq_norms = np.sum(data * data, axis=1, keepdims=True, dtype=float)
    centroid_sq_norms = np.sum(centroids * centroids, axis=1, dtype=float)
    distances = data_sq_norms + centroid_sq_norms[np.newaxis, :] - (2.0 * (data @ centroids.T))
    np.maximum(distances, 0.0, out=distances)
    labels = np.argmin(distances, axis=1)
    min_distances = distances[np.arange(data.shape[0]), labels]
    return labels, min_distances


def _recompute_centroids(
    data: np.ndarray,
    labels: np.ndarray,
    *,
    n_clusters: int,
    rng: np.random.Generator,
) -> np.ndarray:
    centroids = np.zeros((n_clusters, data.shape[1]), dtype=float)
    counts = np.bincount(labels, minlength=n_clusters)
    np.add.at(centroids, labels, data)

    non_empty = counts > 0
    if np.any(non_empty):
        centroids[non_empty] /= counts[non_empty, np.newaxis]

    empty_clusters = np.flatnonzero(~non_empty)
    if empty_clusters.size:
        random_indices = rng.integers(0, data.shape[0], size=int(empty_clusters.size))
        centroids[empty_clusters] = data[random_indices]
    return centroids


def _run_single_kmeans(
    data: np.ndarray,
    config: BackendConfig,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    rng = np.random.default_rng(seed)
    centroids = _init_centroids_plus_plus(data, config.n_clusters, rng)

    for _ in range(config.max_iter):
        labels, _ = _assign_clusters(data, centroids)
        new_centroids = _recompute_centroids(
            data,
            labels,
            n_clusters=config.n_clusters,
            rng=rng,
        )
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


def _silhouette_score_exact(
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


def _sample_stratified_indices(
    labels: np.ndarray,
    *,
    max_samples: int,
    seed: int,
) -> tuple[np.ndarray | None, str | None]:
    unique_labels, cluster_sizes = np.unique(labels, return_counts=True)
    if np.any(cluster_sizes < 2):
        return None, "silhouette requires at least 2 genes in every cluster"

    n_objects = int(labels.shape[0])
    if max_samples >= n_objects:
        return np.arange(n_objects, dtype=np.int64), None

    min_required = int(2 * unique_labels.size)
    if max_samples < min_required:
        return None, f"sampled silhouette requires at least {min_required} total samples"

    target = np.minimum(cluster_sizes.astype(np.int64, copy=False), 2)
    remaining = int(max_samples - target.sum())
    capacity = cluster_sizes.astype(np.int64, copy=False) - target

    while remaining > 0 and np.any(capacity > 0):
        active = capacity > 0
        weights = capacity[active].astype(np.float64)
        weights /= weights.sum()
        extra = np.zeros_like(target)
        extra_active = np.floor(weights * remaining).astype(np.int64)
        if extra_active.sum() == 0:
            extra[np.flatnonzero(active)[0]] = 1
        else:
            extra[np.flatnonzero(active)] = extra_active
        extra = np.minimum(extra, capacity)
        added = int(extra.sum())
        if added <= 0:
            break
        target += extra
        capacity -= extra
        remaining -= added

    if remaining > 0 and np.any(capacity > 0):
        for idx in np.flatnonzero(capacity > 0):
            if remaining <= 0:
                break
            target[idx] += 1
            remaining -= 1

    rng = np.random.default_rng(seed)
    sampled: list[np.ndarray] = []
    for cluster_id, take in zip(unique_labels.tolist(), target.tolist(), strict=True):
        cluster_indices = np.flatnonzero(labels == cluster_id)
        if take >= cluster_indices.size:
            picked = cluster_indices
        else:
            picked = np.sort(rng.choice(cluster_indices, size=int(take), replace=False))
        sampled.append(picked.astype(np.int64, copy=False))
    return np.sort(np.concatenate(sampled)), None


def _silhouette_score_sampled(
    data: np.ndarray,
    labels: np.ndarray,
    *,
    max_samples: int,
    seed: int,
) -> tuple[float | None, str | None, int]:
    sample_indices, sample_reason = _sample_stratified_indices(
        labels,
        max_samples=max_samples,
        seed=seed,
    )
    if sample_indices is None:
        return None, sample_reason, 0

    score, reason = _silhouette_score_exact(data[sample_indices], labels[sample_indices])
    if reason is not None:
        return None, reason, int(sample_indices.size)
    return float(score), sample_reason, int(sample_indices.size)


@beartype
def run_pca_kmeans(
    matrix: GeneProfileMatrix,
    config: BackendConfig,
    *,
    analysis_values: np.ndarray | None = None,
    silhouette_config: SilhouetteConfig | None = None,
) -> dict[str, object]:
    t0 = perf_counter()
    values = matrix.values if analysis_values is None else analysis_values
    embedding, components, ratio, pca_solver, pca_reason = _fit_pca(values, config)
    t_pca = perf_counter()
    labels, centroids, inertia = _fit_kmeans(embedding, config)
    t_kmeans = perf_counter()
    silhouette_policy = _resolve_silhouette_policy(
        embedding.shape[0],
        SilhouetteConfig(seed=config.seed) if silhouette_config is None else silhouette_config,
    )
    if not silhouette_policy.enabled:
        silhouette = None
        silhouette_reason = silhouette_policy.reason
        silhouette_n_samples = 0
    elif silhouette_policy.mode == SilhouetteMode.EXACT.value:
        silhouette, compute_reason = _silhouette_score_exact(embedding, labels)
        silhouette_reason = compute_reason or silhouette_policy.reason
        silhouette_n_samples = int(embedding.shape[0])
    else:
        assert silhouette_policy.n_samples is not None
        silhouette, compute_reason, silhouette_n_samples = _silhouette_score_sampled(
            embedding,
            labels,
            max_samples=silhouette_policy.n_samples,
            seed=(silhouette_config.seed if silhouette_config is not None else config.seed),
        )
        silhouette_reason = compute_reason or silhouette_policy.reason
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
            "pca_solver": pca_solver,
            "pca_solver_reason": pca_reason,
            "pca_exact_max_matrix_size": config.pca_exact_max_matrix_size,
            "cluster_sizes": {
                int(label): int(size)
                for label, size in zip(unique_labels.tolist(), cluster_sizes.tolist(), strict=True)
            },
            "silhouette_mode": silhouette_policy.mode,
            "silhouette_n_samples": int(silhouette_n_samples),
            "silhouette_exact_max_genes": (
                None if silhouette_config is None else silhouette_config.exact_max_genes
            ),
            "silhouette_reason": silhouette_reason,
            "timings_s": {
                "pca": round(t_pca - t0, 6),
                "kmeans": round(t_kmeans - t_pca, 6),
                "silhouette": round(t_silhouette - t_kmeans, 6),
                "total": round(t_silhouette - t0, 6),
            },
        },
    }


@beartype
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
    pca_kmeans = run_pca_kmeans(
        matrix,
        config.backend,
        analysis_values=analysis_values,
        silhouette_config=config.silhouette,
    )
    t_kmeans = perf_counter()
    hierarchical_policy = _resolve_hierarchical_policy(matrix.n_genes, config)

    if (
        config.cluster_source is ClusterSource.HIERARCHICAL
        and hierarchical_policy.mode is not HierarchicalMode.EXACT
    ):
        raise ValueError(
            "Hierarchical cluster_source requested, but hierarchical clustering was "
            f"not exact ({hierarchical_policy.reason}). Increase hierarchical.max_genes, "
            "switch to hierarchical.mode=exact, or use kmeans cluster_source."
        )

    if hierarchical_policy.mode is HierarchicalMode.EXACT:
        hierarchical = run_hierarchical(
            analysis_values,
            config.hierarchical,
            n_clusters=config.backend.n_clusters,
        )
    elif hierarchical_policy.mode is HierarchicalMode.SUBSAMPLE:
        hierarchical = run_hierarchical_subsample(
            analysis_values,
            config.hierarchical,
            n_clusters=config.backend.n_clusters,
            labels_hint=np.asarray(pca_kmeans["labels"], dtype=np.int64),
            seed=config.backend.seed,
        )
    else:
        hierarchical = {
            "labels": None,
            "linkage_matrix": None,
            "leaf_order": None,
            "distance": config.hierarchical.distance.value,
            "linkage": config.hierarchical.linkage.value,
            "timings_s": {},
            "subsampled": False,
            "effective_n_genes": 0,
        }
    t_hier = perf_counter()

    if config.cluster_source is ClusterSource.KMEANS:
        selected_labels = pca_kmeans["labels"]
    else:
        selected_labels = hierarchical["labels"]

    leaf_gene_ids = None
    if hierarchical["leaf_order"] is not None:
        leaf_gene_ids = [
            matrix.genes[int(idx)].gene_id for idx in np.asarray(hierarchical["leaf_order"]).tolist()
        ]

    return GeneClusterResult(
        feature_matrix=matrix,
        labels=np.asarray(selected_labels, dtype=np.int64),
        cluster_source=config.cluster_source.value,
        kmeans_labels=np.asarray(pca_kmeans["labels"], dtype=np.int64),
        hierarchical_labels=(
            None
            if hierarchical["labels"] is None
            else np.asarray(hierarchical["labels"], dtype=np.int64)
        ),
        embedding=np.asarray(pca_kmeans["embedding"], dtype=float),
        components=np.asarray(pca_kmeans["components"], dtype=float),
        centroids=np.asarray(pca_kmeans["centroids"], dtype=float),
        explained_variance_ratio=np.asarray(pca_kmeans["explained_variance_ratio"], dtype=float),
        inertia=float(pca_kmeans["inertia"]),
        silhouette_score=pca_kmeans["silhouette_score"],
        linkage_matrix=(
            None
            if hierarchical["linkage_matrix"] is None
            else np.asarray(hierarchical["linkage_matrix"], dtype=float)
        ),
        leaf_order=(
            None
            if hierarchical["leaf_order"] is None
            else np.asarray(hierarchical["leaf_order"], dtype=np.int64)
        ),
        metadata={
            "normalization_mode": config.gene_profile.normalization.value,
            "cluster_source": config.cluster_source.value,
            "pca_kmeans": pca_kmeans["metadata"],
            "hierarchical": {
                "enabled": hierarchical_policy.mode is not HierarchicalMode.SKIP,
                "mode": hierarchical_policy.mode.value,
                "max_genes": config.hierarchical.max_genes,
                "subsample_genes": config.hierarchical.subsample_genes,
                "diagnostic_only": hierarchical_policy.diagnostic_only,
                "effective_n_genes": hierarchical.get(
                    "effective_n_genes",
                    hierarchical_policy.effective_n_genes,
                ),
                "subsampled": bool(hierarchical.get("subsampled", False)),
                "skipped_reason": (
                    hierarchical_policy.reason
                    if hierarchical_policy.mode is HierarchicalMode.SKIP
                    else None
                ),
                "policy_reason": hierarchical_policy.reason,
                "distance": hierarchical["distance"],
                "linkage": hierarchical["linkage"],
                "leaf_gene_ids": leaf_gene_ids,
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
