from __future__ import annotations

from time import perf_counter

import numpy as np
from beartype import beartype
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist

from .config import HierarchicalConfig, HierarchicalDistance, HierarchicalLinkage


@beartype
def run_hierarchical(
    values: np.ndarray,
    config: HierarchicalConfig,
    *,
    n_clusters: int,
) -> dict[str, object]:
    if values.ndim != 2:
        raise ValueError("Hierarchical input must be 2D")
    if values.shape[0] < 2:
        raise ValueError("Hierarchical clustering requires at least 2 genes")
    if n_clusters < 1 or n_clusters > values.shape[0]:
        raise ValueError("Hierarchical clustering requires 1 <= n_clusters <= n_genes")
    if config.linkage is HierarchicalLinkage.WARD and config.distance is not HierarchicalDistance.EUCLIDEAN:
        raise ValueError("Ward linkage requires euclidean distance")

    t0 = perf_counter()
    condensed = pdist(values, metric=config.distance.value)
    t_dist = perf_counter()
    if condensed.size == 0:
        raise ValueError("Hierarchical clustering requires at least 2 genes")
    if not np.all(np.isfinite(condensed)):
        raise ValueError(
            "Hierarchical distance computation produced non-finite values; "
            "check constant/empty gene profiles or switch distance metric"
        )

    linkage_matrix = linkage(condensed, method=config.linkage.value)
    t_linkage = perf_counter()
    leaf_order = leaves_list(linkage_matrix).astype(np.int64, copy=False)
    labels = (
        fcluster(linkage_matrix, t=n_clusters, criterion="maxclust").astype(np.int64) - 1
    )
    t_labels = perf_counter()

    return {
        "labels": labels,
        "linkage_matrix": linkage_matrix,
        "leaf_order": leaf_order,
        "distance": config.distance.value,
        "linkage": config.linkage.value,
        "subsampled": False,
        "effective_n_genes": int(values.shape[0]),
        "timings_s": {
            "pairwise_distance": round(t_dist - t0, 6),
            "linkage": round(t_linkage - t_dist, 6),
            "cluster_cut": round(t_labels - t_linkage, 6),
            "total": round(t_labels - t0, 6),
        },
    }


def _subsample_indices(
    n_genes: int,
    *,
    max_genes: int,
    labels_hint: np.ndarray | None,
    seed: int,
) -> np.ndarray:
    if max_genes >= n_genes:
        return np.arange(n_genes, dtype=np.int64)

    rng = np.random.default_rng(seed)
    if labels_hint is None or labels_hint.shape != (n_genes,):
        return np.sort(rng.choice(n_genes, size=max_genes, replace=False)).astype(np.int64)

    unique_labels, cluster_sizes = np.unique(labels_hint, return_counts=True)
    if unique_labels.size >= max_genes:
        chosen_clusters = rng.choice(unique_labels, size=max_genes, replace=False)
        picks = [
            int(rng.choice(np.flatnonzero(labels_hint == cluster_label), size=1, replace=False)[0])
            for cluster_label in chosen_clusters.tolist()
        ]
        return np.sort(np.asarray(picks, dtype=np.int64))

    target = np.ones(unique_labels.size, dtype=np.int64)
    remaining = int(max_genes - target.sum())
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

    sampled: list[np.ndarray] = []
    for cluster_label, take in zip(unique_labels.tolist(), target.tolist(), strict=True):
        cluster_indices = np.flatnonzero(labels_hint == cluster_label)
        if take >= cluster_indices.size:
            picked = cluster_indices
        else:
            picked = np.sort(rng.choice(cluster_indices, size=int(take), replace=False))
        sampled.append(picked.astype(np.int64, copy=False))
    return np.sort(np.concatenate(sampled))


@beartype
def run_hierarchical_subsample(
    values: np.ndarray,
    config: HierarchicalConfig,
    *,
    n_clusters: int,
    labels_hint: np.ndarray | None = None,
    seed: int = 0,
) -> dict[str, object]:
    if values.ndim != 2:
        raise ValueError("Hierarchical input must be 2D")
    if values.shape[0] < 2:
        raise ValueError("Hierarchical clustering requires at least 2 genes")
    if config.subsample_genes is None or config.subsample_genes < 2:
        raise ValueError("Hierarchical subsample mode requires subsample_genes >= 2")

    t0 = perf_counter()
    sample_indices = _subsample_indices(
        values.shape[0],
        max_genes=min(config.subsample_genes, values.shape[0]),
        labels_hint=labels_hint,
        seed=seed,
    )
    t_sample = perf_counter()
    result = run_hierarchical(
        values[sample_indices],
        config,
        n_clusters=min(n_clusters, int(sample_indices.size)),
    )
    timings = dict(result["timings_s"])
    timings["subsample_select"] = round(t_sample - t0, 6)
    timings["total"] = round((t_sample - t0) + float(result["timings_s"]["total"]), 6)
    leaf_order_local = np.asarray(result["leaf_order"], dtype=np.int64)

    return {
        "labels": None,
        "linkage_matrix": result["linkage_matrix"],
        "leaf_order": sample_indices[leaf_order_local].astype(np.int64, copy=False),
        "distance": result["distance"],
        "linkage": result["linkage"],
        "subsampled": True,
        "effective_n_genes": int(sample_indices.size),
        "timings_s": timings,
    }
