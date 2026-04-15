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
) -> dict[str, np.ndarray | dict[str, float] | str]:
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
        "timings_s": {
            "pairwise_distance": round(t_dist - t0, 6),
            "linkage": round(t_linkage - t_dist, 6),
            "cluster_cut": round(t_labels - t_linkage, 6),
            "total": round(t_labels - t0, 6),
        },
    }
