from __future__ import annotations

import math

import numpy as np
import polars as pl

from .models import GeneClusterResult


def cluster_metagene_summary(
    result: GeneClusterResult,
) -> pl.DataFrame:
    rows: list[dict[str, int | float | str]] = []
    values = result.feature_matrix.values
    labels = result.labels

    for cluster in sorted(np.unique(labels).tolist()):
        mask = labels == cluster
        cluster_values = values[mask]
        n_genes = int(mask.sum())
        means = cluster_values.mean(axis=0)
        stds = cluster_values.std(axis=0, ddof=0)
        sems = np.divide(
            stds,
            math.sqrt(n_genes),
            out=np.zeros_like(stds, dtype=float),
            where=n_genes > 0,
        )
        for feature_bin, mean, std, sem in zip(
            result.feature_matrix.feature_bins,
            means.tolist(),
            stds.tolist(),
            sems.tolist(),
            strict=True,
        ):
            rows.append(
                {
                    "cluster": int(cluster),
                    "feature_name": feature_bin.feature_name,
                    "segment": feature_bin.segment,
                    "local_bin_index": int(feature_bin.local_bin_index),
                    "global_bin_index": int(feature_bin.global_bin_index),
                    "mean": float(mean),
                    "std": float(std),
                    "sem": float(sem),
                    "n_genes": n_genes,
                }
            )

    return pl.DataFrame(rows)
