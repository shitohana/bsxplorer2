from __future__ import annotations

import warnings

import numpy as np
from beartype import beartype

from .config import NormalizationMode
from .models import FeatureBin, GeneAnnotation, GeneProfileMatrix


@beartype
def coverage_weighted_ratio(
    count_m_sum: int,
    count_total_sum: int,
) -> float:
    if count_total_sum <= 0:
        return float("nan")
    return float(count_m_sum) / float(count_total_sum)


@beartype
def finalize_gene_profile_matrix(
    genes: list[GeneAnnotation],
    feature_bins: list[FeatureBin],
    values: np.ndarray,
    *,
    max_gene_missing_rate: float,
    max_feature_missing_rate: float,
    min_gene_profile_variance: float,
    min_feature_variance: float,
) -> GeneProfileMatrix:
    if values.ndim != 2:
        raise ValueError("values must be 2D")
    if values.shape != (len(genes), len(feature_bins)):
        raise ValueError("values shape does not match genes/feature bins")
    if values.shape[0] == 0:
        raise ValueError("No genes remain for clustering")
    if values.shape[1] == 0:
        raise ValueError("No metagene bins were collected")

    nan_mask = np.isnan(values)
    gene_missing_rate = nan_mask.mean(axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        gene_variance = np.nanvar(values, axis=1)
    gene_variance_finite = np.nan_to_num(gene_variance, nan=-1.0)
    keep_genes = gene_missing_rate <= max_gene_missing_rate
    keep_genes &= gene_variance_finite > max(min_gene_profile_variance, 0.0)
    if not np.any(keep_genes):
        raise ValueError(
            "All genes were filtered out by missing-rate/profile-variance constraints"
        )

    kept_genes = [gene for gene, keep in zip(genes, keep_genes, strict=True) if keep]
    kept_values = values[keep_genes]
    kept_nan_mask = nan_mask[keep_genes]
    kept_gene_missing = gene_missing_rate[keep_genes]
    kept_gene_variance = gene_variance[keep_genes]

    feature_missing_rate = kept_nan_mask.mean(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        feature_variance = np.nanvar(kept_values, axis=0)
    feature_variance_finite = np.nan_to_num(feature_variance, nan=-1.0)
    keep_features = feature_missing_rate <= max_feature_missing_rate
    if kept_values.shape[0] > 1:
        keep_features &= feature_variance_finite >= max(min_feature_variance, 0.0)
    if not np.any(keep_features):
        raise ValueError(
            "All metagene bins were filtered out by missing-rate/variance constraints"
        )

    kept_feature_bins = [
        feature_bin
        for feature_bin, keep in zip(feature_bins, keep_features, strict=True)
        if keep
    ]
    kept_values = kept_values[:, keep_features]
    kept_nan_mask = kept_nan_mask[:, keep_features]
    kept_feature_missing = feature_missing_rate[keep_features]
    kept_feature_variance = feature_variance[keep_features]

    if np.any(kept_nan_mask):
        col_means = np.nanmean(kept_values, axis=0)
        if np.isnan(col_means).any():
            raise ValueError("At least one retained metagene bin contains only missing values")
        missing_rows, missing_cols = np.nonzero(kept_nan_mask)
        kept_values[missing_rows, missing_cols] = col_means[missing_cols]

    return GeneProfileMatrix(
        genes=kept_genes,
        feature_bins=kept_feature_bins,
        values=kept_values,
        gene_missing_rate=kept_gene_missing,
        gene_variance=kept_gene_variance,
        feature_missing_rate=kept_feature_missing,
        feature_variance=kept_feature_variance,
        metadata={
            "dropped_genes": int((~keep_genes).sum()),
            "dropped_features": int((~keep_features).sum()),
            "max_gene_missing_rate": max_gene_missing_rate,
            "max_feature_missing_rate": max_feature_missing_rate,
            "min_gene_profile_variance": min_gene_profile_variance,
            "min_feature_variance": min_feature_variance,
        },
    )


@beartype
def zscore_columns(values: np.ndarray) -> np.ndarray:
    if values.ndim != 2:
        raise ValueError("values must be 2D")
    col_means = values.mean(axis=0, keepdims=True)
    col_stds = values.std(axis=0, keepdims=True)
    if np.any(col_stds <= 0):
        raise ValueError("Analysis matrix contains constant metagene bins after filtering")
    return (values - col_means) / col_stds


@beartype
def zscore_rows(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if values.ndim != 2:
        raise ValueError("values must be 2D")
    row_means = values.mean(axis=1, keepdims=True)
    row_stds = values.std(axis=1, keepdims=True)
    constant_mask = (row_stds[:, 0] <= 0)
    if np.any(constant_mask):
        raise ValueError(
            "Row-wise z-score requested but constant gene profiles remain after filtering"
        )
    return (values - row_means) / row_stds, constant_mask


@beartype
def normalize_matrix(
    values: np.ndarray,
    mode: NormalizationMode,
) -> np.ndarray:
    if mode is NormalizationMode.NONE:
        return values.copy()
    if mode is NormalizationMode.COLUMN_ZSCORE:
        return zscore_columns(values)
    if mode is NormalizationMode.ROW_ZSCORE:
        normalized, _ = zscore_rows(values)
        return normalized
    raise NotImplementedError(f"Unsupported normalization mode: {mode.value}")
