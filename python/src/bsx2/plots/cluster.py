from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np

from bsx2.plots.data import DiscreteRegionData
from bsx2 import io as _io
from bsx2.plots.metagene import Segment, compute_discrete_regions


@dataclass(frozen=True)
class MethylationMatrix:
    """Prepared methylation matrix (regions x bins) with metadata."""

    matrix: np.ndarray  # shape: (n_regions, n_bins)
    region_ids: list[str]
    bins: np.ndarray  # shape: (n_bins,)


@dataclass(frozen=True)
class KMeansResult:
    labels: np.ndarray  # shape: (n_regions,)
    centroids: np.ndarray  # shape: (n_clusters, n_bins)
    inertia: float
    n_clusters: int


@dataclass(frozen=True)
class PCAResult:
    scores: np.ndarray  # shape: (n_regions, n_components)
    loadings: np.ndarray  # shape: (n_bins, n_components)
    explained_variance: np.ndarray  # shape: (n_components,)
    components: int


@dataclass(frozen=True)
class LinkageResult:
    linkage: np.ndarray  # scipy linkage matrix, shape (n_regions-1, 4)
    order: np.ndarray  # leaf order, shape (n_regions,)


def _apply_nan_policy_row(
    row: np.ndarray,
    *,
    drop_nan_rows: bool = False,
    nan_fill: float | None = 0.0,
) -> np.ndarray | None:
    """
    Apply NaN policy to a single row.
    - default: drop row if any NaN/inf;
    - if nan_fill is provided: replace NaN/inf with that value and keep the row.
    """
    row = row.astype(np.float64, copy=True)
    finite = np.isfinite(row)
    if nan_fill is not None:
        return np.where(finite, row, nan_fill)
    if drop_nan_rows and not finite.all():
        return None
    return row


def filter_constant_rows(
    mat: MethylationMatrix,
    *,
    eps: float = 1e-12,
    drop_nan_rows: bool = False,
    nan_fill: float | None = 0.0,
) -> MethylationMatrix:
    """Drop rows без сигнала/почти константные; NaN обрабатываются по флагам."""
    if mat.matrix.size == 0:
        return mat

    new_rows: list[np.ndarray] = []
    new_ids: list[str] = []
    for row, rid in zip(mat.matrix, mat.region_ids):
        row = _apply_nan_policy_row(row, drop_nan_rows=drop_nan_rows, nan_fill=nan_fill)
        if row is None:
            continue
        finite = np.isfinite(row)
        if not finite.any():
            continue
        r = row[finite]
        rng = r.max() - r.min()
        std = r.std()
        if (rng < eps) or (std < eps):
            continue
        new_rows.append(row)
        new_ids.append(rid)

    if not new_rows:
        return MethylationMatrix(matrix=np.empty((0, mat.matrix.shape[1])), region_ids=[], bins=mat.bins)
    return MethylationMatrix(matrix=np.vstack(new_rows), region_ids=new_ids, bins=mat.bins)


def prepare_matrix(
    drd: DiscreteRegionData,
    *,
    norm: Literal["none", "zscore", "minmax"] = "none",
    drop_nan_rows: bool = False,
    nan_fill: float | None = 0.0,
) -> MethylationMatrix:
    """Convert DiscreteRegionData to a dense matrix (regions x bins) with optional per-region normalization."""
    mat, labels, grid = drd._stack_to_common_grid()
    rows: list[np.ndarray] = []
    region_ids: list[str] = []
    for row, rid in zip(mat, labels):
        row = _apply_nan_policy_row(row, drop_nan_rows=drop_nan_rows, nan_fill=nan_fill)
        if row is None:
            continue
        rows.append(row)
        region_ids.append(rid)
    data = np.vstack(rows) if rows else np.empty((0, mat.shape[1]))

    if norm == "zscore":
        mean = np.nanmean(data, axis=1, keepdims=True)
        std = np.nanstd(data, axis=1, keepdims=True)
        std[std == 0.0] = 1.0
        data = (data - mean) / std
    elif norm == "minmax":
        minv = np.nanmin(data, axis=1, keepdims=True)
        maxv = np.nanmax(data, axis=1, keepdims=True)
        span = maxv - minv
        span[span == 0.0] = 1.0
        data = (data - minv) / span
    elif norm == "none":
        pass
    else:
        raise ValueError(f"Unsupported norm: {norm}")

    return MethylationMatrix(matrix=data, region_ids=region_ids, bins=grid)


def run_kmeans(
    mat: MethylationMatrix,
    n_clusters: int,
    *,
    max_iter: int = 300,
    tol: float = 1e-4,
    random_state: int | None = None,
    drop_nan_rows: bool = False,
    nan_fill: float | None = 0.0,
) -> KMeansResult:
    """Simple NumPy-based k-means (Lloyd) on prepared matrix."""
    mat = filter_constant_rows(mat, drop_nan_rows=drop_nan_rows, nan_fill=nan_fill)
    data = np.nan_to_num(mat.matrix, nan=0.0, posinf=0.0, neginf=0.0)
    n_samples, n_features = data.shape
    if n_samples == 0 or n_features == 0:
        raise ValueError("Empty matrix")
    k_eff = min(max(1, n_clusters), n_samples)
    if k_eff < 2:
        raise ValueError("Not enough samples for k-means")

    rng = np.random.default_rng(random_state)
    centroids = np.empty((k_eff, n_features), dtype=np.float64)
    centroids[0] = data[rng.integers(0, n_samples)]
    closest_dist_sq = np.full(n_samples, np.inf, dtype=np.float64)
    for c in range(1, k_eff):
        dist_sq = np.sum((data[:, None, :] - centroids[None, :c, :]) ** 2, axis=2).min(axis=1)
        closest_dist_sq = np.minimum(closest_dist_sq, dist_sq)
        safe = np.nan_to_num(closest_dist_sq, nan=0.0, posinf=0.0, neginf=0.0)
        total = safe.sum()
        probs = np.full(n_samples, 1.0 / n_samples) if total <= 0.0 else safe / total
        centroids[c] = data[rng.choice(n_samples, p=probs)]

    labels = np.zeros(n_samples, dtype=np.int32)
    for _ in range(max_iter):
        distances = np.sum((data[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
        new_labels = distances.argmin(axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for k in range(n_clusters):
            mask = labels == k
            if not np.any(mask):
                centroids[k] = data[rng.integers(0, n_samples)]
            else:
                centroids[k] = data[mask].mean(axis=0)
        shift = np.max(np.linalg.norm(centroids[:, None, :] - centroids[None, :, :], axis=2))
        if shift < tol:
            break

    inertia = float(np.sum((data - centroids[labels]) ** 2))
    return KMeansResult(labels=labels, centroids=centroids, inertia=inertia, n_clusters=k_eff)


def run_pca(
    mat: MethylationMatrix,
    n_components: int = 2,
    *,
    center: bool = True,
    scale: bool = False,
    drop_nan_rows: bool = False,
    nan_fill: float | None = 0.0,
) -> PCAResult:
    """PCA via SVD; returns scores/loadings/explained variance."""
    mat = filter_constant_rows(mat, drop_nan_rows=drop_nan_rows, nan_fill=nan_fill)
    data = np.nan_to_num(mat.matrix, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float64, copy=False)
    if data.size == 0:
        raise ValueError("Empty matrix")

    if center:
        data -= data.mean(axis=0, keepdims=True)
    if scale:
        std = data.std(axis=0, keepdims=True)
        std[std == 0.0] = 1.0
        data /= std

    n_samples, n_features = data.shape
    k = min(n_components, n_samples, n_features)
    U, S, Vt = np.linalg.svd(data, full_matrices=False)
    exp_var = (S**2) / max(n_samples - 1, 1)
    scores = U[:, :k] * S[:k]
    loadings = Vt[:k].T
    return PCAResult(
        scores=scores[:, :k],
        loadings=loadings[:, :k],
        explained_variance=exp_var[:k],
        components=k,
    )


def run_linkage(
    mat: MethylationMatrix,
    *,
    method: str = "ward",
    metric: str = "euclidean",
    drop_nan_rows: bool = False,
    nan_fill: float | None = 0.0,
) -> LinkageResult:
    """Hierarchical clustering via scipy; raises ImportError if scipy is unavailable."""
    try:
        import scipy.cluster.hierarchy as sch  # type: ignore
        from scipy.spatial.distance import pdist  # type: ignore
    except ImportError as e:  # pragma: no cover - optional dependency
        raise ImportError("scipy is required for linkage") from e

    mat = filter_constant_rows(mat, drop_nan_rows=drop_nan_rows, nan_fill=nan_fill)
    data = np.nan_to_num(mat.matrix, nan=0.0, posinf=0.0, neginf=0.0)
    if data.shape[0] < 2:
        raise ValueError("Not enough samples for linkage")
    d = pdist(data, metric=metric)
    if not np.isfinite(d).all():
        d = np.nan_to_num(d, nan=0.0, posinf=0.0, neginf=0.0)
    linkage = sch.linkage(d, method=method)
    order = sch.leaves_list(linkage)
    return LinkageResult(linkage=linkage, order=order)


def reorder_matrix(mat: MethylationMatrix, order: Sequence[int]) -> MethylationMatrix:
    """Reorder rows (regions) by given order."""
    idx = np.asarray(order, dtype=int)
    return MethylationMatrix(matrix=mat.matrix[idx], region_ids=[mat.region_ids[i] for i in idx], bins=mat.bins)


def cluster_subset(mat: MethylationMatrix, labels: np.ndarray, cluster_id: int) -> MethylationMatrix:
    """Select a subset of regions belonging to a specific cluster."""
    if labels.shape[0] != mat.matrix.shape[0]:
        raise ValueError("labels length must match number of regions")
    mask = labels == cluster_id
    return MethylationMatrix(matrix=mat.matrix[mask], region_ids=[r for r, m in zip(mat.region_ids, mask) if m], bins=mat.bins)


def metagene_for_cluster(
    reader: _io.RegionReader,
    contigs: Sequence,
    labels: np.ndarray,
    *,
    cluster_id: int,
    segments: Sequence[Segment],
    agg_method=None,
    reverse_negative: bool = True,
    use_labels: Sequence[str] | None = None,
) -> DiscreteRegionData:
    """Build DiscreteRegionData for contigs that belong to a given cluster."""
    if len(contigs) != labels.shape[0]:
        raise ValueError("contigs length must match labels length")
    contigs_sel = [c for c, lab in zip(contigs, labels) if lab == cluster_id]
    labels_sel = None
    if use_labels is not None:
        labels_sel = [lab for lab, l in zip(use_labels, labels) if l == cluster_id]
    return compute_discrete_regions(
        reader,
        contigs_sel,
        segments=segments,
        agg_method=agg_method,
        reverse_negative=reverse_negative,
        labels=labels_sel,
    )
