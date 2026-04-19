from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from bsx2.clustering.models import FeatureBin, GeneAnnotation, GeneProfileMatrix
from bsx2.viz.data import DiscreteRegionData


@dataclass(frozen=True)
class SyntheticPerfScenario:
    name: str
    n_genes: int
    n_bins: int
    n_clusters: int
    missing_rate: float
    seed: int


def make_synthetic_gene_profile_matrix(
    scenario: SyntheticPerfScenario,
) -> GeneProfileMatrix:
    rng = np.random.default_rng(scenario.seed)
    genes = [
        GeneAnnotation(
            gene_id=f"{scenario.name}_gene_{idx}",
            chrom="chr1",
            start=idx * 100,
            end=(idx + 1) * 100,
            strand="+" if idx % 2 == 0 else "-",
        )
        for idx in range(scenario.n_genes)
    ]
    feature_bins = [
        FeatureBin(
            feature_name=f"bin_{idx + 1}",
            segment="body",
            local_bin_index=idx,
            global_bin_index=idx,
        )
        for idx in range(scenario.n_bins)
    ]

    centers = rng.uniform(0.05, 0.95, size=(scenario.n_clusters, scenario.n_bins))
    values = np.empty((scenario.n_genes, scenario.n_bins), dtype=float)
    for idx in range(scenario.n_genes):
        cluster = idx % scenario.n_clusters
        values[idx] = np.clip(
            centers[cluster] + rng.normal(0.0, 0.03, size=scenario.n_bins),
            0.0,
            1.0,
        )

    if scenario.missing_rate > 0:
        missing_mask = rng.random(values.shape) < scenario.missing_rate
        values[missing_mask] = np.nan

    finite = np.isfinite(values)
    gene_missing_rate = np.mean(~finite, axis=1)
    feature_missing_rate = np.mean(~finite, axis=0)
    gene_variance = np.nanvar(values, axis=1)
    feature_variance = np.nanvar(values, axis=0)

    return GeneProfileMatrix(
        genes=genes,
        feature_bins=feature_bins,
        values=values,
        gene_missing_rate=gene_missing_rate.astype(float, copy=False),
        gene_variance=np.nan_to_num(gene_variance, nan=0.0).astype(float, copy=False),
        feature_missing_rate=feature_missing_rate.astype(float, copy=False),
        feature_variance=np.nan_to_num(feature_variance, nan=0.0).astype(float, copy=False),
        metadata={"scenario": scenario.name},
    )


def make_synthetic_heatmap_drd(
    scenario: SyntheticPerfScenario,
) -> DiscreteRegionData:
    rng = np.random.default_rng(scenario.seed)
    drd = DiscreteRegionData()
    x_vals = (np.arange(scenario.n_bins, dtype=np.float64) + 0.5) / float(scenario.n_bins)
    for idx in range(scenario.n_genes):
        cluster = idx % scenario.n_clusters
        base = 0.15 + 0.25 * cluster
        profile = np.clip(
            base + 0.15 * np.sin((x_vals * np.pi * (cluster + 1)) + idx * 0.03),
            0.0,
            1.0,
        )
        if scenario.missing_rate > 0:
            missing_mask = rng.random(profile.shape) < scenario.missing_rate
            profile = profile.astype(float, copy=True)
            profile[missing_mask] = np.nan
        drd.insert_unchecked(
            x_vals.copy(),
            np.asarray(profile, dtype=np.float64),
            f"{scenario.name}_gene_{idx}",
        )
    return drd


@pytest.fixture
def perf_scenarios() -> dict[str, SyntheticPerfScenario]:
    return {
        "small": SyntheticPerfScenario(
            name="small",
            n_genes=1_000,
            n_bins=100,
            n_clusters=4,
            missing_rate=0.05,
            seed=7,
        ),
        "medium": SyntheticPerfScenario(
            name="medium",
            n_genes=5_000,
            n_bins=100,
            n_clusters=5,
            missing_rate=0.05,
            seed=11,
        ),
        "large": SyntheticPerfScenario(
            name="large",
            n_genes=10_000,
            n_bins=100,
            n_clusters=6,
            missing_rate=0.05,
            seed=13,
        ),
    }
