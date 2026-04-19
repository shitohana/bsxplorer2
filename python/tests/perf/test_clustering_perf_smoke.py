from __future__ import annotations

import os
from pathlib import Path
from time import perf_counter

import numpy as np
import pytest

from bsx2.clustering.backend import cluster_gene_profiles, run_pca_kmeans
from bsx2.clustering.config import (
    BackendConfig,
    ClusterConfig,
    HierarchicalConfig,
    OutputConfig,
    SilhouetteConfig,
)
from bsx2.clustering.io import build_cluster_metrics

from .conftest import SyntheticPerfScenario, make_synthetic_gene_profile_matrix


def _require_extended_perf() -> None:
    if os.environ.get("BSX2_RUN_PERF") != "1":
        pytest.skip("Set BSX2_RUN_PERF=1 to execute medium/large perf smoke scenarios")


def _policy_contract_config(scenario: SyntheticPerfScenario) -> ClusterConfig:
    return ClusterConfig(
        bsx_path=Path("synthetic.bsx"),
        annotation_path=Path("synthetic.gff"),
        backend=BackendConfig(
            n_components=2,
            n_clusters=min(6, scenario.n_clusters),
            seed=scenario.seed,
            n_init=2,
            pca_exact_max_matrix_size=500_000,
        ),
        silhouette=SilhouetteConfig(
            exact_max_genes=2_000,
            max_samples=400,
            seed=scenario.seed,
        ),
        hierarchical=HierarchicalConfig(max_genes=5_000, subsample_genes=2_000),
        output=OutputConfig(output_dir=Path(".")),
    )


@pytest.mark.perf
@pytest.mark.parametrize("scenario_name", ["small"])
def test_run_pca_kmeans_perf_smoke(
    perf_scenarios: dict[str, SyntheticPerfScenario],
    scenario_name: str,
) -> None:
    scenario = perf_scenarios[scenario_name]
    clean_scenario = SyntheticPerfScenario(
        name=scenario.name,
        n_genes=scenario.n_genes,
        n_bins=scenario.n_bins,
        n_clusters=scenario.n_clusters,
        missing_rate=0.0,
        seed=scenario.seed,
    )
    matrix = make_synthetic_gene_profile_matrix(clean_scenario)
    t0 = perf_counter()
    result = run_pca_kmeans(
        matrix,
        BackendConfig(n_components=2, n_clusters=matrix.n_features // 25, seed=7, n_init=2),
        silhouette_config=SilhouetteConfig(exact_max_genes=500, max_samples=400, seed=7),
    )
    total_s = perf_counter() - t0

    assert np.isfinite(total_s)
    assert np.isfinite(result["metadata"]["timings_s"]["pca"])
    assert np.isfinite(result["metadata"]["timings_s"]["kmeans"])
    assert np.isfinite(result["metadata"]["timings_s"]["silhouette"])
    assert result["metadata"]["silhouette_mode"] in {"exact", "sampled", "disabled"}


@pytest.mark.perf
@pytest.mark.parametrize("scenario_name", ["small", "medium", "large"])
def test_cluster_pipeline_perf_smoke(
    perf_scenarios: dict[str, SyntheticPerfScenario],
    scenario_name: str,
) -> None:
    if scenario_name != "small":
        _require_extended_perf()

    scenario = perf_scenarios[scenario_name]
    clean_scenario = SyntheticPerfScenario(
        name=scenario.name,
        n_genes=scenario.n_genes,
        n_bins=scenario.n_bins,
        n_clusters=scenario.n_clusters,
        missing_rate=0.0,
        seed=scenario.seed,
    )
    matrix = make_synthetic_gene_profile_matrix(clean_scenario)
    config = ClusterConfig(
        bsx_path=Path("synthetic.bsx"),
        annotation_path=Path("synthetic.gff"),
        backend=BackendConfig(
            n_components=2,
            n_clusters=min(6, clean_scenario.n_clusters),
            seed=clean_scenario.seed,
            n_init=2,
        ),
        silhouette=SilhouetteConfig(
            exact_max_genes=500,
            max_samples=400,
            seed=clean_scenario.seed,
        ),
        hierarchical=HierarchicalConfig(max_genes=800, subsample_genes=400),
        output=OutputConfig(output_dir=Path(".")),
    )
    result = cluster_gene_profiles(
        matrix,
        config,
    )
    metrics = build_cluster_metrics(result, config)

    assert result.embedding.shape[0] == clean_scenario.n_genes
    assert np.isfinite(result.metadata["timings_s"]["total"])
    assert metrics["performance_analysis"]["version"] == 3
    assert metrics["scalability_guidance"]["mode"] == "safe_auto_defaults"


@pytest.mark.perf
@pytest.mark.parametrize(
    ("scenario_name", "expected_silhouette_mode", "expected_hierarchical_mode", "expected_pca_solver"),
    [
        ("small", "exact", "exact", "exact"),
        ("medium", "sampled", "exact", "exact"),
        ("large", "sampled", "subsample", "truncated"),
    ],
)
def test_reference_perf_policy_contract(
    perf_scenarios: dict[str, SyntheticPerfScenario],
    scenario_name: str,
    expected_silhouette_mode: str,
    expected_hierarchical_mode: str,
    expected_pca_solver: str,
) -> None:
    if scenario_name != "small":
        _require_extended_perf()

    scenario = perf_scenarios[scenario_name]
    clean_scenario = SyntheticPerfScenario(
        name=scenario.name,
        n_genes=scenario.n_genes,
        n_bins=scenario.n_bins,
        n_clusters=scenario.n_clusters,
        missing_rate=0.0,
        seed=scenario.seed,
    )
    matrix = make_synthetic_gene_profile_matrix(clean_scenario)
    config = _policy_contract_config(clean_scenario)
    result = cluster_gene_profiles(matrix, config)
    metrics = build_cluster_metrics(result, config)

    assert result.metadata["pca_kmeans"]["pca_solver"] == expected_pca_solver
    assert result.metadata["pca_kmeans"]["silhouette_mode"] == expected_silhouette_mode
    assert result.metadata["hierarchical"]["mode"] == expected_hierarchical_mode
    assert metrics["scalability_limits"]["reference_scenarios"]["small"]["n_genes"] == 1_000
    assert metrics["perf_non_regression_contract"]["asserts_policy_modes_not_wall_clock"] is True
