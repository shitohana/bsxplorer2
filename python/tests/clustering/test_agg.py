from __future__ import annotations

import numpy as np
import pytest

from bsx2.clustering.agg import finalize_gene_profile_matrix
from bsx2.clustering.models import FeatureBin, GeneAnnotation


def _genes() -> list[GeneAnnotation]:
    return [
        GeneAnnotation("gene_a", "chr1", 0, 100, "+"),
        GeneAnnotation("gene_b", "chr1", 100, 200, "+"),
        GeneAnnotation("gene_c", "chr1", 200, 300, "-"),
        GeneAnnotation("gene_d", "chr1", 300, 400, "-"),
    ]


def _feature_bins() -> list[FeatureBin]:
    return [
        FeatureBin("up_1", "up", 0, 0),
        FeatureBin("body_1", "body", 0, 1),
        FeatureBin("body_2", "body", 1, 2),
        FeatureBin("down_1", "down", 0, 3),
    ]


def test_finalize_gene_profile_matrix_filters_and_imputes_retained_values() -> None:
    values = np.array(
        [
            [0.1, np.nan, 0.3, 0.2],
            [0.2, 0.5, np.nan, 0.4],
            [0.3, 0.7, 0.1, 0.6],
            [np.nan, np.nan, np.nan, np.nan],
        ],
        dtype=float,
    )
    original = values.copy()

    matrix = finalize_gene_profile_matrix(
        genes=_genes(),
        feature_bins=_feature_bins(),
        values=values,
        max_gene_missing_rate=0.6,
        max_feature_missing_rate=0.4,
        min_gene_profile_variance=0.0,
        min_feature_variance=0.0,
    )

    assert [gene.gene_id for gene in matrix.genes] == ["gene_a", "gene_b", "gene_c"]
    assert [feature_bin.feature_name for feature_bin in matrix.feature_bins] == [
        "up_1",
        "body_1",
        "body_2",
        "down_1",
    ]
    np.testing.assert_allclose(
        matrix.values,
        np.array(
            [
                [0.1, 0.6, 0.3, 0.2],
                [0.2, 0.5, 0.2, 0.4],
                [0.3, 0.7, 0.1, 0.6],
            ],
            dtype=float,
        ),
    )
    assert matrix.metadata["dropped_genes"] == 1
    assert matrix.metadata["dropped_features"] == 0
    assert np.isnan(original[0, 1])
    assert np.isnan(original[1, 2])
    assert np.isnan(values[0, 1])
    assert np.isnan(values[1, 2])


def test_finalize_gene_profile_matrix_rejects_all_features_filtered() -> None:
    values = np.array(
        [
            [0.1, np.nan, 0.3, 0.2],
            [0.2, 0.5, np.nan, 0.4],
            [0.3, 0.7, 0.1, 0.6],
            [np.nan, np.nan, np.nan, np.nan],
        ],
        dtype=float,
    )

    with pytest.raises(ValueError, match="All metagene bins were filtered out"):
        finalize_gene_profile_matrix(
            genes=_genes(),
            feature_bins=_feature_bins(),
            values=values,
            max_gene_missing_rate=0.6,
            max_feature_missing_rate=0.2,
            min_gene_profile_variance=0.0,
            min_feature_variance=10.0,
        )
