from __future__ import annotations

import numpy as np

from bsx2.clustering.config import HierarchicalConfig, HierarchicalDistance
from bsx2.clustering.hierarchical import run_hierarchical


def test_run_hierarchical_returns_linkage_and_leaf_order() -> None:
    values = np.array(
        [
            [0.0, 0.1],
            [0.1, 0.2],
            [1.0, 1.1],
            [1.1, 1.3],
        ],
        dtype=float,
    )

    result = run_hierarchical(
        values,
        HierarchicalConfig(distance=HierarchicalDistance.EUCLIDEAN),
        n_clusters=2,
    )

    assert result["linkage_matrix"].shape == (3, 4)
    assert len(result["leaf_order"]) == 4
    assert set(result["labels"].tolist()) == {0, 1}
