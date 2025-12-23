import numpy as np

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.polars_html import (
    discrete_to_long_pl,
    line_df,
    heatmap_df,
    dist_df,
)


def _make_drd(n_regions: int = 3, n_bins: int = 10) -> DiscreteRegionData:
    drd = DiscreteRegionData()
    x = np.linspace(0, 1, n_bins)
    for i in range(n_regions):
        # keep densities within [0, 1] to satisfy beartype validator
        y = np.linspace(0, 1 - 0.01 * i, n_bins)
        drd.insert(x, y, f"r{i+1}")
    return drd


def test_discrete_to_long_pl() -> None:
    drd = _make_drd()
    arr = discrete_to_long_pl(drd)
    assert arr.shape[0] == 3 * 10


def test_line_df() -> None:
    drd = _make_drd()
    x, y = line_df(drd, agg="mean")
    assert len(x) == 10 and len(y) == 10
    assert x[0] == 0.0 and x[-1] == 1.0


def test_heatmap_df_and_dist_df() -> None:
    drd = _make_drd()
    z, regions, bins = heatmap_df(drd)
    assert z.shape == (3, 10)
    assert len(regions) == 3
    assert len(bins) == 10
    ddf = dist_df(drd)
    assert len(ddf) == 3 * 10
