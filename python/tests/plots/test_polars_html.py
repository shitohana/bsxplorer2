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
    df = discrete_to_long_pl(drd)
    assert df.height == 3 * 10
    assert set(df.columns) == {"region", "bin", "x", "density"}


def test_line_df() -> None:
    drd = _make_drd()
    df = line_df(drd, agg="mean")
    assert df.height == 10
    xs = df.select("x").to_series().to_list()
    assert xs[0] == 0.0 and xs[-1] == 1.0


def test_heatmap_df_and_dist_df() -> None:
    drd = _make_drd()
    hdf, regions, bins = heatmap_df(drd)
    assert not hdf.is_empty()
    assert len(regions) == 3
    assert len(bins) == 10
    ddf = dist_df(drd)
    assert not ddf.is_empty()

