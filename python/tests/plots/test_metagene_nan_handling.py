import numpy as np
import pandas as pd

from bsx2.plots.data import DiscreteRegionData
from bsx2.plots.polars_html import line_df, dist_df


def test_line_mean_ignores_nan():
    drd = DiscreteRegionData()
    drd.insert(np.array([0.0, 0.5, 1.0]), np.array([0.1, np.nan, 0.3]))
    x, y = line_df(drd, agg="mean", as_percent=False)
    # Expect mean over finite values only: (0.1 + 0.3) / 2 = 0.2
    assert np.isclose(y[1], 0.2)


def test_dist_df_drops_nan():
    drd = DiscreteRegionData()
    drd.insert(np.array([0.0, 1.0]), np.array([np.nan, 0.5]))
    out = dist_df(drd, as_percent=False)
    # out: list of (bin, density, region)
    densities = [v for _, v, _ in out]
    assert all(np.isfinite(densities))
    assert all(np.isclose(densities, 0.5))
