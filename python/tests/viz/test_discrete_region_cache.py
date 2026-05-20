import numpy as np

from bsx2.viz.compute.discrete_region_cache import DiscreteRegionData, load_discrete_region_data, save_discrete_region_data


def test_save_load_cache(tmp_path):
    drd = DiscreteRegionData()
    drd.insert(np.array([0.0, 1.0]), np.array([0.2, 0.4]), "x")
    save_discrete_region_data(drd, tmp_path)
    loaded = load_discrete_region_data(tmp_path)
    assert len(loaded.positions) == 1
