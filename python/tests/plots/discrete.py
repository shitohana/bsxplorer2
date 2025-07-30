import numpy as np
import pytest
from pydantic import ValidationError

from bsx2.plots.data import DiscreteRegionData


def _set_seed():
    np.random.seed(42)


@pytest.fixture
def discrete_region_data() -> DiscreteRegionData:
    data = DiscreteRegionData()
    num_regions = np.random.randint(5, 7) # Between 5 and 6 regions

    for i in range(num_regions):
        num_points = np.random.randint(5, 21) # Between 5 and 20 points per region
        positions = np.sort(np.random.rand(num_points))
        densities = np.random.rand(num_points)

        label = f"region {i+1}"

        data.insert(positions, densities, label)

    return data


def test_discrete_region_data_validation() -> None:
    # Test case 1: Different lengths of positions and densities lists
    with pytest.raises(expected_exception=ValidationError):
        DiscreteRegionData(positions=[np.array([1.0])], densities=[])

    # Test case 2: Not all position, density pairs have equal lengths
    with pytest.raises(expected_exception=ValidationError):
        DiscreteRegionData(positions=[np.array([1.0, 2.0]), np.array([3.0])], densities=[np.array([10.0, 20.0]), np.array([30.0, 40.0])])

    # Test case 3: Lengths of labels and positions don't match
    with pytest.raises(expected_exception=ValidationError):
        DiscreteRegionData(
            positions=[np.array([1.0])],
            densities=[np.array([10.0])],
            labels=["label1", "label2"],
        )

    # Test case 4: Position array is not sorted in ascending order
    with pytest.raises(expected_exception=ValidationError):
        DiscreteRegionData(
            positions=[np.array([2.0, 1.0])],
            densities=[np.array([20.0, 10.0])],
        )


def test_insert(discrete_region_data):
    before_addition = len(discrete_region_data)
    discrete_region_data.insert(np.array([.2, .3, .5]), np.array([.2, .4, .5]))
    after_addition = len(discrete_region_data)
    assert after_addition == before_addition + 1


def test_find(discrete_region_data):
    found = discrete_region_data.find("region 1")
    assert found is not None


def test_to_line_plot(discrete_region_data):
    n_fragments = 10
    lp_data = discrete_region_data.to_line_plot(n_fragments)
    assert len(lp_data) == n_fragments

    orig_mean = np.concatenate(discrete_region_data.densities).mean()
    calc_mean = lp_data.y.mean()
    assert abs(orig_mean - calc_mean) < 0.1


def test_to_pandas(discrete_region_data):
    pd_data = discrete_region_data.to_pandas()
    assert len(pd_data) == len(discrete_region_data)