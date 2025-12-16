import numpy as np
import pytest

from bsx2.plots import compute_discrete_regions, segments_total_bins, Segment
from bsx2.plots.data import DiscreteRegionData


def test_reverse_negative_flips_positions_and_values():
    class DummyBatch:
        def __init__(self):
            self._xs = [0.0, 0.25, 1.0]
            self._ys = [0.1, 0.4, 0.9]

        def discretise(self, n_bins, agg):
            return self._xs, self._ys

    class DummyContig:
        strand_str = "-"

    class DummyReader:
        def iter_contigs(self, contigs):
            for _ in contigs:
                yield DummyBatch()

    contigs = [DummyContig()]
    drd = compute_discrete_regions(DummyReader(), contigs, segments=[Segment("region", 3)], agg_method=None, reverse_negative=True)
    assert len(drd) == 1
    pos = drd.positions[0]
    val = drd.densities[0]
    assert np.allclose(pos, np.array([0.0, 0.75, 1.0]))
    assert np.allclose(val, np.array([0.9, 0.4, 0.1]))


def test_segments_total_bins_arbitrary():
    segs = [Segment("left", 40), Segment("right", 60)]
    assert segments_total_bins(segs) == 100


def test_discrete_region_data_bad_density_raises():
    drd = DiscreteRegionData()
    with pytest.raises(Exception):
        drd.insert(np.array([0.0, 1.0]), np.array([-0.1, 1.1]))
