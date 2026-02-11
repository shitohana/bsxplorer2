import numpy as np

from bsx2.plots import metagene as mg
from bsx2.plots.data import DiscreteRegionData


def _make_drd(label: str, xs, ys, weights=None) -> DiscreteRegionData:
    drd = DiscreteRegionData()
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    if weights is None:
        drd.insert(x, y, label=label)
    else:
        w = np.asarray(weights, dtype=float)
        drd.insert(x, y, label=label, weights=w)
    return drd


def test_combine_parts_geometry_constant_bins():
    segments = [mg.Segment("up", 2), mg.Segment("body", 2), mg.Segment("down", 2)]
    parts = ["upstream_gene", "gene", "downstream_gene"]
    xs = [0.25, 0.75]  # two points -> two bins per segment

    drd_map = {
        "upstream_gene": _make_drd("gene1", xs, [0.1, 0.1]),
        "gene": _make_drd("gene1", xs, [0.5, 0.5]),
        "downstream_gene": _make_drd("gene1", xs, [0.9, 0.9]),
    }
    combined = mg.combine_parts_drd(drd_map, segments=segments, parts_order=parts)
    _, y = mg._line_from_points(combined, segments)

    assert np.allclose(y[:2], 0.1)
    assert np.allclose(y[2:4], 0.5)
    assert np.allclose(y[4:], 0.9)


def test_line_from_points_weighted_mean():
    segments = [mg.Segment("region", 1)]
    drd = _make_drd("gene1", [0.1, 0.2], [1.0, 0.0], weights=[1.0, 9.0])
    _, y = mg._line_from_points(drd, segments, value_mode="weighted")
    assert np.allclose(y, [0.1])


class _DummySeries:
    def __init__(self, values):
        self._values = list(values)

    def to_list(self):
        return list(self._values)


class _DummyBatch:
    def __init__(self, positions, densities, weights):
        self._positions = positions
        self._densities = densities
        self._weights = weights

    def position(self):
        return _DummySeries(self._positions)

    def density(self):
        return _DummySeries(self._densities)

    def count_total(self):
        return _DummySeries(self._weights)


class _DummyReader:
    def __init__(self, batches):
        self._batches = list(batches)

    def iter_contigs(self, contigs):
        for batch in self._batches:
            yield batch


class _DummyContig:
    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self._strand = strand

    def strand_str(self):
        return self._strand


def test_reverse_negative_strand_flips_profile():
    positions = [0, 1, 2]
    densities = [0.1, 0.2, 0.3]
    weights = [1, 1, 1]

    reader_plus = _DummyReader([_DummyBatch(positions, densities, weights)])
    contig_plus = _DummyContig(0, 10, "+")
    drd_plus = mg.compute_discrete_regions(
        reader_plus,
        [contig_plus],
        segments=[mg.Segment("region", 3)],
        reverse_negative=True,
        mode="raw",
        x_mode="relative",
    )

    reader_minus = _DummyReader([_DummyBatch(positions, densities, weights)])
    contig_minus = _DummyContig(0, 10, "-")
    drd_minus = mg.compute_discrete_regions(
        reader_minus,
        [contig_minus],
        segments=[mg.Segment("region", 3)],
        reverse_negative=True,
        mode="raw",
        x_mode="relative",
    )

    assert np.allclose(drd_plus.positions[0], [0.0, 0.1, 0.2])
    assert np.allclose(drd_minus.positions[0], [0.8, 0.9, 1.0])
    assert np.allclose(drd_minus.densities[0], [0.3, 0.2, 0.1])
