from __future__ import annotations

import numpy as np
import pytest

bsx2 = pytest.importorskip("bsx2")
if not hasattr(bsx2, "Context"):
    pytest.skip("bsx2 extension is unavailable", allow_module_level=True)

chrline = pytest.importorskip("bsx2.viz.chrline")

Context = bsx2.Context
ChrEmptyPolicy = chrline.ChrEmptyPolicy
ChrLineStat = chrline.ChrLineStat
compute_chr_line_data = chrline.compute_chr_line_data


class FakeSeries:
    def __init__(self, values):
        self._values = np.asarray(values, dtype=object)

    def to_numpy(self):
        return self._values


class FakeBatch:
    def __init__(
        self,
        seqname: str,
        *,
        position,
        count_m,
        count_total,
        context,
    ):
        self._seqname = seqname
        self._position = FakeSeries(position)
        self._count_m = FakeSeries(count_m)
        self._count_total = FakeSeries(count_total)
        self._context = FakeSeries(context)

    def seqname(self):
        return self._seqname

    def position(self):
        return self._position

    def count_m(self):
        return self._count_m

    def count_total(self):
        return self._count_total

    def context(self):
        return self._context


class IterableReader:
    def __init__(self, batches):
        self._batches = list(batches)

    def __iter__(self):
        return iter(self._batches)


class QueryReader:
    def __init__(self, order, batches_by_chr):
        self._order = list(order)
        self._batches = dict(batches_by_chr)
        self.reset_calls = 0

    def chr_order(self):
        return list(self._order)

    def query(self, contig):
        seqname = getattr(contig, "seqname", None)
        if callable(seqname):
            seqname = seqname()
        return self._batches.get(str(seqname))

    def reset(self):
        self.reset_calls += 1


def test_compute_chr_line_data_iterable_reader_weighted_windows():
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1, 20_000, 50_001, 70_000],
                count_m=[2, 3, 4, 1],
                count_total=[4, 6, 8, 2],
                context=[True, True, True, True],
            ),
            FakeBatch(
                "chr1",
                position=[99_999],
                count_m=[5],
                count_total=[10],
                context=[True],
            ),
        ]
    )

    data = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 120_000},
    )

    assert data.chromosomes() == ["chr1"]
    assert len(data.tracks) == 1
    track = data.tracks[0]

    np.testing.assert_allclose(track.value[:2], np.array([0.5, 0.5]))
    assert np.isnan(track.value[2])
    np.testing.assert_array_equal(track.site_count, np.array([2, 3, 0]))
    np.testing.assert_allclose(
        track.x_bp,
        np.array([25_000.5, 75_000.5, 110_000.5], dtype=np.float64),
    )


def test_compute_chr_line_data_chh_context_uses_null_context_values():
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1, 2, 3, 50_001],
                count_m=[1, 9, 9, 3],
                count_total=[2, 10, 10, 6],
                context=[None, True, False, None],
            )
        ]
    )

    data = compute_chr_line_data(
        reader,
        context=Context.CHH,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 50_100},
    )
    track = data.tracks[0]

    np.testing.assert_allclose(track.value, np.array([0.5, 0.5]))
    np.testing.assert_array_equal(track.site_count, np.array([1, 1]))


def test_compute_chr_line_data_drop_policy_keeps_real_chr_offsets():
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1],
                count_m=[1],
                count_total=[2],
                context=[True],
            ),
            FakeBatch(
                "chr2",
                position=[1],
                count_m=[3],
                count_total=[6],
                context=[True],
            ),
        ]
    )

    data = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 100_000, "chr2": 50_000},
        empty_policy=ChrEmptyPolicy.DROP,
    )

    assert data.chromosomes() == ["chr1", "chr2"]
    assert len(data.tracks) == 2
    assert data.tracks[0].chr_length_bp == 100_000
    assert data.tracks[1].offset_bp == 100_000
    np.testing.assert_allclose(data.tracks[1].x_global_bp, np.array([125_000.5]))


def test_chrline_global_curve_inserts_nan_boundary_between_chromosomes():
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1],
                count_m=[1],
                count_total=[2],
                context=[True],
            ),
            FakeBatch(
                "chr2",
                position=[1],
                count_m=[3],
                count_total=[6],
                context=[True],
            ),
        ]
    )

    data = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 100_000, "chr2": 50_000},
        empty_policy=ChrEmptyPolicy.DROP,
    )

    x_vals, y_vals = data.global_curve()

    np.testing.assert_allclose(x_vals, np.array([25_000.5, 100_000.0, 125_000.5]))
    assert np.isnan(y_vals[1])
    np.testing.assert_allclose(y_vals[[0, 2]], np.array([0.5, 0.5]))


def test_compute_chr_line_data_supports_query_chr_order_reader():
    reader = QueryReader(
        ["chrA", "chrB"],
        {
            "chrA": FakeBatch(
                "chrA",
                position=[1, 55_000],
                count_m=[1, 4],
                count_total=[2, 8],
                context=[True, True],
            ),
            "chrB": FakeBatch(
                "chrB",
                position=[50_001],
                count_m=[1],
                count_total=[2],
                context=[True],
            ),
        },
    )

    data = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chrA": 60_000, "chrB": 60_000},
        empty_policy=ChrEmptyPolicy.ZERO,
    )

    assert reader.reset_calls == 1
    assert data.chromosomes() == ["chrA", "chrB"]
    np.testing.assert_allclose(data.tracks[0].value, np.array([0.5, 0.5]))
    np.testing.assert_allclose(data.tracks[1].value, np.array([0.0, 0.5]))


def test_compute_chr_line_data_rejects_chr_lengths_smaller_than_observed():
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[2_000],
                count_m=[1],
                count_total=[2],
                context=[True],
            )
        ]
    )

    with pytest.raises(ValueError, match="smaller than observed max position"):
        compute_chr_line_data(
            reader,
            context=Context.CG,
            bin_size_bp=500,
            chr_lengths={"chr1": 1_000},
        )


def test_compute_chr_line_data_supports_unweighted_mean_stat():
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1, 2],
                count_m=[1, 1],
                count_total=[2, 10],
                context=[True, True],
            )
        ]
    )

    weighted = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 50_000},
        stat=ChrLineStat.WEIGHTED_MEAN,
    )
    mean = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 50_000},
        stat="mean",
    )

    w_track = weighted.tracks[0]
    m_track = mean.tracks[0]

    assert weighted.stat is ChrLineStat.WEIGHTED_MEAN
    assert mean.stat is ChrLineStat.MEAN
    np.testing.assert_allclose(w_track.value, np.array([2.0 / 12.0]))
    np.testing.assert_allclose(m_track.value, np.array([(0.5 + 0.1) / 2.0]))


def test_compute_chr_line_data_smoothing_changes_profile():
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1, 50_001, 100_001, 150_001, 200_001],
                count_m=[0, 1, 0, 1, 0],
                count_total=[1, 1, 1, 1, 1],
                context=[True, True, True, True, True],
            )
        ]
    )

    base = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 250_000},
        stat="mean",
    )
    smoothed = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 250_000},
        stat="mean",
        smooth={"method": "savgol", "window_length": 5, "polyorder": 2},
    )

    base_vals = base.tracks[0].value
    smooth_vals = smoothed.tracks[0].value
    assert not np.allclose(base_vals, smooth_vals)
    assert np.all((0.0 <= smooth_vals) & (smooth_vals <= 1.0))
