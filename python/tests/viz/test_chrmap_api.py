from __future__ import annotations

import numpy as np
import pytest

bsx2 = pytest.importorskip("bsx2")
if not hasattr(bsx2, "Context"):
    pytest.skip("bsx2 extension is unavailable", allow_module_level=True)

viz = pytest.importorskip("bsx2.viz")
chrmap = pytest.importorskip("bsx2.viz.chrmap")

Context = bsx2.Context


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


def _is_holoviews_object(obj) -> bool:
    return obj.__class__.__module__.startswith("holoviews")


def _reader() -> IterableReader:
    return IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1, 25_000, 75_001],
                count_m=[1, 3, 4],
                count_total=[2, 6, 8],
                context=[True, True, True],
            ),
            FakeBatch(
                "chr2",
                position=[10_000, 90_000],
                count_m=[2, 1],
                count_total=[4, 2],
                context=[True, True],
            ),
        ]
    )


def test_chrmap_module_exposes_chrline_aliases() -> None:
    assert chrmap.ChromosomeMethylationTrack is chrmap.ChrLineTrack
    assert chrmap.ChromosomeMethylationMapData is chrmap.ChrLineData
    assert chrmap.ChromosomeMethylationMapComposer is chrmap.ChrLinePlotComposer


def test_compute_chromosome_methylation_map_data_returns_chrline_data() -> None:
    data = chrmap.compute_chromosome_methylation_map_data(
        _reader(),
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 100_000, "chr2": 100_000},
    )

    assert isinstance(data, chrmap.ChrLineData)
    assert data.chromosomes() == ["chr1", "chr2"]
    np.testing.assert_allclose(data.tracks[0].value, np.array([0.5, 0.5]))
    np.testing.assert_allclose(data.tracks[1].value, np.array([0.5, 0.5]))


def test_build_chromosome_methylation_map_returns_holoviews_object() -> None:
    data = chrmap.compute_chromosome_methylation_map_data(
        _reader(),
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 100_000, "chr2": 100_000},
    )

    plot = chrmap.build_chromosome_methylation_map(
        data,
        name="sample",
    )

    assert _is_holoviews_object(plot)


def test_one_step_chromosome_methylation_map_returns_holoviews_object() -> None:
    plot = chrmap.chromosome_methylation_map(
        _reader(),
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 100_000, "chr2": 100_000},
        name="sample",
    )

    assert _is_holoviews_object(plot)


def test_public_plots_surface_exports_chrmap_api() -> None:
    assert viz.ChromosomeMethylationTrack is chrmap.ChromosomeMethylationTrack
    assert viz.ChromosomeMethylationMapData is chrmap.ChromosomeMethylationMapData
    assert viz.ChromosomeMethylationMapComposer is chrmap.ChromosomeMethylationMapComposer
    assert (
        viz.compute_chromosome_methylation_map_data
        is chrmap.compute_chromosome_methylation_map_data
    )
    assert viz.build_chromosome_methylation_map is chrmap.build_chromosome_methylation_map
    assert viz.chromosome_methylation_map is chrmap.chromosome_methylation_map
