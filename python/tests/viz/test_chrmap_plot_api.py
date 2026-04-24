from __future__ import annotations

import numpy as np
import pytest

bsx2 = pytest.importorskip("bsx2")
if not hasattr(bsx2, "Context"):
    pytest.skip("bsx2 extension is unavailable", allow_module_level=True)

viz = pytest.importorskip("bsx2.viz")
chrline = pytest.importorskip("bsx2.viz.chrline")

Context = bsx2.Context
compute_chr_line_data = chrline.compute_chr_line_data
build_chromosome_manhattan_plot = viz.build_chromosome_manhattan_plot


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


def test_public_chromosome_manhattan_plot_wrapper_returns_holoviews_object() -> None:
    reader = IterableReader(
        [
            FakeBatch(
                "chr1",
                position=[1, 25_000, 75_000],
                count_m=[1, 3, 2],
                count_total=[2, 6, 4],
                context=[True, True, True],
            ),
            FakeBatch(
                "chr2",
                position=[1, 60_000],
                count_m=[2, 5],
                count_total=[4, 10],
                context=[True, True],
            ),
        ]
    )
    data = compute_chr_line_data(
        reader,
        context=Context.CG,
        bin_size_bp=50_000,
        chr_lengths={"chr1": 100_000, "chr2": 100_000},
    )

    plot = build_chromosome_manhattan_plot(
        data,
        title="Chromosome methylation Manhattan plot",
        width=900,
        height=360,
        point_size=6,
    )

    assert _is_holoviews_object(plot)
    assert plot.__class__.__name__ in {"Overlay", "Points"}
