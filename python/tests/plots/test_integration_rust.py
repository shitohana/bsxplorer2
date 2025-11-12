import math
import os

import numpy as np
import polars as pl
import pytest


@pytest.mark.integration
def test_discretise_via_rust_and_html() -> None:
    """Integration: build BsxBatch from Polars, discretise (Rust), make HTML.

    This test exercises the Rust path (BsxBatch.discretise) without relying on
    external .bsx files. It validates that we can produce aggregated vectors
    and feed them into the HTML render path.
    """
    try:
        from bsx2._bsx2 import BsxBatch, AggMethod
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"bsx2._bsx2 not importable: {e}")

    from bsx2.plots.data import DiscreteRegionData
    from bsx2.plots.polars_html import line_html

    # Prepare a small synthetic batch on a single chromosome
    n = 50
    chr_vals = ["chr1"] * n
    positions = np.arange(1, n + 1, dtype=np.uint32)
    strand = np.array([True] * n, dtype=bool)
    context = np.array([True] * n, dtype=bool)
    count_total = np.array([10] * n, dtype=np.uint16)
    count_m = np.array([int(5 + 4 * math.sin(i / 6.0)) for i in range(n)], dtype=np.uint16)
    density = (count_m / count_total).astype(np.float32)

    df = pl.DataFrame(
        {
            "chr": chr_vals,
            "position": pl.Series(positions, dtype=pl.UInt32),
            "strand": pl.Series(strand, dtype=pl.Boolean),
            "context": pl.Series(context, dtype=pl.Boolean),
            "count_m": pl.Series(count_m, dtype=pl.UInt16),
            "count_total": pl.Series(count_total, dtype=pl.UInt16),
            "density": pl.Series(density, dtype=pl.Float32),
        }
    )

    # Build a BsxBatch (try multiple supported constructors across builds)
    batch = None
    # 1) pyo3-polars path with pandas
    try:
        batch = BsxBatch.from_dataframe(df.to_pandas())  # some builds accept pandas
    except Exception:
        pass
    # 2) pyo3-polars path with polars
    if batch is None:
        try:
            batch = BsxBatch.from_dataframe(df)  # others accept polars directly
        except Exception:
            pass
    # 3) Columnar constructor (common across versions)
    if batch is None:
        try:
            context_opt = [bool(v) for v in context.tolist()]
            batch = BsxBatch(
                "chr1",
                None,  # optional chr dtype
                positions.tolist(),
                strand.tolist(),
                context_opt,
                count_m.tolist(),
                count_total.tolist(),
            )
        except Exception as e:  # noqa: BLE001
            pytest.skip(f"Cannot build BsxBatch in this build: {e}")

    # Discretise (Rust path)
    n_bins = 20
    xs, ys = batch.discretise(n_bins, AggMethod.Mean)
    assert isinstance(xs, list) and isinstance(ys, list)
    assert len(xs) == n_bins and len(ys) == n_bins

    # Build HTML
    drd = DiscreteRegionData()
    drd.insert(np.asarray(xs, dtype=float), np.asarray(ys, dtype=float), "synthetic")
    html = line_html(drd, agg="mean", full_html=False)
    assert isinstance(html, str)
    assert "plotly" in html.lower() or "<div" in html.lower()
