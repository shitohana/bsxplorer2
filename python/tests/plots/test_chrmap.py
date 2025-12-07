import numpy as np
import polars as pl
import pytest

from bsx2.plots.chrmap import (
    prepare_chr_line_data,
    prepare_chr_box_data,
)


def _make_df(include_ci: bool = False) -> pl.DataFrame:
    rows = []
    for chr_name in ["chr1", "chr2"]:
        for w in range(3):
            rows.append(
                {
                    "chr": chr_name,
                    "window": w,
                    "sum": 10 + w,
                    "count": 20,
                    "upper": 0.7 if include_ci else None,
                    "lower": 0.3 if include_ci else None,
                }
            )
    df = pl.DataFrame(rows)
    return df if include_ci else df.drop(["upper", "lower"])


def test_prepare_chr_line_data_basic():
    df = _make_df()
    data = prepare_chr_line_data(df, smooth=0)
    assert data.y.shape[0] == len(df)
    assert len(data.x_ticks) == len(data.x_labels) == 2
    assert data.borders[0] == 0 and data.borders[-1] == len(df)
    assert data.lower is None and data.upper is None
    # densities should be percent of sum/count
    assert np.allclose(data.y[:2], np.array([10, 11]) / 20 * 100)


def test_prepare_chr_line_data_smooth_and_ci():
    df = _make_df(include_ci=True)
    data = prepare_chr_line_data(df, smooth=3)
    assert data.lower is not None and data.upper is not None
    # smoothing keeps length
    assert data.lower.shape == data.y.shape


def test_prepare_chr_box_data():
    df = _make_df()
    box = prepare_chr_box_data(df)
    assert box.labels == ["chr1", "chr2"]
    assert len(box.values) == 2
    # each chromosome has 3 window densities
    assert box.values[0].shape[0] == 3
