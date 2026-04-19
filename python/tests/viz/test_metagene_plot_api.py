from __future__ import annotations

import importlib

import numpy as np
import pytest

bsx2 = pytest.importorskip("bsx2")
if not hasattr(bsx2, "Context"):
    pytest.skip("bsx2 extension is unavailable", allow_module_level=True)

viz = pytest.importorskip("bsx2.viz")
plot_data = pytest.importorskip("bsx2.viz.data")

DiscreteRegionData = plot_data.DiscreteRegionData
MetageneProfileSegment = viz.MetageneProfileSegment
box_plot = viz.box_plot
heatmap = viz.heatmap
line_plot = viz.line_plot
violin_plot = viz.violin_plot


def _is_holoviews_object(obj) -> bool:
    return obj.__class__.__module__.startswith("holoviews")


def _sample_drd() -> DiscreteRegionData:
    drd = DiscreteRegionData()
    positions = np.array([0.1, 0.3, 0.5, 0.7, 0.9], dtype=float)
    drd.insert(
        positions,
        np.array([0.20, 0.30, 0.40, 0.55, 0.60], dtype=float),
        "gene_a",
    )
    drd.insert(
        positions,
        np.array([0.10, np.nan, 0.35, 0.45, 0.55], dtype=float),
        "gene_b",
    )
    return drd


def _segments() -> list[MetageneProfileSegment]:
    return [
        MetageneProfileSegment("up", 2),
        MetageneProfileSegment("body", 2),
        MetageneProfileSegment("down", 2),
    ]


def test_public_line_plot_wrapper_returns_holoviews_curve() -> None:
    plot = line_plot(
        _sample_drd(),
        name="sample",
        segments=_segments(),
        smooth=None,
    )

    assert _is_holoviews_object(plot)
    assert plot.__class__.__name__ == "Curve"


def test_public_heatmap_wrapper_returns_holoviews_object() -> None:
    plot = heatmap(
        _sample_drd(),
        segments=_segments(),
        rank_rows=4,
    )

    assert _is_holoviews_object(plot)
    assert plot.__class__.__name__ in {"Image", "Overlay"}


def test_public_box_plot_wrapper_returns_boxwhisker() -> None:
    plot = box_plot(
        _sample_drd(),
        segments=_segments(),
        as_percent=False,
    )

    assert _is_holoviews_object(plot)
    assert plot.__class__.__name__ == "BoxWhisker"


def test_public_violin_plot_wrapper_returns_violin() -> None:
    plot = violin_plot(
        _sample_drd(),
        per_region=True,
        as_percent=False,
    )

    assert _is_holoviews_object(plot)
    assert plot.__class__.__name__ == "Violin"


def test_public_surface_does_not_expose_html_helpers() -> None:
    assert not hasattr(viz, "box_html")
    assert not hasattr(viz, "violin_html")

    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("bsx2.viz.compat")

    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("bsx2.viz.polars_html")

    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("bsx2.plots")


def test_public_composers_do_not_expose_to_html() -> None:
    assert not hasattr(viz.HeatmapPlotComposer, "to_html")
    assert not hasattr(viz.ChrLinePlotComposer, "to_html")
    assert not hasattr(viz.GeneEmbeddingPlotComposer, "to_html")
    assert not hasattr(viz.GeneDendrogramPlotComposer, "to_html")
