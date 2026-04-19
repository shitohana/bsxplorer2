from __future__ import annotations

import os
from time import perf_counter

import numpy as np
import pytest

from bsx2.viz import MetageneProfileSegment, heatmap, line_plot

from .conftest import SyntheticPerfScenario, make_synthetic_heatmap_drd


def _is_holoviews_object(obj) -> bool:
    return obj.__class__.__module__.startswith("holoviews")


def _require_extended_perf() -> None:
    if os.environ.get("BSX2_RUN_PERF") != "1":
        pytest.skip("Set BSX2_RUN_PERF=1 to execute medium/large perf smoke scenarios")


@pytest.mark.perf
@pytest.mark.parametrize("scenario_name", ["small"])
def test_viz_perf_smoke(
    perf_scenarios: dict[str, SyntheticPerfScenario],
    scenario_name: str,
) -> None:
    scenario = perf_scenarios[scenario_name]
    drd = make_synthetic_heatmap_drd(scenario)
    segments = [MetageneProfileSegment("body", scenario.n_bins)]

    t_line = perf_counter()
    line = line_plot(drd, name=scenario.name, segments=segments, smooth=None)
    line_s = perf_counter() - t_line

    t_heatmap = perf_counter()
    hm = heatmap(drd, segments=segments, rank_rows=100)
    heatmap_s = perf_counter() - t_heatmap

    assert _is_holoviews_object(line)
    assert _is_holoviews_object(hm)
    assert np.isfinite(line_s)
    assert np.isfinite(heatmap_s)


@pytest.mark.perf
@pytest.mark.parametrize("scenario_name", ["medium", "large"])
def test_viz_perf_smoke_extended(
    perf_scenarios: dict[str, SyntheticPerfScenario],
    scenario_name: str,
) -> None:
    _require_extended_perf()
    scenario = perf_scenarios[scenario_name]
    drd = make_synthetic_heatmap_drd(scenario)
    segments = [MetageneProfileSegment("body", scenario.n_bins)]

    line = line_plot(drd, name=scenario.name, segments=segments, smooth=None)
    hm = heatmap(drd, segments=segments, rank_rows=100)

    assert _is_holoviews_object(line)
    assert _is_holoviews_object(hm)
