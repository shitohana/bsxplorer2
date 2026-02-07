import numpy as np
import pytest

from bsx2.plots import metagene as mg


def _base_cfg(**overrides):
    cfg = {
        "method": "savgol",
        "window_length": 5,
        "polyorder": 2,
        "apply": "post",
        "mode": "interp",
        "nan_policy": "interp",
        "per_segment": False,
    }
    cfg.update(overrides)
    return cfg


def test_legacy_smooth_int_coercion():
    cfg = mg._coerce_smooth_config(5, total_bins=20)
    assert cfg is not None
    assert cfg["method"] == "savgol"
    assert cfg["window_length"] == 5
    assert cfg["polyorder"] == 2
    assert cfg["apply"] == "post"


def test_validate_even_window_length():
    cfg = _base_cfg(window_length=4, polyorder=2)
    with pytest.raises(ValueError):
        mg._validate_savgol_config(cfg, series_len=10)


def test_validate_window_length_interp_too_long():
    cfg = _base_cfg(window_length=9, polyorder=2, mode="interp")
    with pytest.raises(ValueError):
        mg._validate_savgol_config(cfg, series_len=7)


def test_per_segment_length_mismatch():
    cfg = _base_cfg(per_segment=True)
    y = np.arange(6, dtype=float)
    with pytest.raises(ValueError):
        mg._apply_savgol_smoothing(y, cfg, segments=[mg.Segment("a", 2), mg.Segment("b", 5)])


def test_nan_policy_raise():
    pytest.importorskip("scipy")
    cfg = _base_cfg(nan_policy="raise")
    y = np.array([0.1, np.nan, 0.2, 0.3, 0.4], dtype=float)
    with pytest.raises(ValueError):
        mg._apply_savgol_smoothing(y, cfg, segments=None)


def test_nan_policy_mask_keeps_nans():
    pytest.importorskip("scipy")
    cfg = _base_cfg(nan_policy="mask")
    y = np.array([0.1, np.nan, 0.2, 0.3, 0.4], dtype=float)
    out = mg._apply_savgol_smoothing(y, cfg, segments=None)
    assert np.isnan(out[1])


def test_per_segment_no_bleed():
    pytest.importorskip("scipy")
    segments = [mg.Segment("left", 5), mg.Segment("right", 5)]
    y = np.array([0, 0, 0, 0, 0, 10, 10, 10, 10, 10], dtype=float)
    cfg_all = _base_cfg(window_length=5, polyorder=1, per_segment=False)
    cfg_seg = _base_cfg(window_length=5, polyorder=1, per_segment=True)
    out_all = mg._apply_savgol_smoothing(y, cfg_all, segments=segments)
    out_seg = mg._apply_savgol_smoothing(y, cfg_seg, segments=segments)
    # With per_segment, the boundary should stay exact 0/10
    assert out_seg[4] == 0
    assert out_seg[5] == 10
    # Without per_segment, smoothing should bleed across the boundary
    assert out_all[4] > 0
    assert out_all[5] < 10
