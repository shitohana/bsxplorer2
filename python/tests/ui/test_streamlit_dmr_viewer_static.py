from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
APP_PATH = ROOT / "apps" / "bsx2_dmr_viewer.py"


def load_app_module():
    spec = importlib.util.spec_from_file_location("bsx2_dmr_viewer_static", APP_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_app_file_exists():
    assert APP_PATH.exists()


def test_app_has_required_title_and_scope_note():
    text = APP_PATH.read_text(encoding="utf-8")
    assert "BSX2 DMR Evidence Viewer" in text
    assert "does not run DMR calling" in text


def test_app_has_no_local_hardcoded_paths():
    text = APP_PATH.read_text(encoding="utf-8")
    forbidden = ["/mnt/g", "G:\\", "/mnt/c", "C:\\", "/home/"]
    assert not any(marker in text for marker in forbidden)


def test_static_helpers_filter_table():
    app = load_app_module()
    df = pd.DataFrame(
        {
            "dmr_id": ["r1", "r2"],
            "context": ["CG", "CHH"],
            "region_delta": [0.3, 0.05],
            "region_q_value": [0.01, 0.2],
            "evidence_class": ["strong", "weak"],
        }
    )
    normalized = app.normalize_columns(df)
    filtered = app.apply_filters(
        normalized,
        contexts=["CG"],
        evidence_classes=["strong"],
        q_threshold=0.05,
        abs_delta_threshold=0.2,
        top_n=10,
    )
    assert filtered["region_id"].tolist() == ["r1"]
