import json
from pathlib import Path

import pandas as pd

from bsx2.viz.dmr_curve_data import normalize_dmr_curve_columns
from bsx2.viz.dmr_curves import (
    DmrCurveSpec,
    make_dmr_box_curve,
    make_dmr_chromosome_curve,
    make_dmr_heatmap_curve,
    make_dmr_line_curve,
    make_dmr_pca_curve,
    make_dmr_violin_curve,
)


def _dmr_table():
    return pd.DataFrame(
        {
            "chr": ["chr1", "chr1", "chr2"],
            "start_bp": [10, 40, 5],
            "end_bp": [20, 60, 25],
            "context": ["CG", "CHG", "CG"],
            "delta": [0.25, -0.3, 0.1],
            "q": [0.01, 0.2, 0.04],
            "evidence_class": ["strong", "weak", "moderate"],
        }
    )


def _region_counts():
    return pd.DataFrame(
        {
            "region_id": ["r1", "r1", "r1", "r1", "r2", "r2", "r2", "r2"],
            "sample_id": ["s1", "s2", "s3", "s4"] * 2,
            "Y": [20, 22, 5, 6, 8, 9, 20, 23],
            "m": [40, 40, 40, 40, 30, 30, 30, 30],
        }
    )


def _design():
    return pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3", "s4"],
            "condition": ["A", "A", "B", "B"],
        }
    )


def test_curve_spec_json_roundtrip(tmp_path):
    spec = DmrCurveSpec(
        curve_id="c1",
        curve_type="dmr_chromosome",
        title="Title",
        description="Desc",
        input_tables={"dmr": "x.tsv"},
        data_mapping={"x": "chrom"},
        plot_params={},
        filters={},
        renderer="table",
    )
    path = tmp_path / "spec.json"
    spec.to_json(path)
    loaded = DmrCurveSpec.from_json(path)
    assert loaded.to_dict() == spec.to_dict()
    assert loaded.validate() == []


def test_column_normalization_aliases():
    df, warnings = normalize_dmr_curve_columns(_dmr_table())
    assert warnings == ["region_id was synthesized from chrom/start/end"]
    assert {"chrom", "start", "end", "q_value", "delta", "region_id"}.issubset(df.columns)
    assert df.loc[0, "chrom"] == "chr1"


def test_chromosome_curve_from_synthetic_dmr():
    result = make_dmr_chromosome_curve(_dmr_table(), render=False)
    assert result.table["n_dmrs"].sum() == 3
    assert result.spec.curve_type == "dmr_chromosome"


def test_box_and_violin_curve_specs_from_synthetic_dmr():
    box = make_dmr_box_curve(dmr_table=_dmr_table(), render=False)
    violin = make_dmr_violin_curve(dmr_table=_dmr_table(), render=False)
    assert box.spec.curve_type == "dmr_box"
    assert violin.spec.curve_type == "dmr_violin"
    assert len(box.table) == 3
    assert len(violin.table) == 3


def test_pca_curve_from_region_counts_and_design():
    result = make_dmr_pca_curve(_region_counts(), _design(), render=False)
    assert result.spec.curve_type == "dmr_pca"
    assert {"sample_id", "PC1", "PC2"}.issubset(result.table.columns)
    assert result.summary["n_samples"] == 4


def test_heatmap_curve_from_region_counts_and_design():
    result = make_dmr_heatmap_curve(_region_counts(), _design(), render=False)
    assert result.spec.curve_type == "dmr_heatmap"
    assert "region_id" in result.table.columns
    assert result.summary["n_samples"] == 4


def test_missing_positional_signal_warns_not_crash():
    result = make_dmr_line_curve(_dmr_table(), render=False)
    assert result.table is None
    assert result.warnings


def test_dashboard_payload_is_json_serializable():
    result = make_dmr_chromosome_curve(_dmr_table(), render=False)
    payload = result.to_dashboard_payload()
    json.dumps(payload)
    assert payload["curve_type"] == "dmr_chromosome"


def test_no_hardcoded_local_paths_in_dmr_curves():
    path = Path(__file__).resolve().parents[2] / "src" / "bsx2" / "viz" / "dmr_curves.py"
    text = path.read_text(encoding="utf-8")
    for forbidden in ["/mnt/g", "G:\\", "/mnt/c", "C:\\"]:
        assert forbidden not in text
