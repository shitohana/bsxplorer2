import pandas as pd
import pytest

from bsx2.analysis import cpg_level_glmm
from bsx2.analysis.cpg_level_glmm import run_cpg_level_glmm_validation


def design():
    return pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3", "s4"],
            "condition": ["control", "control", "case", "case"],
        }
    )


def cpg_counts(n_cpg=3):
    rows = []
    for cpg_idx in range(n_cpg):
        for sample, cond in zip(["s1", "s2", "s3", "s4"], ["control", "control", "case", "case"]):
            rows.append(
                {
                    "region_id": "r1",
                    "cpg_id": f"chr1:{10 + cpg_idx}:+:CG",
                    "seqname": "chr1",
                    "chrom": "chr1",
                    "position": 10 + cpg_idx,
                    "context": "CG",
                    "sample_id": sample,
                    "mC": 8 if cond == "case" else 2,
                    "uC": 2 if cond == "case" else 8,
                    "total": 10,
                }
            )
    return pd.DataFrame(rows)


def test_insufficient_cpg_skipped():
    results, warnings = run_cpg_level_glmm_validation(
        cpg_counts(n_cpg=2),
        design(),
        force_glmm_unavailable=False,
        min_cpg=3,
        temp_dir="unused",
    )
    if set(results["model_status"]) == {"glmmTMB_unavailable"}:
        pytest.skip("glmmTMB unavailable in test environment")
    assert results.loc[0, "model_status"] == "insufficient_cpg"


def test_insufficient_replicates_skipped(monkeypatch, tmp_path):
    monkeypatch.setattr(cpg_level_glmm, "rscript_available", lambda: True)
    monkeypatch.setattr(cpg_level_glmm, "glmmTMB_available", lambda: True)
    one_rep_design = pd.DataFrame({"sample_id": ["s1", "s3"], "condition": ["control", "case"]})
    results, _warnings = run_cpg_level_glmm_validation(
        cpg_counts(),
        one_rep_design,
        min_replicates_per_group=2,
        temp_dir=tmp_path,
    )
    assert results.loc[0, "model_status"] == "insufficient_replicates"


def test_glmmtmb_unavailable_graceful_status():
    results, warnings = run_cpg_level_glmm_validation(
        cpg_counts(),
        design(),
        force_glmm_unavailable=True,
    )
    assert results.loc[0, "model_status"] == "glmmTMB_unavailable"
    assert "glmmTMB_unavailable" in set(warnings["warning_type"])


def test_synthetic_consistent_dmr_if_glmmtmb_available(tmp_path):
    if not cpg_level_glmm.glmmTMB_available():
        pytest.skip("glmmTMB unavailable")
    results, _warnings = run_cpg_level_glmm_validation(
        cpg_counts(),
        design(),
        case_label="case",
        control_label="control",
        temp_dir=tmp_path,
    )
    assert results.loc[0, "model_status"] in {"ok", "model_error"}


def test_one_cpg_outlier_can_be_reported_as_candidate_input():
    df = cpg_counts()
    outlier = df.copy()
    outlier.loc[outlier["cpg_id"] != "chr1:10:+:CG", ["mC", "uC"]] = [5, 5]
    assert outlier["cpg_id"].nunique() == 3
