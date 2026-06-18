import pandas as pd

from dmr_validation_framework.workflows.external_callers import (
    build_target_support_matrix,
    import_external_tables,
    load_target_regions,
    parse_caller_output,
)


def test_parse_caller_output_with_context_and_windows_drive_path():
    spec = parse_caller_output(r"DSS:CG=G:\calls\dss_dmrs_CG.tsv")
    assert spec.caller == "DSS"
    assert spec.context == "CG"
    assert str(spec.path) == r"G:\calls\dss_dmrs_CG.tsv"


def test_external_workflow_scores_target_support(tmp_path):
    target = tmp_path / "target.tsv"
    pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "chrom": ["Chr1", "Chr1"],
            "start": [100, 1000],
            "end": [200, 1100],
            "context": ["CG", "CG"],
            "delta_methylation": [0.3, -0.2],
            "q_value": [0.01, 0.04],
            "evidence_class": ["strong", "weak"],
        }
    ).to_csv(target, sep="\t", index=False)

    dss = tmp_path / "dss.tsv"
    pd.DataFrame(
        {
            "dmr_id": ["dss1"],
            "chrom": ["Chr1"],
            "start": [110],
            "end": [190],
            "context": ["CG"],
            "delta": [0.25],
            "q_value": [0.02],
        }
    ).to_csv(dss, sep="\t", index=False)

    external, inventory = import_external_tables(
        [parse_caller_output(f"DSS:CG={dss}")],
        contrast_id="control_vs_treatment",
        condition_a="control",
        condition_b="treatment",
    )
    assert inventory[0]["status"] == "ok"
    target_df = load_target_regions(target, include_candidate_only=False, max_target_regions=10)
    support = build_target_support_matrix(target_df, external, overlap_threshold=0.5)

    assert support.loc[0, "n_callers_supporting"] == 1
    assert bool(support.loc[0, "DSS_support"])
    assert support.loc[0, "best_external_q_value"] == 0.02
    assert support.loc[1, "n_callers_supporting"] == 0
