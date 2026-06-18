import pandas as pd

from bsx2.analysis.external_dmr_callers import ADAPTERS, adapter_for_caller
from bsx2.analysis.dmr_harmonization import add_dmr_tiers, build_caller_support_matrix


def test_support_matrix_overlap():
    a = pd.DataFrame({"source_caller": ["DSS"], "chrom": ["A"], "start": [0], "end": [10], "context": ["CG"], "q_value": [0.01], "delta": [0.2]})
    b = pd.DataFrame({"source_caller": ["methylKit"], "chrom": ["A"], "start": [2], "end": [9], "context": ["CG"], "q_value": [0.02], "delta": [0.3]})
    out = build_caller_support_matrix([a, b])
    assert out["n_callers_supporting"].iloc[0] == 2


def test_selected_external_caller_adapters_are_registered():
    callers = {
        "dmrseq",
        "RADMeth",
        "MOABS",
        "methylSig",
        "BSmooth",
        "BiSeq",
        "DMRcate",
        "MethyLasso",
        "HMM-DM",
        "comb-p",
    }
    assert callers.issubset(set(ADAPTERS))


def test_dmrseq_adapter_maps_common_columns(tmp_path):
    path = tmp_path / "dmrseq.tsv"
    pd.DataFrame(
        {
            "chr": ["Chr1"],
            "start": [10],
            "end": [20],
            "beta": [0.4],
            "pval": [0.001],
            "qval": [0.01],
            "L": [5],
        }
    ).to_csv(path, sep="\t", index=False)
    out = adapter_for_caller("dmrseq")(path, context="CG").read()
    assert out.loc[0, "source_caller"] == "dmrseq"
    assert out.loc[0, "chrom"] == "Chr1"
    assert out.loc[0, "delta"] == 0.4
    assert out.loc[0, "q_value"] == 0.01
    assert out.loc[0, "n_cytosines"] == 5
    assert out.loc[0, "context"] == "CG"


def test_new_beta_binomial_adapter_maps_aliases(tmp_path):
    path = tmp_path / "radmeth.tsv"
    pd.DataFrame(
        {
            "chromosome": ["Chr2"],
            "begin": [100],
            "stop": [130],
            "meth_diff": [-0.25],
            "p_value": [0.02],
            "fdr": [0.04],
            "n_cpg": [7],
        }
    ).to_csv(path, sep="\t", index=False)
    out = adapter_for_caller("RADMeth")(path).read()
    assert out.loc[0, "source_caller"] == "RADMeth"
    assert out.loc[0, "chrom"] == "Chr2"
    assert out.loc[0, "delta"] == -0.25
    assert out.loc[0, "q_value"] == 0.04
    assert "model_family=beta_binomial" in out.loc[0, "method_notes"]


def test_one_based_caller_coordinates_normalized_to_half_open(tmp_path):
    path = tmp_path / "dss.tsv"
    pd.DataFrame(
        {"chr": ["Chr1"], "start": [101], "end": [200], "diff": [0.2], "fdr": [0.01]}
    ).to_csv(path, sep="\t", index=False)
    out = adapter_for_caller("DSS")(path, context="CG").read()
    # 1-based inclusive [101, 200] -> 0-based half-open [100, 200)
    assert int(out.loc[0, "start"]) == 100
    assert int(out.loc[0, "end"]) == 200


def test_opposite_delta_sign_conventions_agree_after_harmonization(tmp_path):
    # Same biology: condition_b is hypermethylated relative to condition_a.
    # methylKit reports b - a (positive); DSS reports group1 - group2 (negative).
    # Without sign harmonization the two callers would look discordant.
    methylkit = tmp_path / "methylkit.tsv"
    pd.DataFrame(
        {"chr": ["Chr1"], "start": [101], "end": [200], "meth.diff": [0.3], "qvalue": [0.01]}
    ).to_csv(methylkit, sep="\t", index=False)
    dss = tmp_path / "dss.tsv"
    pd.DataFrame(
        {"chr": ["Chr1"], "start": [101], "end": [200], "diff": [-0.3], "fdr": [0.01]}
    ).to_csv(dss, sep="\t", index=False)

    mk = adapter_for_caller("methylKit")(methylkit, context="CG").read()
    ds = adapter_for_caller("DSS")(dss, context="CG").read()
    assert mk.loc[0, "delta"] > 0
    assert ds.loc[0, "delta"] > 0  # DSS sign flipped onto condition_b - condition_a

    support = build_caller_support_matrix([mk, ds])
    assert support.loc[0, "n_callers_supporting"] == 2
    assert support.loc[0, "direction_consensus"] == "same"
    assert not bool(support.loc[0, "caller_conflict_flag"])


def test_tier_1_consensus_validated():
    tables = [
        pd.DataFrame({"dmr_id": ["r1"], "source_caller": ["dmrseq"], "chrom": ["Chr1"], "start": [100], "end": [200], "context": ["CG"], "q_value": [0.01], "delta": [0.25]}),
        pd.DataFrame({"dmr_id": ["r1_b"], "source_caller": ["RADMeth"], "chrom": ["Chr1"], "start": [110], "end": [205], "context": ["CG"], "q_value": [0.03], "delta": [0.22]}),
        pd.DataFrame({"dmr_id": ["r1_c"], "source_caller": ["DMRcate"], "chrom": ["Chr1"], "start": [105], "end": [195], "context": ["CG"], "q_value": [0.04], "delta": [0.21]}),
    ]
    support = build_caller_support_matrix(tables)
    validation = pd.DataFrame(
        {
            "region_id": ["r1"],
            "final_confidence_class": ["HIGH_CONFIDENCE"],
            "model_agreement_status": ["GLMM_CONFIRMED"],
            "direction_changed": [False],
            "loo_direction_changed": [False],
            "ci_includes_zero": [False],
        }
    )
    out = add_dmr_tiers(support, validation)
    assert out.loc[0, "final_dmr_tier"] == "TIER_1_CONSENSUS_VALIDATED"


def test_tier_2_single_caller_validated():
    support = build_caller_support_matrix(
        [
            pd.DataFrame({"dmr_id": ["single"], "source_caller": ["dmrseq"], "chrom": ["Chr1"], "start": [100], "end": [200], "context": ["CG"], "q_value": [0.01], "delta": [0.25]})
        ]
    )
    validation = pd.DataFrame({"region_id": ["single"], "final_confidence_class": ["HIGH_CONFIDENCE"], "model_agreement_status": ["GLMM_CONFIRMED"]})
    out = add_dmr_tiers(support, validation)
    assert out.loc[0, "final_dmr_tier"] == "TIER_2_SINGLE_CALLER_VALIDATED"


def test_tier_4_caller_discordant():
    support = build_caller_support_matrix(
        [
            pd.DataFrame({"source_caller": ["dmrseq"], "chrom": ["Chr1"], "start": [100], "end": [200], "context": ["CG"], "q_value": [0.01], "delta": [0.25]}),
            pd.DataFrame({"source_caller": ["RADMeth"], "chrom": ["Chr1"], "start": [105], "end": [190], "context": ["CG"], "q_value": [0.02], "delta": [-0.20]}),
        ]
    )
    assert support.loc[0, "final_dmr_tier"] == "TIER_4_CALLER_DISCORDANT"


def test_tier_5_rejected_by_audit():
    support = build_caller_support_matrix(
        [
            pd.DataFrame({"dmr_id": ["unstable"], "source_caller": ["dmrseq"], "chrom": ["Chr1"], "start": [100], "end": [200], "context": ["CG"], "q_value": [0.01], "delta": [0.25]})
        ]
    )
    validation = pd.DataFrame({"region_id": ["unstable"], "final_confidence_class": ["HIGH_CONFIDENCE"], "model_agreement_status": ["GLMM_CONFIRMED"], "loo_direction_changed": [True]})
    out = add_dmr_tiers(support, validation)
    assert out.loc[0, "final_dmr_tier"] == "TIER_5_REJECTED_BY_AUDIT"
