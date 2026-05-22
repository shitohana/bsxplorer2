import pandas as pd

from bsx2.analysis.glmm_comparison import compare_beta_binomial_glm_vs_glmm


def test_glm_glmm_confirmed_rule():
    glm = pd.DataFrame({"region_id": ["r1"], "context": ["CG"], "glm_p_value": [0.001], "glm_q_value": [0.01], "glm_status": ["ok"]})
    glmm = pd.DataFrame({"region_id": ["r1"], "context": ["CG"], "glmm_p_value": [0.002], "glmm_q_value": [0.02], "glmm_status": ["ok"]})
    out = compare_beta_binomial_glm_vs_glmm(glm, glmm)
    assert bool(out.loc[0, "confirmed_by_glmm"])
    assert not bool(out.loc[0, "glm_only_candidate"])


def test_glm_only_candidate_rule():
    glm = pd.DataFrame({"region_id": ["r1"], "context": ["CG"], "glm_p_value": [0.001], "glm_q_value": [0.01], "glm_status": ["ok"]})
    glmm = pd.DataFrame({"region_id": ["r1"], "context": ["CG"], "glmm_p_value": [0.2], "glmm_q_value": [0.2], "glmm_status": ["ok"]})
    out = compare_beta_binomial_glm_vs_glmm(glm, glmm)
    assert bool(out.loc[0, "glm_only_candidate"])
    assert bool(out.loc[0, "glmm_more_conservative"])


def test_insufficient_glmm_data_rule():
    glm = pd.DataFrame({"region_id": ["r1"], "context": ["CG"], "glm_p_value": [0.001], "glm_q_value": [0.01], "glm_status": ["ok"]})
    glmm = pd.DataFrame({"region_id": ["r1"], "context": ["CG"], "glmm_p_value": [pd.NA], "glmm_q_value": [pd.NA], "glmm_status": ["insufficient_cpg"]})
    out = compare_beta_binomial_glm_vs_glmm(glm, glmm)
    assert bool(out.loc[0, "glm_only_candidate"])
    assert "insufficient data" in out.loc[0, "interpretation"]
