import pandas as pd

from bsx2.analysis.dmr_harmonization import build_caller_support_matrix


def test_support_matrix_overlap():
    a = pd.DataFrame({"source_caller": ["DSS"], "chrom": ["A"], "start": [0], "end": [10], "context": ["CG"], "q_value": [0.01], "delta": [0.2]})
    b = pd.DataFrame({"source_caller": ["methylKit"], "chrom": ["A"], "start": [2], "end": [9], "context": ["CG"], "q_value": [0.02], "delta": [0.3]})
    out = build_caller_support_matrix([a, b])
    assert out["n_callers_supporting"].iloc[0] == 2
