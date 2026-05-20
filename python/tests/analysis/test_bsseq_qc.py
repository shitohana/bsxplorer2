import pandas as pd

from bsx2.analysis.bsseq_qc import compute_counts_qc


def test_counts_qc():
    qc = compute_counts_qc(pd.DataFrame({"mC": [1, 2], "uC": [3, 4], "context": ["CG", "CHG"]}))
    assert qc["total_coverage"] == 10
