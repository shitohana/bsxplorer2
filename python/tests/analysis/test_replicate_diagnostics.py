import pandas as pd

from bsx2.analysis.replicate_diagnostics import compute_pairwise_replicate_consistency


def test_low_replicate_count_warning():
    counts = pd.DataFrame({"region_id": ["r1", "r1"], "sample_id": ["s1", "s2"], "Y": [1, 3], "m": [4, 4]})
    design = pd.DataFrame({"sample_id": ["s1", "s2"], "condition": ["a", "b"]})
    out = compute_pairwise_replicate_consistency(counts, design)
    assert "low_replicate_count" in out["dispersion_warning"].iloc[0]
