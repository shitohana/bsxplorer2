import pandas as pd

from bsx2.analysis.region_signal import RegionSignalConfig, aggregate_region_signal


def test_region_sum_and_alias():
    regions = pd.DataFrame({"region_id": ["r1"], "chrom": ["chrA01"], "start": [1], "end": [3]})
    counts = pd.DataFrame({"sample_id": ["s1"], "chrom": ["A01"], "position": [2], "mC": [3], "uC": [1], "context": ["CG"], "strand": ["+"]})
    out = aggregate_region_signal(regions, counts, RegionSignalConfig(seqname_aliases={"chrA01": "A01"}))
    assert out["mC"].iloc[0] == 3


def test_zero_coverage_flag():
    regions = pd.DataFrame({"region_id": ["r1"], "chrom": ["A"], "start": [1], "end": [3]})
    counts = pd.DataFrame({"sample_id": ["s1"], "chrom": ["A"], "position": [2], "mC": [0], "uC": [0]})
    out = aggregate_region_signal(regions, counts)
    assert out["coverage_qc"].iloc[0] == "zero_coverage"
