import pandas as pd
import pytest

from bsx2.analysis.region_cpg_counts import (
    REGION_CPG_COUNTS_COLUMNS,
    extract_region_cpg_counts_pandas,
)


def counts():
    return pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s1"],
            "chrom": ["chr1", "chr1", "chr1"],
            "position": [10, 20, 30],
            "strand": ["+", "-", "+"],
            "context": ["CG", "CHG", "CG"],
            "mC": [4, 2, 7],
            "uC": [6, 3, 3],
            "total": [10, 5, 10],
        }
    )


def regions():
    return pd.DataFrame(
        {
            "region_id": ["r1", "r2"],
            "chrom": ["chr1", "chr1"],
            "start": [1, 100],
            "end": [25, 200],
        }
    )


def test_per_cpg_extraction_schema_stable():
    out = extract_region_cpg_counts_pandas(counts(), regions(), sample_id="s1")
    assert list(out.columns) == REGION_CPG_COUNTS_COLUMNS
    assert len(out) == 2


def test_cpg_id_stable():
    out = extract_region_cpg_counts_pandas(counts(), regions().head(1), sample_id="s1")
    assert out.loc[0, "cpg_id"] == "chr1:10:+:CG"


def test_context_filter_works():
    out = extract_region_cpg_counts_pandas(counts(), regions().head(1), context="CG")
    assert set(out["context"]) == {"CG"}
    assert len(out) == 1


def test_min_total_filter_works():
    out = extract_region_cpg_counts_pandas(counts(), regions().head(1), min_total=6)
    assert len(out) == 1
    assert out["total"].min() >= 6


def test_invalid_context_raises():
    with pytest.raises(ValueError):
        extract_region_cpg_counts_pandas(counts(), regions(), context="BAD")
