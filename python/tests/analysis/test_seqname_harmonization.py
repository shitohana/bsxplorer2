import pandas as pd

from bsx2.analysis.seqname_harmonization import apply_seqname_aliases, compare_seqname_sets, normalize_seqname


def test_chr_alias_maps():
    assert normalize_seqname("chrA01", {"chrA01": "A01"}) == "A01"


def test_arbitrary_scaffold_preserved():
    assert normalize_seqname("scaffold_42") == "scaffold_42"


def test_missing_aliases_warn_not_crash():
    result = compare_seqname_sets(["A", "B"], ["A", "C"])
    assert result["n_common"] == 1
    assert "B" in result["left_only"]


def test_dataframe_not_mutated():
    df = pd.DataFrame({"chrom": ["chr1"]})
    out = apply_seqname_aliases(df, "chrom", {"chr1": "1"})
    assert df["chrom"].iloc[0] == "chr1"
    assert out["chrom"].iloc[0] == "1"
