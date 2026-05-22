import pandas as pd
import pytest

from bsx2.analysis.dmr_functional_prioritization import (
    FUNCTIONAL_PRIORITIZATION_COLUMNS,
    join_expression_evidence,
    link_dmrs_to_genes,
    link_dmrs_to_te,
    prioritize_dmr_linked_genes,
)


def _dmrs():
    return pd.DataFrame({
        "dmr_id": ["prom_hyper", "body", "intergenic", "prom_hypo"],
        "chrom": ["chr1", "chr1", "chr1", "chr2"],
        "start": [900, 1500, 7000, 900],
        "end": [950, 1550, 7100, 950],
        "context": ["CG", "CG", "CHH", "CG"],
        "region_delta": [0.3, 0.2, 0.1, -0.3],
        "region_q_value": [0.01, 0.02, 0.2, 0.01],
        "evidence_class": ["strong", "moderate", "weak", "strong"],
    })


def _genes():
    return pd.DataFrame({
        "gene_id": ["gene1", "gene1", "gene1", "gene2"],
        "chrom": ["chr1", "chr1", "chr1", "chr2"],
        "start": [1000, 1400, 1700, 1000],
        "end": [2000, 1600, 1800, 2000],
        "strand": ["+", "+", "+", "+"],
        "feature_type": ["gene", "exon", "intron", "gene"],
    })


def test_dmr_promoter_gene_body_and_intergenic_links():
    linked = link_dmrs_to_genes(_dmrs(), _genes(), promoter_upstream=1000, promoter_downstream=200, max_distance=1000)
    assert linked.loc[linked["dmr_id"] == "prom_hyper", "link_type"].iloc[0] == "promoter"
    assert linked.loc[linked["dmr_id"] == "body", "link_type"].iloc[0] == "gene_body"
    assert linked.loc[linked["dmr_id"] == "body", "exon_overlap"].iloc[0]
    assert linked.loc[linked["dmr_id"] == "intergenic", "link_type"].iloc[0] == "intergenic"
    assert linked.loc[linked["dmr_id"] == "intergenic", "nearest_gene"].iloc[0] == "gene1"


def test_te_overlap_annotation():
    te = pd.DataFrame({
        "te_id": ["te1"],
        "chrom": ["chr1"],
        "start": [1480],
        "end": [1580],
        "te_family": ["Gypsy"],
        "te_class": ["LTR"],
    })
    linked = link_dmrs_to_te(link_dmrs_to_genes(_dmrs(), _genes()), te)
    row = linked[linked["dmr_id"] == "body"].iloc[0]
    assert row["te_overlap"]
    assert row["te_family"] == "Gypsy"


def test_expression_join_and_direction_consistency_rules():
    expr = pd.DataFrame({
        "gene_id": ["gene1", "gene2"],
        "log2FC": [-1.2, 1.5],
        "p_value": [0.001, 0.001],
        "q_value": [0.01, 0.01],
    })
    result = prioritize_dmr_linked_genes(
        _dmrs(),
        _genes(),
        expression_df=expr,
        promoter_upstream=1000,
        promoter_downstream=200,
        dmr_contrast_label="B_vs_A",
        expression_contrast_label="B_vs_A",
    )
    hyper = result[result["dmr_id"] == "prom_hyper"].iloc[0]
    hypo = result[result["dmr_id"] == "prom_hypo"].iloc[0]
    body = result[result["dmr_id"] == "body"].iloc[0]
    assert hyper["direction_consistency"] == "consistent"
    assert hypo["direction_consistency"] == "consistent"
    assert body["direction_consistency"] == "context_dependent"
    assert hyper["functional_support_class"] == "high_confidence_candidate"


def test_expression_without_contrast_labels_is_not_evaluated():
    expr = pd.DataFrame({"gene_id": ["gene1"], "log2FC": [-1.2], "p_value": [0.001], "q_value": [0.01]})
    result = prioritize_dmr_linked_genes(_dmrs(), _genes(), expression_df=expr, promoter_upstream=1000, promoter_downstream=200)
    row = result[result["dmr_id"] == "prom_hyper"].iloc[0]
    assert row["direction_consistency"] == "not_evaluated_contrast_not_validated"
    assert "expression_contrast_not_validated" in row["missing_evidence_flags"]


def test_expression_mismatch_with_require_flag_raises():
    expr = pd.DataFrame({"gene_id": ["gene1"], "log2FC": [-1.2], "p_value": [0.001], "q_value": [0.01]})
    with pytest.raises(ValueError, match="contrast labels do not match"):
        prioritize_dmr_linked_genes(
            _dmrs(),
            _genes(),
            expression_df=expr,
            promoter_upstream=1000,
            promoter_downstream=200,
            dmr_contrast_label="B_vs_A",
            expression_contrast_label="A_vs_B",
            require_matched_contrast=True,
        )


def test_missing_rna_seq_flag_and_stable_schema():
    result = prioritize_dmr_linked_genes(_dmrs(), _genes(), promoter_upstream=1000, promoter_downstream=200)
    assert list(result.columns) == FUNCTIONAL_PRIORITIZATION_COLUMNS
    assert "no_expression_data" in result.loc[result["dmr_id"] == "prom_hyper", "missing_evidence_flags"].iloc[0]
    assert result["functional_support_score"].between(0, 100).all()


def test_te_missing_and_no_overlap_flags_are_separated():
    missing = prioritize_dmr_linked_genes(_dmrs(), _genes(), promoter_upstream=1000, promoter_downstream=200)
    assert "te_input_missing" in missing.loc[missing["dmr_id"] == "prom_hyper", "missing_evidence_flags"].iloc[0]

    te = pd.DataFrame({"te_id": ["te_far"], "chrom": ["chr1"], "start": [9000], "end": [9100], "te_family": ["Gypsy"], "te_class": ["LTR"]})
    no_overlap = prioritize_dmr_linked_genes(_dmrs(), _genes(), te_df=te, promoter_upstream=1000, promoter_downstream=200)
    assert "te_no_overlap" in no_overlap.loc[no_overlap["dmr_id"] == "prom_hyper", "missing_evidence_flags"].iloc[0]


def test_chromatin_missing_and_no_overlap_flags_are_separated():
    missing = prioritize_dmr_linked_genes(_dmrs(), _genes(), promoter_upstream=1000, promoter_downstream=200)
    assert "chromatin_input_missing" in missing.loc[missing["dmr_id"] == "prom_hyper", "missing_evidence_flags"].iloc[0]

    peaks = pd.DataFrame({"peak_id": ["peak_far"], "chrom": ["chr1"], "start": [9000], "end": [9100], "signal": [5], "peak_type": ["ATAC"]})
    no_overlap = prioritize_dmr_linked_genes(_dmrs(), _genes(), chromatin_df=peaks, promoter_upstream=1000, promoter_downstream=200)
    assert "chromatin_no_overlap" in no_overlap.loc[no_overlap["dmr_id"] == "prom_hyper", "missing_evidence_flags"].iloc[0]


def test_join_expression_evidence_without_expression_is_nonfatal():
    linked = link_dmrs_to_genes(_dmrs(), _genes(), promoter_upstream=1000, promoter_downstream=200)
    joined = join_expression_evidence(linked, None)
    assert joined["expression_direction"].eq("no_expression_data").all()
    assert not joined["expression_supported"].any()
