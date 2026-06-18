"""Public analysis helpers added around the BSX2 thesis evidence workflow.

The package exports lightweight, additive utilities for seqname harmonization,
interval signal aggregation, assembly diagnostics, replicate diagnostics,
BS-seq QC import, and external DMR schema harmonization. These modules do not
run raw pipelines or external DMR callers.
"""

from .assembly_compatibility import check_coordinate_bounds, compare_coordinate_sources, normalize_coordinate_table, read_genome_sizes, write_assembly_compatibility_report
from .bsseq_qc import classify_bsseq_qc, compute_counts_qc, parse_bismark_alignment_report, parse_bismark_dedup_report, parse_bismark_mbias_report, write_bsseq_qc_outputs
from .dmr_harmonization import (
    CALLER_MODEL_FAMILY,
    TIER_DESCRIPTIONS,
    add_dmr_tiers,
    assign_final_dmr_tier,
    assign_missing_dmr_ids,
    build_caller_support_matrix,
    build_tier_reasons,
    caller_model_family,
    canonical_dmr_columns,
    compute_multi_caller_evidence_score,
    compute_validation_robustness_score,
    normalize_dmr_coordinates,
    validate_canonical_dmr_schema,
    write_schema_validation_report,
)
from .dmr_functional_prioritization import add_functional_prioritization_scores, join_expression_evidence, link_dmrs_to_chromatin, link_dmrs_to_genes, link_dmrs_to_te, prioritize_dmr_linked_genes, read_chromatin_peaks, read_dmr_evidence_table, read_expression_table, read_gene_annotation, read_te_annotation, write_dmr_functional_prioritization_table
from .cpg_level_glmm import glmmTMB_available, read_design_table, read_dmr_evidence_for_glmm, read_region_cpg_counts, rscript_available, run_cpg_level_glmm_validation
from .coverage_set_qc import aggregate_counts_common_cpg, aggregate_counts_per_sample, build_cpg_coverage_qc, compare_per_sample_vs_common_sets, normalize_cpg_counts_table
from .external_dmr_callers import (
    ADAPTERS,
    BSmoothAdapter,
    BiSeqAdapter,
    CombPAdapter,
    DMRcateAdapter,
    DMRseqAdapter,
    DSSAdapter,
    GenericBedAdapter,
    HmmdmAdapter,
    MethyLassoAdapter,
    MetileneAdapter,
    MethylKitAdapter,
    MethylSigAdapter,
    MOABSAdapter,
    RADMethAdapter,
    adapter_for_caller,
)
from .glmm_comparison import compare_beta_binomial_glm_vs_glmm, read_cpg_level_glmm_results, read_region_level_glm_results
from .region_cpg_counts import extract_region_cpg_counts, extract_region_cpg_counts_pandas, extract_region_cpg_counts_rust, rust_region_cpg_extractor_available
from .region_signal import RegionSignalConfig, aggregate_region_signal, available_region_signal_backends, normalize_counts_table, normalize_region_table, write_region_signal_qc, write_region_signal_table
from .region_signal_rust import aggregate_region_signal_rust, rust_region_aggregator_available
from .replicate_diagnostics import compute_pairwise_replicate_consistency, compute_region_replicate_summary, write_replicate_diagnostics_outputs
from .seqname_harmonization import apply_seqname_aliases, compare_seqname_sets, normalize_seqname, read_seqname_aliases, validate_seqname_compatibility, write_seqname_compatibility_report

__all__ = [name for name in list(globals()) if not name.startswith("_")]
