"""Public analysis helpers added around the BSX2 thesis evidence workflow.

The package exports lightweight, additive utilities for seqname harmonization,
interval signal aggregation, assembly diagnostics, replicate diagnostics,
BS-seq QC import, and external DMR schema harmonization. These modules do not
run raw pipelines or external DMR callers.
"""

from .assembly_compatibility import check_coordinate_bounds, compare_coordinate_sources, normalize_coordinate_table, read_genome_sizes, write_assembly_compatibility_report
from .bsseq_qc import classify_bsseq_qc, compute_counts_qc, parse_bismark_alignment_report, parse_bismark_dedup_report, parse_bismark_mbias_report, write_bsseq_qc_outputs
from .dmr_harmonization import assign_missing_dmr_ids, build_caller_support_matrix, canonical_dmr_columns, normalize_dmr_coordinates, validate_canonical_dmr_schema, write_schema_validation_report
from .external_dmr_callers import DMRseqAdapter, DSSAdapter, GenericBedAdapter, MetileneAdapter, MethylKitAdapter, adapter_for_caller
from .region_signal import RegionSignalConfig, aggregate_region_signal, available_region_signal_backends, normalize_counts_table, normalize_region_table, write_region_signal_qc, write_region_signal_table
from .region_signal_rust import aggregate_region_signal_rust, rust_region_aggregator_available
from .replicate_diagnostics import compute_pairwise_replicate_consistency, compute_region_replicate_summary, write_replicate_diagnostics_outputs
from .seqname_harmonization import apply_seqname_aliases, compare_seqname_sets, normalize_seqname, read_seqname_aliases, validate_seqname_compatibility, write_seqname_compatibility_report

__all__ = [name for name in list(globals()) if not name.startswith("_")]
