# Thesis Code Overview

## Purpose of Thesis Code Layer

The thesis code layer documents and hardens additive BSX2 components around
regional methylation evidence. It keeps reusable logic path-free and separates
stable package modules from CLI wrappers and runtime artifacts.

## Implemented Stable Modules

- `bsx2.analysis.seqname_harmonization`: arbitrary seqname alias checks for non-model organisms; no liftover.
- `bsx2.analysis.assembly_compatibility`: coordinate bounds diagnostics; no assembly-equivalence inference.
- `bsx2.analysis.region_signal`: interval aggregation for genes, promoters, BED, DMR and external regions with `backend="auto|rust|pandas"`.
- `bsx2.analysis.region_signal_rust`: Rust-backed indexed region count aggregation wrapper for high-throughput RegionSignal use.
- `bsx2.analysis.region_cpg_counts`: Rust-backed per-CpG count extraction for predefined regions used by confirmatory validation.
- `bsx2.viz.compute.discrete_region_cache`: compressed cache for extracted regional point data.

## Implemented Validation / Experimental Modules

- Regional Evidence Model and beta-binomial aggregated GLM validation remain separate evidence layers.
- `bsx2.analysis.cpg_level_glmm`: optional CpG-level GLMM confirmatory validation for selected top-N DMR candidates.
- `bsx2.analysis.glmm_comparison`: compares aggregated GLM support with CpG-level GLMM confirmation status.
- `bsx2.analysis.dmr_harmonization` and `external_dmr_callers` standardize external candidate schemas without running external callers.
- `bsx2.analysis.replicate_diagnostics` is diagnostic only, not a CpG-level GLMM.
- `bsx2.analysis.bsseq_qc` imports available QC reports and count summaries without estimating conversion failure unless spike-in/report evidence exists.

## Runtime Artifacts and Freezes

Runtime thesis artifacts live outside reusable code under project-specific result
directories. Frozen artifacts are preserved and are not recalculated by cleanup.

## What Is Not Claimed

This layer does not implement a genome-wide CpG-level DMR caller, dispersion
shrinkage, production DSS/methylKit backends, TE/pericentromere analysis, or
causal methylation-expression inference.

## Reproduce Smoke Checks

Use the commands recorded in
`code_cleanup_for_thesis/final_validation_commands.sh`. If pytest is not
available, the direct synthetic smoke checks provide a local fallback.

## Thesis Support

The modules support thesis chapters by making the architecture explicit:
seqname-safe non-model workflows, production-oriented Rust aggregation path, reusable display
caches, external candidate harmonization, optional CpG-level GLMM confirmatory
validation for selected DMRs, and documented limitations.
