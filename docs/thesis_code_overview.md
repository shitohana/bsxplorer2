# Thesis Code Overview

## Purpose of Thesis Code Layer

The thesis code layer documents and hardens additive BSX2 components around
regional methylation evidence. It keeps reusable logic path-free and separates
stable package modules from CLI wrappers and runtime artifacts.

## Implemented Stable Modules

- `bsx2.analysis.seqname_harmonization`: arbitrary seqname alias checks for non-model organisms; no liftover.
- `bsx2.analysis.assembly_compatibility`: coordinate bounds diagnostics; no assembly-equivalence inference.
- `bsx2.analysis.region_signal`: pandas MVP interval aggregation for genes, promoters, BED, DMR and external regions.
- `bsx2.viz.compute.discrete_region_cache`: compressed cache for extracted regional point data.

## Implemented Validation / Experimental Modules

- Regional Evidence Model and beta-binomial aggregated GLM validation remain separate evidence layers.
- `bsx2.analysis.dmr_harmonization` and `external_dmr_callers` standardize external candidate schemas without running external callers.
- `bsx2.analysis.replicate_diagnostics` is diagnostic only, not a CpG-level GLMM.
- `bsx2.analysis.bsseq_qc` imports available QC reports and count summaries without estimating conversion failure unless spike-in/report evidence exists.

## Runtime Artifacts and Freezes

Runtime thesis artifacts live outside reusable code under project-specific result
directories. Frozen artifacts are preserved and are not recalculated by cleanup.

## What Is Not Claimed

This layer does not implement full CpG-level beta-binomial GLMM, dispersion
shrinkage, production DSS/methylKit backends, TE/pericentromere analysis, or
methylation-expression integration.

## Reproduce Smoke Checks

Use the commands recorded in
`code_cleanup_for_thesis/final_validation_commands.sh`. If pytest is not
available, the direct synthetic smoke checks provide a local fallback.

## Thesis Support

The modules support thesis chapters by making the architecture explicit:
seqname-safe non-model workflows, interval-based aggregation, reusable display
caches, external candidate harmonization, and documented limitations.
