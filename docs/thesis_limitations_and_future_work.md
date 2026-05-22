# Thesis Limitations and Future Work

## Implemented

- Seqname harmonization for arbitrary contigs and scaffolds.
- Assembly coordinate checks with warning-only reports.
- RegionSignal API for unified interval aggregation.
- Rust-backed `.bsx` RegionSignal count aggregation path with pandas fallback.
- DiscreteRegionData cache for reuse after extraction.
- Regional Evidence Model.
- Beta-binomial aggregated GLM validation.
- External DMR harmonization through canonical schema.
- DSS/methylKit candidate-level benchmark artifacts.
- Replicate diagnostics.
- BS-seq QC importers.

## Partially Implemented

- External caller benchmark is a controlled candidate-level demonstration, not a
  production genome-wide backend.
- Beta-binomial validation is an aggregated GLM layer, not a CpG-level
  random-effect GLMM.
- QC importers do not estimate conversion failure without spike-in or explicit
  conversion reports.

## Future Work

- CpG-level GLMM.
- Dispersion shrinkage.
- Production external caller backend.
- TE/pericentromere analysis.
- Methylation-expression integration.
- Full production-scale Rust RegionSignal benchmark on multiple large `.bsx`
  files.
- Full CLI tutorials.
