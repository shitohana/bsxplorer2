# Thesis Limitations and Future Work

## Implemented

- Seqname harmonization for arbitrary contigs and scaffolds.
- Assembly coordinate checks with warning-only reports.
- RegionSignal API for unified interval aggregation.
- Rust-backed indexed `.bsx` RegionSignal count aggregation path with pandas fallback.
- Rust-backed per-CpG count extraction for predefined regions.
- DiscreteRegionData cache for reuse after extraction.
- Regional Evidence Model.
- Beta-binomial aggregated GLM validation.
- Optional CpG-level GLMM confirmatory validation for selected top-N DMR candidates.
- External DMR harmonization through canonical schema.
- DSS/methylKit candidate-level benchmark artifacts.
- Replicate diagnostics.
- BS-seq QC importers.

## Partially Implemented

- External caller benchmark is a controlled candidate-level demonstration, not a
  production genome-wide backend.
- CpG-level GLMM validation is confirmatory and top-N only; it is not a
  genome-wide DMR caller and depends on R/glmmTMB availability.
- QC importers do not estimate conversion failure without spike-in or explicit
  conversion reports.
- Rust RegionSignal exposes `chunk_size` for future chunked implementations; the
  current `.bsx` path performs indexed per-region queries.

## Future Work

- Dispersion shrinkage.
- Production external caller backend.
- TE/pericentromere analysis.
- Methylation-expression integration.
- Full production-scale Rust RegionSignal benchmark on multiple large `.bsx`
  files.
- Full CLI tutorials.
