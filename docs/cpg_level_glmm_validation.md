# CpG-level GLMM Confirmatory Validation

## Purpose

The CpG-level GLMM layer is an optional confirmatory validation step for
selected DMR candidates. It strengthens the existing aggregated beta-binomial
validation by modelling observations at:

```text
region_id x cpg_id x sample_id
```

It does not replace the main DMR statistical model, and it is not a
genome-wide DMR caller.

## Why Aggregated GLM Is Useful But Limited

The existing aggregated validation works on region-level sample summaries:

```text
region_id x sample_id
```

This is efficient and useful for broad validation, but it collapses CpG-level
heterogeneity. A region can appear strong because of a local high-coverage or
high-effect CpG, even when other CpGs inside the interval are less consistent.

## What CpG-level GLMM Adds

For top-N DMR candidates, the confirmatory GLMM uses:

```text
cbind(mC, uC) ~ condition + covariates + (1 | cpg_id) + (1 | sample_id)
```

and compares it with:

```text
cbind(mC, uC) ~ covariates + (1 | cpg_id) + (1 | sample_id)
```

with `glmmTMB::betabinomial(link = "logit")`.

This adds:

- individual CpG observations;
- CpG-specific baseline variation;
- sample-level variation;
- better resistance to local outlier CpGs;
- a stricter confirmatory check for top candidates.

## Architecture

- Rust extracts per-CpG counts from indexed `.bsx` files for predefined
  regions.
- Python orchestrates selection, QC, reporting, and comparison.
- R/glmmTMB fits the confirmatory GLMM when available.

If `Rscript` or `glmmTMB` is unavailable, the workflow writes
`model_status = glmmTMB_unavailable` and does not crash.

## Outputs

The per-CpG extraction layer writes:

- `region_cpg_counts.tsv`
- `region_cpg_counts_summary.md`
- `region_cpg_counts_manifest.json`
- `warnings.tsv`

The GLMM layer writes:

- `cpg_level_glmm_results.tsv`
- `cpg_level_glmm_summary.md`
- `cpg_level_glmm_manifest.json`
- `cpg_level_glmm_warnings.tsv`

The comparison layer writes:

- `glm_vs_glmm_comparison.tsv`

## Interpretation

Use this wording in thesis text:

> CpG-level GLMM confirmatory validation was applied to selected DMR candidates
> to assess whether region-level methylation differences remained supported
> after accounting for CpG-specific baseline and sample-level variation.

Avoid wording that implies:

- automatic genome-wide DMR calling;
- proof of biological causality;
- replacement of DSS, dmrseq, methylKit, or metilene;
- proof of gene function.

## Limitations

- It is computationally heavier than aggregated validation.
- It requires enough CpGs and at least the configured number of replicates per
  condition.
- Model convergence can fail.
- It is run only for selected top-N candidates.
- It does not run raw FASTQ/Bismark processing.
- It does not change frozen DMR/evidence artifacts.
