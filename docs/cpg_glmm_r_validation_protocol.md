# Reproducible CpG-level GLMM Validation Protocol

This protocol records how the optional CpG-level GLMM confirmatory validation
layer is validated with an actual R/glmmTMB environment. It is a validation
bundle only: it does not change the DMR statistical model, does not run raw
FASTQ/Bismark processing, and does not run DSS, methylKit, dmrseq, or metilene.

## Output Location

Runtime validation artifacts are written under:

```text
/mnt/g/bsx2_raw/bnapus_gse202609_6sample/cpg_level_glmm_r_validation/
```

The reproducibility command file in that directory is:

```text
cpg_glmm_r_validation_commands.sh
```

## Environment Proof

The validation uses an isolated conda-style environment on `G:`:

```text
/mnt/g/bsx2_raw/bnapus_gse202609_6sample/cpg_level_glmm_r_validation/envs/bsx2-glmm-validation
```

The recorded environment specification is:

```yaml
name: bsx2-glmm-validation
channels:
  - conda-forge
dependencies:
  - r-base
  - r-glmmtmb
  - r-data.table
  - r-jsonlite
  - r-optparse
```

The preflight records:

- `Rscript` availability.
- `glmmTMB` availability.
- `glmmTMB::betabinomial(link = "logit")` availability.
- R session information.
- R package versions.
- source git commit hash.

If Rscript or glmmTMB are unavailable, the validation is not considered
successful. The protocol should report the missing dependency and preserve the
installation commands instead of silently falling back to a non-GLMM check.

## Synthetic Truth Cases

Synthetic input files are generated with fixed seeds:

- `synthetic/inputs/synthetic_region_cpg_counts.tsv`
- `synthetic/inputs/synthetic_design.tsv`
- `synthetic/inputs/synthetic_truth.tsv`

The design contains two conditions, `control` and `case`, with six replicates
per group. Synthetic regions cover five expected behaviors:

- `null_region`: no true methylation effect, expected not significant.
- `consistent_dmr`: many CpGs shift consistently, expected to be confirmed.
- `one_cpg_outlier`: aggregated signal is driven by a high-contribution CpG;
  expected to be GLM-only/GLMM-conservative or explicitly flagged by top-CpG
  contribution.
- `low_coverage_region`: expected to be skipped for coverage.
- `insufficient_cpg_region`: expected to be skipped for too few CpGs.

The synthetic GLMM run uses:

```text
scripts/run_cpg_level_glmm_validation.py
```

with `--min-cpg 3`, `--min-replicates-per-group 2`, case label `case`, and
control label `control`.

## GLM vs GLMM Comparison

The synthetic region-level GLM table is generated only as a compact validation
input for comparison. It is not a new DMR caller. Comparison is performed with:

```text
scripts/compare_beta_binomial_glm_vs_glmm.py
```

Expected checks:

- `consistent_dmr` is confirmed by CpG-level GLMM confirmatory validation.
- `null_region` remains non-significant.
- `one_cpg_outlier` is either not confirmed by GLMM or is clearly flagged by
  high top-CpG contribution.
- low-data regions are reported as insufficient rather than forced into a fit.

## Diagnostic Plots

Plots are generated with:

```text
scripts/plot_cpg_glmm_validation.py
```

Plot labels must use the wording "CpG-level GLMM confirmatory validation" and
must not label the method as genome-wide DMR calling.

## Real-data Smoke Policy

The real-data smoke is intentionally small and top-N only. It requires existing
`.bsx` methylation files because the smoke validates Rust/.bsx per-CpG
extraction plus R/glmmTMB fitting.

If no `.bsx` files are present under the B. napus sample root, the real-data
smoke is skipped with an explicit reason. Processed TSV tables must not be used
to fake the real-data extraction validation.

## Acceptance Outputs

The validation bundle writes:

- `r_glmmtmb_preflight.tsv`
- `r_session_info.txt`
- `r_package_versions.tsv`
- `synthetic/synthetic_cpg_level_glmm_results.tsv`
- `synthetic/synthetic_glm_vs_glmm_comparison.tsv`
- `synthetic/synthetic_truth_evaluation.tsv`
- diagnostic plots under `synthetic/`
- `real_bsx_candidate_inventory.tsv`
- `real/real_data_smoke_skipped.md` when no `.bsx` files are available
- `cpg_glmm_r_validation_report.md`
- `cpg_glmm_r_validation_acceptance.tsv`
- `cpg_glmm_r_validation_manifest.json`
- `cpg_glmm_r_validation_commands.sh`

Acceptance requires actual `model_status = ok` GLMM fits on synthetic regions,
expected synthetic truth behavior, explicit real-smoke status, and preserved
reproducibility commands.

## Thesis-safe Wording

Use:

```text
CpG-level GLMM was successfully executed in a reproducible R/glmmTMB
environment on controlled synthetic cases and, when available, small real DMR
candidates.
```

Do not claim that this layer is a genome-wide DMR caller, that it replaces
DSS/dmrseq-style analyses, or that it proves causal gene function.
