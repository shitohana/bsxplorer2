# BSX2 DMR Evidence Viewer

`apps/bsx2_dmr_viewer.py` is a lightweight dark-theme Streamlit frontend for inspecting already generated BSX2 DMR and evidence outputs. It is a viewer only and is intentionally separate from DMR calling and statistical modeling code.

## Purpose

The viewer helps inspect DMR regions, evidence scores, optional beta-binomial validation results, external caller support, and annotation or enrichment outputs from existing TSV files. It is suitable for demonstrations and thesis defense workflows where the statistical outputs have already been produced.

The UI uses one dark theme. Inputs live in the sidebar; the main page contains the guide, empty-state card, summary cards, plots, and tables. The sidebar file recommendations are shown as wrapping filename pills so long TSV names do not overflow the panel.

## How To Launch

```bash
streamlit run apps/bsx2_dmr_viewer.py
```

The app uses upload widgets in the sidebar by default. A local-path mode is available only when running the app locally and must be explicitly enabled in the sidebar with `Use local file paths instead of uploads`.

## What Files To Upload

### Required: Main DMR/evidence TSV

Use this required uploader for the primary region table. Best option:

- `dmr_evidence_scores.tsv`

Also supported:

- `dmr_regions.tsv`
- `dmr_region_count_tests.tsv`

This table drives region browsing, filtering, summary cards, and most plots.

### Optional: Beta-binomial validation TSV

Optional supporting input. Recommended files:

- `dmr_beta_binomial_tests.tsv`
- `dmr_beta_binomial_tests_top1000_adjusted.tsv`
- `dmr_beta_binomial_tests_top100_real.tsv`

This adds complementary beta-binomial confirmation metrics when matching columns are present.

### Optional: Caller support matrix TSV

Optional supporting input. Recommended files:

- `dmr_caller_support_matrix.tsv`
- `external_caller_support_matrix.tsv`

This shows whether external callers or internal BSX2 support a region.

### Optional: Annotation / enrichment TSV

Optional supporting input. Recommended files:

- `differential_methylated_features.tsv`
- `dmr_annotation_enrichment.tsv`
- `plant_region_methylation_summary.tsv`

This provides biological context such as promoter, gene body, intergenic, or enrichment summaries.

## Input Column Expectations

Column names are matched flexibly where possible:

- chromosome: `chrom`, `chr`, `chromosome`, `seqname`
- coordinates: `start`, `start_bp`, `begin`, `end`, `end_bp`, `stop`
- context: `context`, `methylation_context`
- delta: `delta`, `mean_delta`, `region_delta`, `beta_binom_delta`
- p-value: `p_value`, `p`, `pvalue`, `pval`, `region_p_value`
- q-value: `q_value`, `q`, `fdr`, `qvalue`, `qval`, `region_q_value`, `beta_binom_q_value`
- evidence class: `evidence_class`, `class`
- caller support: `n_callers_supporting`, `caller_support_count`

Missing optional columns do not stop the app. Related filters or plots are disabled and shown as not available.

## Dark Theme

The app injects local CSS to keep Streamlit's main container, sidebar, uploaders, text inputs, tabs, metric cards, dataframes, buttons, and alert boxes in one dark visual style. No global `.streamlit/config.toml` is required.

Some internal Streamlit widgets use version-specific markup, so very small native details may vary by Streamlit version, but the app targets stable `data-testid` selectors where possible.

## What The App Does Not Do

- no DMR calling;
- no raw FASTQ/BAM/Bismark processing;
- no DSS, methylKit, dmrseq, or metilene execution;
- no evidence recalculation;
- no uploaded-file persistence in the repository;
- no statistical-model replacement.

External harmonization standardizes schemas, not caller-specific statistical assumptions. Beta-binomial validation is a complementary evidence layer; aggregated mode is a GLM validation layer, not a random-effect GLMM.
