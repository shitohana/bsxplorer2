# DMR Visualization Curves

## Purpose

The DMR visualization layer represents DMRs as ordinary genomic intervals plus optional region/sample methylation counts. It creates reusable curve specifications that can be saved as JSON, rendered into summary tables or figures, and embedded into Streamlit or another dashboard without rewriting visualization logic.

This layer is not a new DMR caller and does not replace BSX2 Regional Evidence or beta-binomial validation.

## Relation To Existing BSX2 Visualization

The curve layer is intentionally lightweight. It reuses existing BSX2 concepts where possible:

- Line/meta-region profiles use the existing `DiscreteRegionData` and `bsx2.viz.render.line.line_plot` path when positional methylation signal is available.
- Heatmap-style DMR methylation matrices are compatible with the existing heatmap layer and can be rendered as matrix tables or simple figures.
- Box and violin summaries mirror the existing `bsx2.viz.render.box` and `bsx2.viz.render.violin` distribution views.
- Chromosome-level DMR distributions are represented as interval counts per chromosome and can be connected to chromosome-level renderers.
- PCA uses the same sample-by-feature matrix idea as dimensional reduction views, but stores a portable result table for dashboard use.

## Input Files

Common inputs include:

- `dmr_evidence_scores.tsv`
- `dmr_regions.tsv`
- `dmr_region_count_tests.tsv`
- `region_counts_top1000_extended.tsv`
- `experimental_design_samples_extended.tsv`
- `dmr_beta_binomial_tests.tsv`
- `dmr_beta_binomial_tests_top1000_adjusted.tsv`
- `external_caller_support_matrix.tsv`
- `dmr_caller_support_matrix.tsv`
- annotation or enrichment TSV files

Column aliases are normalized for common names such as `chr/chrom/seqname`, `start/start_bp`, `end/end_bp`, `q/fdr/q_value`, and `delta/mean_delta`.

## Curve Types

- `dmr_line`: DMR-centered or scaled methylation profile. Requires positional methylation signal or precomputed `DiscreteRegionData`.
- `dmr_chromosome`: DMR count distribution by chromosome.
- `dmr_box`: Box plot table for delta, q-value, or methylation by condition.
- `dmr_violin`: Violin plot table for delta or methylation distributions.
- `dmr_pca`: PCA of samples using a DMR region methylation matrix from `Y / m`.
- `dmr_heatmap`: DMR x sample methylation matrix.
- `dmr_volcano`: delta versus `-log10(q_value)`.
- `dmr_evidence_bar`: evidence class counts.
- `dmr_caller_support`: external caller support distribution.

## Example CLI

```bash
python scripts/build_dmr_curve_bundle.py \
  --dmr-table path/to/dmr_evidence_scores.tsv \
  --region-counts path/to/region_counts_top1000_extended.tsv \
  --design path/to/experimental_design_samples_extended.tsv \
  --caller-support path/to/external_caller_support_matrix.tsv \
  --out-dir dmr_curves/
```

The output directory contains:

- `dmr_curve_bundle.json`
- `dmr_curve_manifest.tsv`
- `dmr_curve_bundle_summary.md`
- `specs/*.json`
- `tables/*.tsv`
- `figures/*.png` when rendering dependencies are available
- `warnings.tsv`

## Limitations

- A DMR table alone supports summary plots, not methylation profiles.
- Line/meta-DMR curves require positional methylation signal or `DiscreteRegionData`.
- PCA and heatmap curves require `region_counts` plus a sample design table.
- Curve specifications are visualization metadata, not statistical results.
- External caller curves display harmonized support but do not standardize caller-specific statistical assumptions.
