# Plant DMR Validation Run Template

This is the reusable terminal template for plant WGBS/RRBS DMR validation runs.
It follows the current framework UX:

1. optional annotation package normalization
2. DMR caller outputs
3. controlled GLM/GLMM confirmation
4. validation checks
5. optional expression/RNA-seq support normalization
6. summary reports and plots
7. optional LaTeX-ready figure folder

The template is dataset-agnostic. It assumes the dataset root contains:

```text
<ROOT>/
  raw_cx/
    <sample>.CX_report.txt.gz
  metadata/
    sample_manifest.tsv
```

`sample_manifest.tsv` must contain at least:

```text
sample_id    file    condition
```

`condition` should contain exactly two groups for the default run:

```text
control
treatment
```

## 0. PowerShell Setup

```powershell
cd C:\Users\sysue\bsx\bsxplorer2

$env:PYTHONPATH="$PWD\python\src;$PWD"

# Change these per plant/dataset.
$ROOT="G:\bsx2_raw\tomato_gse190780_6sample"
$RSCRIPT="C:\Program Files\R\R-4.6.0\bin\Rscript.exe"

# Optional. Use only if metilene is installed and available from Windows.
$METILENE=""

# One run id keeps all output folders synchronized.
$RUN_ID="500k_v1"
$MAX_SITES=500000
$EXTERNAL="$ROOT\external_callers_combined_$RUN_ID"
```

Quick input check:

```powershell
Get-ChildItem "$ROOT\raw_cx"
Get-Content "$ROOT\metadata\sample_manifest.tsv" -TotalCount 10
& $RSCRIPT -e "pkgs <- c('data.table','matrixStats','DSS','methylKit','bsseq','BiocParallel','GenomicRanges','SummarizedExperiment','glmmTMB'); x <- sapply(pkgs, requireNamespace, quietly=TRUE); print(x); if (any(!x)) quit(status=1)"
```

## 0b. Current Tomato GSE190780 Run

Use this concrete setup for the current tomato benchmark where caller outputs
already exist in `external_callers_combined_500k_v3`.

```powershell
cd C:\Users\sysue\bsx\bsxplorer2

$env:PYTHONPATH="$PWD\python\src;$PWD"

$ROOT="G:\bsx2_raw\tomato_gse190780_6sample"
$RSCRIPT="C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
$RUN_ID="tomato_500k_final_no_dmrseq"
$EXTERNAL="$ROOT\external_callers_combined_500k_v3"
```

Then run sections 1, 3, and 6 below. Section 2 is only needed when caller
outputs must be recomputed from raw CX reports.

## 1. Annotation Package

Use this before biological downstream analysis, metagene projection, or gene/TE overlap.
For tomato GSE190780/SL3.0, this downloads SGN ITAG3.2 gene models, RepeatModeler repeat regions, and GO annotations, then writes the common framework format:

```powershell
python -B -m dmr_validation_framework.cli workflow plant-annotation -- `
  --root "$ROOT" `
  --profile tomato_itag3_2
```

Main outputs:

```text
<ROOT>\annotation_itag3_2\
  tomato_itag3_2_genes.normalized.gff3
  metagene_gene_regions.bed
  metagene_gene_regions_qc.tsv
  seqname_mapping.tsv
  te_regions.bed
  te_regions_qc.tsv
  go_annotations.tsv
  annotation_manifest.json
  annotation_preflight_report.md
```

## 1b. Expression / RNA-Seq Support

This layer is optional, but it is required for `fig5_29_expression_support_levels`.
The framework accepts any gene-level expression table, including RNA-seq TPM,
FPKM, normalized counts, or array expression, as long as it has a gene column.

For RNA-seq or any already gene-level expression matrix:

```powershell
python -B -m dmr_validation_framework.cli workflow expression-support -- `
  --root "$ROOT" `
  --expression-table "path\to\gene_expression_tpm_or_counts.tsv" `
  --gene-id-column "gene_id" `
  --out-dir "$ROOT\expression_support_$RUN_ID"
```

Main output:

```text
<ROOT>\expression_support_<RUN_ID>\expression_support.tsv
```

The input can be TSV/CSV/TXT, Parquet/Feather, or Excel (`.xlsx/.xls`) if
`pandas` can read it in the current environment.

For tomato `GSE190782`, the public transcriptome subseries is `GSE190779`.
It is not RNA-seq; it is Affymetrix Tomato Genome Array expression. GEO reports
RMA-normalized expression values in the series matrix and uses platform `GPL4741`.
This can still be normalized as expression evidence, but `GPL4741` uses older
probe/transcript identifiers. If the platform annotation does not contain Solyc
gene IDs, the command will produce a WARN summary and `fig5_29` will remain
skipped until a probe-to-Solyc mapping is provided.

```powershell
New-Item -ItemType Directory -Force "$ROOT\expression_gse190779" | Out-Null

Invoke-WebRequest `
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190779/matrix/GSE190779_series_matrix.txt.gz" `
  -OutFile "$ROOT\expression_gse190779\GSE190779_series_matrix.txt.gz"

Invoke-WebRequest `
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4741/annot/GPL4741.annot.gz" `
  -OutFile "$ROOT\expression_gse190779\GPL4741.annot.gz"

python -B -m dmr_validation_framework.cli workflow expression-support -- `
  --root "$ROOT" `
  --series-matrix "$ROOT\expression_gse190779\GSE190779_series_matrix.txt.gz" `
  --platform-annot "$ROOT\expression_gse190779\GPL4741.annot.gz" `
  --out-dir "$ROOT\expression_gse190779"
```

If you have a probe-to-gene mapping table:

```powershell
python -B -m dmr_validation_framework.cli workflow expression-support -- `
  --root "$ROOT" `
  --series-matrix "$ROOT\expression_gse190779\GSE190779_series_matrix.txt.gz" `
  --probe-gene-map "path\to\probe_to_solyc.tsv" `
  --out-dir "$ROOT\expression_gse190779"
```

For another plant, use the same command shape with local annotation inputs:

```powershell
python -B -m dmr_validation_framework.cli workflow plant-annotation -- `
  --root "$ROOT" `
  --profile custom `
  --species "Species name" `
  --assembly "assembly_name" `
  --annotation-version "annotation_version" `
  --output-prefix "species_annotation" `
  --gene-gff "path\to\genes.gff3" `
  --te-gff "path\to\repeats.gff3" `
  --go-table "path\to\go.tsv"
```

## 2. DMR + GLM/GLMM + Validation + Reports

Use this as the default command for a practical limited benchmark.

```powershell
python -B -m dmr_validation_framework.cli workflow run-pipeline -- `
  --root "$ROOT" `
  --rscript "$RSCRIPT" `
  --external-root "$ROOT\external_callers_combined_$RUN_ID" `
  --glm-glmm-out-dir "$ROOT\glm_glmm_methylkit_top30_$RUN_ID" `
  --validation-out-dir "$ROOT\validation_core_$RUN_ID" `
  --report-dir "$ROOT\dmr_validation_report_$RUN_ID" `
  --contexts "CG,CHG,CHH" `
  --callers "DSS,methylKit,BSmooth,comb-p,metilene" `
  --profile core `
  --q-threshold 0.10 `
  --min-coverage 5 `
  --min-sites 3 `
  --max-gap 1000 `
  --max-sites-per-context $MAX_SITES `
  --read-chunk-size 250000 `
  --glm-glmm-caller methylKit `
  --glm-glmm-top-n-per-context 10 `
  --skip-harmonize
```

If `metilene` is available from Windows, add:

```powershell
  --metilene "$METILENE"
```

If `metilene` is not available, either remove it from `--callers`:

```powershell
--callers "DSS,methylKit,BSmooth,comb-p"
```

or run metilene separately in WSL and then rerun the reporting pipeline from existing outputs.

Do not add `dmrseq` to the default command unless runtime is acceptable. If needed:

```powershell
--callers "DSS,methylKit,dmrseq,BSmooth,comb-p,metilene" `
--dmrseq-permutations 5
```

## 3. Reuse Existing Caller Outputs

Use this when caller outputs already exist and only GLM/GLMM, validation and reports should be rebuilt.

```powershell
python -B -m dmr_validation_framework.cli workflow run-pipeline -- `
  --root "$ROOT" `
  --rscript "$RSCRIPT" `
  --external-root "$ROOT\external_callers_combined_$RUN_ID" `
  --glm-glmm-out-dir "$ROOT\glm_glmm_methylkit_top30_$RUN_ID" `
  --validation-out-dir "$ROOT\validation_core_$RUN_ID" `
  --report-dir "$ROOT\dmr_validation_report_$RUN_ID" `
  --contexts "CG,CHG,CHH" `
  --callers "DSS,methylKit,BSmooth,comb-p,metilene" `
  --profile core `
  --q-threshold 0.10 `
  --min-coverage 5 `
  --glm-glmm-caller methylKit `
  --glm-glmm-top-n-per-context 10 `
  --skip-callers `
  --skip-harmonize
```

This is the same pattern used for the tomato run after the caller files were already present.

## 4. GLM/GLMM Only

Use this when caller summary is done and only the confirmatory GLMM layer needs to be rebuilt.

```powershell
python -B -m dmr_validation_framework.cli workflow glm-glmm-validation -- `
  --root "$ROOT" `
  --external-root "$ROOT\external_callers_combined_$RUN_ID" `
  --out-dir "$ROOT\glm_glmm_methylkit_top30_$RUN_ID" `
  --rscript "$RSCRIPT" `
  --contexts "CG,CHG,CHH" `
  --glm-caller methylKit `
  --top-n-per-context 10 `
  --q-threshold 0.10 `
  --min-total 5 `
  --min-cpg 3
```

For a quick smoke test:

```powershell
--top-n-per-context 1
```

## 5. Reports Only

Caller summary:

```powershell
python -B -m dmr_validation_framework.cli report caller-summary -- `
  --external-root "$ROOT\external_callers_combined_$RUN_ID" `
  --out-dir "$ROOT\dmr_validation_report_$RUN_ID\caller_summary" `
  --contexts "CG,CHG,CHH" `
  --callers "DSS,methylKit,BSmooth,comb-p,metilene"
```

Validation summary:

```powershell
python -B -m dmr_validation_framework.cli report validation-summary -- `
  --validation-dir "$ROOT\validation_core_$RUN_ID" `
  --out-dir "$ROOT\dmr_validation_report_$RUN_ID\validation_summary"
```

GLM/GLMM plots:

```powershell
python -B -m dmr_validation_framework.cli report glmm-plots -- `
  --comparison "$ROOT\glm_glmm_methylkit_top30_$RUN_ID\glm_vs_glmm\glm_vs_glmm_comparison.tsv" `
  --glmm-results "$ROOT\glm_glmm_methylkit_top30_$RUN_ID\glmm\cpg_level_glmm_results.tsv" `
  --region-cpg-counts "$ROOT\glm_glmm_methylkit_top30_$RUN_ID\region_cpg_counts.tsv" `
  --design "$ROOT\glm_glmm_methylkit_top30_$RUN_ID\sample_design.tsv" `
  --out-dir "$ROOT\glm_glmm_methylkit_top30_$RUN_ID\glmm_plots"
```

## 6. LaTeX-Ready Folder

Build a generic `thesis_figures_final/for_latex` bundle from caller reports,
validation outputs, GLM/GLMM plots, optional plant annotation, optional
expression support, and lightweight extended metagene/random-control outputs:

```powershell
python -B -m dmr_validation_framework.cli workflow thesis-figures -- `
  --root "$ROOT" `
  --external-root "$EXTERNAL" `
  --report-dir "$ROOT\dmr_validation_report_$RUN_ID" `
  --validation-dir "$ROOT\validation_core_$RUN_ID" `
  --glm-glmm-dir "$ROOT\glm_glmm_methylkit_top30_$RUN_ID" `
  --annotation-dir "$ROOT\annotation_itag3_2" `
  --out-dir "$ROOT\thesis_figures_final_$RUN_ID" `
  --figure-prefix "fig_dmr" `
  --random-control-iterations 100 `
  --max-random-regions-per-iteration 5000
```

For the current tomato benchmark from section `0b`, use the explicit existing
caller-output path and a tomato-specific figure prefix:

```powershell
python -B -m dmr_validation_framework.cli workflow thesis-figures -- `
  --root "$ROOT" `
  --external-root "$EXTERNAL" `
  --report-dir "$ROOT\dmr_validation_report_$RUN_ID" `
  --validation-dir "$ROOT\validation_core_$RUN_ID" `
  --glm-glmm-dir "$ROOT\glm_glmm_methylkit_top30_$RUN_ID" `
  --annotation-dir "$ROOT\annotation_itag3_2" `
  --out-dir "$ROOT\thesis_figures_final_$RUN_ID" `
  --figure-prefix "fig_tomato_dmr" `
  --random-control-iterations 100 `
  --max-random-regions-per-iteration 5000
```

If expression support was normalized into a non-default folder, pass it
explicitly:

```powershell
--expression-table "$ROOT\expression_support_$RUN_ID\expression_support.tsv"
```

Main outputs:

```text
<ROOT>\thesis_figures_final_<RUN_ID>\
  for_latex\
    *.pdf
  figures_pdf\
    *.pdf
  figures_png\
    *.png
  tables\
    thesis_figure_manifest.tsv
  reports\
    thesis_figure_bundle.html
```

The bundle creates the same reusable Chapter 5 slots when inputs exist:

```text
fig5_08  overlap threshold sensitivity
fig5_09  common-CpG QC
fig5_10  region CpG/coverage distribution
fig5_11  replicate vs pooled delta sensitivity
fig5_12  bootstrap CI forest plot
fig5_13  delta vs coverage
fig5_22  observed vs random occupancy profile
fig5_23  center vs interval-overlap projection sensitivity
fig5_24  top20 DMR-linked loci TE status
fig5_25  top50 TE-like composition
fig5_26  TE enrichment summary
fig5_27  CG/CHG/CHH context distribution
fig5_28  feature class distribution
fig5_29  expression support levels
fig5_30  DMR-linked GO mapping coverage limitation
```

If a required input is unavailable, the figure is not silently omitted: it is
listed as `skipped` in:

```text
<ROOT>\thesis_figures_final_<RUN_ID>\tables\thesis_figure_manifest.tsv
```

`fig5_22` and `fig5_23` are now generated directly by `workflow thesis-figures`
when caller DMR rows and `annotation_itag3_2\metagene_gene_regions_qc.tsv` are
available. The intermediate tables are written to:

```text
<ROOT>\thesis_figures_final_<RUN_ID>\tables\extended_metagene\
  mapped_dmr_centers.tsv
  random_control_bin_empirical_p.tsv
  random_control_density_by_section.tsv
  random_control_occupancy_summary.tsv
  projection_sensitivity_density.tsv
  projection_sensitivity_summary.tsv
```

`fig5_29` requires `expression_support.tsv` with a `gene_id`, `gene`, `id`, or
`locus` column. If no gene-level expression support is available, it is listed
as `skipped` in the figure manifest rather than silently omitted.

This generic bundle can build validation, GLM/GLMM, annotation, TE/repeat, GO
coverage, and seqname-compatibility figures. For a thesis-quality biological
chapter, organism-specific evidence can still be added separately:

- gene annotation
- TE annotation
- GO/enrichment
- RNA-seq/expression support
- literature interpretation

These are not part of the generic DMR validation framework.

## 7. Acceptance Checks

```powershell
Get-ChildItem "$ROOT\external_callers_combined_$RUN_ID"
Get-ChildItem "$ROOT\glm_glmm_methylkit_top30_$RUN_ID\glm_vs_glmm"
Get-ChildItem "$ROOT\validation_core_$RUN_ID"
Get-ChildItem "$ROOT\dmr_validation_report_$RUN_ID"
Get-ChildItem "$ROOT\thesis_figures_final_$RUN_ID\for_latex" -ErrorAction SilentlyContinue
```

Minimal expected files:

```text
external_callers_combined_<RUN_ID>/external_caller_run_status.tsv
external_callers_combined_<RUN_ID>/methylkit/methylkit_dmrs_CG.tsv
glm_glmm_methylkit_top30_<RUN_ID>/glm_vs_glmm/glm_vs_glmm_comparison.tsv
validation_core_<RUN_ID>/glm_glmm_status_summary.tsv
dmr_validation_report_<RUN_ID>/caller_summary/caller_summary.html
dmr_validation_report_<RUN_ID>/validation_summary/validation_summary.html
```

## 8. Multi-Plant Batch Template

```powershell
$DATASETS = @(
  @{
    Name="tomato_gse190780"
    Root="G:\bsx2_raw\tomato_gse190780_6sample"
    MaxSites=500000
  },
  @{
    Name="oryza_gse202715"
    Root="G:\bsx2_raw\oryza_gse202715_full_final_workflow"
    MaxSites=500000
  }
)

foreach ($D in $DATASETS) {
  $ROOT=$D.Root
  $RUN_ID="500k_v1"
  $MAX_SITES=$D.MaxSites

  python -B -m dmr_validation_framework.cli workflow run-pipeline -- `
    --root "$ROOT" `
    --rscript "$RSCRIPT" `
    --external-root "$ROOT\external_callers_combined_$RUN_ID" `
    --glm-glmm-out-dir "$ROOT\glm_glmm_methylkit_top30_$RUN_ID" `
    --validation-out-dir "$ROOT\validation_core_$RUN_ID" `
    --report-dir "$ROOT\dmr_validation_report_$RUN_ID" `
    --contexts "CG,CHG,CHH" `
    --callers "DSS,methylKit,BSmooth,comb-p,metilene" `
    --profile core `
    --q-threshold 0.10 `
    --min-coverage 5 `
    --min-sites 3 `
    --max-gap 1000 `
    --max-sites-per-context $MAX_SITES `
    --read-chunk-size 250000 `
    --glm-glmm-caller methylKit `
    --glm-glmm-top-n-per-context 10 `
    --skip-harmonize
}
```

Use `--skip-callers` in the batch command if caller outputs already exist.
