# DMR Validation Framework Terminal UX

This document shows the intended command-line workflow for running external DMR callers,
validating their outputs, and building summary reports.

## 1. End-to-End Pipeline

```powershell
cd C:\Users\sysue\bsx\bsxplorer2

$env:PYTHONPATH="$PWD\python\src;$PWD"
$ROOT="G:\bsx2_raw\tomato_gse190780_6sample"
$RSCRIPT="C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
# Optional. Required only if metilene should be executed, not just prepared.
$METILENE="C:\tools\metilene\metilene.exe"

python -B -m dmr_validation_framework.cli workflow run-pipeline -- `
  --root "$ROOT" `
  --rscript "$RSCRIPT" `
  --external-root "$ROOT\external_callers_limited_500k_v2" `
  --glm-glmm-out-dir "$ROOT\glm_glmm_limited_500k_v2" `
  --validation-out-dir "$ROOT\validation_core_500k_v2" `
  --report-dir "$ROOT\dmr_validation_report_500k_v2" `
  --contexts "CG,CHG,CHH" `
  --callers "DSS,methylKit,dmrseq,BSmooth,comb-p,metilene" `
  --profile core `
  --min-coverage 5 `
  --min-sites 3 `
  --max-gap 1000 `
  --q-threshold 0.10 `
  --max-sites-per-context 500000 `
  --read-chunk-size 250000 `
  --glm-glmm-caller methylKit `
  --glm-glmm-top-n-per-context 10 `
  --dmrseq-permutations 5 `
  --metilene "$METILENE"
```

Main outputs:

- `external_callers_limited_500k_v2/`: caller-native DMR outputs.
- `glm_glmm_limited_500k_v2/`: controlled top-N GLM-like vs CpG-level GLMM confirmation layer.
- `validation_core_500k_v2/`: framework validation checks.
- `dmr_validation_report_500k_v2/pipeline_manifest.json`: run manifest.
- `dmr_validation_report_500k_v2/caller_summary/caller_summary.html`: caller/context report.
- `dmr_validation_report_500k_v2/validation_summary/validation_summary.html`: validation report.

The default pipeline order is:

1. external callers
2. optional caller harmonization
3. controlled GLM/GLMM confirmation on top-N methylKit candidates
4. validation checks
5. HTML/TSV/PNG reports

Use `--skip-glm-glmm` when `glmmTMB` is unavailable or when you only need caller summaries.

## 2. Run Caller Stage Only

```powershell
python -B -m dmr_validation_framework.cli workflow run-external-callers -- `
  --root "$ROOT" `
  --external-root "$ROOT\external_callers_limited_500k_v2" `
  --rscript "$RSCRIPT" `
  --contexts "CG,CHG,CHH" `
  --callers "DSS,methylKit,dmrseq,BSmooth,comb-p,metilene" `
  --min-coverage 5 `
  --min-sites 3 `
  --max-gap 1000 `
  --q-threshold 0.10 `
  --max-sites-per-context 500000 `
  --read-chunk-size 250000 `
  --dmrseq-permutations 5 `
  --skip-harmonize
```

## 3. Build Reports From Existing Outputs

```powershell
python -B -m dmr_validation_framework.cli report caller-summary -- `
  --external-root "$ROOT\external_callers_limited_500k_v2" `
  --out-dir "$ROOT\dmr_validation_report_500k_v2\caller_summary" `
  --contexts "CG,CHG,CHH" `
  --callers "DSS,methylKit,dmrseq,BSmooth,comb-p,metilene"

python -B -m dmr_validation_framework.cli report validation-summary -- `
  --validation-dir "$ROOT\validation_core_500k_v2" `
  --out-dir "$ROOT\dmr_validation_report_500k_v2\validation_summary"
```

## 4. Run GLM/GLMM Confirmation From Existing Caller Outputs

This stage assumes the caller stage already produced methylKit DMRs and
`methylkit_input/<context>/*.bismark_coverage.tsv` count files.

```powershell
python -B -m dmr_validation_framework.cli workflow glm-glmm-validation -- `
  --root "$ROOT" `
  --external-root "$ROOT\external_callers_limited_500k_v2" `
  --out-dir "$ROOT\glm_glmm_limited_500k_v2" `
  --rscript "$RSCRIPT" `
  --contexts "CG,CHG,CHH" `
  --glm-caller methylKit `
  --top-n-per-context 10 `
  --q-threshold 0.10 `
  --min-total 5 `
  --min-cpg 3
```

Main outputs:

- `selected_methylkit_top10_per_context_regions.tsv`
- `methylkit_as_glm_results.tsv`
- `region_cpg_counts.tsv`
- `glmm/cpg_level_glmm_results.tsv`
- `glm_vs_glmm/glm_vs_glmm_comparison.tsv`
- `glmm_plots/*.png`

For a quick smoke test, use `--top-n-per-context 1`.

## 5. Interpretation Note

If `--max-sites-per-context` is non-zero, the run is a limited high-coverage benchmark.
It is suitable for framework validation, caller agreement, and robustness diagnostics,
but it is not a full genome-wide DMR catalogue.

`dmrseq` requires loci to have coverage in at least one sample from each condition.
The runner filters loci with zero coverage in all samples of either condition before
calling `dmrseq`.
Use `--dmrseq-permutations 5` or `10` for a controlled benchmark run. Use `0` to keep
the package default.

`metilene` requires an external binary. If `--metilene` is omitted and `metilene` is
not available on `PATH`, the runner writes `metilene_input/metilene_<context>.tsv`
and records `input_ready`; it does not produce final metilene DMRs.

The GLM/GLMM layer is also intentionally limited. It does not run a genome-wide GLMM
caller. It validates a controlled top-N set of predefined DMR candidates with
CpG-level count models and reports whether the aggregated caller signal remains
supported after accounting for CpG/sample-level variability.
