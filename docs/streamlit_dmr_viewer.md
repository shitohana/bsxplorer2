# BSX2 DMR Evidence Viewer

`apps/bsx2_dmr_viewer.py` is a lightweight Streamlit demo frontend for inspecting existing BSX2 DMR and evidence outputs. It is intentionally separate from the statistical pipeline.

Run:

```bash
streamlit run apps/bsx2_dmr_viewer.py
```

Expected input files are TSV tables such as:

- `dmr_evidence_scores.tsv`
- `dmr_region_count_tests.tsv`
- `dmr_caller_support_matrix.tsv`
- `dmr_beta_binomial_tests.tsv`

The app supports file upload and optional local path input. It does not save uploaded files into the repository and it does not create a project cache directory.

## Features

- dark-theme viewer for DMR/evidence tables;
- filters for context, evidence class, q-value, absolute delta, caller support, and top N rows;
- summary cards for total regions, q-value support, evidence classes, beta-binomial confirmation, and external caller support when those columns are available;
- overview, table, distribution, caller-support, and method-note tabs;
- download of the filtered TSV table.

## Limitations

- viewer only; no raw processing, DMR calling, or evidence recalculation;
- no Bismark, DSS, methylKit, dmrseq, or metilene execution;
- q-values can originate from different statistical layers and should not be mixed without method context;
- external harmonization standardizes columns, not caller-specific assumptions;
- beta-binomial aggregated output is a GLM validation layer, not a random-effect GLMM;
- large TSV uploads are previewed by reading only the first rows.
