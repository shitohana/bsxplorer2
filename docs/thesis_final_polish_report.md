# Thesis Final Polish Report

## Scope

This polish applied only safety, wording, and edge-case fixes around:

- Rust-backed region count aggregation.
- DMR-linked candidate gene prioritization.

No DMR statistical model was changed. No frozen artifacts were modified. No
raw/Bismark/DSS/methylKit/dmrseq/metilene workflow was run.

## Changes

### RNA-seq Contrast Safety

- Added optional CLI arguments:
  - `--dmr-contrast-label`
  - `--expression-contrast-label`
  - `--require-matched-contrast`
- Direction consistency is now evaluated only when DMR and expression contrast
  labels are both provided and match.
- If expression evidence is provided without contrast labels,
  `direction_consistency` becomes
  `not_evaluated_contrast_not_validated`.
- If labels are provided but differ, `direction_consistency` becomes
  `not_evaluated_contrast_mismatch`, unless `--require-matched-contrast` is
  used, in which case a clear `ValueError` is raised.
- Summary and manifest now record contrast labels and
  `expression_direction_consistency_evaluated`.

### Strand Backend Parity

- `backend="auto"` with `strand_policy="opposite"` now routes to pandas for
  pandas-compatible count tables.
- `backend="rust"` with `strand_policy="opposite"` now raises:

```text
strand_policy='opposite' is not supported by Rust backend; use backend='pandas' or backend='auto'.
```

- Backend discovery notes now document Rust/pandas strand-policy support.

### Missing Evidence Flags

Functional prioritization row-level flags now distinguish:

- `te_input_missing`
- `te_no_overlap`
- `chromatin_input_missing`
- `chromatin_no_overlap`

This avoids conflating missing optional inputs with provided-but-nonoverlapping
evidence.

### Wording

Rust aggregation docs now use conservative wording:

- "Rust-backed indexed region count aggregation"
- "production-oriented Rust aggregation path"
- "`chunk_size` is exposed for future chunked implementations; the current
  `.bsx` path performs indexed per-region queries"

## Validation

Commands run:

- `python3 -m py_compile` on changed Python files.
- `python3 -m compileall -q python/src/bsx2/analysis scripts`.
- `python3 scripts/prioritize_dmr_linked_genes.py --help`.
- `python3 scripts/aggregate_region_signal.py --help`.
- Direct synthetic smoke checks for:
  - expression present without contrast labels;
  - matched expression/DMR labels;
  - mismatched labels with required match;
  - separated TE missing/no-overlap flags;
  - separated chromatin missing/no-overlap flags;
  - `backend="auto"` plus `strand_policy="opposite"`;
  - `backend="rust"` plus `strand_policy="opposite"`.
- CLI synthetic smoke checks for prioritization and aggregation.
- `git diff --check` on changed files.

`pytest` was attempted but is not available in this environment:

```text
/usr/bin/python3: No module named pytest
```

## Runtime Smoke Output Location

Temporary validation outputs were written under:

```text
/mnt/g/bsx2_raw/bnapus_gse202609_6sample/final_polish_validation/
```

No frozen artifacts or processed B. napus proof-run outputs were modified.
