# DMR Functional Prioritization

This module adds an optional candidate gene prioritization layer on top of
existing DMR evidence outputs.

It does not prove gene function, infer causality, run RNA-seq alignment, run
DESeq2/edgeR, call ATAC-seq or ChIP-seq peaks, or change the DMR statistical
model. It only combines already prepared tables and reports missing evidence
explicitly.

## Inputs

Required:

- Existing DMR evidence table with chromosome, start, end, delta, q-value, and
  evidence class columns.
- Gene annotation in GFF3, GTF, BED, or TSV-like format.

Optional:

- TE annotation in GFF3, GTF, BED, or TSV-like format.
- Ready-made RNA-seq differential expression table.
- Ready-made chromatin peak table in BED/narrowPeak-like format.

## Python API

```python
from bsx2.analysis import (
    prioritize_dmr_linked_genes,
    read_dmr_evidence_table,
    read_gene_annotation,
    read_expression_table,
)

dmrs = read_dmr_evidence_table("dmr_evidence_scores.tsv")
genes = read_gene_annotation("annotation.gff3.gz")
expression = read_expression_table("differential_expression.tsv")

prioritized = prioritize_dmr_linked_genes(
    dmrs,
    genes,
    expression_df=expression,
    promoter_upstream=2000,
    promoter_downstream=200,
    max_distance=10000,
)
```

## CLI

```bash
python scripts/prioritize_dmr_linked_genes.py \
  --dmr-evidence dmr_evidence_scores.tsv \
  --genes annotation.gff3.gz \
  --expression differential_expression.tsv \
  --te-annotation te_annotation.bed \
  --chromatin-peaks peaks.narrowPeak \
  --out-dir functional_prioritization
```

Outputs:

- `dmr_functional_prioritization.tsv`
- `dmr_functional_prioritization_summary.md`
- `dmr_functional_prioritization_manifest.json`
- `warnings.tsv`

## Scoring

The score is rule-based:

```text
functional_support_score =
  dmr_evidence_component
  + annotation_component
  + expression_component
  + direction_consistency_component
  + te_component
  + chromatin_component
  - qc_penalty
```

The score prioritizes DMR-linked candidates for review. It is not a formal
functional validation statistic.

Direction consistency rules are conservative:

- promoter hypermethylation plus gene downregulation is marked consistent;
- promoter hypomethylation plus gene upregulation is marked consistent;
- gene-body methylation is marked context-dependent;
- TE-associated DMRs are interpreted separately as repeat/TE methylation
  candidates;
- missing RNA-seq evidence is reported as `no_expression_data`, not treated as
  strong negative evidence.

## Limitations

- Candidate prioritization is not proof of gene function.
- The module does not infer causality between methylation and expression.
- Stronger interpretation requires external RNA-seq, TE annotation, chromatin
  evidence, and biological validation.
- Missing evidence is preserved in `missing_evidence_flags` and `warnings.tsv`.
- The DMR statistical model and Regional Evidence Model are unchanged.
