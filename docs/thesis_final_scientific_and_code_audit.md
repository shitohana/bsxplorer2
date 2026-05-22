# Thesis Final Scientific And Code Audit

Scope: read-only audit of two recently added layers:

1. Rust-backed region count aggregation.
2. DMR-linked gene prioritization / functional evidence prioritization.

This report is an audit artifact only. It does not change implementation,
statistical methods, frozen outputs, or runtime results.

## Executive Summary

Overall status: acceptable for thesis use with caveats.

The Rust-backed aggregation layer is framed as an engineering acceleration path
for region-level count aggregation, not as a DMR caller. The documentation states
that it does not change the DMR caller, Regional Evidence Model, p-values,
confidence intervals, evidence scoring, QC interpretation, or visualization
logic (`docs/rust_region_aggregation_backend.md:5-11`). It also states that
statistical testing remains outside Rust (`docs/rust_region_aggregation_backend.md:53-64`).

The functional prioritization layer is also framed conservatively. It is called
candidate gene prioritization, says it does not prove gene function, does not
infer causality, does not run RNA-seq/chromatin processing, and does not change
the DMR statistical model (`docs/dmr_functional_prioritization.md:3-9`,
`docs/dmr_functional_prioritization.md:96-103`).

Primary risks are not overclaiming in code, but interpretation and edge-case
behavior:

- Rust `chunk_size` is exposed but is not clearly active in the `.bsx`
  `aggregate_counts_for_regions` path.
- Rust and pandas schemas match, but strand policy parity is imperfect:
  pandas accepts `opposite`; Rust currently supports `both`, `plus`, `minus`,
  and `region_strand`.
- Functional prioritization assumes the RNA-seq contrast direction matches the
  DMR contrast direction; this is not validated.
- Missing TE/chromatin flags conflate "annotation not supplied" with "supplied
  but no overlap".
- Class names such as `high_confidence_candidate` are acceptable because they
  retain "candidate", but thesis wording should clarify that confidence refers
  to prioritization support, not validated gene function.

## Audit Checklist

| Check | Status | Evidence | Notes |
| --- | --- | --- | --- |
| Rust aggregation not named DMR-caller | PASS | Rust docs say the backend does not change the DMR caller and is not a new DMR calling model (`docs/rust_region_aggregation_backend.md:9-11`, `docs/rust_region_aggregation_backend.md:138-144`). | Safe wording. |
| Functional prioritization not named proof of gene function | PASS | Docs and CLI summary explicitly say it does not prove function or causality (`docs/dmr_functional_prioritization.md:6-9`, `scripts/prioritize_dmr_linked_genes.py:48`). | Good for defense wording. |
| RNA-seq/TE/chromatin accepted as ready tables | PASS | CLI args and docs describe optional TE annotation, ready-made RNA-seq table, and BED/narrowPeak-like chromatin peaks (`docs/dmr_functional_prioritization.md:19-23`, `scripts/prioritize_dmr_linked_genes.py:80-82`). | No alignment/DE pipeline is invoked. |
| Missing evidence flags exist | PASS | Output schema includes `missing_evidence_flags` (`python/src/bsx2/analysis/dmr_functional_prioritization.py:18-41`); flags are populated for no gene, no expression, no TE overlap/annotation, no chromatin overlap/data (`python/src/bsx2/analysis/dmr_functional_prioritization.py:566-578`). | Needs wording caveat because no-overlap and no-input are conflated. |
| DMR statistical model unchanged | PASS | Rust docs state p-value/CI/evidence scoring remain unchanged (`docs/rust_region_aggregation_backend.md:9-11`); functional docs state DMR model unchanged (`docs/dmr_functional_prioritization.md:96-103`); CLI manifest notes no DMR model step was run (`scripts/prioritize_dmr_linked_genes.py:149-152`). | Good. |
| Rust/pandas output schema identical | PASS | Pandas schema is `REGION_SIGNAL_OUTPUT_COLUMNS` (`python/src/bsx2/analysis/region_signal.py:47-62`); Rust wrapper schema is `REGION_SIGNAL_SCHEMA` with the same columns (`python/src/bsx2/analysis/region_signal_rust.py:22-37`); Rust binding returns the same keys (`python/rust/region_aggregation.rs:183-200`). | Schema parity is strong. |
| `backend_used` in summary | PASS for CLI summary | CLI summary writes `backend_used` (`scripts/aggregate_region_signal.py:80-93`). | API return does not include backend metadata. |
| Clear limitations in docs | PASS | Rust limitations are explicit (`docs/rust_region_aggregation_backend.md:127-136`); functional limitations are explicit (`docs/dmr_functional_prioritization.md:96-103`). | Adequate. |

## Rust-Backed Region Count Aggregation

### What Is Scientifically Safe

- The docs correctly frame Rust as acceleration for count aggregation, not a new
  statistical method (`docs/rust_region_aggregation_backend.md:35-51`).
- Downstream Regional Evidence statistics remain in Python and are explicitly
  listed outside Rust (`docs/rust_region_aggregation_backend.md:53-64`).
- The public Python API retains pandas fallback and backend discovery
  (`python/src/bsx2/analysis/region_signal.py:264-297`).
- The result schema is stable across Rust and pandas:
  `region_id, seqname, chrom, start, end, strand, context, sample_id, mC, uC,
  total, n_cytosines, mean_methylation, coverage_qc`
  (`python/src/bsx2/analysis/region_signal.py:47-62`,
  `python/src/bsx2/analysis/region_signal_rust.py:22-37`).

### Weak Spots / Risks

1. `chunk_size` is exposed but not clearly used in the `.bsx` reader path.

   `AggregationOptions` includes `chunk_size` (`bsxplorer2/src/tools/region_aggregation.rs:41-48`), and PyO3 accepts it
   (`python/rust/region_aggregation.rs:21-34`). However, the production `.bsx`
   path `aggregate_counts_for_regions` iterates over sorted regions and performs
   one `reader.query` per region (`bsxplorer2/src/tools/region_aggregation.rs:158-219`).
   The separate in-memory `aggregate_counts_for_region_chunk` path exists
   (`bsxplorer2/src/tools/region_aggregation.rs:222-253`), but the audit did not
   find evidence that `chunk_size` drives `.bsx` query chunking. Thesis wording
   should say "Rust-backed indexed aggregation", not "fully streaming chunked
   genome-wide aggregation" unless later improved.

2. Strand policy parity is incomplete.

   Python accepts `opposite` in `_normalize_strand_policy` and CLI choices
   (`python/src/bsx2/analysis/region_signal.py:142-156`,
   `scripts/aggregate_region_signal.py:32-38`). Rust `StrandPolicy::parse`
   supports `both`, `plus`, `minus`, and `region_strand`, but not `opposite`
   (`bsxplorer2/src/tools/region_aggregation.rs:21-39`). If a user requests
   `backend="rust"` with `strand_policy="opposite"`, the Rust binding should
   reject it. That is a clear error, not silent corruption, but it is a parity
   gap.

3. Coordinate convention should be stated in docs.

   Pandas and Rust filtering both use inclusive region membership for positions
   (`python/src/bsx2/analysis/region_signal.py:240-253`,
   `bsxplorer2/src/tools/region_aggregation.rs:293-319`). Rust builds the
   `Contig` query with `end.saturating_add(1)` before final filtering
   (`bsxplorer2/src/tools/region_aggregation.rs:200-207`). This appears
   internally consistent, but the public docs should explicitly say whether input
   regions are interpreted as 1-based inclusive, 0-based closed, or generic
   point-inclusive intervals.

4. `backend_used` is visible in CLI summaries, not in API DataFrame metadata.

   This satisfies CLI reporting (`scripts/aggregate_region_signal.py:80-93`),
   but downstream programmatic workflows cannot inspect backend choice from the
   returned DataFrame alone.

### Thesis-Safe Wording

Use:

> Rust-backed RegionSignal accelerates indexed aggregation of methylated and
> unmethylated counts over predefined regions while preserving the existing DMR
> calling and Regional Evidence statistical model.

Avoid:

> Rust DMR caller, Rust statistical DMR model, Rust proof of DMR significance,
> fully streaming genome-wide aggregator.

## DMR-Linked Gene Prioritization

### What Is Scientifically Safe

- Module docstring says it consumes existing evidence and does not call DMRs,
  prove function, infer causality, or run RNA-seq/chromatin pipelines
  (`python/src/bsx2/analysis/dmr_functional_prioritization.py:1-7`).
- Documentation repeats the same limitations and states that stronger
  interpretation requires external RNA-seq, TE, chromatin evidence, and
  biological validation (`docs/dmr_functional_prioritization.md:96-103`).
- The CLI only accepts tables; it does not invoke RNA-seq alignment, DESeq2,
  edgeR, ATAC-seq, ChIP-seq, or DMR recalculation
  (`scripts/prioritize_dmr_linked_genes.py:76-86`, `scripts/prioritize_dmr_linked_genes.py:149-152`).
- Missing optional inputs are converted into warnings, not crashes
  (`scripts/prioritize_dmr_linked_genes.py:99-110`).
- The output schema includes stable evidence and interpretation columns
  (`python/src/bsx2/analysis/dmr_functional_prioritization.py:18-41`).

### Weak Spots / Risks

1. RNA-seq contrast direction is assumed, not validated.

   Direction consistency uses DMR `region_delta` and RNA-seq `log2FC`
   (`python/src/bsx2/analysis/dmr_functional_prioritization.py:486-502`).
   The module does not verify that the expression table contrast has the same
   condition ordering as the DMR contrast. This is the most important
   methodological caveat: promoter hyper + down is only meaningful if both
   deltas are aligned to the same condition comparison.

2. Functional score is heuristic.

   The score is a rule-based weighted sum (`python/src/bsx2/analysis/dmr_functional_prioritization.py:505-520`), and classes are heuristic
   (`python/src/bsx2/analysis/dmr_functional_prioritization.py:522-533`).
   Documentation correctly says it is not a formal functional validation
   statistic (`docs/dmr_functional_prioritization.md:68-84`). In thesis text,
   treat it as prioritization for review, not a validated predictor.

3. Missing evidence flags conflate missing annotation with no overlap.

   Row-level flags use `no_te_overlap_or_no_te_annotation` and
   `no_chromatin_overlap_or_no_chromatin_data`
   (`python/src/bsx2/analysis/dmr_functional_prioritization.py:573-576`).
   This is conservative but less precise. The CLI warnings capture whether
   optional inputs were missing globally (`scripts/prioritize_dmr_linked_genes.py:99-110`).

4. Gene/exon/intron linking may be limited by annotation ID conventions.

   GFF parsing chooses `gene_id`, `ID`, `Parent`, or `Name`
   (`python/src/bsx2/analysis/dmr_functional_prioritization.py:117-124`).
   Exon/intron overlap is checked only for features with the same normalized
   `gene_id` as the chosen linked gene (`python/src/bsx2/analysis/dmr_functional_prioritization.py:372-379`).
   Some GFFs use transcript IDs or parent chains that may require more explicit
   parent-child resolution.

5. `high_confidence_candidate` can be misunderstood.

   The class retains "candidate", so it is acceptable, but thesis wording should
   define it as "high support for prioritization", not confirmed gene function.

### Thesis-Safe Wording

Use:

> The functional prioritization layer links existing DMR evidence to gene,
> repeat, expression, and chromatin annotations and assigns a conservative
> candidate support score for hypothesis generation.

Avoid:

> Functional gene discovery, proof of regulatory function, causal methylation
> effect, validated expression regulation.

## Final Risk Ranking

| Risk | Severity | Why It Matters | Recommended Handling |
| --- | --- | --- | --- |
| RNA-seq contrast direction not validated | Medium | Direction consistency can be wrong if expression log2FC is oriented differently from DMR delta. | State as limitation; require contrast metadata before using expression-based conclusions. |
| Rust `chunk_size` not active in `.bsx` path | Medium | Could overpromise chunked behavior/performance. | Describe as Rust-backed indexed aggregation; reserve "full chunked engine" for future work. |
| Rust strand policy gap for `opposite` | Low/Medium | Backend parity issue for one policy. | Document as unsupported in Rust or route `opposite` to pandas until implemented. |
| TE/chromatin missing flags are coarse | Low | Conservative but less diagnostic. | Use CLI warnings plus row flags together. |
| GFF parent-child resolution is shallow | Low/Medium | Exon/intron labels may be incomplete for complex annotations. | Use gene/promoter/gene_body as primary claims; treat exon/intron as best effort. |
| Programmatic `backend_used` missing from returned DataFrame | Low | CLI summary has it, API does not. | Use CLI summary for thesis runs; add API metadata later if needed. |

## Final Verdict

Both layers are defensible for thesis support if described narrowly.

Rust-backed aggregation supports an engineering contribution: high-throughput
count aggregation over predefined regions with a pandas fallback and unchanged
statistical model.

DMR-linked functional prioritization supports a biological interpretation aid:
candidate ranking from existing DMR, gene, TE, expression, and chromatin tables.
It should not be described as functional proof or causality inference.

No audited text was found that directly claims Rust aggregation is a DMR caller
or that functional prioritization proves gene function. The remaining risks are
mostly around edge-case behavior and wording precision.
