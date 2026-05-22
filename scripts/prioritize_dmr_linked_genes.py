#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis import (  # noqa: E402
    prioritize_dmr_linked_genes,
    read_chromatin_peaks,
    read_dmr_evidence_table,
    read_expression_table,
    read_gene_annotation,
    read_te_annotation,
    write_dmr_functional_prioritization_table,
)


def _sha256_if_small(path: Path, max_bytes: int = 500_000_000) -> str:
    if not path.exists() or path.stat().st_size > max_bytes:
        return "not_computed"
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_warnings(path: Path, warnings: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("warning_type\tmessage\tseverity\n")
        for warning in warnings:
            handle.write(f"{warning['warning_type']}\t{warning['message']}\t{warning['severity']}\n")


def _summary_text(result, warnings: list[dict[str, str]], args: argparse.Namespace) -> str:
    class_counts = result["functional_support_class"].value_counts(dropna=False).to_dict() if not result.empty else {}
    link_counts = result["link_type"].value_counts(dropna=False).to_dict() if not result.empty else {}
    lines = [
        "# DMR Functional Prioritization Summary",
        "",
        "This run prioritizes candidate genes linked to DMR evidence. It does not prove gene function or infer causality.",
        "",
        f"- DMR evidence table: {args.dmr_evidence}",
        f"- gene annotation: {args.genes}",
        f"- TE annotation: {args.te_annotation or 'not provided'}",
        f"- RNA-seq differential expression: {args.expression or 'not provided'}",
        f"- chromatin peaks: {args.chromatin_peaks or 'not provided'}",
        f"- promoter_upstream: {args.promoter_upstream}",
        f"- promoter_downstream: {args.promoter_downstream}",
        f"- max_distance: {args.max_distance}",
        f"- output rows: {len(result)}",
        "",
        "## Functional Support Classes",
    ]
    lines.extend(f"- {key}: {value}" for key, value in sorted(class_counts.items()))
    lines.append("")
    lines.append("## Link Types")
    lines.extend(f"- {key}: {value}" for key, value in sorted(link_counts.items()))
    lines.append("")
    lines.append("## Warnings")
    if warnings:
        lines.extend(f"- {w['warning_type']}: {w['message']}" for w in warnings)
    else:
        lines.append("- none")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Prioritize candidate genes linked to existing DMR evidence.")
    parser.add_argument("--dmr-evidence", required=True, help="Existing DMR evidence table.")
    parser.add_argument("--genes", required=True, help="Gene annotation GFF/GTF/BED/TSV.")
    parser.add_argument("--te-annotation", help="Optional TE annotation GFF/GTF/BED/TSV.")
    parser.add_argument("--expression", help="Optional ready-made RNA-seq differential expression table.")
    parser.add_argument("--chromatin-peaks", help="Optional BED/narrowPeak-like chromatin peak table.")
    parser.add_argument("--out-dir", required=True, help="Output directory.")
    parser.add_argument("--promoter-upstream", type=int, default=2000)
    parser.add_argument("--promoter-downstream", type=int, default=200)
    parser.add_argument("--max-distance", type=int, default=10000)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[dict[str, str]] = []

    dmr_df = read_dmr_evidence_table(args.dmr_evidence)
    gene_df = read_gene_annotation(args.genes)
    te_df = None
    expression_df = None
    chromatin_df = None

    if args.te_annotation:
        te_df = read_te_annotation(args.te_annotation)
    else:
        warnings.append({"warning_type": "te_annotation_missing", "message": "TE evidence was not provided; TE-linked support is unavailable.", "severity": "info"})
    if args.expression:
        expression_df = read_expression_table(args.expression)
    else:
        warnings.append({"warning_type": "expression_missing", "message": "RNA-seq differential expression table was not provided; expression support is marked as no_expression_data.", "severity": "info"})
    if args.chromatin_peaks:
        chromatin_df = read_chromatin_peaks(args.chromatin_peaks)
    else:
        warnings.append({"warning_type": "chromatin_missing", "message": "Chromatin peak evidence was not provided; chromatin support is unavailable.", "severity": "info"})

    result = prioritize_dmr_linked_genes(
        dmr_df,
        gene_df,
        te_df,
        expression_df,
        chromatin_df,
        promoter_upstream=args.promoter_upstream,
        promoter_downstream=args.promoter_downstream,
        max_distance=args.max_distance,
    )

    result_path = out_dir / "dmr_functional_prioritization.tsv"
    summary_path = out_dir / "dmr_functional_prioritization_summary.md"
    manifest_path = out_dir / "dmr_functional_prioritization_manifest.json"
    warnings_path = out_dir / "warnings.tsv"

    write_dmr_functional_prioritization_table(result, result_path)
    _write_warnings(warnings_path, warnings)
    summary_path.write_text(_summary_text(result, warnings, args), encoding="utf-8")

    outputs = [result_path, summary_path, warnings_path]
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "dmr_evidence": args.dmr_evidence,
            "genes": args.genes,
            "te_annotation": args.te_annotation,
            "expression": args.expression,
            "chromatin_peaks": args.chromatin_peaks,
        },
        "parameters": {
            "promoter_upstream": args.promoter_upstream,
            "promoter_downstream": args.promoter_downstream,
            "max_distance": args.max_distance,
        },
        "outputs": {str(path): {"size_bytes": path.stat().st_size, "sha256": _sha256_if_small(path)} for path in outputs},
        "warnings": warnings,
        "notes": [
            "Candidate gene prioritization only; this module does not prove gene function.",
            "No DMR statistical model, frozen artifact, RNA-seq alignment, DESeq2, edgeR, ATAC-seq, or ChIP-seq step was run.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
