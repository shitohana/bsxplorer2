#!/usr/bin/env python3
"""Build reusable BSX2 DMR visualization curve bundle."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _ensure_package_path() -> None:
    root = Path(__file__).resolve().parents[1]
    src = root / "python" / "src"
    if src.exists() and str(src) not in sys.path:
        sys.path.insert(0, str(src))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build reusable DMR curve specs/tables from BSX2 DMR outputs.",
    )
    parser.add_argument("--dmr-table", required=True, help="DMR/evidence TSV, e.g. dmr_evidence_scores.tsv")
    parser.add_argument("--region-counts", help="Optional region_counts TSV with region_id/sample_id/Y/m")
    parser.add_argument("--design", help="Optional design TSV with sample_id/condition")
    parser.add_argument("--beta-binom", help="Optional beta-binomial validation TSV")
    parser.add_argument("--caller-support", help="Optional caller support matrix TSV")
    parser.add_argument("--annotation", help="Optional annotation/enrichment TSV")
    parser.add_argument("--out-dir", required=True, help="Output directory for bundle specs/tables/figures")
    parser.add_argument("--top-n", type=int, default=1000, help="Maximum DMRs/variable regions for matrix curves")
    parser.add_argument("--min-total", type=int, default=5, help="Minimum region coverage for methylation matrix values")
    parser.add_argument(
        "--formats",
        default="json,tsv,png",
        help="Comma-separated output formats for generated artifacts; specs/tables are always saved when possible.",
    )
    parser.add_argument(
        "--no-render",
        action="store_true",
        help="Save specs/tables only; skip optional figure rendering.",
    )
    parser.add_argument(
        "--thesis-ready-only",
        action="store_true",
        help="Write tables/figures only for curves with quality_status=thesis_ready; manifest still records all curves.",
    )
    parser.add_argument(
        "--prefer-full-dmr-table",
        action="store_true",
        help="Search near --dmr-table for a fuller region-level DMR table and use it when better.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    _ensure_package_path()
    from bsx2.viz.dmr_curves import build_default_dmr_curve_bundle, save_dmr_curve_bundle_outputs
    from bsx2.viz.dmr_curve_data import find_best_dmr_table

    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dmr_table = args.dmr_table
    selection_lines = [
        "# DMR Table Selection Report",
        "",
        f"Requested DMR table: `{args.dmr_table}`",
    ]
    best, reason = find_best_dmr_table(args.dmr_table)
    if best is not None:
        selection_lines.append(f"Best nearby DMR table: `{best}`")
        selection_lines.append(f"Reason: {reason}")
        if args.prefer_full_dmr_table and Path(args.dmr_table).resolve() != best.resolve():
            dmr_table = str(best)
            selection_lines.append("Action: using best nearby full DMR table because --prefer-full-dmr-table was set.")
        else:
            selection_lines.append("Action: using requested DMR table.")
    else:
        selection_lines.append(f"Best nearby DMR table: not found ({reason})")
        selection_lines.append("Action: using requested DMR table.")
    (out_dir / "dmr_table_selection_report.md").write_text("\n".join(selection_lines) + "\n", encoding="utf-8")
    bundle, results = build_default_dmr_curve_bundle(
        dmr_table=dmr_table,
        region_counts=args.region_counts,
        design=args.design,
        beta_binom=args.beta_binom,
        caller_support=args.caller_support,
        annotation=args.annotation,
        out_dir=None,
        top_n=args.top_n,
        min_total=args.min_total,
        render=not args.no_render,
        thesis_ready_only=args.thesis_ready_only,
    )
    formats = tuple(item.strip() for item in args.formats.split(",") if item.strip())
    save_dmr_curve_bundle_outputs(bundle, results, out_dir, formats=formats, thesis_ready_only=args.thesis_ready_only)
    _write_input_semantic_qc(args, dmr_table, out_dir)
    n_warnings = sum(len(result.warnings) for result in results)
    print(f"wrote DMR curve bundle: {out_dir}")
    print(f"curves: {len(results)}")
    print(f"warnings: {n_warnings}")
    return 0


def _write_input_semantic_qc(args: argparse.Namespace, dmr_table: str, out_dir: Path) -> None:
    import csv
    import pandas as pd

    from bsx2.viz.dmr_curve_data import (
        classify_annotation_table,
        classify_design_table,
        classify_dmr_table,
        classify_region_counts_table,
        read_annotation_table,
        read_beta_binom_table,
        read_caller_support_table,
        read_design_table,
        read_dmr_table,
        read_region_counts,
    )

    specs = [
        ("dmr_table", dmr_table, read_dmr_table, classify_dmr_table, "chromosome/evidence/delta/volcano"),
        ("region_counts", args.region_counts, read_region_counts, classify_region_counts_table, "pca/heatmap/methylation distributions"),
        ("design", args.design, read_design_table, classify_design_table, "sample labels and condition groups"),
        ("beta_binom", args.beta_binom, read_beta_binom_table, None, "beta-binomial summary"),
        ("caller_support", args.caller_support, read_caller_support_table, None, "caller support summary"),
        ("annotation", args.annotation, read_annotation_table, classify_annotation_table, "annotation composition/enrichment"),
    ]
    rows = []
    region_counts_df = None
    for role, path, reader, classifier, usable in specs:
        exists = bool(path and Path(path).exists())
        columns = ""
        n_rows = "not_checked"
        semantic_type = "not_provided"
        warnings: list[str] = []
        if exists:
            df, read_warnings = reader(path)
            warnings.extend(read_warnings)
            columns = ",".join(str(col) for col in df.columns)
            n_rows = len(df)
            if role == "region_counts":
                region_counts_df = df
            if role == "design":
                semantic_type, class_warnings = classify_design_table(df, region_counts_df)
            elif classifier:
                semantic_type, class_warnings = classifier(df)
            else:
                semantic_type = role + "_table"
                class_warnings = []
            warnings.extend(class_warnings)
        rows.append(
            {
                "input_role": role,
                "path": path or "",
                "exists": exists,
                "n_rows_estimate": n_rows,
                "columns": columns,
                "semantic_type": semantic_type,
                "usable_for": usable if exists else "",
                "warnings": ";".join(dict.fromkeys(warnings)),
            }
        )
    with (out_dir / "dmr_input_semantic_qc.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
