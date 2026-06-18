#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


from dmr_validation_framework.core.harmonization import ADAPTERS, TIER_DESCRIPTIONS, add_dmr_tiers, adapter_for_caller, build_caller_support_matrix, write_schema_validation_report


TIER_COLUMNS = [
    "region_id",
    "chrom",
    "start",
    "end",
    "context",
    "n_callers_supporting",
    "supporting_callers",
    "n_model_families_supporting",
    "supporting_model_families",
    "direction_consensus",
    "fraction_same_direction",
    "n_opposite_direction",
    "mean_reciprocal_overlap",
    "min_reciprocal_overlap",
    "best_q_value",
    "n_callers_q05",
    "n_callers_q10",
    "max_abs_delta",
    "median_abs_delta",
    "mean_abs_delta",
    "caller_conflict_flag",
    "multi_caller_evidence_score",
    "multi_caller_evidence_class",
    "validation_robustness_score",
    "final_confidence_class",
    "model_agreement_status",
    "final_dmr_tier",
    "tier_reasons",
]


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else "," if path.suffix.lower() == ".csv" else None
    return pd.read_csv(path, sep=sep, engine="python")


def write_tier_summary_outputs(tiers: pd.DataFrame, out_dir: Path) -> None:
    summary = (
        tiers.groupby("final_dmr_tier", dropna=False)
        .size()
        .reset_index(name="n_regions")
        .sort_values(["final_dmr_tier"])
    )
    summary["description"] = summary["final_dmr_tier"].map(TIER_DESCRIPTIONS).fillna("")
    summary.to_csv(out_dir / "dmr_consensus_tiers_summary.tsv", sep="\t", index=False)
    table = tiers[[column for column in TIER_COLUMNS if column in tiers.columns]].head(200)
    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\"><title>DMR consensus tiers</title>",
            "<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse;margin:16px 0}td,th{border:1px solid #ddd;padding:5px 7px;font-size:12px}th{background:#f4f4f4}.note{max-width:980px}</style>",
            "</head><body>",
            "<h1>DMR consensus tiers</h1>",
            "<p class=\"note\"><strong>Tier is not a ranking of DMR callers.</strong> Tier is a confidence ranking of a DMR candidate based on independent evidence consistency. Caller-specific q-values are used only as threshold evidence (q&lt;=0.05 or q&lt;=0.10), not as direct cross-caller rankings.</p>",
            "<h2>Tier summary</h2>",
            summary.to_html(index=False, escape=True),
            "<h2>Top rows</h2>",
            table.to_html(index=False, escape=True),
            "</body></html>",
        ]
    )
    (out_dir / "dmr_consensus_tiers_summary.html").write_text(body + "\n", encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Import external DMR candidates into the canonical BSX2 schema.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--caller", required=True, choices=sorted(ADAPTERS))
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--contrast-id", default="")
    parser.add_argument("--condition-a", default="")
    parser.add_argument("--condition-b", default="")
    parser.add_argument("--context")
    parser.add_argument("--counts-manifest")
    parser.add_argument("--dmr-region-count-tests")
    parser.add_argument("--dmr-evidence-scores")
    parser.add_argument("--annotation")
    parser.add_argument("--validation-audit", help="Optional enhanced validation audit TSV, e.g. delta_weighting_sensitivity.tsv.")
    parser.add_argument("--overlap-threshold", type=float, default=0.5)
    parser.add_argument("--enable-evidence-join", default="false")
    parser.add_argument("--enable-beta-binomial-join", default="false")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adapter_cls = adapter_for_caller(args.caller)
    adapter = adapter_cls(args.input, contrast_id=args.contrast_id, condition_a=args.condition_a, condition_b=args.condition_b, context=args.context)
    canonical = adapter.read()
    canonical.to_csv(out / "external_dmr_candidates_canonical.tsv", sep="\t", index=False)
    write_schema_validation_report(canonical, out / "external_dmr_schema_validation.tsv")
    pd.DataFrame({"warning": adapter.warnings or [""]}).to_csv(out / "external_dmr_import_warnings.tsv", sep="\t", index=False)
    canonical.to_csv(out / "dmr_evidence_harmonized.tsv", sep="\t", index=False)
    summary = pd.DataFrame([{"source_caller": args.caller, "n_candidates": len(canonical)}])
    summary.to_csv(out / "dmr_caller_comparison_summary.tsv", sep="\t", index=False)
    support = build_caller_support_matrix([canonical], overlap_threshold=args.overlap_threshold)
    if args.validation_audit:
        support = add_dmr_tiers(support, read_table(args.validation_audit))
    support.to_csv(out / "dmr_caller_support_matrix.tsv", sep="\t", index=False)
    support[[column for column in TIER_COLUMNS if column in support.columns]].to_csv(out / "dmr_consensus_tiers.tsv", sep="\t", index=False)
    write_tier_summary_outputs(support, out)
    (out / "dmr_external_method_audit.md").write_text(
        "# External DMR Method Audit\n\n"
        "Schema harmonization only; external callers were not run.\n\n"
        "Tier is not a ranking of DMR callers. Tier is a confidence ranking of a DMR candidate based on independent evidence consistency.\n",
        encoding="utf-8",
    )
    (out / "external_dmr_harmonization_manifest.json").write_text(json.dumps({"created_at": datetime.now(timezone.utc).isoformat(), "input": args.input, "caller": args.caller, "outputs": sorted(p.name for p in out.iterdir())}, indent=2), encoding="utf-8")
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
