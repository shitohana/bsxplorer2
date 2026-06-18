#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from dmr_validation_framework.core.io import read_table as read_table_auto

from dmr_validation_framework.core.harmonization import TIER_DESCRIPTIONS, add_dmr_tiers, build_caller_support_matrix


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
    return read_table_auto(path)


def read_canonical_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    df = read_table(path)
    if "source_caller" not in df.columns:
        df["source_caller"] = path.stem
    return df


def write_outputs(tiers: pd.DataFrame, out_dir: Path, inputs: list[str], validation_audit: str | None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    tiers.to_csv(out_dir / "dmr_caller_support_matrix.tsv", sep="\t", index=False)
    tiers[[column for column in TIER_COLUMNS if column in tiers.columns]].to_csv(out_dir / "dmr_consensus_tiers.tsv", sep="\t", index=False)

    summary = (
        tiers.groupby("final_dmr_tier", dropna=False)
        .size()
        .reset_index(name="n_regions")
        .sort_values(["final_dmr_tier"])
    )
    summary["description"] = summary["final_dmr_tier"].map(TIER_DESCRIPTIONS).fillna("")
    summary.to_csv(out_dir / "dmr_consensus_tiers_summary.tsv", sep="\t", index=False)

    display = tiers[[column for column in TIER_COLUMNS if column in tiers.columns]].head(200)
    html = "\n".join(
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
            display.to_html(index=False, escape=True),
            "</body></html>",
        ]
    )
    (out_dir / "dmr_consensus_tiers_summary.html").write_text(html + "\n", encoding="utf-8")

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "canonical_inputs": inputs,
        "validation_audit": validation_audit,
        "outputs": [
            "dmr_caller_support_matrix.tsv",
            "dmr_consensus_tiers.tsv",
            "dmr_consensus_tiers_summary.tsv",
            "dmr_consensus_tiers_summary.html",
        ],
        "note": "Tier is not a ranking of DMR callers. Tier is a confidence ranking of a DMR candidate based on independent evidence consistency.",
    }
    (out_dir / "dmr_consensus_tiers_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build tier-based DMR consensus from harmonized canonical caller tables.")
    parser.add_argument("--canonical-dmr", action="append", required=True, help="Canonical DMR TSV/CSV. Repeat once per caller.")
    parser.add_argument("--validation-audit", help="Optional enhanced validation audit TSV, e.g. delta_weighting_sensitivity.tsv.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--overlap-threshold", type=float, default=0.5)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    canonical_tables = [read_canonical_table(path) for path in args.canonical_dmr]
    tiers = build_caller_support_matrix(canonical_tables, overlap_threshold=args.overlap_threshold)
    if args.validation_audit:
        tiers = add_dmr_tiers(tiers, read_table(args.validation_audit))
    write_outputs(tiers, Path(args.out_dir), args.canonical_dmr, args.validation_audit)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
