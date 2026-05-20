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
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    _ensure_package_path()
    from bsx2.viz.dmr_curves import build_default_dmr_curve_bundle, save_dmr_curve_bundle_outputs

    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bundle, results = build_default_dmr_curve_bundle(
        dmr_table=args.dmr_table,
        region_counts=args.region_counts,
        design=args.design,
        beta_binom=args.beta_binom,
        caller_support=args.caller_support,
        annotation=args.annotation,
        out_dir=None,
        top_n=args.top_n,
        min_total=args.min_total,
        render=not args.no_render,
    )
    formats = tuple(item.strip() for item in args.formats.split(",") if item.strip())
    save_dmr_curve_bundle_outputs(bundle, results, out_dir, formats=formats)
    n_warnings = sum(len(result.warnings) for result in results)
    print(f"wrote DMR curve bundle: {out_dir}")
    print(f"curves: {len(results)}")
    print(f"warnings: {n_warnings}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
