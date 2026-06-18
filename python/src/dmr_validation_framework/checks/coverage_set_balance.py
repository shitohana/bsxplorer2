#!/usr/bin/env python3
"""Audit coverage-set balance for DMR/metagene downstream validation.

This check quantifies whether regional methylation summaries use different
covered CpG/cytosine sets across samples and how a strict/common-CpG set changes
regional methylation deltas. It does not call DMRs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


from dmr_validation_framework.core.coverage import (  # noqa: E402
    aggregate_counts_per_sample,
    aggregate_counts_common_cpg,
    build_cpg_coverage_qc,
    normalize_cpg_counts_table,
)
from dmr_validation_framework.core.io import ensure_out_dir, find_preferred_file, read_table, write_tsv  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-cpg-counts", type=Path)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--condition-column", default="condition")
    parser.add_argument("--min-coverage", type=int, default=1)
    parser.add_argument("--min-covered-per-group", type=int)
    parser.add_argument("--min-covered-A", type=int)
    parser.add_argument("--min-covered-B", type=int)
    parser.add_argument("--min-common-cpgs", type=int, default=3)
    parser.add_argument("--min-region-total-coverage", type=int)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def find_file(name: str) -> Path | None:
    return find_preferred_file([name])


def skipped(out_dir: Path, reason: str) -> int:
    row = {"status": "SKIPPED", "notes": reason}
    for name in [
        "coverage_set_region_qc.tsv",
        "coverage_set_sample_qc.tsv",
        "coverage_set_delta_shift.tsv",
        "coverage_set_summary.tsv",
    ]:
        write_tsv(out_dir / name, [row])
    return 0


def direction(value: float | int | None) -> str:
    if value is None or pd.isna(value):
        return "NA"
    value = float(value)
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def region_delta(df: pd.DataFrame, value_col: str, condition_order: list[str]) -> pd.Series:
    means = df.groupby(["region_id", "condition"])[value_col].mean().unstack("condition")
    if len(condition_order) < 2 or condition_order[0] not in means.columns or condition_order[1] not in means.columns:
        return pd.Series(dtype=float)
    return means[condition_order[1]] - means[condition_order[0]]


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    counts_path = args.region_cpg_counts or find_file("region_cpg_counts.tsv")
    design_path = args.design or find_file("design.tsv")
    if not counts_path or not design_path:
        return skipped(out_dir, "per-CpG counts table or design.tsv not found")

    counts = read_table(counts_path)
    design = read_table(design_path)
    if "sample_id" not in design.columns or args.condition_column not in design.columns:
        return skipped(out_dir, f"design table lacks sample_id or {args.condition_column}")
    required = {"region_id", "cpg_id", "sample_id", "mC", "uC", "total"}
    missing = sorted(required - set(counts.columns))
    if missing:
        return skipped(out_dir, f"region CpG counts missing columns: {', '.join(missing)}")

    min_covered_per_group = args.min_covered_per_group
    if args.min_covered_A is not None or args.min_covered_B is not None:
        min_covered_per_group = {
            "A": args.min_covered_A if args.min_covered_A is not None else args.min_covered_per_group or 1,
            "B": args.min_covered_B if args.min_covered_B is not None else args.min_covered_per_group or 1,
        }

    common_cpg, region_qc = build_cpg_coverage_qc(
        counts,
        design,
        condition_col=args.condition_column,
        min_coverage=args.min_coverage,
        min_covered_per_group=min_covered_per_group,
        min_common_cpgs=args.min_common_cpgs,
        min_region_total_coverage=args.min_region_total_coverage,
    )
    normalized_counts = normalize_cpg_counts_table(counts).drop(columns=["condition"]).merge(
        design[["sample_id", args.condition_column]].rename(columns={args.condition_column: "condition"}),
        on="sample_id",
        how="left",
    )
    common_region = aggregate_counts_common_cpg(
        normalized_counts,
        common_cpg,
        design,
        min_common_cpgs=args.min_common_cpgs,
        min_region_total_coverage=args.min_region_total_coverage,
    )
    per_sample = aggregate_counts_per_sample(normalized_counts, min_coverage=args.min_coverage)
    condition_order = list(design[args.condition_column].astype(str).drop_duplicates())[:2]

    per_counts = per_sample[["region_id", "sample_id", "condition", "n_cpg_per_sample", "N_per_sample"]].rename(
        columns={"N_per_sample": "total_coverage_per_sample"}
    )
    sample_qc = per_counts.merge(
        common_region[
            [
                "region_id",
                "sample_id",
                "n_common_cpgs",
                "N_common",
                "coverage_qc_status",
                "min_region_total_coverage",
                "region_sample_coverage_status",
            ]
        ],
        on=["region_id", "sample_id"],
        how="left",
    )
    sample_qc["n_cpg_common"] = sample_qc["n_common_cpgs"].fillna(0)
    sample_qc["fraction_sample_cpg_in_common"] = np.where(
        sample_qc["n_cpg_per_sample"] > 0,
        sample_qc["n_cpg_common"] / sample_qc["n_cpg_per_sample"],
        np.nan,
    )
    sample_qc["total_coverage_common"] = sample_qc["N_common"].fillna(0)
    sample_qc["min_region_total_coverage"] = sample_qc["min_region_total_coverage"].where(
        sample_qc["min_region_total_coverage"].notna(),
        args.min_region_total_coverage,
    )
    sample_qc["region_sample_coverage_status"] = sample_qc["region_sample_coverage_status"].fillna("no_region_coverage")
    sample_qc["status"] = sample_qc["coverage_qc_status"].fillna("no_common_cpgs")
    sample_qc["notes"] = (
        "min_region_total_coverage not required; only N_common > 0 checked"
        if args.min_region_total_coverage is None
        else f"requires N_common >= {args.min_region_total_coverage}"
    )
    sample_qc = sample_qc[
        [
            "region_id",
            "sample_id",
            "condition",
            "n_cpg_per_sample",
            "n_cpg_common",
            "fraction_sample_cpg_in_common",
            "total_coverage_per_sample",
            "total_coverage_common",
            "min_region_total_coverage",
            "region_sample_coverage_status",
            "status",
            "notes",
        ]
    ]

    sample_status = (
        sample_qc.groupby("region_id")["region_sample_coverage_status"]
        .agg(
            n_samples_pass_region_coverage=lambda s: int((s == "pass_region_coverage").sum()),
            n_samples_fail_region_coverage=lambda s: int((s != "pass_region_coverage").sum()),
        )
        .reset_index()
    )
    region_qc = region_qc.merge(sample_status, on="region_id", how="left")
    region_qc["n_samples_pass_region_coverage"] = region_qc["n_samples_pass_region_coverage"].fillna(0).astype(int)
    region_qc["n_samples_fail_region_coverage"] = region_qc["n_samples_fail_region_coverage"].fillna(0).astype(int)

    def region_coverage_status(row: pd.Series) -> str:
        if row["coverage_qc_status"] == "no_common_cpgs":
            return "no_common_cpgs"
        if row["coverage_qc_status"] == "insufficient_common_cpgs":
            return "insufficient_common_cpgs"
        if args.min_region_total_coverage is None:
            return "pass_common_cpg_region_coverage_not_required"
        if int(row["n_samples_fail_region_coverage"]) == 0:
            return "pass_common_cpg_and_region_coverage"
        return "pass_common_cpg_but_low_region_coverage"

    region_qc["region_coverage_status"] = region_qc.apply(region_coverage_status, axis=1)

    delta_per = region_delta(per_sample, "mu_hat_per_sample", condition_order)
    delta_common = region_delta(common_region.rename(columns={"mu_hat_common": "mu_hat_common"}), "mu_hat_common", condition_order)
    delta = pd.DataFrame({"delta_per_sample": delta_per, "delta_common": delta_common}).reset_index()
    delta = delta.merge(region_qc[["region_id", "n_cpg_union", "n_cpg_common", "coverage_qc_status"]], on="region_id", how="left")
    delta["abs_delta_shift"] = (pd.to_numeric(delta["delta_common"], errors="coerce") - pd.to_numeric(delta["delta_per_sample"], errors="coerce")).abs()
    delta["direction_per_sample"] = delta["delta_per_sample"].map(direction)
    delta["direction_common"] = delta["delta_common"].map(direction)
    delta["direction_changed"] = (delta["direction_per_sample"] != delta["direction_common"]).astype(object)
    delta.loc[delta["direction_per_sample"].eq("NA") | delta["direction_common"].eq("NA"), "direction_changed"] = pd.NA
    delta = delta[
        [
            "region_id",
            "delta_per_sample",
            "delta_common",
            "abs_delta_shift",
            "direction_per_sample",
            "direction_common",
            "direction_changed",
            "n_cpg_union",
            "n_cpg_common",
            "coverage_qc_status",
        ]
    ]

    n_total = len(region_qc)
    n_pass = int((region_qc["coverage_qc_status"] == "pass_common_cpg").sum())
    n_insufficient = int((region_qc["coverage_qc_status"] == "insufficient_common_cpgs").sum())
    n_no_common = int((region_qc["coverage_qc_status"] == "no_common_cpgs").sum())
    n_sample_pass = int((sample_qc["region_sample_coverage_status"] == "pass_region_coverage").sum())
    n_sample_insufficient = int((sample_qc["region_sample_coverage_status"] == "insufficient_region_coverage").sum())
    n_sample_no = int((sample_qc["region_sample_coverage_status"] == "no_region_coverage").sum())
    n_region_pass_coverage = int((region_qc["region_coverage_status"] == "pass_common_cpg_and_region_coverage").sum())
    n_region_low_coverage = int((region_qc["region_coverage_status"] == "pass_common_cpg_but_low_region_coverage").sum())
    n_direction_changed = int(delta["direction_changed"].fillna(False).sum()) if not delta.empty else 0
    warn_for_region_coverage = args.min_region_total_coverage is not None and (
        n_sample_insufficient > 0 or n_sample_no > 0 or n_region_low_coverage > 0
    )
    summary = [
        {
            "n_regions_total": n_total,
            "n_regions_pass_common_cpg": n_pass,
            "n_regions_insufficient_common_cpg": n_insufficient,
            "n_regions_no_common_cpgs": n_no_common,
            "min_region_total_coverage": args.min_region_total_coverage,
            "n_region_sample_pass_region_coverage": n_sample_pass,
            "n_region_sample_insufficient_region_coverage": n_sample_insufficient,
            "n_region_sample_no_region_coverage": n_sample_no,
            "n_regions_pass_common_cpg_and_region_coverage": n_region_pass_coverage,
            "n_regions_pass_common_cpg_but_low_region_coverage": n_region_low_coverage,
            "mean_fraction_common": float(region_qc["fraction_common"].mean()) if n_total else np.nan,
            "median_fraction_common": float(region_qc["fraction_common"].median()) if n_total else np.nan,
            "n_direction_changed": n_direction_changed,
            "max_abs_delta_shift": float(delta["abs_delta_shift"].max()) if not delta.empty else np.nan,
            "status": "WARN" if n_direction_changed > 0 or (n_insufficient + n_no_common) > 0 or warn_for_region_coverage else "PASS",
            "notes": (
                f"delta_common is condition {condition_order[1]} minus {condition_order[0]}; "
                "coverage-set QC controls covered-position set differences but does not eliminate all biases; "
                + (
                    "min_region_total_coverage not required; only N_common > 0 checked"
                    if args.min_region_total_coverage is None
                    else f"requires N_common >= {args.min_region_total_coverage}"
                )
            )
            if len(condition_order) >= 2
            else "could not infer two-condition order",
        }
    ]

    region_qc.to_csv(out_dir / "coverage_set_region_qc.tsv", sep="\t", index=False)
    sample_qc.to_csv(out_dir / "coverage_set_sample_qc.tsv", sep="\t", index=False)
    delta.to_csv(out_dir / "coverage_set_delta_shift.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(out_dir / "coverage_set_summary.tsv", sep="\t", index=False)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
