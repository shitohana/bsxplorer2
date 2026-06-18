#!/usr/bin/env python3
"""Bootstrap confidence intervals for replicate-level delta methylation.

This is a lightweight QC/robustness layer for the descriptive effect size
Delta_R^{rep}. It resamples biological replicates within each condition and
does not perform DMR testing, compute p-values, or control genome-wide FDR.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


from dmr_validation_framework.core.coverage import aggregate_counts_per_sample_chunked  # noqa: E402
from dmr_validation_framework.core.io import ensure_out_dir, find_preferred_file, read_table, write_tsv  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--input", type=Path, help="Optional region x sample count/proportion table.")
    parser.add_argument("--design", type=Path, help="Optional design table with sample_id and condition.")
    parser.add_argument("--min-coverage", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return str(lower[candidate.lower()])
    return None


def find_file(names: list[str]) -> Path | None:
    return find_preferred_file(names)


def skipped(out_dir: Path, reason: str) -> int:
    rows = [{"status": "SKIPPED", "notes": reason}]
    write_tsv(out_dir / "delta_bootstrap_ci.tsv", rows)
    write_tsv(out_dir / "delta_bootstrap_summary.tsv", rows)
    return 0


def load_design(path: Path | None) -> pd.DataFrame | None:
    design_path = path or find_file(["design.tsv", "design.csv", "sample_metadata.tsv", "metadata.tsv"])
    if not design_path or not design_path.exists():
        return None
    design = read_table(design_path)
    sample_col = find_col(design, ["sample_id", "sample", "sample_name"])
    cond_col = find_col(design, ["condition", "group", "treatment"])
    if not sample_col or not cond_col:
        return None
    return design.rename(columns={sample_col: "sample_id", cond_col: "condition"})[["sample_id", "condition"]]


def normalize_region_sample_table(df: pd.DataFrame, design: pd.DataFrame | None) -> tuple[pd.DataFrame, str]:
    region_col = find_col(df, ["region_id", "dmr_id", "harmonized_region_id"])
    sample_col = find_col(df, ["sample_id", "sample", "sample_name"])
    cond_col = find_col(df, ["condition", "group", "treatment"])
    m_col = find_col(df, ["M", "mC", "Y", "methylated", "M_common", "M_per_sample", "sum_mC"])
    n_col = find_col(df, ["N", "total", "m", "coverage", "N_common", "N_per_sample", "sum_total"])
    mu_col = find_col(
        df,
        [
            "mu_hat",
            "mu_hat_common",
            "mu_hat_per_sample",
            "mean_methylation",
            "weighted_methylation",
            "methylation",
        ],
    )
    if not region_col or not sample_col:
        raise ValueError("region-sample table lacks region_id/sample_id columns")
    if not ((m_col and n_col) or mu_col):
        raise ValueError("region-sample table lacks M/N counts and methylation proportion columns")
    out = pd.DataFrame(
        {
            "region_id": df[region_col].astype(str),
            "sample_id": df[sample_col].astype(str),
            "condition": df[cond_col].astype(str) if cond_col else pd.NA,
        }
    )
    notes = "using count-based methylation proportions"
    if m_col and n_col:
        out["M"] = pd.to_numeric(df[m_col], errors="coerce")
        out["N"] = pd.to_numeric(df[n_col], errors="coerce")
        out["mu_hat"] = np.where(out["N"] > 0, out["M"] / out["N"], np.nan)
    else:
        out["M"] = np.nan
        out["N"] = np.nan
        out["mu_hat"] = pd.to_numeric(df[mu_col], errors="coerce")
        notes = "using precomputed methylation proportions; counts unavailable"
    if out["condition"].isna().all():
        if design is None:
            raise ValueError("condition missing and design table not found")
        out = out.drop(columns=["condition"]).merge(design, on="sample_id", how="left")
    return out, notes


def load_region_sample_values(args: argparse.Namespace) -> tuple[pd.DataFrame, str, str]:
    design = load_design(args.design)
    if args.input:
        if not args.input.exists():
            raise FileNotFoundError(f"input table does not exist: {args.input}")
        values, notes = normalize_region_sample_table(read_table(args.input), design)
        return values, str(args.input), notes

    region_sample_path = find_file(
        [
            "aggregated_region_counts.tsv",
            "region_sample_counts.tsv",
            "aggregated_glm_input.tsv",
            "region_counts.tsv",
            "coverage_set_common_region_counts.tsv",
            "coverage_set_per_sample_region_counts.tsv",
            "region_level_counts.tsv",
        ]
    )
    if region_sample_path:
        values, notes = normalize_region_sample_table(read_table(region_sample_path), design)
        return values, str(region_sample_path), notes

    cpg_path = find_file(["region_cpg_counts.tsv", "selected_region_cpg_counts.tsv", "extracted_cpg_counts.tsv"])
    if cpg_path:
        if design is None:
            raise ValueError("region CpG counts found but design table was not found")
        per = aggregate_counts_per_sample_chunked(
            cpg_path,
            design=design,
            condition_column="condition",
            min_coverage=args.min_coverage,
            chunk_size=args.chunk_size,
        )
        values = per.rename(columns={"M_per_sample": "M", "N_per_sample": "N"})[
            ["region_id", "sample_id", "condition", "M", "N", "mu_hat_per_sample"]
        ].rename(columns={"mu_hat_per_sample": "mu_hat"})
        return values, str(cpg_path), "aggregated per-CpG counts to region x sample proportions"

    raise FileNotFoundError("no usable region x sample count/proportion table found")


def direction(value: float) -> str:
    if pd.isna(value) or value == 0:
        return "zero"
    return "positive" if value > 0 else "negative"


def bootstrap_delta(a_values: np.ndarray, b_values: np.ndarray, n_bootstrap: int, rng: np.random.Generator) -> np.ndarray:
    a_idx = rng.integers(0, len(a_values), size=(n_bootstrap, len(a_values)))
    b_idx = rng.integers(0, len(b_values), size=(n_bootstrap, len(b_values)))
    return b_values[b_idx].mean(axis=1) - a_values[a_idx].mean(axis=1)


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    if args.n_bootstrap <= 0:
        return skipped(out_dir, "--n-bootstrap must be positive")
    if not 0 < args.alpha < 1:
        return skipped(out_dir, "--alpha must be between 0 and 1")
    try:
        values, source, source_notes = load_region_sample_values(args)
    except Exception as exc:
        return skipped(out_dir, str(exc))

    values = values.dropna(subset=["condition", "mu_hat"]).copy()
    if values.empty:
        return skipped(out_dir, "no rows with condition and methylation proportion")
    groups = list(values["condition"].astype(str).drop_duplicates())
    group_note = ""
    if len(groups) < 2:
        return skipped(out_dir, "fewer than two conditions in input")
    if len(groups) > 2:
        groups = sorted(groups)[:2]
        values = values[values["condition"].astype(str).isin(groups)].copy()
        group_note = f"more than two groups found; using first two sorted groups {groups}"
    else:
        group_note = f"delta is {groups[1]} minus {groups[0]}"

    rng = np.random.default_rng(args.seed)
    rows: list[dict] = []
    n_skipped = 0
    small_n_warnings = 0
    for region_id, sub in values.groupby("region_id", sort=False):
        a = pd.to_numeric(sub[sub["condition"].astype(str) == groups[0]]["mu_hat"], errors="coerce").dropna().to_numpy()
        b = pd.to_numeric(sub[sub["condition"].astype(str) == groups[1]]["mu_hat"], errors="coerce").dropna().to_numpy()
        if len(a) == 0 or len(b) == 0:
            n_skipped += 1
            rows.append(
                {
                    "region_id": region_id,
                    "delta_rep": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "ci_width": np.nan,
                    "ci_includes_zero": pd.NA,
                    "direction_observed": "zero",
                    "direction_stable": pd.NA,
                    "n_samples_A": len(a),
                    "n_samples_B": len(b),
                    "n_bootstrap": args.n_bootstrap,
                    "alpha": args.alpha,
                    "status": "SKIPPED",
                    "notes": "missing one condition for region",
                }
            )
            continue
        delta_rep = float(np.mean(b) - np.mean(a))
        note_parts = [group_note, source_notes, f"source={source}"]
        status = "PASS"
        if len(a) < 2 or len(b) < 2:
            n_skipped += 1
            status = "WARN_TOO_FEW_REPLICATES"
            note_parts.append("fewer than two biological replicates in at least one condition")
            ci_low = ci_high = ci_width = np.nan
            ci_includes_zero = pd.NA
            direction_stable = pd.NA
        else:
            if len(a) <= 2 or len(b) <= 2:
                small_n_warnings += 1
                status = "WARN_SMALL_N"
                note_parts.append("bootstrap limited by small number of biological replicates")
            deltas = bootstrap_delta(a, b, args.n_bootstrap, rng)
            ci_low = float(np.quantile(deltas, args.alpha / 2))
            ci_high = float(np.quantile(deltas, 1 - args.alpha / 2))
            ci_width = ci_high - ci_low
            ci_includes_zero = bool(ci_low <= 0 <= ci_high)
            direction_stable = bool(not ci_includes_zero)
        rows.append(
            {
                "region_id": region_id,
                "delta_rep": delta_rep,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "ci_width": ci_width,
                "ci_includes_zero": ci_includes_zero,
                "direction_observed": direction(delta_rep),
                "direction_stable": direction_stable,
                "n_samples_A": len(a),
                "n_samples_B": len(b),
                "n_bootstrap": args.n_bootstrap,
                "alpha": args.alpha,
                "status": status,
                "notes": "; ".join(note_parts),
            }
        )

    result = pd.DataFrame(rows)
    ok = result[~result["ci_width"].isna()].copy()
    if ok.empty:
        status = "SKIPPED"
        notes = "no regions had enough replicates for bootstrap CI"
    else:
        max_width = float(ok["ci_width"].max())
        frac_stable = float(ok["direction_stable"].fillna(False).mean())
        status = "WARN" if small_n_warnings > 0 or max_width > 0.5 or frac_stable < 0.5 else "PASS"
        notes = "Bootstrap CI is diagnostic for effect-size uncertainty only; it is not a p-value, DMR test, or FDR-controlled significance procedure."
        if small_n_warnings:
            notes += " bootstrap limited by small number of biological replicates."
    summary = [
        {
            "n_regions": int(len(result)),
            "n_regions_bootstrap_ok": int(len(ok)),
            "n_regions_direction_stable": int(ok["direction_stable"].fillna(False).sum()) if not ok.empty else 0,
            "fraction_direction_stable": float(ok["direction_stable"].fillna(False).mean()) if not ok.empty else np.nan,
            "mean_ci_width": float(ok["ci_width"].mean()) if not ok.empty else np.nan,
            "median_ci_width": float(ok["ci_width"].median()) if not ok.empty else np.nan,
            "max_ci_width": float(ok["ci_width"].max()) if not ok.empty else np.nan,
            "n_warn_too_few_replicates": int(small_n_warnings + (result["status"] == "WARN_TOO_FEW_REPLICATES").sum()),
            "n_skipped": int(n_skipped),
            "status": status,
            "notes": notes,
        }
    ]
    result.to_csv(out_dir / "delta_bootstrap_ci.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(out_dir / "delta_bootstrap_summary.tsv", sep="\t", index=False)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
