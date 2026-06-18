#!/usr/bin/env python3
"""Check replicate-level versus pooled read-level delta sensitivity."""

from __future__ import annotations

import argparse
import html
from pathlib import Path

import numpy as np
import pandas as pd


from dmr_validation_framework.core.io import ensure_out_dir, write_tsv  # noqa: E402
from dmr_validation_framework.core.delta_weighting import (
    bootstrap_delta_stats,
    classify_confidence,
    compute_loo,
    direction_changed,
    dispersion_score,
    load_cpg_coverage_metrics,
    load_glm_glmm_agreement,
    load_region_sample_counts,
    safe_ratio,
    sign,
)
from dmr_validation_framework.reports.delta_weighting import write_html_report


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-sample-counts", type=Path)
    parser.add_argument("--region-cpg-counts", type=Path)
    parser.add_argument("--design", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--condition-column", default="condition")
    parser.add_argument("--min-coverage", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--bootstrap-alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--glm-glmm-comparison", type=Path)
    parser.add_argument("--glmm-results", type=Path)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--max-abs-delta-shift", type=float, default=0.05)
    parser.add_argument("--max-relative-delta-shift", type=float, default=1.0)
    parser.add_argument("--max-coverage-imbalance", type=float, default=3.0)
    parser.add_argument("--no-html-report", action="store_true")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def skipped(out_dir: Path, reason: str) -> int:
    row = {"status": "SKIPPED", "notes": reason}
    write_tsv(out_dir / "delta_weighting_sensitivity.tsv", [row])
    write_tsv(out_dir / "delta_weighting_loo.tsv", [row])
    write_tsv(out_dir / "delta_weighting_summary.tsv", [row])
    (out_dir / "delta_weighting_summary.html").write_text(
        f"<html><body><h1>Delta weighting sensitivity</h1><p>SKIPPED: {html.escape(reason)}</p></body></html>\n",
        encoding="utf-8",
    )
    return 0


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    if args.n_bootstrap <= 0:
        return skipped(out_dir, "--n-bootstrap must be positive")
    if not 0 < args.bootstrap_alpha < 1:
        return skipped(out_dir, "--bootstrap-alpha must be between 0 and 1")
    try:
        counts, source = load_region_sample_counts(args)
    except Exception as exc:
        return skipped(out_dir, str(exc))
    counts = counts.dropna(subset=["condition"]).copy()
    groups = list(counts["condition"].astype(str).drop_duplicates())
    if len(groups) < 2:
        return skipped(out_dir, "fewer than two conditions in region-sample counts")
    if len(groups) > 2:
        groups = sorted(groups)[:2]
        counts = counts[counts["condition"].astype(str).isin(groups)].copy()
        group_note = f"more than two groups found; using first two sorted groups {groups}"
    else:
        group_note = f"delta is {groups[1]} minus {groups[0]}"
    counts["mu_hat"] = np.where(counts["N"] > 0, counts["M"] / counts["N"], np.nan)
    cpg_metrics = load_cpg_coverage_metrics(args, groups, source)
    cpg_by_region = cpg_metrics.set_index("region_id").to_dict(orient="index") if not cpg_metrics.empty else {}
    glm_glmm, glm_glmm_source = load_glm_glmm_agreement(args)
    glm_by_region = glm_glmm.set_index("region_id").to_dict(orient="index") if not glm_glmm.empty else {}
    rng = np.random.default_rng(args.seed)
    rows: list[dict] = []
    loo_rows: list[dict] = []
    for region_id, sub in counts.groupby("region_id"):
        a = sub[sub["condition"].astype(str) == groups[0]]
        b = sub[sub["condition"].astype(str) == groups[1]]
        if a.empty or b.empty:
            continue
        delta_rep = b["mu_hat"].mean() - a["mu_hat"].mean()
        mu_pool_a = a["M"].sum() / a["N"].sum() if a["N"].sum() > 0 else np.nan
        mu_pool_b = b["M"].sum() / b["N"].sum() if b["N"].sum() > 0 else np.nan
        delta_pool = mu_pool_b - mu_pool_a
        direction_rep = sign(delta_rep)
        direction_pool = sign(delta_pool)
        delta_shift = delta_pool - delta_rep
        abs_delta_shift = abs(delta_shift)
        direction_changed_flag = direction_changed(delta_rep, delta_pool)
        loo_stats, region_loo_rows = compute_loo(str(region_id), sub, groups, delta_rep)
        loo_rows.extend(region_loo_rows)
        a_values = pd.to_numeric(a["mu_hat"], errors="coerce").dropna().to_numpy(dtype=float)
        b_values = pd.to_numeric(b["mu_hat"], errors="coerce").dropna().to_numpy(dtype=float)
        bootstrap = bootstrap_delta_stats(
            a_values,
            b_values,
            n_bootstrap=args.n_bootstrap,
            alpha=args.bootstrap_alpha,
            rng=rng,
        )
        mean_coverage_a = float(pd.to_numeric(a["N"], errors="coerce").mean())
        mean_coverage_b = float(pd.to_numeric(b["N"], errors="coerce").mean())
        coverage_ratio = safe_ratio(mean_coverage_b, mean_coverage_a)
        coverage_imbalance_score = max(coverage_ratio, 1 / coverage_ratio) if pd.notna(coverage_ratio) and coverage_ratio > 0 else np.nan
        cpg_row = cpg_by_region.get(str(region_id), {})
        if not cpg_row and "n_cpg" in sub.columns:
            a_cpg = pd.to_numeric(a["n_cpg"], errors="coerce").dropna()
            b_cpg = pd.to_numeric(b["n_cpg"], errors="coerce").dropna()
            cpg_row = {
                "n_cpg_A": int(a_cpg.max()) if len(a_cpg) else pd.NA,
                "n_cpg_B": int(b_cpg.max()) if len(b_cpg) else pd.NA,
                "n_cpg_common": pd.NA,
                "fraction_common_cpg": pd.NA,
                "coverage_metric_source": "region_sample_counts_n_cpg",
            }
        dispersion_a = dispersion_score(a)
        dispersion_b = dispersion_score(b)
        agreement = {
            "glm_q_value": pd.NA,
            "glmm_q_value": pd.NA,
            "glm_significant": pd.NA,
            "glmm_significant": pd.NA,
            "same_direction": pd.NA,
            "model_agreement_status": "NOT_AVAILABLE",
            "glm_status": "not_available",
            "glmm_status": "not_available",
            "glm_glmm_source": glm_glmm_source,
        }
        agreement.update(glm_by_region.get(str(region_id), {}))
        row = {
            "region_id": region_id,
            "delta_rep": delta_rep,
            "delta_pool": delta_pool,
            "delta_shift": delta_shift,
            "abs_delta_shift": abs_delta_shift,
            "relative_delta_shift": safe_ratio(abs_delta_shift, abs(delta_rep)),
            "direction_rep": direction_rep,
            "direction_pool": direction_pool,
            "direction_changed": direction_changed_flag,
            "loo_min_delta": loo_stats["loo_min_delta"],
            "loo_max_delta": loo_stats["loo_max_delta"],
            "loo_direction_changed": loo_stats["loo_direction_changed"],
            "most_influential_sample": loo_stats["most_influential_sample"],
            "most_influential_sample_abs_shift": loo_stats["most_influential_sample_abs_shift"],
            **bootstrap,
            "n_samples_A": a["sample_id"].nunique(),
            "n_samples_B": b["sample_id"].nunique(),
            "total_coverage_A": a["N"].sum(),
            "total_coverage_B": b["N"].sum(),
            "mean_coverage_A": mean_coverage_a,
            "mean_coverage_B": mean_coverage_b,
            "coverage_ratio": coverage_ratio,
            "coverage_ratio_B_over_A": b["N"].sum() / a["N"].sum() if a["N"].sum() > 0 else np.nan,
            "coverage_imbalance_score": coverage_imbalance_score,
            "n_cpg_A": cpg_row.get("n_cpg_A", pd.NA),
            "n_cpg_B": cpg_row.get("n_cpg_B", pd.NA),
            "n_cpg_common": cpg_row.get("n_cpg_common", pd.NA),
            "fraction_common_cpg": cpg_row.get("fraction_common_cpg", pd.NA),
            "dispersion_A": dispersion_a,
            "dispersion_B": dispersion_b,
            "overdispersion_score": np.nanmax([dispersion_a, dispersion_b]) if not pd.isna(dispersion_a) or not pd.isna(dispersion_b) else np.nan,
            **agreement,
            "status": "PASS",
            "notes": (
                f"{group_note}; source={source}; "
                f"coverage_metrics={cpg_row.get('coverage_metric_source', 'not_available')}; "
                f"glm_glmm={glm_glmm_source}; "
                "replicate-level delta is primary; pooled, LOO and bootstrap are robustness diagnostics"
            ),
        }
        final_class, reasons = classify_confidence(row, args)
        row["final_confidence_class"] = final_class
        row["confidence_reasons"] = reasons
        rows.append(row)
    if not rows:
        return skipped(out_dir, "no evaluable region-sample groups")
    out = pd.DataFrame(rows)
    max_shift = float(out["abs_delta_shift"].max())
    max_relative_shift = float(pd.to_numeric(out["relative_delta_shift"], errors="coerce").max())
    n_changed = int(out["direction_changed"].sum())
    n_loo_changed = int(out["loo_direction_changed"].fillna(False).astype(bool).sum())
    n_ci_includes_zero = int(out["ci_includes_zero"].fillna(False).astype(bool).sum())
    n_high_coverage_imbalance = int((pd.to_numeric(out["coverage_imbalance_score"], errors="coerce") > args.max_coverage_imbalance).sum())
    n_low_or_weak = int(out["final_confidence_class"].isin(["LOW_CONFIDENCE", "WEAK_CONFIDENCE"]).sum())
    status = (
        "WARN"
        if n_changed > 0
        or max_shift > args.max_abs_delta_shift
        or max_relative_shift > args.max_relative_delta_shift
        or n_loo_changed > 0
        or n_ci_includes_zero > 0
        or n_high_coverage_imbalance > 0
        or n_low_or_weak > 0
        else "PASS"
    )
    summary = pd.DataFrame(
        [
            {
                "n_regions": len(out),
                "mean_abs_delta_shift": float(out["abs_delta_shift"].mean()),
                "median_abs_delta_shift": float(out["abs_delta_shift"].median()),
                "max_abs_delta_shift": max_shift,
                "mean_relative_delta_shift": float(pd.to_numeric(out["relative_delta_shift"], errors="coerce").mean()),
                "median_relative_delta_shift": float(pd.to_numeric(out["relative_delta_shift"], errors="coerce").median()),
                "max_relative_delta_shift": max_relative_shift,
                "n_direction_changed": n_changed,
                "fraction_direction_changed": n_changed / len(out),
                "n_loo_direction_changed": n_loo_changed,
                "n_ci_includes_zero": n_ci_includes_zero,
                "n_high_coverage_imbalance": n_high_coverage_imbalance,
                "n_glmm_confirmed": int((out["model_agreement_status"] == "GLMM_CONFIRMED").sum()),
                "n_glm_only": int((out["model_agreement_status"] == "GLM_ONLY").sum()),
                "n_direction_discordant": int((out["model_agreement_status"] == "DIRECTION_DISCORDANT").sum()),
                "n_low_confidence": int((out["final_confidence_class"] == "LOW_CONFIDENCE").sum()),
                "n_weak_confidence": int((out["final_confidence_class"] == "WEAK_CONFIDENCE").sum()),
                "n_moderate_confidence": int((out["final_confidence_class"] == "MODERATE_CONFIDENCE").sum()),
                "n_high_confidence": int((out["final_confidence_class"] == "HIGH_CONFIDENCE").sum()),
                "status": status,
                "notes": "replicate-level delta gives biological replicates equal weight; pooled, LOO, bootstrap, coverage and GLM/GLMM fields are diagnostic robustness checks",
            }
        ]
    )
    out.to_csv(out_dir / "delta_weighting_sensitivity.tsv", sep="\t", index=False)
    write_tsv(out_dir / "delta_weighting_loo.tsv", loo_rows or [{"status": "SKIPPED", "notes": "no LOO rows"}])
    summary.to_csv(out_dir / "delta_weighting_summary.tsv", sep="\t", index=False)
    if not args.no_html_report:
        write_html_report(out_dir, out, summary)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
