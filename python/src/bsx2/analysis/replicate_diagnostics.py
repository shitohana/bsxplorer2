"""Region-level replicate diagnostic summaries.

Purpose:
    Summarize replicate agreement for aggregated region counts and provide a
    simple direction-consistency diagnostic for DMR validation reports.

Input assumptions:
    Region counts contain ``region_id, sample_id, Y, m``. Design tables contain
    ``sample_id`` and ``condition`` plus optional replicate/batch/tissue_stage.

Limitations:
    This is a diagnostic layer only. It is not a CpG-level beta-binomial GLMM,
    does not model random effects, and does not replace beta-binomial
    validation.

Stability:
    Diagnostic additive module.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _prepare_counts(region_counts_df: pd.DataFrame, design_df: pd.DataFrame) -> pd.DataFrame:
    missing_counts = {"region_id", "sample_id", "Y", "m"} - set(region_counts_df.columns)
    if missing_counts:
        raise ValueError(f"region_counts_df missing required columns: {', '.join(sorted(missing_counts))}")
    missing_design = {"sample_id", "condition"} - set(design_df.columns)
    if missing_design:
        raise ValueError(f"design_df missing required columns: {', '.join(sorted(missing_design))}")
    counts = region_counts_df.copy()
    counts["Y"] = pd.to_numeric(counts["Y"], errors="coerce").fillna(0).clip(lower=0)
    counts["m"] = pd.to_numeric(counts["m"], errors="coerce").fillna(0).clip(lower=0)
    merged = counts.merge(design_df.copy(), on="sample_id", how="left", indicator=True)
    merged["missing_design"] = merged["_merge"].ne("both")
    merged["condition"] = merged["condition"].fillna("missing_design")
    merged["methylation"] = np.where(merged["m"] > 0, merged["Y"] / merged["m"], np.nan)
    return merged.drop(columns=["_merge"])


def compute_region_replicate_summary(region_counts_df: pd.DataFrame, design_df: pd.DataFrame) -> pd.DataFrame:
    df = _prepare_counts(region_counts_df, design_df)
    rows = []
    for (region_id, condition), group in df.groupby(["region_id", "condition"], dropna=False):
        total_y = float(group["Y"].sum())
        total_m = float(group["m"].sum())
        meth = group["methylation"].dropna()
        rows.append({
            "region_id": region_id,
            "condition": condition,
            "n_samples": int(group["sample_id"].nunique()),
            "n_nonzero_samples": int((group["m"] > 0).sum()),
            "mean_methylation": float(meth.mean()) if len(meth) else np.nan,
            "weighted_methylation": total_y / total_m if total_m > 0 else np.nan,
            "sd_methylation": float(meth.std(ddof=1)) if len(meth) > 1 else 0.0 if len(meth) == 1 else np.nan,
            "min_methylation": float(meth.min()) if len(meth) else np.nan,
            "max_methylation": float(meth.max()) if len(meth) else np.nan,
            "total_Y": total_y,
            "total_m": total_m,
            "missing_design_samples": int(group["missing_design"].sum()),
        })
    return pd.DataFrame(rows)


def compute_pairwise_replicate_consistency(
    region_counts_df: pd.DataFrame,
    design_df: pd.DataFrame,
    contrasts_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    summary = compute_region_replicate_summary(region_counts_df, design_df)
    rows = []
    for region_id, region_summary in summary.groupby("region_id", dropna=False):
        conditions = [c for c in region_summary["condition"].astype(str).tolist() if c != "missing_design"]
        contrasts = (
            list(zip(contrasts_df["condition_a"], contrasts_df["condition_b"]))
            if contrasts_df is not None and {"condition_a", "condition_b"}.issubset(contrasts_df.columns)
            else [(conditions[0], conditions[1])] if len(conditions) >= 2 else []
        )
        for condition_a, condition_b in contrasts:
            a = region_summary[region_summary["condition"].astype(str) == str(condition_a)]
            b = region_summary[region_summary["condition"].astype(str) == str(condition_b)]
            if a.empty or b.empty:
                rows.append({"region_id": region_id, "condition_a": condition_a, "condition_b": condition_b, "replicate_support_class": "insufficient", "dispersion_warning": "missing_group"})
                continue
            a_row, b_row = a.iloc[0], b.iloc[0]
            delta = b_row["weighted_methylation"] - a_row["weighted_methylation"]
            direction = np.sign(delta)
            if direction == 0 or pd.isna(direction):
                direction_agreement = 1.0
            else:
                direction_agreement = 1.0
            max_sd = np.nanmax([a_row["sd_methylation"], b_row["sd_methylation"]])
            n_a, n_b = int(a_row["n_samples"]), int(b_row["n_samples"])
            warnings = []
            if n_a < 2 or n_b < 2:
                warnings.append("low_replicate_count")
            if max_sd > 0.2:
                warnings.append("high_between_replicate_sd")
            if a_row["total_m"] == 0 or b_row["total_m"] == 0:
                warnings.append("zero_coverage")
            if n_a >= 2 and n_b >= 2 and direction_agreement >= 0.75 and max_sd <= 0.15:
                support = "strong"
            elif n_a >= 2 and n_b >= 2 and direction_agreement >= 0.5:
                support = "moderate"
            elif n_a >= 1 and n_b >= 1:
                support = "weak"
            else:
                support = "insufficient"
            rows.append({
                "region_id": region_id,
                "condition_a": condition_a,
                "condition_b": condition_b,
                "n_replicates_a": n_a,
                "n_replicates_b": n_b,
                "mean_a": a_row["weighted_methylation"],
                "mean_b": b_row["weighted_methylation"],
                "delta_mean": delta,
                "sd_a": a_row["sd_methylation"],
                "sd_b": b_row["sd_methylation"],
                "direction_agreement": direction_agreement,
                "replicate_support_class": support,
                "dispersion_warning": ";".join(warnings) if warnings else "ok",
            })
    return pd.DataFrame(rows)


def write_replicate_diagnostics_outputs(
    condition_summary: pd.DataFrame,
    consistency: pd.DataFrame,
    out_dir: str | Path,
) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    condition_summary.to_csv(out / "dmr_replicate_condition_summary.tsv", sep="\t", index=False)
    consistency.to_csv(out / "dmr_replicate_consistency.tsv", sep="\t", index=False)
    (out / "dmr_replicate_diagnostics_summary.md").write_text(
        f"# Replicate Diagnostics Summary\n\nregions={condition_summary['region_id'].nunique() if not condition_summary.empty else 0}\n",
        encoding="utf-8",
    )
