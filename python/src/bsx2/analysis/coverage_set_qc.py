"""Coverage-set QC for DMR/metagene downstream validation.

This module makes the CpG/cytosine set used for region summaries explicit.
Per-sample aggregation is preserved for descriptive summaries, while
strict/common-CpG aggregation can be used as an optional QC/statistical
validation input.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd


CANONICAL_CPG_COLUMNS = [
    "region_id",
    "cpg_id",
    "sample_id",
    "condition",
    "chrom",
    "position",
    "context",
    "mC",
    "uC",
    "total",
]


def _find_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def normalize_cpg_counts_table(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize long per-CpG counts to a stable schema.

    Output columns:
    region_id, cpg_id, sample_id, condition, chrom, position, context, mC, uC, total.
    """
    raw = df.copy()
    aliases = {
        "region_id": ("region_id", "dmr_id", "harmonized_region_id"),
        "cpg_id": ("cpg_id", "cytosine_id", "site_id"),
        "sample_id": ("sample_id", "sample", "sample_name"),
        "condition": ("condition", "group", "treatment"),
        "mC": ("mC", "mc", "count_m", "methylated", "Y"),
        "uC": ("uC", "uc", "count_u", "unmethylated"),
        "total": ("total", "count_total", "coverage", "N", "n"),
        "chrom": ("chrom", "chr", "seqname"),
        "position": ("position", "pos", "start"),
        "context": ("context", "ctx"),
    }
    out = pd.DataFrame(index=raw.index)
    for target, candidates in aliases.items():
        col = _find_column(raw, candidates)
        if col is not None:
            out[target] = raw[col]
    missing = sorted({"region_id", "sample_id", "mC", "total"} - set(out.columns))
    if missing:
        raise ValueError(f"CpG counts table is missing required columns or aliases: {', '.join(missing)}")
    if "chrom" not in out.columns:
        out["chrom"] = "NA"
    if "position" not in out.columns:
        out["position"] = pd.NA
    if "context" not in out.columns:
        out["context"] = "NA"
    if "condition" not in out.columns:
        out["condition"] = pd.NA
    if "cpg_id" not in out.columns:
        out["cpg_id"] = (
            out["chrom"].astype(str)
            + ":"
            + pd.to_numeric(out["position"], errors="coerce").fillna(-1).astype(int).astype(str)
            + ":"
            + out["context"].astype(str)
        )
    out["mC"] = pd.to_numeric(out["mC"], errors="coerce").fillna(0)
    out["total"] = pd.to_numeric(out["total"], errors="coerce").fillna(0)
    if "uC" not in out.columns:
        out["uC"] = out["total"] - out["mC"]
    else:
        out["uC"] = pd.to_numeric(out["uC"], errors="coerce").fillna(out["total"] - out["mC"])
    out["region_id"] = out["region_id"].astype(str)
    out["cpg_id"] = out["cpg_id"].astype(str)
    out["sample_id"] = out["sample_id"].astype(str)
    out["condition"] = out["condition"].astype("string")
    out["chrom"] = out["chrom"].astype(str)
    out["context"] = out["context"].astype(str)
    out["position"] = pd.to_numeric(out["position"], errors="coerce")
    return out[CANONICAL_CPG_COLUMNS]


def _with_design(
    cpg_counts_df: pd.DataFrame,
    design_df: pd.DataFrame | None,
    *,
    sample_col: str = "sample_id",
    condition_col: str = "condition",
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], str]:
    counts = normalize_cpg_counts_table(cpg_counts_df)
    notes: list[str] = []
    if design_df is None:
        if counts["condition"].isna().all():
            raise ValueError("design_df is required when counts table has no condition column")
        design = counts[[sample_col, condition_col]].dropna().drop_duplicates().copy()
    else:
        design = design_df.copy()
        sample_alias = _find_column(design, (sample_col, "sample_id", "sample", "sample_name"))
        condition_alias = _find_column(design, (condition_col, "condition", "group", "treatment"))
        if sample_alias is None or condition_alias is None:
            raise ValueError("design_df must contain sample_id/sample and condition/group columns")
        design = design.rename(columns={sample_alias: sample_col, condition_alias: condition_col})
        design = design[[sample_col, condition_col]].dropna().drop_duplicates().copy()
    design[sample_col] = design[sample_col].astype(str)
    design[condition_col] = design[condition_col].astype(str)
    counts = counts.drop(columns=["condition"]).merge(design, on=sample_col, how="left")
    if counts[condition_col].isna().any():
        notes.append("some count rows did not match design samples and were excluded")
        counts = counts.dropna(subset=[condition_col]).copy()
    groups = sorted(design[condition_col].dropna().astype(str).unique())
    if len(groups) < 2:
        raise ValueError("coverage-set QC requires at least two conditions")
    if len(groups) > 2:
        notes.append(f"more than two groups found; using first two sorted groups: {groups[:2]}")
        groups = groups[:2]
        design = design[design[condition_col].isin(groups)].copy()
        counts = counts[counts[condition_col].isin(groups)].copy()
    return counts, design, groups, "; ".join(notes)


def _resolve_group_thresholds(
    design: pd.DataFrame,
    groups: list[str],
    *,
    condition_col: str,
    min_covered_per_group: int | dict[str, int] | None,
) -> tuple[int, int]:
    if isinstance(min_covered_per_group, dict):
        return (
            int(min_covered_per_group.get(groups[0], min_covered_per_group.get("A", 1))),
            int(min_covered_per_group.get(groups[1], min_covered_per_group.get("B", 1))),
        )
    if min_covered_per_group is not None:
        k = int(min_covered_per_group)
        return k, k
    counts = design.groupby(condition_col).size()
    return int(counts.loc[groups[0]]), int(counts.loc[groups[1]])


def build_cpg_coverage_qc(
    cpg_counts_df: pd.DataFrame,
    design_df: pd.DataFrame | None = None,
    region_col: str = "region_id",
    cpg_col: str = "cpg_id",
    sample_col: str = "sample_id",
    condition_col: str = "condition",
    total_col: str = "total",
    min_coverage: int = 1,
    min_covered_per_group: int | dict[str, int] | None = None,
    min_common_cpgs: int = 3,
    min_region_total_coverage: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build common-CpG flags and region-level coverage QC tables."""
    counts, design, groups, design_note = _with_design(
        cpg_counts_df,
        design_df,
        sample_col=sample_col,
        condition_col=condition_col,
    )
    k_a, k_b = _resolve_group_thresholds(
        design,
        groups,
        condition_col=condition_col,
        min_covered_per_group=min_covered_per_group,
    )
    counts = counts.rename(columns={region_col: "region_id", cpg_col: "cpg_id", sample_col: "sample_id", total_col: "total"})
    counts["_covered"] = pd.to_numeric(counts["total"], errors="coerce").fillna(0) >= int(min_coverage)
    covered = (
        counts[counts["_covered"]]
        .drop_duplicates(["region_id", "cpg_id", "sample_id", condition_col])
        .groupby(["region_id", "cpg_id", condition_col])
        .size()
        .unstack(condition_col, fill_value=0)
    )
    union_index = counts[["region_id", "cpg_id"]].drop_duplicates().set_index(["region_id", "cpg_id"]).index
    covered = covered.reindex(union_index, fill_value=0)
    for group in groups:
        if group not in covered.columns:
            covered[group] = 0
    common = covered.reset_index()
    common["covered_A"] = common[groups[0]].astype(int)
    common["covered_B"] = common[groups[1]].astype(int)
    common["is_common"] = (common["covered_A"] >= k_a) & (common["covered_B"] >= k_b)

    def reason(row: pd.Series) -> str:
        if bool(row["is_common"]):
            return "common"
        parts = []
        if int(row["covered_A"]) < k_a:
            parts.append(f"{groups[0]}_coverage_below_k")
        if int(row["covered_B"]) < k_b:
            parts.append(f"{groups[1]}_coverage_below_k")
        return ";".join(parts) or "not_common"

    common["reason"] = common.apply(reason, axis=1)
    common_cpg_df = common[["region_id", "cpg_id", "covered_A", "covered_B", "is_common", "reason"]].copy()

    sample_counts = design.groupby(condition_col)["sample_id"].nunique()
    rows = []
    for region_id, sub in common_cpg_df.groupby("region_id"):
        n_union = int(len(sub))
        n_common = int(sub["is_common"].sum())
        if n_common >= int(min_common_cpgs):
            status = "pass_common_cpg"
        elif n_common > 0:
            status = "insufficient_common_cpgs"
        else:
            status = "no_common_cpgs"
        rows.append(
            {
                "region_id": region_id,
                "n_cpg_union": n_union,
                "n_cpg_common": n_common,
                "fraction_common": n_common / n_union if n_union else 0.0,
                "n_samples_A": int(sample_counts.loc[groups[0]]),
                "n_samples_B": int(sample_counts.loc[groups[1]]),
                "k_A": k_a,
                "k_B": k_b,
                "min_common_cpgs": int(min_common_cpgs),
                "min_region_total_coverage": min_region_total_coverage,
                "coverage_qc_status": status,
                "region_coverage_status": (
                    "pass_common_cpg_region_coverage_not_required"
                    if status == "pass_common_cpg"
                    else status
                ),
                "notes": f"A={groups[0]}; B={groups[1]}; min_coverage={min_coverage}; {design_note}".strip("; "),
            }
        )
    return common_cpg_df, pd.DataFrame(rows)


def aggregate_counts_per_sample(
    cpg_counts_df: pd.DataFrame,
    min_coverage: int = 1,
) -> pd.DataFrame:
    """Aggregate region/sample counts over each sample-specific covered set I_s(R)."""
    counts = normalize_cpg_counts_table(cpg_counts_df)
    covered = counts[pd.to_numeric(counts["total"], errors="coerce").fillna(0) >= int(min_coverage)].copy()
    out = (
        covered.groupby(["region_id", "sample_id", "condition"], dropna=False)
        .agg(
            M_per_sample=("mC", "sum"),
            U_per_sample=("uC", "sum"),
            N_per_sample=("total", "sum"),
            n_cpg_per_sample=("cpg_id", "nunique"),
        )
        .reset_index()
    )
    out["mu_hat_per_sample"] = np.where(out["N_per_sample"] > 0, out["M_per_sample"] / out["N_per_sample"], np.nan)
    return out[
        [
            "region_id",
            "sample_id",
            "condition",
            "M_per_sample",
            "U_per_sample",
            "N_per_sample",
            "mu_hat_per_sample",
            "n_cpg_per_sample",
        ]
    ]


def aggregate_counts_common_cpg(
    cpg_counts_df: pd.DataFrame,
    common_cpg_df: pd.DataFrame,
    design_df: pd.DataFrame | None = None,
    sample_col: str = "sample_id",
    condition_col: str = "condition",
    min_common_cpgs: int = 3,
    min_region_total_coverage: int | None = None,
) -> pd.DataFrame:
    """Aggregate region/sample counts over the shared I_common(R) set."""
    counts, design, _groups, _note = _with_design(cpg_counts_df, design_df, sample_col=sample_col, condition_col=condition_col)
    common = common_cpg_df.copy()
    common = common[common["is_common"].astype(bool)][["region_id", "cpg_id"]].drop_duplicates()
    n_common = common.groupby("region_id")["cpg_id"].nunique().rename("n_common_cpgs")
    filtered = counts.merge(common, on=["region_id", "cpg_id"], how="inner")
    grouped = (
        filtered.groupby(["region_id", "sample_id"], dropna=False)
        .agg(M_common=("mC", "sum"), U_common=("uC", "sum"), N_common=("total", "sum"))
        .reset_index()
    )
    regions = pd.Series(common["region_id"].drop_duplicates(), name="region_id")
    samples = design[["sample_id", "condition"]].drop_duplicates()
    complete = regions.to_frame().merge(samples, how="cross") if not regions.empty else pd.DataFrame(columns=["region_id", "sample_id", "condition"])
    out = complete.merge(grouped, on=["region_id", "sample_id"], how="left")
    for col in ("M_common", "U_common", "N_common"):
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0)
    out = out.merge(n_common.reset_index(), on="region_id", how="left")
    out["n_common_cpgs"] = out["n_common_cpgs"].fillna(0).astype(int)
    out["mu_hat_common"] = np.where(out["N_common"] > 0, out["M_common"] / out["N_common"], np.nan)
    out["coverage_qc_status"] = np.select(
        [out["n_common_cpgs"] >= int(min_common_cpgs), out["n_common_cpgs"] > 0],
        ["pass_common_cpg", "insufficient_common_cpgs"],
        default="no_common_cpgs",
    )
    if min_region_total_coverage is None:
        out["region_sample_coverage_status"] = np.where(
            out["N_common"] > 0,
            "pass_region_coverage",
            "no_region_coverage",
        )
    else:
        threshold = int(min_region_total_coverage)
        out["region_sample_coverage_status"] = np.select(
            [out["N_common"] >= threshold, out["N_common"] > 0],
            ["pass_region_coverage", "insufficient_region_coverage"],
            default="no_region_coverage",
        )
    out["min_region_total_coverage"] = min_region_total_coverage
    return out[
        [
            "region_id",
            "sample_id",
            "condition",
            "M_common",
            "U_common",
            "N_common",
            "mu_hat_common",
            "n_common_cpgs",
            "coverage_qc_status",
            "min_region_total_coverage",
            "region_sample_coverage_status",
        ]
    ]


def compare_per_sample_vs_common_sets(
    per_sample_region_df: pd.DataFrame,
    common_region_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compare per-sample regional methylation with common-CpG aggregation."""
    merged = per_sample_region_df.merge(
        common_region_df[["region_id", "sample_id", "condition", "mu_hat_common", "n_common_cpgs", "coverage_qc_status"]],
        on=["region_id", "sample_id", "condition"],
        how="outer",
    )
    merged["delta_mu_common_minus_per_sample"] = pd.to_numeric(merged["mu_hat_common"], errors="coerce") - pd.to_numeric(
        merged["mu_hat_per_sample"], errors="coerce"
    )
    shift = merged["delta_mu_common_minus_per_sample"].abs()
    merged["coverage_set_shift_status"] = np.select(
        [shift.isna(), shift > 0.1, shift > 0.05],
        ["not_evaluable", "large_shift", "moderate_shift"],
        default="small_or_no_shift",
    )
    return merged[
        [
            "region_id",
            "sample_id",
            "condition",
            "mu_hat_per_sample",
            "mu_hat_common",
            "delta_mu_common_minus_per_sample",
            "n_cpg_per_sample",
            "n_common_cpgs",
            "coverage_set_shift_status",
        ]
    ]
