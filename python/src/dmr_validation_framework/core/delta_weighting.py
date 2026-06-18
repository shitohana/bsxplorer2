"""Delta-weighting robustness computations."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from dmr_validation_framework.core.coverage import aggregate_counts_per_sample_chunked, normalize_cpg_counts_table
from dmr_validation_framework.core.io import find_preferred_file, read_table

def find_file(names: list[str]) -> Path | None:
    return find_preferred_file(names)

def find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None

def load_region_sample_counts(args: argparse.Namespace) -> tuple[pd.DataFrame, str]:
    if args.region_sample_counts and args.region_sample_counts.exists():
        df = read_table(args.region_sample_counts)
        source = str(args.region_sample_counts)
    else:
        cpg_path = args.region_cpg_counts or find_file(["region_cpg_counts.tsv"])
        design_path = args.design or find_file(["design.tsv"])
        if cpg_path and design_path:
            design = read_table(design_path)
            sample_col = find_col(design, ["sample_id", "sample"])
            cond_col = find_col(design, [args.condition_column, "condition", "group", "treatment"])
            if not sample_col or not cond_col:
                raise ValueError("design table lacks sample/condition columns")
            design = design.rename(columns={sample_col: "sample_id", cond_col: "condition"})[["sample_id", "condition"]]
            per = aggregate_counts_per_sample_chunked(
                cpg_path,
                design=design,
                condition_column="condition",
                min_coverage=args.min_coverage,
                chunk_size=getattr(args, "chunk_size", 1_000_000),
            )
            per = per.rename(
                columns={
                    "M_per_sample": "M",
                    "N_per_sample": "N",
                    "U_per_sample": "U",
                    "n_cpg_per_sample": "n_cpg",
                }
            )
            return per[["region_id", "sample_id", "condition", "M", "U", "N", "n_cpg"]], str(cpg_path)
        raise FileNotFoundError("no region-sample counts or region CpG counts found")

    region_col = find_col(df, ["region_id", "dmr_id", "harmonized_region_id"])
    sample_col = find_col(df, ["sample_id", "sample", "sample_name"])
    cond_col = find_col(df, [args.condition_column, "condition", "group", "treatment"])
    m_col = find_col(df, ["M", "mC", "Y", "methylated", "M_common", "M_per_sample", "sum_mC"])
    u_col = find_col(df, ["U", "uC", "unmethylated", "U_common", "U_per_sample", "sum_uC"])
    n_col = find_col(df, ["N", "total", "m", "coverage", "N_common", "N_per_sample", "sum_total"])
    n_cpg_col = find_col(df, ["n_cpg", "n_cpg_per_sample", "n_sites", "n_cytosines", "n_common_cpgs"])
    if not region_col or not sample_col or not m_col or not n_col:
        raise ValueError("region-sample table lacks region/sample/M/N columns")
    out = pd.DataFrame(
        {
            "region_id": df[region_col].astype(str),
            "sample_id": df[sample_col].astype(str),
            "condition": df[cond_col].astype(str) if cond_col else pd.NA,
            "M": pd.to_numeric(df[m_col], errors="coerce").fillna(0),
            "N": pd.to_numeric(df[n_col], errors="coerce").fillna(0),
        }
    )
    out["U"] = pd.to_numeric(df[u_col], errors="coerce").fillna(0) if u_col else out["N"] - out["M"]
    out["n_cpg"] = pd.to_numeric(df[n_cpg_col], errors="coerce") if n_cpg_col else pd.NA
    if out["condition"].isna().all():
        design_path = args.design or find_file(["design.tsv"])
        if not design_path:
            raise ValueError("condition missing and design.tsv not found")
        design = read_table(design_path)
        sample_alias = find_col(design, ["sample_id", "sample"])
        cond_alias = find_col(design, [args.condition_column, "condition", "group", "treatment"])
        out = out.drop(columns=["condition"]).merge(
            design.rename(columns={sample_alias: "sample_id", cond_alias: "condition"})[["sample_id", "condition"]],
            on="sample_id",
            how="left",
        )
    return out, source

def sign(value: float) -> str:
    if pd.isna(value) or value == 0:
        return "zero"
    return "positive" if value > 0 else "negative"

def safe_ratio(numerator: float | int, denominator: float | int) -> float:
    if pd.isna(numerator) or pd.isna(denominator) or float(denominator) == 0:
        return np.nan
    return float(numerator) / float(denominator)

def direction_changed(base: float, alternative: float) -> bool:
    base_direction = sign(base)
    alt_direction = sign(alternative)
    return base_direction != alt_direction and "zero" not in {base_direction, alt_direction}

def region_delta(sub: pd.DataFrame, groups: list[str]) -> float:
    a = sub[sub["condition"].astype(str) == groups[0]]
    b = sub[sub["condition"].astype(str) == groups[1]]
    if a.empty or b.empty:
        return np.nan
    return float(b["mu_hat"].mean() - a["mu_hat"].mean())

def compute_loo(
    region_id: str,
    sub: pd.DataFrame,
    groups: list[str],
    delta_rep: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    valid: list[tuple[str, float, float, bool]] = []
    for sample_id in sub["sample_id"].astype(str).drop_duplicates():
        removed = sub[sub["sample_id"].astype(str) == sample_id]
        kept = sub[sub["sample_id"].astype(str) != sample_id]
        delta_without = region_delta(kept, groups)
        shift = delta_without - delta_rep if pd.notna(delta_without) and pd.notna(delta_rep) else np.nan
        changed = direction_changed(delta_rep, delta_without) if pd.notna(delta_without) else pd.NA
        removed_condition = ";".join(sorted(removed["condition"].astype(str).dropna().unique()))
        row = {
            "region_id": region_id,
            "removed_sample_id": sample_id,
            "removed_condition": removed_condition,
            "delta_without_sample_i": delta_without,
            "delta_shift_without_sample_i": shift,
            "abs_delta_shift_without_sample_i": abs(shift) if pd.notna(shift) else np.nan,
            "direction_without_sample_i": sign(delta_without),
            "loo_direction_changed": changed,
        }
        rows.append(row)
        if pd.notna(delta_without) and pd.notna(shift):
            valid.append((sample_id, float(delta_without), float(abs(shift)), bool(changed)))
    if not valid:
        return (
            {
                "loo_min_delta": np.nan,
                "loo_max_delta": np.nan,
                "loo_direction_changed": pd.NA,
                "most_influential_sample": pd.NA,
                "most_influential_sample_abs_shift": np.nan,
            },
            rows,
        )
    most = max(valid, key=lambda item: item[2])
    deltas = [item[1] for item in valid]
    return (
        {
            "loo_min_delta": min(deltas),
            "loo_max_delta": max(deltas),
            "loo_direction_changed": any(item[3] for item in valid),
            "most_influential_sample": most[0],
            "most_influential_sample_abs_shift": most[2],
        },
        rows,
    )

def bootstrap_delta_stats(
    a_values: np.ndarray,
    b_values: np.ndarray,
    *,
    n_bootstrap: int,
    alpha: float,
    rng: np.random.Generator,
) -> dict[str, object]:
    if len(a_values) < 2 or len(b_values) < 2:
        return {
            "delta_bootstrap_mean": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "ci_width": np.nan,
            "ci_includes_zero": pd.NA,
            "bootstrap_status": "WARN_TOO_FEW_REPLICATES",
        }
    a_idx = rng.integers(0, len(a_values), size=(n_bootstrap, len(a_values)))
    b_idx = rng.integers(0, len(b_values), size=(n_bootstrap, len(b_values)))
    deltas = b_values[b_idx].mean(axis=1) - a_values[a_idx].mean(axis=1)
    ci_low = float(np.quantile(deltas, alpha / 2))
    ci_high = float(np.quantile(deltas, 1 - alpha / 2))
    return {
        "delta_bootstrap_mean": float(deltas.mean()),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "ci_width": ci_high - ci_low,
        "ci_includes_zero": bool(ci_low <= 0 <= ci_high),
        "bootstrap_status": "PASS",
    }

def dispersion_score(group: pd.DataFrame) -> float:
    values = pd.to_numeric(group["mu_hat"], errors="coerce").dropna()
    totals = pd.to_numeric(group.loc[values.index, "N"], errors="coerce")
    if len(values) < 2 or (totals > 0).sum() < 2:
        return np.nan
    pooled = safe_ratio(group["M"].sum(), group["N"].sum())
    if pd.isna(pooled):
        return np.nan
    observed_var = float(values.var(ddof=1))
    binomial_var = float((pooled * (1 - pooled) / totals[totals > 0]).mean())
    if binomial_var <= 0:
        return np.nan
    return max(0.0, observed_var / binomial_var - 1.0)

def load_cpg_coverage_metrics(args: argparse.Namespace, groups: list[str], count_source: str) -> pd.DataFrame:
    if args.region_cpg_counts:
        cpg_path = args.region_cpg_counts
    elif args.region_sample_counts:
        return pd.DataFrame()
    else:
        source_path = Path(count_source)
        cpg_path = source_path if source_path.name == "region_cpg_counts.tsv" else find_file(["region_cpg_counts.tsv"])
    if not cpg_path or not cpg_path.exists():
        return pd.DataFrame()
    cpg = normalize_cpg_counts_table(read_table(cpg_path))
    design_path = args.design or find_file(["design.tsv"])
    if design_path and design_path.exists():
        design = read_table(design_path)
        sample_col = find_col(design, ["sample_id", "sample"])
        cond_col = find_col(design, [args.condition_column, "condition", "group", "treatment"])
        if sample_col and cond_col:
            design = design.rename(columns={sample_col: "sample_id", cond_col: "condition"})[["sample_id", "condition"]]
            cpg = cpg.drop(columns=["condition"]).merge(design, on="sample_id", how="left")
    cpg = cpg.dropna(subset=["condition"]).copy()
    cpg = cpg[cpg["condition"].astype(str).isin(groups)].copy()
    cpg = cpg[pd.to_numeric(cpg["total"], errors="coerce").fillna(0) >= int(args.min_coverage)].copy()
    if cpg.empty:
        return pd.DataFrame()
    rows: list[dict[str, object]] = []
    for region_id, sub in cpg.groupby("region_id"):
        sets = {
            group: set(sub[sub["condition"].astype(str) == group]["cpg_id"].astype(str))
            for group in groups[:2]
        }
        common = sets[groups[0]] & sets[groups[1]]
        union = sets[groups[0]] | sets[groups[1]]
        rows.append(
            {
                "region_id": str(region_id),
                "n_cpg_A": len(sets[groups[0]]),
                "n_cpg_B": len(sets[groups[1]]),
                "n_cpg_common": len(common),
                "fraction_common_cpg": safe_ratio(len(common), len(union)),
                "coverage_metric_source": str(cpg_path),
            }
        )
    return pd.DataFrame(rows)

def status_failed(value: object) -> bool:
    text = str(value).strip().lower()
    return any(token in text for token in ["fail", "error", "unavailable"])

def status_insufficient(value: object) -> bool:
    text = str(value).strip().lower()
    return any(token in text for token in ["insufficient", "skipped", "not_tested", "low_coverage"])

def scalar_bool(value: object) -> bool | None:
    if pd.isna(value):
        return None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"true", "1", "yes", "pass"}:
            return True
        if text in {"false", "0", "no"}:
            return False
    return bool(value)

def load_glm_glmm_agreement(args: argparse.Namespace) -> tuple[pd.DataFrame, str]:
    path = args.glm_glmm_comparison or find_file(["glm_vs_glmm_comparison.tsv"])
    if not path or not path.exists():
        return pd.DataFrame(), "not_found"
    df = read_table(path)
    region_col = find_col(df, ["region_id", "dmr_id", "harmonized_region_id"])
    if not region_col:
        return pd.DataFrame(), f"region_id missing in {path}"
    glm_q_col = find_col(df, ["glm_q_value", "glm_q", "aggregated_glm_q", "region_q_value", "q_value"])
    glmm_q_col = find_col(df, ["glmm_q_value", "glmm_q"])
    glm_status_col = find_col(df, ["glm_status", "aggregated_glm_status"])
    glmm_status_col = find_col(df, ["glmm_status", "model_status"])
    glm_delta_col = find_col(df, ["glm_delta", "region_delta", "region_delta_callus_minus_seedling"])
    glmm_delta_col = find_col(df, ["glmm_delta", "delta_methylation"])
    out = pd.DataFrame({"region_id": df[region_col].astype(str)})
    out["glm_q_value"] = pd.to_numeric(df[glm_q_col], errors="coerce") if glm_q_col else pd.NA
    out["glmm_q_value"] = pd.to_numeric(df[glmm_q_col], errors="coerce") if glmm_q_col else pd.NA
    out["glm_status"] = df[glm_status_col].astype(str) if glm_status_col else "unknown"
    out["glmm_status"] = df[glmm_status_col].astype(str) if glmm_status_col else "unknown"
    glm_delta = pd.to_numeric(df[glm_delta_col], errors="coerce") if glm_delta_col else pd.Series(pd.NA, index=df.index)
    glmm_delta = pd.to_numeric(df[glmm_delta_col], errors="coerce") if glmm_delta_col else pd.Series(pd.NA, index=df.index)

    glmm_path = args.glmm_results or find_file(["cpg_level_glmm_results.tsv", "real_cpg_level_glmm_results.tsv"])
    if glmm_path and glmm_path.exists():
        try:
            glmm_raw = read_table(glmm_path)
            glmm_region_col = find_col(glmm_raw, ["region_id", "dmr_id", "harmonized_region_id"])
            glmm_delta_extra_col = find_col(glmm_raw, ["delta_methylation", "glmm_delta"])
            glmm_q_extra_col = find_col(glmm_raw, ["glmm_q_value", "q_value", "qvalue", "fdr"])
            glmm_status_extra_col = find_col(glmm_raw, ["glmm_status", "model_status", "status"])
            if glmm_region_col:
                glmm_extra = pd.DataFrame({"region_id": glmm_raw[glmm_region_col].astype(str)})
                if glmm_delta_extra_col:
                    glmm_extra["_glmm_delta_extra"] = pd.to_numeric(glmm_raw[glmm_delta_extra_col], errors="coerce")
                if glmm_q_extra_col:
                    glmm_extra["_glmm_q_extra"] = pd.to_numeric(glmm_raw[glmm_q_extra_col], errors="coerce")
                if glmm_status_extra_col:
                    glmm_extra["_glmm_status_extra"] = glmm_raw[glmm_status_extra_col].astype(str)
                glmm_extra = glmm_extra.drop_duplicates("region_id", keep="first")
                out = out.merge(glmm_extra, on="region_id", how="left")
                if "_glmm_delta_extra" in out:
                    glmm_delta = glmm_delta.fillna(out["_glmm_delta_extra"])
                if "_glmm_q_extra" in out:
                    out["glmm_q_value"] = pd.to_numeric(out["glmm_q_value"], errors="coerce").fillna(out["_glmm_q_extra"])
                if "_glmm_status_extra" in out:
                    missing_status = out["glmm_status"].astype(str).isin({"", "nan", "None", "unknown"})
                    out.loc[missing_status, "glmm_status"] = out.loc[missing_status, "_glmm_status_extra"]
                out = out.drop(columns=[c for c in out.columns if c.startswith("_glmm_")])
        except Exception:
            pass

    out["glm_significant"] = pd.to_numeric(out["glm_q_value"], errors="coerce") <= args.q_threshold
    out["glmm_significant"] = (
        (pd.to_numeric(out["glmm_q_value"], errors="coerce") <= args.q_threshold)
        & ~out["glmm_status"].map(status_failed)
        & ~out["glmm_status"].map(status_insufficient)
    )
    if glm_delta.notna().any() and glmm_delta.notna().any():
        out["same_direction"] = [
            pd.NA if sign(a) == "zero" or sign(b) == "zero" else sign(a) == sign(b)
            for a, b in zip(glm_delta, glmm_delta)
        ]
    else:
        out["same_direction"] = pd.NA

    glm_sig = out["glm_significant"].map(lambda value: scalar_bool(value) is True)
    glmm_sig = out["glmm_significant"].map(lambda value: scalar_bool(value) is True)
    same_dir = out["same_direction"].map(scalar_bool)
    glm_failed = out["glm_status"].map(status_failed)
    glmm_failed = out["glmm_status"].map(status_failed)
    glmm_insufficient = out["glmm_status"].map(status_insufficient)
    out["model_agreement_status"] = "NOT_SIGNIFICANT"
    out.loc[glmm_sig & ~glm_sig, "model_agreement_status"] = "GLMM_ONLY"
    out.loc[glm_sig & ~glmm_sig, "model_agreement_status"] = "GLM_ONLY"
    out.loc[glm_sig & glmm_sig, "model_agreement_status"] = "GLMM_CONFIRMED"
    out.loc[(same_dir == False) & (glm_sig | glmm_sig), "model_agreement_status"] = "DIRECTION_DISCORDANT"  # noqa: E712
    out.loc[glmm_insufficient, "model_agreement_status"] = "INSUFFICIENT_DATA"
    out.loc[glm_failed | glmm_failed, "model_agreement_status"] = "MODEL_FAILED"
    out["_sort_q"] = pd.to_numeric(out["glmm_q_value"], errors="coerce").fillna(
        pd.to_numeric(out["glm_q_value"], errors="coerce")
    ).fillna(1.0)
    out = out.sort_values("_sort_q").drop_duplicates("region_id", keep="first").drop(columns=["_sort_q"])
    out["glm_glmm_source"] = str(path)
    return out, str(path)

def bool_true(value: object) -> bool:
    return scalar_bool(value) is True

def greater_than(value: object, threshold: float) -> bool:
    return pd.notna(value) and float(value) > threshold

def classify_confidence(row: dict[str, object], args: argparse.Namespace) -> tuple[str, str]:
    severe: list[str] = []
    moderate: list[str] = []
    if int(row.get("n_samples_A", 0)) < 2 or int(row.get("n_samples_B", 0)) < 2:
        severe.append("too_few_replicates")
    if bool_true(row.get("direction_changed")):
        severe.append("pooled_delta_direction_changed")
    if bool_true(row.get("loo_direction_changed")):
        severe.append("leave_one_out_direction_changed")
    if str(row.get("model_agreement_status", "")) == "DIRECTION_DISCORDANT":
        severe.append("glm_glmm_direction_discordant")
    if greater_than(row.get("abs_delta_shift"), args.max_abs_delta_shift):
        moderate.append("large_absolute_delta_shift")
    if greater_than(row.get("relative_delta_shift"), args.max_relative_delta_shift):
        moderate.append("large_relative_delta_shift")
    if bool_true(row.get("ci_includes_zero")):
        moderate.append("bootstrap_ci_includes_zero")
    if str(row.get("bootstrap_status", "")) == "WARN_TOO_FEW_REPLICATES":
        moderate.append("bootstrap_too_few_replicates")
    if greater_than(row.get("coverage_imbalance_score"), args.max_coverage_imbalance):
        moderate.append("coverage_imbalance")
    if str(row.get("model_agreement_status", "")) in {"MODEL_FAILED", "INSUFFICIENT_DATA", "GLM_ONLY"}:
        moderate.append(str(row.get("model_agreement_status")).lower())
    if greater_than(row.get("overdispersion_score"), 10.0):
        moderate.append("high_overdispersion_diagnostic")
    if severe:
        return "LOW_CONFIDENCE", ";".join(severe + moderate)
    if len(moderate) >= 2:
        return "WEAK_CONFIDENCE", ";".join(moderate)
    if moderate:
        return "MODERATE_CONFIDENCE", ";".join(moderate)
    return "HIGH_CONFIDENCE", "all_primary_checks_stable"
