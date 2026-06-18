"""Compare aggregated beta-binomial GLM results with CpG-level GLMM validation."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


GLM_VS_GLMM_COLUMNS = [
    "region_id",
    "context",
    "glm_p_value",
    "glm_q_value",
    "glmm_p_value",
    "glmm_q_value",
    "glm_status",
    "glmm_status",
    "confirmed_by_glmm",
    "glm_only_candidate",
    "glmm_more_conservative",
    "insufficient_glmm_data",
    "interpretation",
]


def _find_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def read_region_level_glm_results(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    region_col = _find_column(df, ("region_id", "dmr_id"))
    if region_col is None:
        raise ValueError("region-level GLM table must contain region_id or dmr_id")
    out = pd.DataFrame({"region_id": df[region_col].astype(str)})
    context_col = _find_column(df, ("context",))
    p_col = _find_column(df, ("glm_p_value", "region_p_value", "p_value"))
    q_col = _find_column(df, ("glm_q_value", "region_q_value", "q_value", "padj", "fdr"))
    status_col = _find_column(df, ("glm_status", "model_status", "status"))
    out["context"] = df[context_col].astype(str) if context_col else "NA"
    out["glm_p_value"] = pd.to_numeric(df[p_col], errors="coerce") if p_col else pd.NA
    out["glm_q_value"] = pd.to_numeric(df[q_col], errors="coerce") if q_col else pd.NA
    out["glm_status"] = df[status_col].astype(str) if status_col else "ok"
    return out


def read_cpg_level_glmm_results(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "region_id" not in df.columns:
        raise ValueError("CpG-level GLMM table must contain region_id")
    out = pd.DataFrame({"region_id": df["region_id"].astype(str)})
    out["context"] = df["context"].astype(str) if "context" in df.columns else "NA"
    out["glmm_p_value"] = pd.to_numeric(df["p_value"], errors="coerce") if "p_value" in df.columns else pd.NA
    out["glmm_q_value"] = pd.to_numeric(df["q_value"], errors="coerce") if "q_value" in df.columns else pd.NA
    out["glmm_status"] = df["model_status"].astype(str) if "model_status" in df.columns else "unknown"
    return out


def compare_beta_binomial_glm_vs_glmm(
    glm_df: pd.DataFrame,
    glmm_df: pd.DataFrame,
    *,
    q_threshold: float = 0.05,
) -> pd.DataFrame:
    merged = glm_df.merge(
        glmm_df,
        on=["region_id", "context"],
        how="outer",
        suffixes=("_glm", "_glmm"),
    )
    if "glm_status" not in merged.columns:
        merged["glm_status"] = "unknown"
    if "glmm_status" not in merged.columns:
        merged["glmm_status"] = "missing"
    glm_sig = pd.to_numeric(merged["glm_q_value"], errors="coerce") <= q_threshold
    glmm_sig = (
        (pd.to_numeric(merged["glmm_q_value"], errors="coerce") <= q_threshold)
        & (merged["glmm_status"].astype(str) == "ok")
    )
    skipped_data = merged["glmm_status"].astype(str).isin({
        "insufficient_cpg",
        "insufficient_replicates",
        "too_many_zero_coverage",
        "missing_condition",
        "low_coverage",
    })
    unavailable_or_error = merged["glmm_status"].astype(str).isin({
        "glmmTMB_unavailable",
        "glmm_temp_dir_missing",
        "model_error",
        "model_no_lrt",
        "model_not_converged",
        "model_singular",
    })
    merged["confirmed_by_glmm"] = glm_sig & glmm_sig
    merged["glm_only_candidate"] = glm_sig & ~glmm_sig
    merged["glmm_more_conservative"] = glm_sig & (~glmm_sig | skipped_data | unavailable_or_error)
    merged["insufficient_glmm_data"] = skipped_data

    interpretations = []
    for idx, row in merged.iterrows():
        if row["confirmed_by_glmm"]:
            interpretations.append("Aggregated GLM signal is confirmed by confirmatory CpG-level GLMM.")
        elif bool(skipped_data.loc[idx]):
            interpretations.append("CpG-level GLMM has insufficient data for confirmatory interpretation.")
        elif bool(unavailable_or_error.loc[idx]):
            interpretations.append("CpG-level GLMM was not available or did not fit; keep GLM result as unconfirmed.")
        elif row["glm_only_candidate"]:
            interpretations.append("Aggregated GLM candidate was not confirmed by CpG-level GLMM.")
        else:
            interpretations.append("No GLM/GLMM confirmation pattern detected.")
    merged["interpretation"] = interpretations

    for column in GLM_VS_GLMM_COLUMNS:
        if column not in merged.columns:
            merged[column] = pd.NA
    return merged[GLM_VS_GLMM_COLUMNS]
