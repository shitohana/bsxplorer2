"""Diagnostic plots for confirmatory CpG-level GLMM validation."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


STATUS_COLORS = {
    "confirmed_by_glmm": "#2CA02C",
    "glm_only_candidate": "#F28E2B",
    "insufficient_glmm_data": "#D62728",
    "other_or_missing": "#7F7FBA",
}


def _plt():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def _truthy(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def _neglog10_q(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=numeric.index, dtype=float)
    positive = numeric > 0
    out.loc[positive] = -np.log10(numeric.loc[positive])
    zero_mask = numeric == 0
    if zero_mask.any():
        finite = out.replace([np.inf, -np.inf], np.nan).dropna()
        zero_position = 50.0 if finite.empty else max(float(finite.max()) * 1.05, float(finite.max()) + 1.0)
        out.loc[zero_mask] = zero_position
    return out.replace([np.inf, -np.inf], np.nan)


def plot_glm_vs_glmm_scatter(
    comparison_df: pd.DataFrame,
    out_path: str | Path,
    *,
    q_threshold: float = 0.10,
) -> None:
    plt = _plt()
    df = comparison_df.copy()
    df["glm_q_value"] = pd.to_numeric(df["glm_q_value"], errors="coerce")
    df["glmm_q_value"] = pd.to_numeric(df["glmm_q_value"], errors="coerce")
    df["_x"] = _neglog10_q(df["glm_q_value"])
    df["_y"] = _neglog10_q(df["glmm_q_value"])
    df = df.dropna(subset=["_x", "_y"]).copy()

    if "confirmed_by_glmm" in df.columns:
        confirmed = _truthy(df["confirmed_by_glmm"])
    else:
        confirmed = pd.Series(False, index=df.index)
    if "glm_only_candidate" in df.columns:
        glm_only = _truthy(df["glm_only_candidate"])
    else:
        glm_only = pd.Series(False, index=df.index)
    if "insufficient_glmm_data" in df.columns:
        insufficient = _truthy(df["insufficient_glmm_data"])
    else:
        insufficient = pd.Series(False, index=df.index)

    df["_status"] = "other_or_missing"
    df.loc[insufficient, "_status"] = "insufficient_glmm_data"
    df.loc[glm_only, "_status"] = "glm_only_candidate"
    df.loc[confirmed, "_status"] = "confirmed_by_glmm"

    fig, ax = plt.subplots(figsize=(7, 5.2))
    for status, group in df.groupby("_status", sort=False):
        ax.scatter(
            group["_x"],
            group["_y"],
            s=42,
            alpha=0.82,
            color=STATUS_COLORS.get(status, "#666666"),
            label=status.replace("_", " "),
            edgecolor="white",
            linewidth=0.35,
        )

    if q_threshold > 0:
        threshold_line = float(-np.log10(max(q_threshold, 1e-300)))
        ax.axvline(threshold_line, color="#4D4D4D", linestyle="--", linewidth=0.9)
        ax.axhline(threshold_line, color="#4D4D4D", linestyle="--", linewidth=0.9)

    if df.empty:
        x_upper = 1.0
        y_upper = 1.0
    else:
        x_upper = max(float(df["_x"].max()) * 1.08, threshold_line * 1.25 if q_threshold > 0 else 1.0)
        y_upper = max(float(df["_y"].max()) * 1.08, threshold_line * 1.25 if q_threshold > 0 else 1.0)
    ax.set_xlim(0, x_upper)
    ax.set_ylim(0, y_upper)
    ax.set_xlabel("-log10 aggregated GLM q-value")
    ax.set_ylabel("-log10 CpG-level GLMM q-value")
    ax.set_title("Aggregated GLM vs CpG-level GLMM confirmatory validation")
    if not df.empty:
        ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_glmm_status_summary(glmm_df: pd.DataFrame, out_path: str | Path) -> None:
    plt = _plt()
    counts = glmm_df["model_status"].fillna("NA").value_counts()
    fig, ax = plt.subplots(figsize=(7, 4))
    counts.plot(kind="bar", ax=ax)
    ax.set_ylabel("Regions")
    ax.set_title("CpG-level GLMM confirmatory validation model status")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_glmm_confirmation_by_context(comparison_df: pd.DataFrame, out_path: str | Path) -> None:
    plt = _plt()
    df = comparison_df.copy()
    df["category"] = df["confirmed_by_glmm"].map({True: "confirmed", False: "not_confirmed"})
    table = df.groupby(["context", "category"]).size().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(7, 4))
    table.plot(kind="bar", stacked=True, ax=ax)
    ax.set_ylabel("Regions")
    ax.set_title("Confirmed vs GLM-only by context")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_top_cpg_contribution(cpg_counts_df: pd.DataFrame, out_path: str | Path) -> None:
    plt = _plt()
    df = cpg_counts_df.copy()
    df["total"] = pd.to_numeric(df["total"], errors="coerce").fillna(0)
    by_cpg = df.groupby(["region_id", "cpg_id"], as_index=False)["total"].sum()
    by_region = by_cpg.groupby("region_id")["total"].transform("sum")
    by_cpg["coverage_fraction"] = by_cpg["total"] / by_region.replace(0, pd.NA)
    top = by_cpg.groupby("region_id")["coverage_fraction"].max().dropna()
    fig, ax = plt.subplots(figsize=(6, 4))
    top.hist(bins=30, ax=ax)
    ax.set_xlabel("Top CpG coverage fraction per region")
    ax.set_ylabel("Regions")
    ax.set_title("Top CpG contribution in CpG-level GLMM confirmatory validation input")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_example_cpg_consistency(
    cpg_counts_df: pd.DataFrame,
    *,
    region_id: str,
    out_path: str | Path,
    design_df: Optional[pd.DataFrame] = None,
    condition_column: str = "condition",
    case_label: Optional[str] = None,
    control_label: Optional[str] = None,
    region_delta: Optional[float] = None,
) -> None:
    plt = _plt()
    df = cpg_counts_df[cpg_counts_df["region_id"].astype(str) == str(region_id)].copy()
    if "position" not in df.columns and "pos" in df.columns:
        df["position"] = df["pos"]
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["methylation"] = pd.to_numeric(df["mC"], errors="coerce") / pd.to_numeric(df["total"], errors="coerce").replace(0, pd.NA)
    if design_df is not None and condition_column in design_df.columns:
        df = df.merge(design_df[["sample_id", condition_column]], on="sample_id", how="left")
        groups = sorted(df[condition_column].dropna().astype(str).unique())
        if case_label is None or control_label is None:
            if len(groups) >= 2:
                control_label, case_label = groups[:2]
        if case_label and control_label:
            summary = df.groupby(["position", condition_column])["methylation"].mean().unstack()
            y = summary.get(str(case_label), pd.Series(dtype=float)) - summary.get(str(control_label), pd.Series(dtype=float))
        else:
            y = df.groupby("position")["methylation"].mean()
    else:
        y = df.groupby("position")["methylation"].mean()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(y.index, y.values, alpha=0.8)
    if region_delta is not None:
        ax.axhline(float(region_delta), color="red", linestyle="--", label="region-level delta")
        ax.legend()
    ax.set_xlabel("CpG genomic position")
    ax.set_ylabel("Methylation difference" if design_df is not None else "Mean methylation")
    ax.set_title(f"CpG-level GLMM confirmatory validation consistency: {region_id}")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
