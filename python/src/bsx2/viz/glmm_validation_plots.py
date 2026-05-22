"""Diagnostic plots for confirmatory CpG-level GLMM validation."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


def _plt():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def plot_glm_vs_glmm_scatter(comparison_df: pd.DataFrame, out_path: str | Path) -> None:
    plt = _plt()
    df = comparison_df.copy()
    df["glm_q_value"] = pd.to_numeric(df["glm_q_value"], errors="coerce")
    df["glmm_q_value"] = pd.to_numeric(df["glmm_q_value"], errors="coerce")
    x = -np.log10(df["glm_q_value"].clip(lower=1e-300).astype(float))
    y = -np.log10(df["glmm_q_value"].clip(lower=1e-300).astype(float))
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x, y, alpha=0.7)
    ax.set_xlabel("-log10 aggregated GLM q-value")
    ax.set_ylabel("-log10 confirmatory CpG-level GLMM q-value")
    ax.set_title("Aggregated GLM vs confirmatory CpG-level GLMM")
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
    ax.set_title("confirmatory CpG-level GLMM model status")
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
    ax.set_title("Top CpG contribution in confirmatory CpG-level GLMM input")
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
    ax.set_title(f"confirmatory CpG-level GLMM consistency: {region_id}")
    fig.tight_layout()
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
