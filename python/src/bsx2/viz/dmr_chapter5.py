"""Chapter 5 figure helpers for DMR/metagene validation outputs.

The functions in this module are intentionally thin rendering wrappers around
prepared BSX2 downstream tables. They do not run DMR calling, do not recompute
statistical results, and do not read raw sequencing data. The goal is to keep
chapter/report figures in the same BSX2 visualization namespace while reusing
existing validation artifacts.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


SEEDLING_COLOR = "#1f77b4"
CALLUS_COLOR = "#ff7f0e"
CONDITION_COLORS = {"Seedling": SEEDLING_COLOR, "Callus": CALLUS_COLOR}
CONTEXT_COLORS = {"CG": "#4C78A8", "CHG": "#54A24B", "CHH": "#E45756"}
STATUS_COLORS = {
    "PASS": "#2CA02C",
    "WARN": "#F28E2B",
    "SKIPPED": "#9E9E9E",
    "FAIL": "#D62728",
    "confirmed_by_glmm": "#2CA02C",
    "glm_only_candidate": "#F28E2B",
    "glmm_more_conservative": "#F28E2B",
    "neither": "#9E9E9E",
    "insufficient_glmm_data": "#D62728",
    "direction_stable": "#2CA02C",
    "direction_not_stable": "#D62728",
}


@dataclass
class SavedFigure:
    png: str
    pdf: str
    svg: str


def setup_style() -> None:
    """Apply a restrained Matplotlib style compatible with BSX2 reports."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 300,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.18,
            "grid.linewidth": 0.7,
        }
    )


def save_figure(fig, stem: str, out_dirs: Mapping[str, Path], formats: Sequence[str]) -> SavedFigure:
    """Save a Matplotlib figure to the requested format directories."""

    setup_style()
    outputs: dict[str, str] = {"png": "", "pdf": "", "svg": ""}
    for fmt in formats:
        if fmt not in outputs:
            continue
        out_dir = out_dirs[fmt]
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"{stem}.{fmt}"
        fig.savefig(path, bbox_inches="tight")
        outputs[fmt] = str(path)
    try:
        import matplotlib.pyplot as plt

        plt.close(fig)
    except Exception:
        pass
    return SavedFigure(outputs["png"], outputs["pdf"], outputs["svg"])


def read_tsv(path: str | Path, *, nrows: int | None = None) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", nrows=nrows)


def numeric_series(values) -> pd.Series:
    return pd.to_numeric(values, errors="coerce")


def cap_neglog10_q(q_values: pd.Series, cap: float = 50.0) -> pd.Series:
    q = pd.to_numeric(q_values, errors="coerce")
    q = q.where(q > 0, np.nan)
    y = -np.log10(q)
    return y.replace([np.inf, -np.inf], np.nan).clip(upper=cap)


def add_metagene_boundaries(ax, n_bins: int, upstream_bins: int = 20, body_bins: int = 100) -> None:
    for x in (upstream_bins - 0.5, upstream_bins + body_bins - 0.5):
        ax.axvline(x, color="#4D4D4D", linestyle="--", linewidth=0.8)
    ax.set_xticks([upstream_bins / 2, upstream_bins + body_bins / 2, upstream_bins + body_bins + (n_bins - upstream_bins - body_bins) / 2])
    ax.set_xticklabels(["upstream 2 kb", "gene body", "downstream 2 kb"])


def compressed_heatmap_matrix(df: pd.DataFrame, *, max_rows: int = 900) -> pd.DataFrame:
    """Return a row-ranked/compressed matrix from a gene/bin table."""

    bin_cols = [col for col in df.columns if str(col).startswith("bin_")]
    mat = df[bin_cols].apply(pd.to_numeric, errors="coerce")
    if mat.empty:
        return mat
    body_cols = [col for col in bin_cols if 20 <= int(str(col).split("_")[-1]) < 120]
    score_cols = body_cols or bin_cols
    order = mat[score_cols].mean(axis=1, skipna=True).sort_values(ascending=False).index
    mat = mat.loc[order]
    if len(mat) <= max_rows:
        return mat
    chunks = np.array_split(np.arange(len(mat)), max_rows)
    rows = []
    for chunk in chunks:
        rows.append(mat.iloc[chunk].mean(axis=0, skipna=True))
    return pd.DataFrame(rows, columns=bin_cols)


def render_image_file(src: str | Path, *, title: str | None = None, figsize: tuple[float, float] = (8, 5)):
    """Wrap an existing PNG/JPEG into a Matplotlib figure for uniform export."""

    setup_style()
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(mpimg.imread(src))
    ax.set_axis_off()
    if title:
        ax.set_title(title)
    return fig


def plot_image_grid(
    image_paths: Sequence[str | Path],
    labels: Sequence[str],
    *,
    title: str,
    ncols: int = 3,
    figsize: tuple[float, float] = (13, 4.8),
):
    setup_style()
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    n = len(image_paths)
    ncols = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for idx, path in enumerate(image_paths):
        ax = axes.ravel()[idx]
        ax.imshow(mpimg.imread(path))
        ax.set_title(labels[idx])
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    return fig


def plot_runtime_by_stage_group(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    if "stage_group" not in df.columns:
        df = df.copy()
        stage = df.get("stage_name", pd.Series(dtype=str)).astype(str)
        df["stage_group"] = np.select(
            [
                stage.str.contains("visual", case=False, na=False),
                stage.str.contains("regional|evidence", case=False, na=False),
                stage.str.contains("dss|methyl|external", case=False, na=False),
                stage.str.contains("glmm|glm", case=False, na=False),
            ],
            ["visualization", "regional_evidence", "external_callers", "GLM/GLMM"],
            default="final_assembly",
        )
        df = df.groupby("stage_group", as_index=False)["wall_seconds"].sum()
    else:
        df = df.copy()

    df["wall_minutes"] = pd.to_numeric(df.get("wall_minutes", df["wall_seconds"] / 60.0), errors="coerce")
    order = df.sort_values("wall_minutes", ascending=False)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(order["stage_group"].astype(str), order["wall_minutes"], color="#4C78A8")
    ax.set_ylabel("Wall time, min")
    ax.set_xlabel("Workflow stage group")
    ax.set_title("Workflow runtime by stage group")
    ax.tick_params(axis="x", rotation=30)
    for i, v in enumerate(order["wall_minutes"]):
        if np.isfinite(v):
            ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    return fig


def plot_metagene_profiles(profile_by_context: Mapping[str, pd.DataFrame]):
    setup_style()
    import matplotlib.pyplot as plt

    contexts = ["CG", "CHG", "CHH"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), sharey=False)
    for ax, context in zip(axes, contexts, strict=True):
        df = profile_by_context.get(context, pd.DataFrame()).copy()
        if df.empty:
            ax.text(0.5, 0.5, f"{context}: missing input", ha="center", va="center")
            ax.set_axis_off()
            continue
        value_col = "weighted_methylation" if "weighted_methylation" in df else "gene_profile_mean"
        x_col = "bin" if "bin" in df else "relative_position"
        if "condition" in df:
            for condition, grp in df.groupby("condition", sort=False):
                color = CONDITION_COLORS.get(str(condition), "#666666")
                ax.plot(grp[x_col], grp[value_col], label=str(condition), color=color, linewidth=1.8)
        else:
            ax.plot(df[x_col], df[value_col], color=CONTEXT_COLORS.get(context, "#4C78A8"), linewidth=1.8)
            if {"band_lower", "band_upper"}.issubset(df.columns):
                ax.fill_between(df[x_col], df["band_lower"], df["band_upper"], color=CONTEXT_COLORS.get(context, "#4C78A8"), alpha=0.18)
        n_bins = int(pd.to_numeric(df[x_col], errors="coerce").max()) + 1 if x_col == "bin" else len(df)
        add_metagene_boundaries(ax, n_bins)
        ax.set_title(context)
        ax.set_ylabel("Weighted methylation")
        ax.set_ylim(0, min(1.0, max(0.15, np.nanmax(pd.to_numeric(df[value_col], errors="coerce")) * 1.15)))
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle("Metagene methylation profiles by context", y=1.03)
    fig.tight_layout()
    return fig


def plot_metagene_heatmaps(matrix_by_context: Mapping[str, pd.DataFrame]):
    setup_style()
    import matplotlib.pyplot as plt

    contexts = ["CG", "CHG", "CHH"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=False)
    for ax, context in zip(axes, contexts, strict=True):
        df = matrix_by_context.get(context, pd.DataFrame())
        if df.empty:
            ax.text(0.5, 0.5, f"{context}: missing input", ha="center", va="center")
            ax.set_axis_off()
            continue
        mat = compressed_heatmap_matrix(df)
        image = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", interpolation="nearest", cmap="viridis", vmin=0, vmax=1)
        add_metagene_boundaries(ax, len(mat.columns))
        ax.set_title(context)
        ax.set_ylabel("Genes ranked by gene-body methylation")
        ax.set_yticks([])
    cbar = fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02)
    cbar.set_label("Weighted methylation")
    fig.suptitle("Metagene methylation heatmaps by context", y=1.02)
    fig.tight_layout()
    return fig


def plot_volcano(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 5.1))
    work = df.copy()
    delta_col = "delta_methylation" if "delta_methylation" in work else "region_delta"
    q_col = "q_value" if "q_value" in work else "region_q_value"
    class_col = "evidence_class" if "evidence_class" in work else None
    work["neglog10q"] = work["neg_log10_q_capped"] if "neg_log10_q_capped" in work else cap_neglog10_q(work[q_col])
    for label, grp in work.groupby(class_col, dropna=False) if class_col else [("candidate", work)]:
        ax.scatter(grp[delta_col], grp["neglog10q"], s=5, alpha=0.45, label=str(label))
    ax.axhline(-np.log10(0.05), color="#444444", linestyle="--", linewidth=0.9)
    ax.axvline(0, color="#444444", linewidth=0.8)
    ax.set_xlabel("Delta methylation (Callus - Seedling)")
    ax.set_ylabel("-log10(q), capped at 50")
    ax.set_title("Regional Evidence volcano with MSU7 annotation")
    if class_col:
        ax.legend(markerscale=2, frameon=False, ncols=2)
    fig.tight_layout()
    return fig


def plot_bar(df: pd.DataFrame, *, x: str, y: str, title: str, ylabel: str, color: str = "#4C78A8", rotate: int = 30):
    setup_style()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    work = df.copy()
    ax.bar(work[x].astype(str), pd.to_numeric(work[y], errors="coerce"), color=color)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=rotate)
    for i, v in enumerate(pd.to_numeric(work[y], errors="coerce")):
        if np.isfinite(v):
            ax.text(i, v, f"{int(v):,}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    return fig


def plot_stacked_percent(df: pd.DataFrame, *, x: str, stack: str, value: str, title: str, ylabel: str = "Percent"):
    setup_style()
    import matplotlib.pyplot as plt

    pivot = df.pivot_table(index=x, columns=stack, values=value, aggfunc="sum", fill_value=0)
    perc = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0) * 100.0
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    bottom = np.zeros(len(perc))
    cmap = plt.get_cmap("tab20")
    for idx, col in enumerate(perc.columns):
        vals = perc[col].to_numpy(dtype=float)
        ax.bar(perc.index.astype(str), vals, bottom=bottom, label=str(col), color=cmap(idx % 20))
        bottom += np.nan_to_num(vals)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=25)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return fig


def plot_chrom_distribution(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    work = df.copy()
    y_col = "regions_per_mb" if "regions_per_mb" in work else "dmr_per_mb" if "dmr_per_mb" in work else None
    if y_col is None:
        numeric = [col for col in work.columns if col not in {"chrom", "context"} and pd.api.types.is_numeric_dtype(work[col])]
        y_col = numeric[0]
    pivot = work.pivot_table(index="chrom", columns="context", values=y_col, aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(sorted(pivot.index, key=lambda x: int(str(x).replace("Chr", "")) if str(x).replace("Chr", "").isdigit() else 999))
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    x = np.arange(len(pivot))
    width = 0.25
    for idx, context in enumerate(["CG", "CHG", "CHH"]):
        if context in pivot:
            ax.bar(x + (idx - 1) * width, pivot[context], width, label=context, color=CONTEXT_COLORS.get(context))
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.astype(str), rotation=45)
    ax.set_ylabel("Significant regions per Mb")
    ax.set_title("Normalized DMR chromosome distribution")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_overlap_threshold_sensitivity(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    work = df.copy()
    if "threshold" not in work.columns:
        return None
    work["threshold"] = pd.to_numeric(work["threshold"], errors="coerce")
    metric_col = "jaccard_like" if "jaccard_like" in work.columns else "n_pairs"
    if metric_col not in work.columns:
        return None
    work[metric_col] = pd.to_numeric(work[metric_col], errors="coerce")
    work = work.dropna(subset=["threshold", metric_col])
    if work.empty:
        return None

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    if "matching_policy" in work.columns:
        preferred = ["best_reciprocal", "one_to_one_greedy", "many_to_many"]
        seen = [str(value) for value in work["matching_policy"].dropna().unique()]
        policies = [policy for policy in preferred if policy in seen]
        policies.extend(sorted(policy for policy in seen if policy not in set(policies)))
        grouped = [(policy, work[work["matching_policy"].astype(str) == policy]) for policy in policies]
    else:
        policies = ["all"]
        grouped = [("all", work)]

    traces: dict[tuple[tuple[float, float, float, float], ...], dict[str, object]] = {}
    for label, grp in grouped:
        agg = (
            grp.groupby("threshold", as_index=False)[metric_col]
            .agg(median="median", min="min", max="max", count="count")
            .sort_values("threshold")
        )
        if agg.empty:
            continue
        signature = tuple(
            (
                round(float(row["threshold"]), 12),
                round(float(row["median"]), 12),
                round(float(row["min"]), 12),
                round(float(row["max"]), 12),
            )
            for _, row in agg.iterrows()
        )
        if signature not in traces:
            traces[signature] = {"labels": [], "agg": agg}
        traces[signature]["labels"].append(label)

    colors = ["#4C78A8", "#F28E2B", "#54A24B", "#E45756", "#7F7F7F"]
    for idx, trace in enumerate(traces.values()):
        agg = trace["agg"]
        labels = trace["labels"]
        if len(labels) > 1 and set(labels) == set(policies):
            label = "all matching policies"
        else:
            label = " / ".join(str(item) for item in labels)
        color = colors[idx % len(colors)]
        x = agg["threshold"].to_numpy(dtype=float)
        y = agg["median"].to_numpy(dtype=float)
        lo = agg["min"].to_numpy(dtype=float)
        hi = agg["max"].to_numpy(dtype=float)
        ax.plot(x, y, marker="o", linewidth=2.0, color=color, label=label)
        if bool(np.nanmax(hi - lo) > 0):
            ax.fill_between(x, lo, hi, color=color, alpha=0.16, linewidth=0)

    if work["threshold"].min() <= 0.5 <= work["threshold"].max():
        ax.axvline(0.5, color="#4D4D4D", linestyle="--", linewidth=0.9, alpha=0.65)
    metric_label = "Jaccard-like agreement" if metric_col == "jaccard_like" else "matched DMR pairs"
    medians = work.groupby("threshold")[metric_col].median()
    note = "median; ribbon = min-max across caller pairs"
    if medians.nunique(dropna=True) == 1:
        note += "; no change across tested thresholds"
    ax.text(0.01, 0.98, note, transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#4D4D4D")
    ax.set_xlabel("Strict reciprocal overlap threshold (tau)")
    ax.set_ylabel(metric_label)
    ax.set_title("Consensus sensitivity to overlap threshold")
    if metric_col == "jaccard_like":
        ymax = float(work[metric_col].max())
        ax.set_ylim(0, min(1.0, max(0.05, ymax + 0.08)))
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    return fig


def plot_coverage_qc(region_qc: pd.DataFrame, sample_qc: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4))
    axes[0].hist(pd.to_numeric(region_qc["fraction_common"], errors="coerce").dropna(), bins=12, color="#4C78A8")
    axes[0].set_title("Fraction common CpG")
    axes[0].set_xlabel("fraction_common")
    axes[0].set_ylabel("Regions")

    status_counts = sample_qc["region_sample_coverage_status"].fillna("missing").value_counts()
    axes[1].bar(status_counts.index.astype(str), status_counts.values, color="#72B7B2")
    axes[1].set_title("Region x sample coverage status")
    axes[1].tick_params(axis="x", rotation=30)

    axes[2].hist(pd.to_numeric(region_qc["n_cpg_common"], errors="coerce").dropna(), bins=12, color="#F28E2B")
    axes[2].set_title("Common CpG per region")
    axes[2].set_xlabel("n_cpg_common")
    fig.suptitle("Coverage-set common-CpG QC", y=1.02)
    fig.tight_layout()
    return fig


def plot_region_cpg_coverage(sample_qc: pd.DataFrame, top_cpg: pd.DataFrame | None = None):
    setup_style()
    import matplotlib.pyplot as plt

    n_panels = 3 if top_cpg is not None and not top_cpg.empty else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(5.0 * n_panels, 4.0))
    axes = np.atleast_1d(axes)
    axes[0].hist(pd.to_numeric(sample_qc["n_cpg_common"], errors="coerce").dropna(), bins=12, color="#4C78A8")
    axes[0].set_title("n common CpG")
    axes[1].hist(pd.to_numeric(sample_qc["total_coverage_common"], errors="coerce").dropna(), bins=12, color="#F28E2B")
    axes[1].set_title("total common coverage")
    if n_panels == 3:
        col = "top_cpg_coverage_fraction"
        axes[2].hist(pd.to_numeric(top_cpg[col], errors="coerce").dropna(), bins=10, color="#59A14F")
        axes[2].set_title("Top CpG contribution")
        axes[2].set_xlabel("fraction")
    fig.suptitle("Region-level CpG count and coverage distribution", y=1.02)
    fig.tight_layout()
    return fig


def plot_top_cpg_contribution(top_cpg: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    col = "top_cpg_coverage_fraction"
    ax.hist(pd.to_numeric(top_cpg[col], errors="coerce").dropna(), bins=12, color="#59A14F")
    ax.set_xlabel("Top CpG coverage fraction")
    ax.set_ylabel("Regions")
    ax.set_title("Top CpG contribution fraction distribution")
    fig.tight_layout()
    return fig


def plot_feature_classes_by_context_panels(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    sets = [item for item in ["significant", "high_confidence"] if item in set(df["candidate_set"].astype(str))]
    if not sets:
        sets = list(df["candidate_set"].astype(str).dropna().unique()[:2])
    fig, axes = plt.subplots(1, len(sets), figsize=(6.8 * len(sets), 4.3), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, candidate_set in zip(axes, sets, strict=True):
        sub = df[df["candidate_set"].astype(str) == candidate_set]
        pivot = sub.pivot_table(index="context", columns="feature_class", values="n_regions", aggfunc="sum", fill_value=0)
        pivot = pivot.reindex(["CG", "CHG", "CHH"]).fillna(0)
        perc = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0) * 100.0
        bottom = np.zeros(len(perc))
        cmap = plt.get_cmap("tab20")
        for idx, col in enumerate(perc.columns):
            vals = perc[col].to_numpy(dtype=float)
            ax.bar(perc.index.astype(str), vals, bottom=bottom, label=str(col), color=cmap(idx % 20))
            bottom += np.nan_to_num(vals)
        ax.set_title(candidate_set.replace("_", " "))
        ax.set_ylabel("Percent within context")
        ax.set_ylim(0, 100)
    axes[-1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.suptitle("DMR candidate genomic feature classes by context", y=1.02)
    fig.tight_layout()
    return fig


def plot_delta_rep_vs_pooled(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    x = pd.to_numeric(df["delta_rep"], errors="coerce")
    y = pd.to_numeric(df["delta_pool"], errors="coerce")
    colors = np.where(df.get("direction_changed", False).astype(bool), "#D62728", "#4C78A8")
    axes[0].scatter(x, y, c=colors, s=38, alpha=0.85)
    lim = np.nanmax(np.abs(np.r_[x, y]))
    axes[0].plot([-lim, lim], [-lim, lim], color="#444444", linestyle="--", linewidth=0.9)
    axes[0].set_xlabel("Delta replicate-level")
    axes[0].set_ylabel("Delta pooled")
    axes[0].set_title("Replicate-level vs pooled delta")
    axes[1].hist(pd.to_numeric(df["abs_delta_shift"], errors="coerce").dropna(), bins=12, color="#F28E2B")
    axes[1].set_title("Absolute delta shift")
    axes[1].set_xlabel("|pooled - replicate|")
    fig.tight_layout()
    return fig


def plot_bootstrap_forest(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    work = df.copy()
    work["abs_delta"] = pd.to_numeric(work["delta_rep"], errors="coerce").abs()
    work = work.sort_values("abs_delta", ascending=True).tail(30)
    y = np.arange(len(work))
    fig, ax = plt.subplots(figsize=(8, max(4.5, 0.26 * len(work))))
    stable = work["direction_stable"].astype(bool)
    colors = np.where(stable, "#2CA02C", "#D62728")
    ax.hlines(y, work["ci_low"], work["ci_high"], color=colors, linewidth=1.5)
    ax.scatter(work["delta_rep"], y, c=colors, s=24, zorder=3)
    ax.axvline(0, color="#444444", linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(work["region_id"].astype(str), fontsize=7)
    ax.set_xlabel("Delta replicate-level with bootstrap CI")
    ax.set_title("Bootstrap CI for replicate-level Delta")
    fig.tight_layout()
    return fig


def plot_delta_vs_coverage(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    work = df.copy()
    delta_col = "delta_rep" if "delta_rep" in work else "delta_methylation"
    cov_col = "total_coverage_A" if "total_coverage_A" in work else None
    if cov_col and "total_coverage_B" in work:
        coverage = pd.to_numeric(work["total_coverage_A"], errors="coerce") + pd.to_numeric(work["total_coverage_B"], errors="coerce")
    else:
        coverage = pd.Series(np.arange(len(work)) + 1)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(np.log10(coverage.replace(0, np.nan)), pd.to_numeric(work[delta_col], errors="coerce"), s=36, alpha=0.8, color="#4C78A8")
    ax.axhline(0, color="#444444", linewidth=0.8)
    ax.set_xlabel("log10 total coverage")
    ax.set_ylabel("Delta methylation")
    ax.set_title("Delta methylation vs coverage")
    fig.tight_layout()
    return fig


def plot_glm_vs_glmm(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    work = df.copy()
    gx = cap_neglog10_q(work["glm_q_value"])
    gy = cap_neglog10_q(work["glmm_q_value"])
    classes = np.select(
        [
            work.get("confirmed_by_glmm", False).astype(bool),
            work.get("glm_only_candidate", False).astype(bool),
            work.get("insufficient_glmm_data", False).astype(bool),
        ],
        ["confirmed_by_glmm", "glm_only_candidate", "insufficient_glmm_data"],
        default="neither",
    )
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for cls in pd.unique(classes):
        mask = classes == cls
        ax.scatter(gx[mask], gy[mask], s=44, color=STATUS_COLORS.get(cls, "#666666"), label=cls, alpha=0.9)
    thr = -np.log10(0.05)
    ax.axvline(thr, color="#444444", linestyle="--", linewidth=0.9)
    ax.axhline(thr, color="#444444", linestyle="--", linewidth=0.9)
    ax.set_xlabel("-log10 aggregated GLM q")
    ax.set_ylabel("-log10 CpG-level GLMM q")
    ax.set_title("Aggregated GLM vs CpG-level GLMM q-values")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def plot_model_status_summary(glm_glmm: pd.DataFrame, glmm: pd.DataFrame | None = None):
    setup_style()
    import matplotlib.pyplot as plt

    rows = []
    for col in ("confirmed_by_glmm", "glm_only_candidate", "glmm_more_conservative", "insufficient_glmm_data"):
        if col in glm_glmm:
            rows.append((col, int(glm_glmm[col].astype(bool).sum())))
    if glmm is not None and "model_status" in glmm:
        for status, n in glmm["model_status"].fillna("missing").value_counts().items():
            rows.append((f"GLMM:{status}", int(n)))
    table = pd.DataFrame(rows, columns=["status", "n"])
    return plot_bar(table, x="status", y="n", title="GLM/GLMM model status summary", ylabel="Regions", color="#72B7B2", rotate=35)


def plot_cpg_consistency(cpg_counts: pd.DataFrame, region_id: str):
    setup_style()
    import matplotlib.pyplot as plt

    df = cpg_counts[cpg_counts["region_id"].astype(str) == region_id].copy()
    if df.empty:
        raise ValueError(f"region not found: {region_id}")
    df["meth"] = pd.to_numeric(df["mC"], errors="coerce") / pd.to_numeric(df["total"], errors="coerce")
    pivot = df.pivot_table(index="position", columns="condition", values="meth", aggfunc="mean")
    if not {"Callus", "Seedling"}.issubset(pivot.columns):
        raise ValueError("condition columns Callus/Seedling not available")
    pivot["delta"] = pivot["Callus"] - pivot["Seedling"]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.scatter(pivot.index, pivot["delta"], color="#4C78A8", s=28)
    ax.axhline(float(pivot["delta"].mean()), color="#D62728", linestyle="--", label="mean CpG delta")
    ax.axhline(0, color="#444444", linewidth=0.8)
    ax.set_xlabel("CpG genomic position")
    ax.set_ylabel("Callus - Seedling methylation")
    ax.set_title(f"CpG-level consistency: {region_id}")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_dmr_metagene_density(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    contexts = ["CG", "CHG", "CHH"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4), sharey=False)
    for ax, context in zip(axes, contexts, strict=True):
        sub = df[df["context"].astype(str) == context]
        if sub.empty:
            ax.text(0.5, 0.5, f"{context}: missing", ha="center", va="center")
            ax.set_axis_off()
            continue
        for candidate_set, grp in sub.groupby("candidate_set", sort=False):
            if str(candidate_set) not in {"significant", "high_confidence", "significant_dmr_candidates", "high_confidence_candidate"}:
                continue
            ax.plot(grp["metagene_bin"], grp["fraction_of_candidate_set_context"], marker=None, linewidth=1.6, label=str(candidate_set))
        add_metagene_boundaries(ax, int(sub["metagene_bin"].max()) + 1)
        ax.set_title(context)
        ax.set_ylabel("Fraction of candidate set")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("DMR candidate density across MSU7 gene model", y=1.02)
    fig.tight_layout()
    return fig


def plot_observed_vs_random(bin_p: pd.DataFrame, density: pd.DataFrame | None = None):
    setup_style()
    import matplotlib.pyplot as plt

    df = bin_p.copy()
    y = "observed_density" if "observed_density" in df else "density_center"
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(df["bin_index"], df[y], color="#4C78A8", linewidth=1.8, label="Observed")
    if "empirical_q_BH" in df:
        sig = df[pd.to_numeric(df["empirical_q_BH"], errors="coerce") < 0.05]
        ax.scatter(sig["bin_index"], sig[y], color="#D62728", s=16, label="BH q < 0.05")
    add_metagene_boundaries(ax, int(df["bin_index"].max()) + 1)
    ax.set_ylabel("Fraction of genes with DMR center")
    ax.set_title("Observed vs random DMR occupancy profile")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig


def plot_projection_sensitivity(df: pd.DataFrame):
    setup_style()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(df["bin_index"], df["density_center"], label="center projection", color="#4C78A8", linewidth=1.8)
    ax.plot(df["bin_index"], df["density_interval_overlap"], label="interval-overlap projection", color="#F28E2B", linewidth=1.8)
    add_metagene_boundaries(ax, int(df["bin_index"].max()) + 1)
    ax.set_ylabel("Density")
    ax.set_title("Projection sensitivity: center vs interval-overlap")
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig
