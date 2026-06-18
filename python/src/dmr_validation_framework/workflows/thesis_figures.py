#!/usr/bin/env python
"""Build a generic LaTeX-ready DMR validation figure bundle."""

from __future__ import annotations

import argparse
import html
import json
import random
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from dmr_validation_framework.core.io import read_table
from dmr_validation_framework.core.stats import bh_qvalues
from dmr_validation_framework.reports.concordance import build_concordance_figure, compute_caller_concordance
from dmr_validation_framework.reports.funnel import build_funnel_figure, compute_funnel_stages
from dmr_validation_framework.reports.genome import build_manhattan_figure
from dmr_validation_framework.reports.palette import PALETTE as COLORS
from dmr_validation_framework.reports.palette import caller_cmap, legend_outside


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a generic thesis_figures_final/for_latex bundle from DMR validation outputs. "
            "The workflow is plant-agnostic and uses optional annotation inputs when available."
        )
    )
    parser.add_argument("--root", help="Dataset root.")
    parser.add_argument("--run-id", default="", help="Optional run id used only for default output naming.")
    parser.add_argument("--external-root", help="External caller output directory.")
    parser.add_argument("--report-dir", help="DMR validation report directory.")
    parser.add_argument("--validation-dir", help="Validation output directory.")
    parser.add_argument("--glm-glmm-dir", help="GLM/GLMM workflow output directory.")
    parser.add_argument("--annotation-dir", help="Plant annotation package directory.")
    parser.add_argument("--out-dir", help="Defaults to <root>/thesis_figures_final_<run-id> or <root>/thesis_figures_final.")
    parser.add_argument("--figure-prefix", default="fig_dmr", help="Prefix for generated LaTeX figure filenames.")
    parser.add_argument("--start-index", type=int, default=1)
    parser.add_argument("--max-forest-rows", type=int, default=40)
    parser.add_argument("--max-table-copy-rows", type=int, default=200_000)
    parser.add_argument("--expression-table", help="Optional gene-level expression evidence table for fig5_29.")
    parser.add_argument("--promoter-window", type=int, default=2000)
    parser.add_argument(
        "--skip-extended-metagene",
        action="store_true",
        help="Do not auto-build lightweight metagene/random-control outputs for fig5_22 and fig5_23.",
    )
    parser.add_argument("--random-control-iterations", type=int, default=100)
    parser.add_argument("--max-random-regions-per-iteration", type=int, default=5000)
    parser.add_argument("--max-metagene-dmrs", type=int, default=100000)
    parser.add_argument("--metagene-seed", type=int, default=202715)
    parser.add_argument("--upstream-len", type=int, default=2000)
    parser.add_argument("--downstream-len", type=int, default=2000)
    parser.add_argument("--upstream-bins", type=int, default=20)
    parser.add_argument("--body-bins", type=int, default=100)
    parser.add_argument("--downstream-bins", type=int, default=20)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def exists(path: Path | None) -> bool:
    return bool(path and path.exists() and path.stat().st_size > 0)


def latest_dir(root: Path, patterns: list[str]) -> Path | None:
    candidates: list[Path] = []
    for pattern in patterns:
        candidates.extend([p for p in root.glob(pattern) if p.is_dir()])
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def resolve_inputs(args: argparse.Namespace) -> dict[str, Path | None]:
    root = Path(args.root) if args.root else None
    run_id = args.run_id.strip()
    out_dir = Path(args.out_dir) if args.out_dir else None
    if out_dir is None:
        if root is None:
            out_dir = Path("outputs") / "thesis_figures_final"
        elif run_id:
            out_dir = root / f"thesis_figures_final_{run_id}"
        else:
            out_dir = root / "thesis_figures_final"

    def explicit_or_latest(value: str | None, patterns: list[str]) -> Path | None:
        if value:
            return Path(value)
        if root is None:
            return None
        if run_id:
            for pattern in patterns:
                candidate = root / pattern.format(run_id=run_id)
                if candidate.exists():
                    return candidate
        return latest_dir(root, [pattern.format(run_id="*") for pattern in patterns])

    report_dir = explicit_or_latest(args.report_dir, ["dmr_validation_report_{run_id}", "dmr_validation_report*"])
    validation_dir = explicit_or_latest(args.validation_dir, ["validation_core_{run_id}", "validation_{run_id}", "validation*"])
    glm_glmm_dir = explicit_or_latest(args.glm_glmm_dir, ["glm_glmm_methylkit_top30_{run_id}", "glm_glmm*"])
    external_root = explicit_or_latest(args.external_root, ["external_callers_combined_{run_id}", "external_callers*"])
    annotation_dir = explicit_or_latest(args.annotation_dir, ["annotation_itag3_2", "annotation*"])
    expression_table = Path(args.expression_table) if args.expression_table else None
    if expression_table is None and root is not None:
        candidates = [
            p
            for pattern in ["expression*/expression_support.tsv", "expression*/*expression_support*.tsv", "*expression_support*.tsv"]
            for p in root.glob(pattern)
            if p.is_file()
        ]
        if candidates:
            expression_table = max(candidates, key=lambda p: p.stat().st_mtime)

    return {
        "root": root,
        "out_dir": out_dir,
        "report_dir": report_dir,
        "validation_dir": validation_dir,
        "glm_glmm_dir": glm_glmm_dir,
        "external_root": external_root,
        "annotation_dir": annotation_dir,
        "expression_table": expression_table,
    }


def safe_read(path: Path | None) -> pd.DataFrame:
    if not exists(path):
        return pd.DataFrame()
    try:
        return read_table(path)
    except Exception:
        try:
            return pd.read_csv(path, sep=None, engine="python")
        except Exception:
            return pd.DataFrame()


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from dmr_validation_framework.reports.theme import apply_matplotlib_theme

    apply_matplotlib_theme()
    return plt


def style_ax(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=COLORS["lightgray"], linewidth=0.8)
    ax.set_axisbelow(True)


class ThesisBundle:
    def __init__(self, out_dir: Path, *, figure_prefix: str, start_index: int) -> None:
        self.out_dir = out_dir
        self.figure_prefix = figure_prefix
        self.counter = start_index
        self.rows: list[dict[str, object]] = []
        for sub in ("figures_pdf", "figures_png", "tables", "reports", "logs", "for_latex"):
            (out_dir / sub).mkdir(parents=True, exist_ok=True)

    def next_id(self, slug: str) -> str:
        figure_id = f"{self.figure_prefix}_{self.counter:02d}_{slug}"
        self.counter += 1
        return figure_id

    def record(
        self,
        *,
        figure_id: str,
        title: str,
        status: str,
        sources: list[Path],
        pdf: Path | None = None,
        png: Path | None = None,
        notes: str = "",
    ) -> None:
        latex_path = ""
        if exists(pdf):
            target = self.out_dir / "for_latex" / pdf.name
            shutil.copy2(pdf, target)
            latex_path = str(target)
        self.rows.append(
            {
                "figure_id": figure_id,
                "title": title,
                "status": status,
                "pdf": str(pdf) if pdf else "",
                "png": str(png) if png else "",
                "for_latex": latex_path,
                "source_files": ";".join(str(p) for p in sources if p),
                "notes": notes,
            }
        )

    def skip(self, slug: str, title: str, sources: list[Path], reason: str) -> None:
        figure_id = self.next_id(slug)
        self.record(figure_id=figure_id, title=title, status="skipped", sources=sources, notes=reason)

    def skip_named(self, figure_id: str, title: str, sources: list[Path], reason: str) -> None:
        self.record(figure_id=figure_id, title=title, status="skipped", sources=sources, notes=reason)

    def save_fig(self, fig, slug: str, title: str, sources: list[Path], notes: str = "") -> None:
        figure_id = self.next_id(slug)
        self.save_named_fig(fig, figure_id, title, sources, notes=notes)

    def save_named_fig(self, fig, figure_id: str, title: str, sources: list[Path], notes: str = "") -> None:
        pdf = self.out_dir / "figures_pdf" / f"{figure_id}.pdf"
        png = self.out_dir / "figures_png" / f"{figure_id}.png"
        fig.savefig(pdf, bbox_inches="tight")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        setup_matplotlib().close(fig)
        self.record(figure_id=figure_id, title=title, status="built", sources=sources, pdf=pdf, png=png, notes=notes)

    def wrap_png(self, source: Path, slug: str, title: str, notes: str = "converted from existing PNG report") -> None:
        if not exists(source):
            self.skip(slug, title, [source], "source PNG missing")
            return
        plt = setup_matplotlib()
        image = plt.imread(source)
        height, width = image.shape[:2]
        fig_w = min(12.0, max(6.0, width / 160.0))
        fig_h = min(9.0, max(4.0, height / 160.0))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.imshow(image)
        ax.set_axis_off()
        self.save_fig(fig, slug, title, [source], notes=notes)

    def copy_table(self, source: Path | None, name: str, max_rows: int) -> None:
        if not exists(source):
            return
        target = self.out_dir / "tables" / name
        try:
            df = read_table(source)
            if len(df) > max_rows:
                df = df.head(max_rows)
            df.to_csv(target, sep="\t", index=False)
        except Exception:
            shutil.copy2(source, target)

    def finalize(self, inputs: dict[str, Path | None]) -> None:
        manifest = pd.DataFrame(self.rows)
        manifest_path = self.out_dir / "tables" / "thesis_figure_manifest.tsv"
        manifest.to_csv(manifest_path, sep="\t", index=False)

        html_rows = manifest.to_html(index=False, escape=True) if not manifest.empty else "<p>No figures.</p>"
        image_tags = []
        for row in self.rows:
            pdf = Path(str(row.get("for_latex", "")))
            if pdf.name and row.get("status") == "built":
                image_tags.append(f"<li><code>{html.escape(pdf.name)}</code> - {html.escape(str(row.get('title', '')))}</li>")
        html_report = "\n".join(
            [
                "<!doctype html>",
                "<html><head><meta charset=\"utf-8\"><title>Thesis figure bundle</title>",
                "<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse;margin:12px 0}td,th{border:1px solid #ddd;padding:5px 7px;font-size:12px}th{background:#f4f4f4}code{background:#f6f6f6;padding:1px 3px}</style>",
                "</head><body>",
                "<h1>Generic DMR Thesis Figure Bundle</h1>",
                "<p>This bundle is generated from framework outputs. Figures are descriptive validation/downstream diagnostics, not a replacement for caller-native statistical inference.</p>",
                "<h2>Input directories</h2>",
                "<pre>" + html.escape(json.dumps({k: str(v) if v else "" for k, v in inputs.items()}, indent=2)) + "</pre>",
                "<h2>LaTeX PDFs</h2>",
                "<ul>" + "\n".join(image_tags) + "</ul>",
                "<h2>Manifest</h2>",
                html_rows,
                "</body></html>",
            ]
        )
        (self.out_dir / "reports" / "thesis_figure_bundle.html").write_text(html_report + "\n", encoding="utf-8")
        (self.out_dir / "for_latex" / "README.md").write_text(
            "# LaTeX-ready DMR figures\n\n"
            "Use the PDF files in this directory for LaTeX inclusion. See `../tables/thesis_figure_manifest.tsv` for sources and status.\n",
            encoding="utf-8",
        )


def plot_candidate_funnel(bundle: ThesisBundle, report_dir: Path | None, glm_glmm_dir: Path | None) -> None:
    caller_rows_path = report_dir / "caller_summary" / "caller_dmr_rows_for_plots.tsv" if report_dir else None
    consensus_candidates = [
        report_dir / "upset" / "dmr_caller_support_matrix.tsv" if report_dir else None,
        report_dir / "dmr_caller_support_matrix.tsv" if report_dir else None,
    ]
    consensus_path = next((p for p in consensus_candidates if exists(p)), None)
    glm_glmm_path = glm_glmm_dir / "glm_vs_glmm" / "glm_vs_glmm_comparison.tsv" if glm_glmm_dir else None
    sources = [p for p in [caller_rows_path, consensus_path, glm_glmm_path] if p]
    stages = compute_funnel_stages(
        consensus=safe_read(consensus_path) if consensus_path else None,
        glm_glmm=safe_read(glm_glmm_path) if glm_glmm_path else None,
        caller_rows=safe_read(caller_rows_path) if caller_rows_path else None,
    )
    if not stages:
        bundle.skip_named("fig5_00_candidate_funnel", "Candidate selection funnel", sources, "no funnel inputs available")
        return
    fig = build_funnel_figure(stages)
    if fig is None:
        bundle.skip_named("fig5_00_candidate_funnel", "Candidate selection funnel", sources, "funnel could not be built")
        return
    bundle.save_named_fig(fig, "fig5_00_candidate_funnel", "Candidate selection funnel", sources)


def plot_volcano(bundle: ThesisBundle, report_dir: Path | None, q_threshold: float = 0.1) -> None:
    path = report_dir / "caller_summary" / "caller_dmr_rows_for_plots.tsv" if report_dir else None
    df = safe_read(path)
    if df.empty or "delta" not in df.columns or "q_value" not in df.columns:
        bundle.skip_named("fig5_02_volcano", "Volcano (delta vs -log10 q)", [path] if path else [], "caller rows with delta/q_value missing")
        return
    df = df.copy()
    df["_x"] = pd.to_numeric(df["delta"], errors="coerce")
    df["_y"] = -np.log10(pd.to_numeric(df["q_value"], errors="coerce").clip(lower=1e-300))
    df = df.dropna(subset=["_x", "_y"])
    caller_col = "_caller" if "_caller" in df.columns else "caller" if "caller" in df.columns else None
    ctx_col = "_context" if "_context" in df.columns else "context" if "context" in df.columns else None
    if df.empty or caller_col is None:
        bundle.skip_named("fig5_02_volcano", "Volcano (delta vs -log10 q)", [path], "no plottable delta/q rows")
        return
    contexts = sorted(df[ctx_col].astype(str).unique()) if ctx_col else ["all"]
    cmap = caller_cmap(df[caller_col].astype(str).tolist())
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, len(contexts), figsize=(max(6, len(contexts) * 4.5), 4.6), squeeze=False)
    sig = float(-np.log10(max(q_threshold, 1e-300)))
    for ax, ctx in zip(axes[0], contexts):
        sub = df if ctx == "all" else df[df[ctx_col].astype(str) == ctx]
        for caller, group in sub.groupby(caller_col):
            ax.scatter(group["_x"], group["_y"], s=18, alpha=0.75, color=cmap.get(str(caller), COLORS["gray"]), label=str(caller), edgecolor="white", linewidth=0.3)
        ax.axhline(sig, color=COLORS["gray"], linestyle="--", linewidth=0.9)
        ax.axvline(0, color=COLORS["gray"], linestyle=":", linewidth=0.8)
        ax.set_title(str(ctx))
        ax.set_xlabel("delta methylation")
        ax.set_ylabel("-log10 q")
        style_ax(ax)
    legend_outside(axes[0][-1])
    fig.suptitle("DMR effect size vs significance", y=1.02)
    bundle.save_named_fig(fig, "fig5_02_volcano", "Volcano (delta vs -log10 q)", [path])


def plot_manhattan(bundle: ThesisBundle, report_dir: Path | None) -> None:
    path = report_dir / "caller_summary" / "caller_dmr_rows_for_plots.tsv" if report_dir else None
    df = safe_read(path)
    if df.empty or not {"chrom", "start", "end"}.issubset(df.columns):
        bundle.skip_named("fig5_03_genome_manhattan", "Genome-wide DMR distribution", [path] if path else [], "caller rows with coordinates missing")
        return
    if "_context" not in df.columns and "context" not in df.columns:
        df = df.copy()
        df["context"] = "NA"
    fig = build_manhattan_figure(df)
    if fig is None:
        bundle.skip_named("fig5_03_genome_manhattan", "Genome-wide DMR distribution", [path], "no plottable rows")
        return
    bundle.save_named_fig(fig, "fig5_03_genome_manhattan", "Genome-wide DMR distribution", [path])


def plot_caller_concordance(bundle: ThesisBundle, report_dir: Path | None) -> None:
    membership_path = report_dir / "upset" / "upset_membership.tsv" if report_dir else None
    membership = safe_read(membership_path)
    if membership.empty:
        bundle.skip_named("fig5_01_caller_concordance", "Caller concordance", [membership_path] if membership_path else [], "upset membership table missing")
        return
    long_df = compute_caller_concordance(membership)
    fig = build_concordance_figure(long_df)
    if fig is None:
        bundle.skip_named("fig5_01_caller_concordance", "Caller concordance", [membership_path], "no caller columns in membership table")
        return
    bundle.save_named_fig(fig, "fig5_01_caller_concordance", "Caller concordance (Jaccard)", [membership_path])


def add_existing_report_figures(bundle: ThesisBundle, report_dir: Path | None, glm_glmm_dir: Path | None) -> None:
    # Only the informative caller distributions and one GLMM-input QC histogram
    # are wrapped here. The low-information validation status-count bars and the
    # GLM/GLMM status plots are intentionally excluded from the PDF bundle:
    # they are covered better by the funnel/volcano/concordance/Manhattan
    # figures, the validation_summary.html, and the glm_vs_glmm tables.
    specs = [
        (report_dir / "caller_summary" / "dmr_count_by_caller_context.png" if report_dir else None, "caller_dmr_count", "DMR count by caller and context"),
        (report_dir / "caller_summary" / "abs_delta_distribution_by_caller_context.png" if report_dir else None, "caller_abs_delta", "Absolute delta distribution by context and caller"),
        (report_dir / "caller_summary" / "dmr_length_distribution_by_caller_context.png" if report_dir else None, "caller_dmr_length", "DMR length distribution by context and caller"),
        (glm_glmm_dir / "glmm_plots" / "top_cpg_contribution.png" if glm_glmm_dir else None, "top_cpg_contribution", "Top CpG contribution (GLMM input QC)"),
    ]
    for path, slug, title in specs:
        if path is not None:
            bundle.wrap_png(path, slug, title)


def plot_overlap_threshold(bundle: ThesisBundle, validation_dir: Path | None) -> None:
    path = validation_dir / "overlap_threshold_sensitivity.tsv" if validation_dir else None
    df = safe_read(path)
    metric_col = "jaccard_like" if "jaccard_like" in df.columns else "n_pairs"
    if df.empty or "threshold" not in df.columns or metric_col not in df.columns:
        bundle.skip_named("fig5_08_overlap_threshold_sensitivity", "Overlap threshold sensitivity", [path] if path else [], "table missing or incompatible")
        return
    df = df.copy()
    df["threshold"] = pd.to_numeric(df["threshold"], errors="coerce")
    df[metric_col] = pd.to_numeric(df[metric_col], errors="coerce")
    df = df.dropna(subset=["threshold", metric_col])
    if df.empty:
        bundle.skip_named("fig5_08_overlap_threshold_sensitivity", "Overlap threshold sensitivity", [path] if path else [], "no numeric threshold sensitivity rows")
        return
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    if "matching_policy" in df.columns:
        preferred = ["best_reciprocal", "one_to_one_greedy", "many_to_many"]
        seen = [str(value) for value in df["matching_policy"].dropna().unique()]
        policies = [policy for policy in preferred if policy in seen]
        policies.extend(sorted(policy for policy in seen if policy not in set(policies)))
        grouped = [(policy, df[df["matching_policy"].astype(str) == policy]) for policy in policies]
    else:
        policies = ["all"]
        grouped = [("all", df)]

    traces: dict[tuple[tuple[float, float, float, float], ...], dict[str, object]] = {}
    for label, sub in grouped:
        agg = (
            sub.groupby("threshold", as_index=False)[metric_col]
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

    colors = [COLORS["blue"], COLORS["orange"], COLORS["green"], COLORS["red"], COLORS["gray"]]
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
    if df["threshold"].min() <= 0.5 <= df["threshold"].max():
        ax.axvline(0.5, color=COLORS["gray"], linestyle="--", linewidth=0.9, alpha=0.7)
    metric_label = "Jaccard-like agreement" if metric_col == "jaccard_like" else "Matched DMR pairs"
    note = "median; ribbon = min-max across caller pairs"
    if df.groupby("threshold")[metric_col].median().nunique(dropna=True) == 1:
        note += "; no change across tested thresholds"
    ax.text(0.01, 0.98, note, transform=ax.transAxes, ha="left", va="top", fontsize=8, color=COLORS["gray"])
    ax.set_xlabel("Strict reciprocal overlap threshold (tau)")
    ax.set_ylabel(metric_label)
    ax.set_title("Consensus sensitivity to overlap threshold")
    if metric_col == "jaccard_like":
        ymax = float(df[metric_col].max())
        ax.set_ylim(0, min(1.0, max(0.05, ymax + 0.08)))
    style_ax(ax)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    bundle.save_named_fig(fig, "fig5_08_overlap_threshold_sensitivity", "Overlap threshold sensitivity", [path])


def plot_delta_weighting(bundle: ThesisBundle, validation_dir: Path | None) -> None:
    path = validation_dir / "delta_weighting_sensitivity.tsv" if validation_dir else None
    df = safe_read(path)
    required = {"delta_rep", "delta_pool", "abs_delta_shift"}
    if df.empty or not required.issubset(df.columns):
        bundle.skip_named("fig5_11_delta_rep_vs_pooled_sensitivity", "Replicate vs pooled delta", [path] if path else [], "table missing or incompatible")
        return
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    delta_rep = pd.to_numeric(df["delta_rep"], errors="coerce")
    delta_pool = pd.to_numeric(df["delta_pool"], errors="coerce")
    changed = df.get("direction_changed", pd.Series(False, index=df.index)).astype(str).str.lower().isin(["true", "1", "yes"])
    colors = np.where(changed, COLORS["red"], COLORS["blue"])
    axes[0].scatter(delta_rep, delta_pool, c=colors, s=40, edgecolor="white", linewidth=0.4)
    lo = np.nanmin([delta_rep.min(), delta_pool.min()])
    hi = np.nanmax([delta_rep.max(), delta_pool.max()])
    if np.isfinite(lo) and np.isfinite(hi):
        axes[0].plot([lo, hi], [lo, hi], color=COLORS["gray"], linestyle="--")
    axes[0].set_xlabel("Replicate-level delta")
    axes[0].set_ylabel("Pooled read-level delta")
    axes[0].set_title("Delta estimates")
    style_ax(axes[0])
    axes[1].hist(pd.to_numeric(df["abs_delta_shift"], errors="coerce").dropna(), bins=20, color=COLORS["orange"], edgecolor="white")
    axes[1].set_xlabel("|delta_pool - delta_rep|")
    axes[1].set_ylabel("Regions")
    axes[1].set_title("Absolute shift")
    style_ax(axes[1])
    fig.suptitle("Replicate-level vs pooled delta sensitivity", y=1.02)
    bundle.save_named_fig(fig, "fig5_11_delta_rep_vs_pooled_sensitivity", "Replicate-level vs pooled delta sensitivity", [path])


def plot_bootstrap_ci(bundle: ThesisBundle, validation_dir: Path | None, max_rows: int) -> None:
    path = validation_dir / "delta_bootstrap_ci.tsv" if validation_dir else None
    df = safe_read(path)
    required = {"region_id", "delta_rep", "ci_low", "ci_high"}
    if df.empty or not required.issubset(df.columns):
        bundle.skip_named("fig5_12_delta_bootstrap_ci_forest", "Bootstrap CI for delta", [path] if path else [], "table missing or incompatible")
        return
    work = df.copy()
    work["delta_rep"] = pd.to_numeric(work["delta_rep"], errors="coerce")
    work["ci_low"] = pd.to_numeric(work["ci_low"], errors="coerce")
    work["ci_high"] = pd.to_numeric(work["ci_high"], errors="coerce")
    work["rank_abs"] = work["delta_rep"].abs()
    work = work.dropna(subset=["delta_rep", "ci_low", "ci_high"]).sort_values("rank_abs", ascending=False).head(max_rows)
    if work.empty:
        bundle.skip_named("fig5_12_delta_bootstrap_ci_forest", "Bootstrap CI for delta", [path], "no plottable rows")
        return
    work = work.sort_values("delta_rep")
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8, max(4.5, len(work) * 0.22)))
    y = np.arange(len(work))
    xerr = np.vstack([work["delta_rep"] - work["ci_low"], work["ci_high"] - work["delta_rep"]])
    stable = work.get("direction_stable", pd.Series(True, index=work.index)).astype(str).str.lower().isin(["true", "1", "yes"])
    ax.errorbar(work["delta_rep"], y, xerr=xerr, fmt="none", ecolor="#9CA3AF", elinewidth=1.2, capsize=2)
    ax.scatter(work["delta_rep"], y, c=np.where(stable, COLORS["blue"], COLORS["red"]), s=24, zorder=3)
    ax.axvline(0, color=COLORS["gray"], linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(work["region_id"].astype(str), fontsize=7)
    ax.set_xlabel("Replicate-level delta with bootstrap CI")
    ax.set_title("Bootstrap CI for DMR effect size")
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_12_delta_bootstrap_ci_forest", "Bootstrap CI for delta", [path])


def plot_coverage_qc(bundle: ThesisBundle, validation_dir: Path | None) -> None:
    summary_path = validation_dir / "coverage_set_summary.tsv" if validation_dir else None
    region_path = validation_dir / "coverage_set_region_qc.tsv" if validation_dir else None
    summary = safe_read(summary_path)
    region = safe_read(region_path)
    if summary.empty and region.empty:
        bundle.skip_named("fig5_09_coverage_set_common_cpg_qc", "Coverage common-CpG QC", [p for p in [summary_path, region_path] if p], "coverage tables missing")
        return
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    if not summary.empty:
        row = summary.iloc[0]
        labels = ["pass", "insufficient", "none"]
        values = [
            int(pd.to_numeric(pd.Series([row.get("n_regions_pass_common_cpg", 0)]), errors="coerce").fillna(0).iloc[0]),
            int(pd.to_numeric(pd.Series([row.get("n_regions_insufficient_common_cpg", 0)]), errors="coerce").fillna(0).iloc[0]),
            int(pd.to_numeric(pd.Series([row.get("n_regions_no_common_cpgs", 0)]), errors="coerce").fillna(0).iloc[0]),
        ]
        axes[0].bar(labels, values, color=[COLORS["green"], COLORS["orange"], COLORS["red"]])
        axes[0].set_title("Common-CpG region status")
        axes[0].set_ylabel("Regions")
    else:
        axes[0].set_axis_off()
    style_ax(axes[0])
    if not region.empty and "fraction_common" in region.columns:
        vals = pd.to_numeric(region["fraction_common"], errors="coerce").dropna()
        axes[1].hist(vals, bins=20, color=COLORS["blue"], edgecolor="white")
        axes[1].set_xlabel("Fraction common CpGs")
        axes[1].set_ylabel("Regions")
        axes[1].set_title("Common-CpG fraction")
        style_ax(axes[1])
    else:
        axes[1].set_axis_off()
    bundle.save_named_fig(fig, "fig5_09_coverage_set_common_cpg_qc", "Coverage common-CpG QC", [p for p in [summary_path, region_path] if p])


def plot_region_cpg_coverage_distribution(bundle: ThesisBundle, validation_dir: Path | None, glm_glmm_dir: Path | None) -> None:
    region_path = validation_dir / "coverage_set_region_qc.tsv" if validation_dir else None
    sample_path = validation_dir / "coverage_set_sample_qc.tsv" if validation_dir else None
    counts_path = glm_glmm_dir / "region_cpg_counts_region_summary.tsv" if glm_glmm_dir else None
    region = safe_read(region_path)
    sample = safe_read(sample_path)
    counts = safe_read(counts_path)
    if region.empty and sample.empty and counts.empty:
        bundle.skip_named(
            "fig5_10_region_cpg_coverage_distribution",
            "Region CpG and coverage distribution",
            [p for p in [region_path, sample_path, counts_path] if p],
            "coverage distribution tables missing",
        )
        return
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    if not region.empty and "n_cpg_common" in region.columns:
        vals = pd.to_numeric(region["n_cpg_common"], errors="coerce").dropna()
        axes[0].hist(vals, bins=min(30, max(5, len(vals))), color=COLORS["blue"], edgecolor="white")
        axes[0].set_xlabel("Common CpGs per region")
        axes[0].set_ylabel("Regions")
        axes[0].set_title("Common CpGs")
        style_ax(axes[0])
    else:
        axes[0].set_axis_off()
    if not sample.empty and "total_coverage_common" in sample.columns:
        vals = pd.to_numeric(sample["total_coverage_common"], errors="coerce").dropna()
        axes[1].hist(np.log10(vals + 1), bins=30, color=COLORS["orange"], edgecolor="white")
        axes[1].set_xlabel("log10(common coverage + 1)")
        axes[1].set_ylabel("Region x sample rows")
        axes[1].set_title("Common coverage")
        style_ax(axes[1])
    elif not counts.empty and "mean_total" in counts.columns:
        vals = pd.to_numeric(counts["mean_total"], errors="coerce").dropna()
        axes[1].hist(np.log10(vals + 1), bins=30, color=COLORS["orange"], edgecolor="white")
        axes[1].set_xlabel("log10(mean total coverage + 1)")
        axes[1].set_ylabel("Rows")
        axes[1].set_title("Mean coverage")
        style_ax(axes[1])
    else:
        axes[1].set_axis_off()
    if not region.empty and "fraction_common" in region.columns:
        vals = pd.to_numeric(region["fraction_common"], errors="coerce").dropna()
        axes[2].hist(vals, bins=20, color=COLORS["green"], edgecolor="white")
        axes[2].set_xlabel("Fraction common CpGs")
        axes[2].set_ylabel("Regions")
        axes[2].set_title("Common fraction")
        style_ax(axes[2])
    else:
        axes[2].set_axis_off()
    fig.suptitle("Region-level CpG and coverage diagnostics", y=1.02)
    bundle.save_named_fig(
        fig,
        "fig5_10_region_cpg_coverage_distribution",
        "Region CpG and coverage distribution",
        [p for p in [region_path, sample_path, counts_path] if p],
    )


def plot_delta_vs_coverage(bundle: ThesisBundle, validation_dir: Path | None) -> None:
    delta_path = validation_dir / "delta_weighting_sensitivity.tsv" if validation_dir else None
    sample_path = validation_dir / "coverage_set_sample_qc.tsv" if validation_dir else None
    delta = safe_read(delta_path)
    sample = safe_read(sample_path)
    if delta.empty or sample.empty or "region_id" not in delta.columns or "region_id" not in sample.columns:
        bundle.skip_named(
            "fig5_13_delta_vs_mean_or_coverage",
            "Delta vs coverage",
            [p for p in [delta_path, sample_path] if p],
            "delta or sample coverage table missing",
        )
        return
    if "delta_rep" not in delta.columns or "total_coverage_common" not in sample.columns:
        bundle.skip_named(
            "fig5_13_delta_vs_mean_or_coverage",
            "Delta vs coverage",
            [p for p in [delta_path, sample_path] if p],
            "required columns missing",
        )
        return
    cov = (
        sample.assign(total_coverage_common=pd.to_numeric(sample["total_coverage_common"], errors="coerce"))
        .groupby("region_id", as_index=False)["total_coverage_common"]
        .sum()
        .rename(columns={"total_coverage_common": "total_common_coverage"})
    )
    work = delta.merge(cov, on="region_id", how="left")
    work["delta_abs"] = pd.to_numeric(work["delta_rep"], errors="coerce").abs()
    work["total_common_coverage"] = pd.to_numeric(work["total_common_coverage"], errors="coerce")
    if "context" not in work.columns:
        work["context"] = work["region_id"].astype(str).str.extract(r"_(CG|CHG|CHH)")[0].fillna("NA")
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.8))
    colors = {"CG": COLORS["blue"], "CHG": COLORS["orange"], "CHH": COLORS["green"], "NA": COLORS["gray"]}
    for context, sub in work.dropna(subset=["delta_abs", "total_common_coverage"]).groupby("context", dropna=False):
        ax.scatter(
            np.log10(sub["total_common_coverage"] + 1),
            sub["delta_abs"],
            s=45,
            alpha=0.85,
            color=colors.get(str(context), COLORS["gray"]),
            label=str(context),
            edgecolor="white",
            linewidth=0.4,
        )
    ax.set_xlabel("log10(total common coverage + 1)")
    ax.set_ylabel("|replicate-level delta|")
    ax.set_title("Effect size vs common coverage")
    legend_outside(ax)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_13_delta_vs_mean_or_coverage", "Delta vs coverage", [delta_path, sample_path])


def _section_for_bin(bin_index: int, args: argparse.Namespace) -> str:
    if bin_index < args.upstream_bins:
        return "upstream"
    if bin_index < args.upstream_bins + args.body_bins:
        return "gene_body"
    return "downstream"


def _map_pos_to_metagene_bin(pos: float, gene, args: argparse.Namespace) -> int | None:
    start = float(gene.start)
    end = float(gene.end)
    strand = str(gene.strand)
    if strand == "+":
        if start - args.upstream_len <= pos < start:
            rel = (pos - (start - args.upstream_len)) / max(1, args.upstream_len)
            return min(args.upstream_bins - 1, max(0, int(rel * args.upstream_bins)))
        if start <= pos <= end:
            rel = (pos - start) / max(1, end - start + 1)
            return args.upstream_bins + min(args.body_bins - 1, max(0, int(rel * args.body_bins)))
        if end < pos <= end + args.downstream_len:
            rel = (pos - end) / max(1, args.downstream_len)
            return args.upstream_bins + args.body_bins + min(args.downstream_bins - 1, max(0, int(rel * args.downstream_bins)))
    elif strand == "-":
        if end < pos <= end + args.upstream_len:
            rel = ((end + args.upstream_len) - pos) / max(1, args.upstream_len)
            return min(args.upstream_bins - 1, max(0, int(rel * args.upstream_bins)))
        if start <= pos <= end:
            rel = (end - pos) / max(1, end - start + 1)
            return args.upstream_bins + min(args.body_bins - 1, max(0, int(rel * args.body_bins)))
        if start - args.downstream_len <= pos < start:
            rel = (start - pos) / max(1, args.downstream_len)
            return args.upstream_bins + args.body_bins + min(args.downstream_bins - 1, max(0, int(rel * args.downstream_bins)))
    return None


def _metagene_bin_interval(gene, bin_index: int, args: argparse.Namespace) -> tuple[float, float]:
    start = float(gene.start)
    end = float(gene.end)
    strand = str(gene.strand)
    if bin_index < args.upstream_bins:
        frac0 = bin_index / args.upstream_bins
        frac1 = (bin_index + 1) / args.upstream_bins
        if strand == "+":
            return start - args.upstream_len + frac0 * args.upstream_len, start - args.upstream_len + frac1 * args.upstream_len
        return end + (1 - frac1) * args.upstream_len, end + (1 - frac0) * args.upstream_len
    if bin_index < args.upstream_bins + args.body_bins:
        local = bin_index - args.upstream_bins
        frac0 = local / args.body_bins
        frac1 = (local + 1) / args.body_bins
        span = end - start + 1
        if strand == "+":
            return start + frac0 * span, start + frac1 * span
        return end - frac1 * span, end - frac0 * span
    local = bin_index - args.upstream_bins - args.body_bins
    frac0 = local / args.downstream_bins
    frac1 = (local + 1) / args.downstream_bins
    if strand == "+":
        return end + frac0 * args.downstream_len, end + frac1 * args.downstream_len
    return start - frac1 * args.downstream_len, start - frac0 * args.downstream_len


def _build_gene_window_index(genes: pd.DataFrame, args: argparse.Namespace) -> tuple[dict[str, dict[str, object]], int]:
    index: dict[str, dict[str, object]] = {}
    max_window = 0
    for chrom, sub in genes.groupby("chrom", sort=False):
        work = sub.copy()
        work["start"] = pd.to_numeric(work["start"], errors="coerce")
        work["end"] = pd.to_numeric(work["end"], errors="coerce")
        work = work.dropna(subset=["start", "end"]).copy()
        if work.empty:
            continue
        work["start"] = work["start"].astype(int)
        work["end"] = work["end"].astype(int)
        work["window_start"] = np.maximum(1, work["start"] - args.downstream_len)
        work["window_end"] = work["end"] + args.upstream_len
        work = work.sort_values("window_start").reset_index(drop=True)
        max_window = max(max_window, int((work["window_end"] - work["window_start"]).max()))
        index[str(chrom)] = {"data": work, "window_starts": work["window_start"].to_numpy(dtype=int)}
    return index, max_window


def _project_dmr_center(chrom: str, center: float, gene_index: dict[str, dict[str, object]], max_window: int, args: argparse.Namespace):
    chrom_index = gene_index.get(chrom)
    if chrom_index is None:
        return None
    genes = chrom_index["data"]
    starts = chrom_index["window_starts"]
    lo = max(0, int(np.searchsorted(starts, center - max_window - 1, side="left")))
    hi = int(np.searchsorted(starts, center, side="right"))
    best = None
    best_dist = None
    for gene in genes.iloc[lo:hi].itertuples(index=False):
        if not (gene.window_start <= center <= gene.window_end):
            continue
        bin_index = _map_pos_to_metagene_bin(center, gene, args)
        if bin_index is None:
            continue
        tss = float(gene.start if str(gene.strand) == "+" else gene.end)
        dist = abs(center - tss)
        if best is None or best_dist is None or dist < best_dist:
            best = gene, int(bin_index), float(dist)
            best_dist = dist
    return best


def _occupancy_density(mapped: pd.DataFrame, n_genes: int, n_bins: int) -> np.ndarray:
    density = np.zeros(n_bins, dtype=float)
    if mapped.empty or "gene_id" not in mapped.columns or "bin_index" not in mapped.columns:
        return density
    pairs = mapped[["gene_id", "bin_index"]].dropna().drop_duplicates()
    for bin_index, count in pairs.groupby("bin_index").size().items():
        try:
            idx = int(bin_index)
        except Exception:
            continue
        if 0 <= idx < n_bins:
            density[idx] = float(count) / max(1, n_genes)
    return density


def build_extended_metagene_outputs(
    report_dir: Path | None,
    annotation_dir: Path | None,
    out_dir: Path,
    args: argparse.Namespace,
) -> Path | None:
    dmr_path = report_dir / "caller_summary" / "caller_dmr_rows_for_plots.tsv" if report_dir else None
    gene_path = annotation_dir / "metagene_gene_regions_qc.tsv" if annotation_dir else None
    dmrs = safe_read(dmr_path)
    genes = safe_read(gene_path)
    extended_dir = out_dir / "tables" / "extended_metagene"
    extended_dir.mkdir(parents=True, exist_ok=True)
    n_bins = args.upstream_bins + args.body_bins + args.downstream_bins
    if dmrs.empty or genes.empty or not {"chrom", "start", "end"}.issubset(dmrs.columns) or not {"chrom", "start", "end", "gene_id", "strand"}.issubset(genes.columns):
        pd.DataFrame(
            [
                {
                    "status": "SKIPPED",
                    "notes": "DMR rows or metagene gene annotation missing",
                    "dmr_rows": str(dmr_path) if dmr_path else "missing",
                    "genes": str(gene_path) if gene_path else "missing",
                }
            ]
        ).to_csv(extended_dir / "extended_metagene_summary.tsv", sep="\t", index=False)
        return extended_dir

    work_dmrs = dmrs.copy()
    work_dmrs["chrom"] = work_dmrs["chrom"].astype(str)
    work_dmrs["start"] = pd.to_numeric(work_dmrs["start"], errors="coerce")
    work_dmrs["end"] = pd.to_numeric(work_dmrs["end"], errors="coerce")
    work_dmrs = work_dmrs.dropna(subset=["start", "end"]).copy()
    work_dmrs["start"] = work_dmrs["start"].astype(int)
    work_dmrs["end"] = work_dmrs["end"].astype(int)
    work_dmrs = work_dmrs[work_dmrs["end"] >= work_dmrs["start"]].copy()
    if work_dmrs.empty:
        return extended_dir
    if len(work_dmrs) > args.max_metagene_dmrs:
        work_dmrs = work_dmrs.sample(args.max_metagene_dmrs, random_state=args.metagene_seed).copy()
        sampling_note = f"DMRs sampled to {args.max_metagene_dmrs} for lightweight metagene audit"
    else:
        sampling_note = "all DMR rows used"
    if "region_id" not in work_dmrs.columns:
        work_dmrs["region_id"] = [f"dmr_{i + 1}" for i in range(len(work_dmrs))]
    if "caller" not in work_dmrs.columns and "_caller" in work_dmrs.columns:
        work_dmrs["caller"] = work_dmrs["_caller"]
    if "context" not in work_dmrs.columns and "_context" in work_dmrs.columns:
        work_dmrs["context"] = work_dmrs["_context"]

    genes = genes.copy()
    genes["chrom"] = genes["chrom"].astype(str)
    gene_index, max_window = _build_gene_window_index(genes, args)
    mapped_rows: list[dict[str, object]] = []
    for row in work_dmrs.itertuples(index=False):
        chrom = str(getattr(row, "chrom"))
        start = int(getattr(row, "start"))
        end = int(getattr(row, "end"))
        center = (start + end) / 2.0
        projected = _project_dmr_center(chrom, center, gene_index, max_window, args)
        if projected is None:
            continue
        gene, bin_index, distance = projected
        mapped_rows.append(
            {
                "region_id": getattr(row, "region_id", ""),
                "caller": getattr(row, "caller", ""),
                "context": getattr(row, "context", ""),
                "chrom": chrom,
                "start": start,
                "end": end,
                "length_bp": end - start + 1,
                "center": center,
                "gene_id": str(gene.gene_id),
                "strand": str(gene.strand),
                "bin_index": bin_index,
                "bin": f"bin_{bin_index}",
                "section": _section_for_bin(bin_index, args),
                "distance_to_tss": distance,
            }
        )
    mapped = pd.DataFrame(mapped_rows)
    mapped.to_csv(extended_dir / "mapped_dmr_centers.tsv", sep="\t", index=False)
    if mapped.empty:
        pd.DataFrame(
            [
                {
                    "status": "SKIPPED",
                    "n_dmr_input": len(work_dmrs),
                    "n_mapped_centers": 0,
                    "notes": "No DMR centers mapped into gene metagene windows",
                }
            ]
        ).to_csv(extended_dir / "extended_metagene_summary.tsv", sep="\t", index=False)
        return extended_dir

    n_genes = int(genes["gene_id"].nunique())
    center_density_all = _occupancy_density(mapped, n_genes, n_bins)

    genes_by_id = genes.drop_duplicates("gene_id").set_index("gene_id")
    interval_density = np.zeros(n_bins, dtype=float)
    changed_rows: list[dict[str, object]] = []
    interval_pairs: set[tuple[str, int]] = set()
    for row in mapped.itertuples(index=False):
        gene_id = str(row.gene_id)
        if gene_id not in genes_by_id.index:
            continue
        gene = genes_by_id.loc[gene_id]
        overlap_by_bin: list[tuple[int, float]] = []
        for bin_index in range(n_bins):
            b0, b1 = _metagene_bin_interval(gene, bin_index, args)
            lo = max(min(b0, b1), float(row.start))
            hi = min(max(b0, b1), float(row.end))
            overlap = max(0.0, hi - lo + 1)
            if overlap > 0:
                interval_pairs.add((gene_id, bin_index))
                overlap_by_bin.append((bin_index, overlap))
        dominant_bin = max(overlap_by_bin, key=lambda x: x[1])[0] if overlap_by_bin else None
        center_bin = int(row.bin_index)
        if dominant_bin is None or dominant_bin != center_bin:
            changed_rows.append(
                {
                    "region_id": row.region_id,
                    "gene_id": gene_id,
                    "chrom": row.chrom,
                    "start": row.start,
                    "end": row.end,
                    "center_bin": center_bin,
                    "dominant_interval_bin": dominant_bin,
                    "center_section": _section_for_bin(center_bin, args),
                    "dominant_interval_section": _section_for_bin(dominant_bin, args) if dominant_bin is not None else "unmapped",
                    "n_interval_bins": len({b for b, _ in overlap_by_bin}),
                }
            )
    for _, bin_index in interval_pairs:
        if 0 <= bin_index < n_bins:
            interval_density[bin_index] += 1.0 / max(1, n_genes)
    density_rows = [
        {
            "bin": f"bin_{i}",
            "bin_index": i,
            "section": _section_for_bin(i, args),
            "density_center": center_density_all[i],
            "density_interval_overlap": interval_density[i],
            "density_difference_interval_minus_center": interval_density[i] - center_density_all[i],
        }
        for i in range(n_bins)
    ]
    pd.DataFrame(density_rows).to_csv(extended_dir / "projection_sensitivity_density.tsv", sep="\t", index=False)
    corr = float(np.corrcoef(center_density_all, interval_density)[0, 1]) if n_bins > 1 and np.isfinite(interval_density).all() else np.nan
    pd.DataFrame(
        [
            {
                "status": "PASS",
                "n_dmr_input": len(work_dmrs),
                "n_center_mapped": len(mapped),
                "n_changed_regions": len(changed_rows),
                "density_profile_correlation": corr,
                "notes": sampling_note,
            }
        ]
    ).to_csv(extended_dir / "projection_sensitivity_summary.tsv", sep="\t", index=False)
    pd.DataFrame(changed_rows or [{"status": "PASS", "notes": "no changed regions"}]).to_csv(
        extended_dir / "projection_changed_regions.tsv", sep="\t", index=False
    )

    random_n = len(mapped) if args.max_random_regions_per_iteration <= 0 else min(len(mapped), args.max_random_regions_per_iteration)
    observed_for_random = mapped.sample(random_n, random_state=args.metagene_seed).copy() if len(mapped) > random_n else mapped.copy()
    observed_density = _occupancy_density(observed_for_random, n_genes, n_bins)
    chrom_sizes = (
        pd.concat(
            [
                genes.groupby("chrom")["end"].max().rename("gene_end"),
                work_dmrs.groupby("chrom")["end"].max().rename("dmr_end"),
            ],
            axis=1,
        )
        .max(axis=1)
        .dropna()
        .astype(int)
        .to_dict()
    )
    chrom_weights = observed_for_random["chrom"].astype(str).value_counts(normalize=True).to_dict()
    chroms = [chrom for chrom in chrom_weights if chrom in chrom_sizes and chrom in gene_index]
    weights = [chrom_weights[chrom] for chrom in chroms]
    rng = random.Random(args.metagene_seed)
    ge_counts = np.zeros(n_bins, dtype=int)
    random_sum = np.zeros(n_bins, dtype=float)
    section_rows: list[dict[str, object]] = []
    random_global_max: list[float] = []
    lengths = observed_for_random["length_bp"].astype(int).clip(lower=1).tolist()
    for iteration in range(max(0, args.random_control_iterations)):
        random_rows = []
        for length in lengths:
            if not chroms:
                break
            chrom = rng.choices(chroms, weights=weights, k=1)[0]
            chrom_size = int(chrom_sizes[chrom])
            if chrom_size <= length + 2:
                continue
            start = rng.randint(1, max(1, chrom_size - length))
            end = start + length - 1
            center = (start + end) / 2.0
            projected = _project_dmr_center(chrom, center, gene_index, max_window, args)
            if projected is None:
                continue
            gene, bin_index, _ = projected
            random_rows.append({"gene_id": str(gene.gene_id), "bin_index": int(bin_index)})
        random_mapped = pd.DataFrame(random_rows)
        random_density = _occupancy_density(random_mapped, n_genes, n_bins)
        random_sum += random_density
        ge_counts += random_density >= observed_density
        random_global_max.append(float(random_density.max()) if len(random_density) else 0.0)
        for section in ["upstream", "gene_body", "downstream"]:
            idx = [i for i in range(n_bins) if _section_for_bin(i, args) == section]
            section_rows.append(
                {
                    "iteration": iteration,
                    "section": section,
                    "random_mean_density": float(random_density[idx].mean()) if idx else 0.0,
                    "n_projected_gene_bins": int(len(random_mapped.drop_duplicates())) if not random_mapped.empty else 0,
                }
            )
    iterations = max(1, args.random_control_iterations)
    random_mean = random_sum / iterations
    empirical_p = (1 + ge_counts) / (iterations + 1)
    empirical_q = bh_qvalues(empirical_p)
    observed_global_max_density = float(observed_density.max()) if len(observed_density) else 0.0
    global_profile_p = (1 + sum(x >= observed_global_max_density for x in random_global_max)) / (iterations + 1)
    random_bin_rows = [
        {
            "bin": f"bin_{i}",
            "bin_index": i,
            "section": _section_for_bin(i, args),
            "observed_density": float(observed_density[i]),
            "random_mean_density": float(random_mean[i]),
            "empirical_p_random_ge_observed": float(empirical_p[i]),
            "empirical_q_BH": float(empirical_q[i]) if not pd.isna(empirical_q[i]) else np.nan,
            "global_profile_p": global_profile_p,
            "observed_global_max_density": observed_global_max_density,
        }
        for i in range(n_bins)
    ]
    pd.DataFrame(random_bin_rows).to_csv(extended_dir / "random_control_bin_empirical_p.tsv", sep="\t", index=False)
    pd.DataFrame(section_rows).to_csv(extended_dir / "random_control_density_by_section.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "status": "PASS",
                "n_observed_mapped_available": len(mapped),
                "n_observed_used_for_random": len(observed_for_random),
                "iterations_requested": args.random_control_iterations,
                "n_genes_denominator": n_genes,
                "observed_global_max_density": observed_global_max_density,
                "global_profile_p": global_profile_p,
                "notes": (
                    f"{sampling_note}; random controls use matched DMR lengths and chromosome distribution; "
                    "chromosome sizes are inferred from annotation and DMR rows when no FASTA index is supplied"
                ),
            }
        ]
    ).to_csv(extended_dir / "random_control_occupancy_summary.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "status": "PASS",
                "n_dmr_input": len(work_dmrs),
                "n_center_mapped": len(mapped),
                "n_genes_denominator": n_genes,
                "n_random_iterations": args.random_control_iterations,
                "n_random_regions_per_iteration": len(observed_for_random),
                "notes": sampling_note,
            }
        ]
    ).to_csv(extended_dir / "extended_metagene_summary.tsv", sep="\t", index=False)
    return extended_dir


def _extended_path(primary_dir: Path | None, fallback_dir: Path | None, filename: str) -> Path | None:
    primary = primary_dir / filename if primary_dir else None
    if exists(primary):
        return primary
    fallback = fallback_dir / filename if fallback_dir else None
    if exists(fallback):
        return fallback
    return primary or fallback


def plot_random_occupancy(bundle: ThesisBundle, validation_dir: Path | None, extended_dir: Path | None = None) -> None:
    bin_path = _extended_path(validation_dir, extended_dir, "random_control_bin_empirical_p.tsv")
    section_path = _extended_path(validation_dir, extended_dir, "random_control_density_by_section.tsv")
    summary_path = _extended_path(validation_dir, extended_dir, "random_control_occupancy_summary.tsv")
    df = safe_read(bin_path)
    if df.empty:
        bundle.skip_named(
            "fig5_22_observed_vs_random_occupancy_profile",
            "Observed vs random occupancy profile",
            [p for p in [bin_path, section_path, summary_path] if p],
            "random-control occupancy tables missing; provide annotation + DMR rows or run extended profile/random_control_occupancy first",
        )
        return
    required = {"bin_index", "observed_density"}
    if not required.issubset(df.columns):
        bundle.skip_named("fig5_22_observed_vs_random_occupancy_profile", "Observed vs random occupancy profile", [bin_path], "required columns missing")
        return
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(9, 4.5))
    df = df.sort_values("bin_index")
    ax.plot(df["bin_index"], df["observed_density"], color=COLORS["blue"], linewidth=2, label="observed")
    if "random_mean_density" in df.columns:
        ax.plot(df["bin_index"], df["random_mean_density"], color=COLORS["orange"], linewidth=1.5, label="random mean")
    if "empirical_q_BH" in df.columns:
        sig = df[pd.to_numeric(df["empirical_q_BH"], errors="coerce") < 0.05]
        if not sig.empty:
            ax.scatter(sig["bin_index"], sig["observed_density"], color=COLORS["red"], s=16, label="BH q<0.05")
    ax.set_xlabel("Metagene bin")
    ax.set_ylabel("Occupancy density")
    ax.set_title("Observed DMR occupancy vs random controls")
    legend_outside(ax)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_22_observed_vs_random_occupancy_profile", "Observed vs random occupancy profile", [p for p in [bin_path, section_path, summary_path] if p])


def plot_projection_sensitivity(bundle: ThesisBundle, validation_dir: Path | None, extended_dir: Path | None = None) -> None:
    density_path = _extended_path(validation_dir, extended_dir, "projection_sensitivity_density.tsv")
    summary_path = _extended_path(validation_dir, extended_dir, "projection_sensitivity_summary.tsv")
    df = safe_read(density_path)
    if df.empty:
        bundle.skip_named(
            "fig5_23_projection_sensitivity_center_vs_interval",
            "Projection sensitivity",
            [p for p in [density_path, summary_path] if p],
            "projection sensitivity tables missing; provide annotation + DMR rows or run extended profile/projection_sensitivity first",
        )
        return
    required = {"bin_index", "density_center", "density_interval_overlap"}
    if not required.issubset(df.columns):
        bundle.skip_named("fig5_23_projection_sensitivity_center_vs_interval", "Projection sensitivity", [density_path], "required columns missing")
        return
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(9, 4.5))
    df = df.sort_values("bin_index")
    ax.plot(df["bin_index"], df["density_center"], color=COLORS["blue"], linewidth=2, label="center")
    ax.plot(df["bin_index"], df["density_interval_overlap"], color=COLORS["orange"], linewidth=2, label="interval overlap")
    ax.set_xlabel("Metagene bin")
    ax.set_ylabel("Density")
    ax.set_title("Projection sensitivity: center vs interval-overlap")
    legend_outside(ax)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_23_projection_sensitivity_center_vs_interval", "Projection sensitivity", [p for p in [density_path, summary_path] if p])


def plot_annotation_gene(bundle: ThesisBundle, annotation_dir: Path | None) -> None:
    path = annotation_dir / "metagene_gene_regions_qc.tsv" if annotation_dir else None
    df = safe_read(path)
    if df.empty or not {"chrom", "length_bp"}.issubset(df.columns):
        bundle.skip("annotation_gene_overview", "Gene annotation overview", [path] if path else [], "gene QC table missing")
        return
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    lengths = pd.to_numeric(df["length_bp"], errors="coerce").dropna()
    axes[0].hist(np.log10(lengths + 1), bins=40, color=COLORS["blue"], edgecolor="white")
    axes[0].set_xlabel("log10(gene length bp + 1)")
    axes[0].set_ylabel("Genes")
    axes[0].set_title("Gene length distribution")
    style_ax(axes[0])
    chrom_counts = df["chrom"].astype(str).value_counts().sort_index()
    axes[1].bar(chrom_counts.index, chrom_counts.values, color=COLORS["green"])
    axes[1].set_ylabel("Genes")
    axes[1].set_title("Genes by seqname")
    axes[1].tick_params(axis="x", rotation=60)
    style_ax(axes[1])
    bundle.save_fig(fig, "annotation_gene_overview", "Gene annotation overview", [path])


def plot_annotation_te(bundle: ThesisBundle, annotation_dir: Path | None) -> None:
    path = annotation_dir / "te_regions_qc.tsv" if annotation_dir else None
    df = safe_read(path)
    if df.empty or "length_bp" not in df.columns:
        bundle.skip("annotation_te_overview", "TE/repeat annotation overview", [path] if path else [], "TE QC table missing")
        return
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    lengths = pd.to_numeric(df["length_bp"], errors="coerce").dropna()
    axes[0].hist(np.log10(lengths + 1), bins=50, color=COLORS["orange"], edgecolor="white")
    axes[0].set_xlabel("log10(TE/repeat length bp + 1)")
    axes[0].set_ylabel("Regions")
    axes[0].set_title("TE/repeat length distribution")
    style_ax(axes[0])
    class_col = "te_class" if "te_class" in df.columns else "feature_type" if "feature_type" in df.columns else None
    if class_col:
        counts = df[class_col].astype(str).value_counts().head(12).sort_values()
        axes[1].barh(counts.index, counts.values, color=COLORS["purple"])
        axes[1].set_xlabel("Regions")
        axes[1].set_title("Top TE/repeat classes")
        style_ax(axes[1])
    else:
        axes[1].set_axis_off()
    bundle.save_fig(fig, "annotation_te_overview", "TE/repeat annotation overview", [path])


def plot_go_coverage(bundle: ThesisBundle, annotation_dir: Path | None) -> None:
    gene_path = annotation_dir / "metagene_gene_regions_qc.tsv" if annotation_dir else None
    go_path = annotation_dir / "go_gene_summary.tsv" if annotation_dir else None
    genes = safe_read(gene_path)
    go = safe_read(go_path)
    if genes.empty or go.empty or "gene_id" not in genes.columns or "gene_id" not in go.columns:
        bundle.skip("go_mapping_coverage", "GO mapping coverage", [p for p in [gene_path, go_path] if p], "gene or GO table missing")
        return
    mapped = int(genes["gene_id"].astype(str).isin(set(go["gene_id"].astype(str))).sum())
    unmapped = int(len(genes) - mapped)
    plt = setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    axes[0].bar(["GO mapped", "No GO"], [mapped, unmapped], color=[COLORS["green"], COLORS["gray"]])
    axes[0].set_ylabel("Genes")
    axes[0].set_title("GO mapping coverage")
    style_ax(axes[0])
    if "n_go_terms" in go.columns:
        vals = pd.to_numeric(go["n_go_terms"], errors="coerce").dropna()
        axes[1].hist(vals, bins=30, color=COLORS["blue"], edgecolor="white")
        axes[1].set_xlabel("GO terms per mapped gene")
        axes[1].set_ylabel("Genes")
        axes[1].set_title("GO term count distribution")
        style_ax(axes[1])
    else:
        axes[1].set_axis_off()
    notes = f"mapped_genes={mapped}; unmapped_genes={unmapped}"
    bundle.save_fig(fig, "go_mapping_coverage", "GO mapping coverage", [gene_path, go_path], notes=notes)


def plot_seqname_compatibility(bundle: ThesisBundle, annotation_dir: Path | None, report_dir: Path | None) -> None:
    seq_path = annotation_dir / "seqname_mapping.tsv" if annotation_dir else None
    dmr_path = report_dir / "caller_summary" / "caller_dmr_rows_for_plots.tsv" if report_dir else None
    seq = safe_read(seq_path)
    dmrs = safe_read(dmr_path)
    if seq.empty or dmrs.empty or "chrom" not in dmrs.columns:
        bundle.skip("seqname_compatibility", "Annotation seqname compatibility", [p for p in [seq_path, dmr_path] if p], "seqname or DMR row table missing")
        return
    ann_col = "normalized_seqname" if "normalized_seqname" in seq.columns else "raw_seqname"
    annotation_seqnames = set(seq[ann_col].astype(str))
    work = dmrs.copy()
    caller_col = "caller" if "caller" in work.columns else "_caller" if "_caller" in work.columns else None
    context_col = "context" if "context" in work.columns else "_context" if "_context" in work.columns else None
    if caller_col is None:
        work["caller"] = "DMR"
        caller_col = "caller"
    if context_col is None:
        work["context"] = "NA"
        context_col = "context"
    work["annotation_seqname_status"] = np.where(work["chrom"].astype(str).isin(annotation_seqnames), "matched", "unmatched")
    summary = work.groupby([caller_col, context_col, "annotation_seqname_status"], dropna=False).size().reset_index(name="n")
    pivot = summary.pivot_table(index=[caller_col, context_col], columns="annotation_seqname_status", values="n", fill_value=0).reset_index()
    if "matched" not in pivot.columns:
        pivot["matched"] = 0
    if "unmatched" not in pivot.columns:
        pivot["unmatched"] = 0
    pivot["label"] = pivot[caller_col].astype(str) + " " + pivot[context_col].astype(str)
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(max(9, len(pivot) * 0.55), 4.8))
    x = np.arange(len(pivot))
    ax.bar(x, pivot["matched"], color=COLORS["green"], label="matched")
    ax.bar(x, pivot["unmatched"], bottom=pivot["matched"], color=COLORS["red"], label="unmatched")
    ax.set_xticks(x)
    ax.set_xticklabels(pivot["label"], rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("DMR rows")
    ax.set_title("DMR seqname compatibility with annotation")
    legend_outside(ax)
    style_ax(ax)
    bundle.save_fig(fig, "seqname_compatibility", "Annotation seqname compatibility", [seq_path, dmr_path])


def _interval_overlap_mask(query: pd.DataFrame, intervals: pd.DataFrame) -> pd.Series:
    if query.empty or intervals.empty:
        return pd.Series(False, index=query.index)
    result = pd.Series(False, index=query.index)
    for chrom, qsub in query.groupby("chrom", sort=False):
        isub = intervals[intervals["chrom"].astype(str).eq(str(chrom))].copy()
        if isub.empty:
            continue
        isub["start"] = pd.to_numeric(isub["start"], errors="coerce")
        isub["end"] = pd.to_numeric(isub["end"], errors="coerce")
        isub = isub.dropna(subset=["start", "end"]).sort_values("start")
        if isub.empty:
            continue
        starts = isub["start"].to_numpy(dtype=float)
        ends = isub["end"].to_numpy(dtype=float)
        cummax_ends = np.maximum.accumulate(ends)
        qstarts = pd.to_numeric(qsub["start"], errors="coerce").to_numpy(dtype=float)
        qends = pd.to_numeric(qsub["end"], errors="coerce").to_numpy(dtype=float)
        idx = np.searchsorted(starts, qends, side="right") - 1
        ok = (idx >= 0) & (cummax_ends[np.clip(idx, 0, len(cummax_ends) - 1)] >= qstarts)
        result.loc[qsub.index] = ok
    return result


def _assign_gene_links(dmrs: pd.DataFrame, genes: pd.DataFrame, promoter_window: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if dmrs.empty or genes.empty:
        return pd.DataFrame()
    gene_groups: dict[str, dict[str, object]] = {}
    for chrom, sub in genes.groupby("chrom", sort=False):
        work = sub.copy()
        work["start"] = pd.to_numeric(work["start"], errors="coerce")
        work["end"] = pd.to_numeric(work["end"], errors="coerce")
        work = work.dropna(subset=["start", "end"]).sort_values("start").reset_index(drop=True)
        starts = work["start"].to_numpy(dtype=int)
        ends = work["end"].to_numpy(dtype=int)
        cummax_ends = np.maximum.accumulate(ends) if len(ends) else np.array([], dtype=int)
        cummax_idx = np.zeros(len(ends), dtype=int)
        if len(ends):
            best_idx = 0
            best_end = ends[0]
            for i, value in enumerate(ends):
                if value >= best_end:
                    best_idx = i
                    best_end = value
                cummax_idx[i] = best_idx
        gene_groups[str(chrom)] = {
            "data": work,
            "starts": starts,
            "ends": ends,
            "cummax_ends": cummax_ends,
            "cummax_idx": cummax_idx,
        }
    for rec in dmrs.itertuples(index=False):
        chrom = str(getattr(rec, "chrom"))
        start = int(getattr(rec, "start"))
        end = int(getattr(rec, "end"))
        center = (start + end) // 2
        genes_chrom = gene_groups.get(chrom)
        gene_id = ""
        feature_class = "intergenic"
        distance = pd.NA
        if genes_chrom is not None:
            genes_data = genes_chrom["data"]
            starts = genes_chrom["starts"]
            ends = genes_chrom["ends"]
            cummax_ends = genes_chrom["cummax_ends"]
            cummax_idx = genes_chrom["cummax_idx"]
            idx_end = np.searchsorted(starts, end, side="right")
            prefix_idx = idx_end - 1
            if prefix_idx >= 0 and len(cummax_ends) and cummax_ends[prefix_idx] >= start:
                idx = int(cummax_idx[prefix_idx])
                gene_id = str(genes_data.at[idx, "gene_id"])
                feature_class = "gene_body"
                distance = 0
            else:
                window_start = max(0, center - promoter_window)
                window_end = center + promoter_window
                idx0 = np.searchsorted(starts, window_start, side="left")
                idx1 = np.searchsorted(starts, window_end, side="right")
                candidates = genes_data.iloc[idx0:idx1].copy()
                if not candidates.empty:
                    candidates["center_distance"] = np.minimum(
                        (center - candidates["start"]).abs(),
                        (center - candidates["end"]).abs(),
                    )
                    best = candidates.sort_values("center_distance").iloc[0]
                    gene_id = str(best["gene_id"])
                    distance = int(best["center_distance"])
                    strand = str(best.get("strand", "."))
                    gstart = int(best["start"])
                    gend = int(best["end"])
                    if strand == "-":
                        feature_class = "promoter" if center >= gend else "downstream"
                    else:
                        feature_class = "promoter" if center <= gstart else "downstream"
        rows.append(
            {
                "region_id": getattr(rec, "region_id", ""),
                "caller": getattr(rec, "caller", getattr(rec, "_caller", "")),
                "chrom": chrom,
                "start": start,
                "end": end,
                "context": getattr(rec, "context", getattr(rec, "_context", "")),
                "delta": getattr(rec, "delta", pd.NA),
                "q_value": getattr(rec, "q_value", pd.NA),
                "gene_id": gene_id,
                "feature_class": feature_class,
                "distance_to_gene": distance,
            }
        )
    return pd.DataFrame(rows)


def build_dmr_gene_summary(report_dir: Path | None, annotation_dir: Path | None, promoter_window: int) -> tuple[pd.DataFrame, pd.DataFrame, list[Path], str]:
    dmr_path = report_dir / "caller_summary" / "caller_dmr_rows_for_plots.tsv" if report_dir else None
    gene_path = annotation_dir / "metagene_gene_regions_qc.tsv" if annotation_dir else None
    te_path = annotation_dir / "te_regions_qc.tsv" if annotation_dir else None
    dmrs = safe_read(dmr_path)
    genes = safe_read(gene_path)
    tes = safe_read(te_path)
    sources = [p for p in [dmr_path, gene_path, te_path] if p]
    if dmrs.empty or genes.empty or not {"chrom", "start", "end"}.issubset(dmrs.columns) or not {"chrom", "start", "end", "gene_id"}.issubset(genes.columns):
        return pd.DataFrame(), pd.DataFrame(), sources, "DMR or gene annotation table missing"
    dmrs = dmrs.copy()
    dmrs["chrom"] = dmrs["chrom"].astype(str)
    dmrs["start"] = pd.to_numeric(dmrs["start"], errors="coerce")
    dmrs["end"] = pd.to_numeric(dmrs["end"], errors="coerce")
    dmrs = dmrs.dropna(subset=["start", "end"]).copy()
    dmrs["start"] = dmrs["start"].astype(int)
    dmrs["end"] = dmrs["end"].astype(int)
    if "region_id" not in dmrs.columns:
        dmrs["region_id"] = [f"dmr_{i + 1}" for i in range(len(dmrs))]
    if "caller" not in dmrs.columns and "_caller" in dmrs.columns:
        dmrs["caller"] = dmrs["_caller"]
    if "context" not in dmrs.columns and "_context" in dmrs.columns:
        dmrs["context"] = dmrs["_context"]
    linked = _assign_gene_links(dmrs, genes, promoter_window)
    if not tes.empty and {"chrom", "start", "end"}.issubset(tes.columns):
        linked["te_overlap"] = _interval_overlap_mask(linked[["chrom", "start", "end"]], tes[["chrom", "start", "end"]])
    else:
        linked["te_overlap"] = False
    if linked.empty:
        return linked, pd.DataFrame(), sources, "No DMR-gene links could be built"
    linked["linked_locus"] = np.where(linked["gene_id"].astype(str).ne(""), linked["gene_id"].astype(str), "intergenic")
    summary = (
        linked.groupby("linked_locus", dropna=False)
        .agg(
            n_dmr=("region_id", "count"),
            n_te_overlap=("te_overlap", "sum"),
            n_cg=("context", lambda x: int((x.astype(str).str.upper() == "CG").sum())),
            n_chg=("context", lambda x: int((x.astype(str).str.upper() == "CHG").sum())),
            n_chh=("context", lambda x: int((x.astype(str).str.upper() == "CHH").sum())),
        )
        .reset_index()
        .sort_values(["n_dmr", "linked_locus"], ascending=[False, True])
    )
    summary["te_like_status"] = np.where(summary["n_te_overlap"] > 0, "TE-overlap", "no TE overlap")
    return linked, summary, sources, ""


def _group_locus_sets(summary: pd.DataFrame) -> dict[str, pd.DataFrame]:
    no_intergenic = summary[summary["linked_locus"].astype(str).ne("intergenic")].copy()
    if no_intergenic.empty:
        no_intergenic = summary.copy()
    return {
        "top20": no_intergenic.head(20),
        "top50": no_intergenic.head(50),
        "all linked": no_intergenic,
    }


def plot_top_loci_te_status(bundle: ThesisBundle, report_dir: Path | None, annotation_dir: Path | None, promoter_window: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    linked, summary, sources, reason = build_dmr_gene_summary(report_dir, annotation_dir, promoter_window)
    if summary.empty:
        bundle.skip_named("fig5_24_top20_dmr_loci_te_status", "Top20 DMR-linked loci TE status", sources, reason)
        return linked, summary
    top = summary[summary["linked_locus"].astype(str).ne("intergenic")].head(20).sort_values("n_dmr")
    if top.empty:
        bundle.skip_named("fig5_24_top20_dmr_loci_te_status", "Top20 DMR-linked loci TE status", sources, "no non-intergenic linked loci")
        return linked, summary
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8, max(5, len(top) * 0.28)))
    colors = np.where(top["te_like_status"].eq("TE-overlap"), COLORS["red"], COLORS["blue"])
    ax.barh(top["linked_locus"], top["n_dmr"], color=colors)
    ax.set_xlabel("DMR links")
    ax.set_title("Top20 DMR-linked loci by DMR count and TE-overlap status")
    style_ax(ax)
    bundle.save_named_fig(
        fig,
        "fig5_24_top20_dmr_loci_te_status",
        "Top20 DMR-linked loci TE status",
        sources,
        notes="Generic status is based on DMR interval overlap with TE/repeat annotation, not organism-specific gene description text.",
    )
    return linked, summary


def plot_te_like_composition(bundle: ThesisBundle, summary: pd.DataFrame, sources: list[Path]) -> None:
    if summary.empty:
        bundle.skip_named("fig5_25_top50_te_like_composition", "TE-like composition", sources, "DMR-locus summary missing")
        return
    groups = _group_locus_sets(summary)
    rows = []
    for group, data in groups.items():
        counts = data["te_like_status"].value_counts().to_dict()
        total = max(1, len(data))
        for status in ["TE-overlap", "no TE overlap"]:
            rows.append({"group": group, "status": status, "fraction": counts.get(status, 0) / total, "n": counts.get(status, 0), "total": total})
    comp = pd.DataFrame(rows)
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bottom = np.zeros(len(groups))
    labels = list(groups)
    for status, color in [("TE-overlap", COLORS["red"]), ("no TE overlap", COLORS["blue"])]:
        vals = np.array([comp[(comp["group"].eq(label)) & (comp["status"].eq(status))]["fraction"].sum() for label in labels]) * 100
        ax.bar(labels, vals, bottom=bottom, color=color, label=status)
        bottom += vals
    ax.set_ylabel("Percent of loci")
    ax.set_title("TE-overlap composition in top loci vs background")
    legend_outside(ax)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_25_top50_te_like_composition", "TE-like composition", sources)


def plot_te_enrichment(bundle: ThesisBundle, summary: pd.DataFrame, sources: list[Path]) -> None:
    if summary.empty:
        bundle.skip_named("fig5_26_te_enrichment_summary", "TE enrichment summary", sources, "DMR-locus summary missing")
        return
    try:
        from scipy.stats import fisher_exact
    except Exception:
        fisher_exact = None
    data_all = summary[summary["linked_locus"].astype(str).ne("intergenic")].copy()
    if data_all.empty:
        bundle.skip_named("fig5_26_te_enrichment_summary", "TE enrichment summary", sources, "no non-intergenic linked loci")
        return
    rows = []
    for label, subset in {"top20": data_all.head(20), "top50": data_all.head(50)}.items():
        a = int(subset["te_like_status"].eq("TE-overlap").sum())
        b = int(len(subset) - a)
        rest = data_all.drop(subset.index)
        c = int(rest["te_like_status"].eq("TE-overlap").sum())
        d = int(len(rest) - c)
        odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        p_value = fisher_exact([[a, b], [c, d]], alternative="greater").pvalue if fisher_exact and len(rest) else np.nan
        rows.append({"test": label, "odds_ratio": odds, "p_value": p_value, "te_overlap": a, "non_te": b})
    tests = pd.DataFrame(rows)
    if tests.empty:
        bundle.skip_named("fig5_26_te_enrichment_summary", "TE enrichment summary", sources, "no enrichment tests available")
        return
    tests["minus_log10_p"] = -np.log10(pd.to_numeric(tests["p_value"], errors="coerce").replace(0, np.nextafter(0, 1)))
    tests["minus_log10_p"] = tests["minus_log10_p"].fillna(0)
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.barh(tests["test"], tests["minus_log10_p"], color=COLORS["red"])
    ax.set_xlabel("-log10(Fisher p)")
    ax.set_title("TE-overlap enrichment in top DMR-linked loci")
    for y, val, odds in zip(tests["test"], tests["minus_log10_p"], tests["odds_ratio"]):
        ax.text(val + 0.03, y, f"OR={odds:.2g}", va="center", fontsize=8)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_26_te_enrichment_summary", "TE enrichment summary", sources, notes="Generic Fisher test over TE-overlap status of DMR-linked loci.")


def plot_context_distribution(bundle: ThesisBundle, summary: pd.DataFrame, sources: list[Path]) -> None:
    if summary.empty:
        bundle.skip_named("fig5_27_top50_context_distribution", "Context distribution", sources, "DMR-locus summary missing")
        return
    groups = _group_locus_sets(summary)
    labels = list(groups)
    contexts = [("CG", "n_cg", COLORS["blue"]), ("CHG", "n_chg", COLORS["orange"]), ("CHH", "n_chh", COLORS["green"])]
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bottom = np.zeros(len(labels))
    for name, col, color in contexts:
        vals = np.array([pd.to_numeric(groups[label][col], errors="coerce").fillna(0).sum() if col in groups[label].columns else 0 for label in labels])
        ax.bar(labels, vals, bottom=bottom, color=color, label=name)
        bottom += vals
    ax.set_ylabel("DMR links")
    ax.set_title("DMR context distribution in linked loci")
    legend_outside(ax)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_27_top50_context_distribution", "Context distribution", sources)


def plot_feature_class_distribution(bundle: ThesisBundle, linked: pd.DataFrame, sources: list[Path]) -> None:
    if linked.empty or "feature_class" not in linked.columns:
        bundle.skip_named("fig5_28_feature_class_distribution", "Feature class distribution", sources, "DMR-gene links missing")
        return
    counts = linked["feature_class"].astype(str).value_counts().reindex(["promoter", "gene_body", "downstream", "intergenic"]).fillna(0)
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(counts.index, counts.values, color=[COLORS["orange"], COLORS["blue"], COLORS["green"], COLORS["gray"]])
    ax.set_ylabel("DMR links")
    ax.set_title("DMR-gene link feature class distribution")
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_28_feature_class_distribution", "Feature class distribution", sources)


def _canonical_gene_id(value: object) -> str:
    text = str(value).strip()
    match = re.search(r"(Solyc\d{2}g\d{6})(?:\.\d+)?", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).lower()
    return text.lower()


def plot_expression_support(bundle: ThesisBundle, summary: pd.DataFrame, expression_table: Path | None) -> None:
    sources = [expression_table] if expression_table else []
    if summary.empty:
        bundle.skip_named(
            "fig5_29_expression_support_levels",
            "Expression support levels",
            sources,
            "DMR-locus summary missing; cannot audit expression-evidence coverage",
        )
        return
    # Render this as an availability audit: even when no expression evidence maps to the
    # DMR-linked loci the figure still documents that limitation (0% coverage) rather than
    # being silently skipped. Expression genes are matched on the canonical Solyc locus.
    expr = safe_read(expression_table)
    gene_col = next((col for col in expr.columns if str(col).lower() in {"gene_id", "gene", "id", "locus"}), None) if not expr.empty else None
    expression_genes = {_canonical_gene_id(value) for value in expr[gene_col].dropna()} if gene_col else set()
    groups = _group_locus_sets(summary)
    labels = list(groups)
    totals = [len(groups[label]) for label in labels]
    available = [int(groups[label]["linked_locus"].map(_canonical_gene_id).isin(expression_genes).sum()) for label in labels]
    available_pct = [100.0 * available[i] / totals[i] if totals[i] else 0.0 for i in range(len(labels))]
    unavailable_pct = [100.0 - available_pct[i] for i in range(len(labels))]
    xticklabels = [f"{label}\n(n={totals[i]})" for i, label in enumerate(labels)]
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bars_av = ax.bar(labels, available_pct, color=COLORS["green"], label="expression evidence")
    ax.bar(labels, unavailable_pct, bottom=available_pct, color=COLORS["gray"], label="no evidence")
    for i, rect in enumerate(bars_av):
        ax.text(
            rect.get_x() + rect.get_width() / 2,
            min(available_pct[i] + 2.5, 96),
            f"{available[i]}/{totals[i]}\n{available_pct[i]:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(xticklabels)
    ax.set_ylim(0, 100)
    ax.set_ylabel("DMR-linked loci (%)")
    title = "Expression evidence coverage for DMR-linked loci"
    if not expression_genes:
        title += "\n(no Solyc-mappable expression source available)"
    ax.set_title(title)
    legend_outside(ax)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_29_expression_support_levels", "Expression support levels", sources)


def plot_dmr_linked_go_coverage(bundle: ThesisBundle, summary: pd.DataFrame, annotation_dir: Path | None) -> None:
    go_path = annotation_dir / "go_gene_summary.tsv" if annotation_dir else None
    go = safe_read(go_path)
    sources = [go_path] if go_path else []
    if summary.empty or go.empty or "gene_id" not in go.columns:
        bundle.skip_named("fig5_30_go_mapping_coverage_limitation", "DMR-linked GO mapping coverage", sources, "DMR-locus summary or GO table missing")
        return
    go_genes = {_canonical_gene_id(value) for value in go["gene_id"].dropna()}
    groups = _group_locus_sets(summary)
    labels = list(groups)
    totals = [len(groups[label]) for label in labels]
    mapped = [int(groups[label]["linked_locus"].map(_canonical_gene_id).isin(go_genes).sum()) for label in labels]
    # Group totals differ by orders of magnitude (top20 vs all-linked), so a shared
    # absolute axis hides the small groups. Plot coverage as a percentage instead and
    # annotate each bar with the absolute counts.
    mapped_pct = [100.0 * mapped[i] / totals[i] if totals[i] else 0.0 for i in range(len(labels))]
    unmapped_pct = [100.0 - mapped_pct[i] for i in range(len(labels))]
    xticklabels = [f"{label}\n(n={totals[i]})" for i, label in enumerate(labels)]
    plt = setup_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bars_mapped = ax.bar(labels, mapped_pct, color=COLORS["green"], label="GO mapped")
    ax.bar(labels, unmapped_pct, bottom=mapped_pct, color=COLORS["gray"], label="unmapped")
    for i, rect in enumerate(bars_mapped):
        ax.text(
            rect.get_x() + rect.get_width() / 2,
            min(mapped_pct[i] + 2.5, 96),
            f"{mapped[i]}/{totals[i]}\n{mapped_pct[i]:.0f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(xticklabels)
    ax.set_ylim(0, 100)
    ax.set_ylabel("DMR-linked loci (%)")
    ax.set_title("GO mapping coverage for DMR-linked loci")
    legend_outside(ax)
    style_ax(ax)
    bundle.save_named_fig(fig, "fig5_30_go_mapping_coverage_limitation", "DMR-linked GO mapping coverage", sources)


def copy_key_tables(bundle: ThesisBundle, inputs: dict[str, Path | None], max_rows: int) -> None:
    report_dir = inputs["report_dir"]
    validation_dir = inputs["validation_dir"]
    glm_glmm_dir = inputs["glm_glmm_dir"]
    annotation_dir = inputs["annotation_dir"]
    expression_table = inputs["expression_table"]
    table_specs = [
        (report_dir / "caller_summary" / "caller_context_summary.tsv" if report_dir else None, "caller_context_summary.tsv"),
        (validation_dir / "run_summary.tsv" if validation_dir else None, "validation_run_summary.tsv"),
        (validation_dir / "delta_weighting_summary.tsv" if validation_dir else None, "delta_weighting_summary.tsv"),
        (validation_dir / "delta_bootstrap_summary.tsv" if validation_dir else None, "delta_bootstrap_summary.tsv"),
        (validation_dir / "coverage_set_summary.tsv" if validation_dir else None, "coverage_set_summary.tsv"),
        (validation_dir / "glm_glmm_status_summary.tsv" if validation_dir else None, "glm_glmm_status_summary.tsv"),
        (glm_glmm_dir / "glm_vs_glmm" / "glm_vs_glmm_comparison.tsv" if glm_glmm_dir else None, "glm_vs_glmm_comparison.tsv"),
        (glm_glmm_dir / "region_cpg_counts_region_summary.tsv" if glm_glmm_dir else None, "region_cpg_counts_region_summary.tsv"),
        (annotation_dir / "annotation_manifest.json" if annotation_dir else None, "annotation_manifest.json"),
        (annotation_dir / "gff3_feature_summary.tsv" if annotation_dir else None, "annotation_gff3_feature_summary.tsv"),
        (annotation_dir / "te_class_summary.tsv" if annotation_dir else None, "annotation_te_class_summary.tsv"),
        (expression_table, "expression_support.tsv"),
    ]
    for source, name in table_specs:
        bundle.copy_table(source, name, max_rows)


def run(args: argparse.Namespace) -> int:
    inputs = resolve_inputs(args)
    out_dir = inputs["out_dir"]
    assert out_dir is not None
    bundle = ThesisBundle(out_dir, figure_prefix=args.figure_prefix, start_index=args.start_index)

    plot_candidate_funnel(bundle, inputs["report_dir"], inputs["glm_glmm_dir"])
    plot_caller_concordance(bundle, inputs["report_dir"])
    plot_volcano(bundle, inputs["report_dir"])
    plot_manhattan(bundle, inputs["report_dir"])
    add_existing_report_figures(bundle, inputs["report_dir"], inputs["glm_glmm_dir"])
    plot_overlap_threshold(bundle, inputs["validation_dir"])
    plot_coverage_qc(bundle, inputs["validation_dir"])
    plot_region_cpg_coverage_distribution(bundle, inputs["validation_dir"], inputs["glm_glmm_dir"])
    plot_delta_weighting(bundle, inputs["validation_dir"])
    plot_bootstrap_ci(bundle, inputs["validation_dir"], args.max_forest_rows)
    plot_delta_vs_coverage(bundle, inputs["validation_dir"])
    extended_metagene_dir = None
    if not args.skip_extended_metagene:
        extended_metagene_dir = build_extended_metagene_outputs(
            inputs["report_dir"],
            inputs["annotation_dir"],
            out_dir,
            args,
        )
    plot_random_occupancy(bundle, inputs["validation_dir"], extended_metagene_dir)
    plot_projection_sensitivity(bundle, inputs["validation_dir"], extended_metagene_dir)
    linked, locus_summary = plot_top_loci_te_status(bundle, inputs["report_dir"], inputs["annotation_dir"], args.promoter_window)
    if not linked.empty:
        linked.to_csv(out_dir / "tables" / "generic_dmr_gene_links.tsv", sep="\t", index=False)
    if not locus_summary.empty:
        locus_summary.to_csv(out_dir / "tables" / "generic_dmr_locus_summary.tsv", sep="\t", index=False)
    dmr_link_sources = [
        inputs["report_dir"] / "caller_summary" / "caller_dmr_rows_for_plots.tsv" if inputs["report_dir"] else None,
        inputs["annotation_dir"] / "metagene_gene_regions_qc.tsv" if inputs["annotation_dir"] else None,
        inputs["annotation_dir"] / "te_regions_qc.tsv" if inputs["annotation_dir"] else None,
    ]
    dmr_link_sources = [p for p in dmr_link_sources if p]
    plot_te_like_composition(bundle, locus_summary, dmr_link_sources)
    plot_te_enrichment(bundle, locus_summary, dmr_link_sources)
    plot_context_distribution(bundle, locus_summary, dmr_link_sources)
    plot_feature_class_distribution(bundle, linked, dmr_link_sources)
    plot_expression_support(bundle, locus_summary, inputs["expression_table"])
    plot_dmr_linked_go_coverage(bundle, locus_summary, inputs["annotation_dir"])
    plot_annotation_gene(bundle, inputs["annotation_dir"])
    plot_annotation_te(bundle, inputs["annotation_dir"])
    plot_go_coverage(bundle, inputs["annotation_dir"])
    plot_seqname_compatibility(bundle, inputs["annotation_dir"], inputs["report_dir"])
    copy_key_tables(bundle, inputs, args.max_table_copy_rows)

    manifest_json = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "inputs": {key: str(value) if value else "" for key, value in inputs.items()},
        "out_dir": str(out_dir),
        "figure_prefix": args.figure_prefix,
        "start_index": args.start_index,
        "extended_metagene_dir": str(extended_metagene_dir) if extended_metagene_dir else "",
    }
    (out_dir / "reports" / "thesis_figure_bundle_manifest.json").write_text(
        json.dumps(manifest_json, indent=2) + "\n",
        encoding="utf-8",
    )
    bundle.finalize(inputs)
    built = sum(1 for row in bundle.rows if row["status"] == "built")
    skipped = sum(1 for row in bundle.rows if row["status"] != "built")
    print(f"thesis figure bundle: {out_dir}")
    print(f"for_latex: {out_dir / 'for_latex'}")
    print(f"figures built: {built}; skipped: {skipped}")
    print(f"manifest: {out_dir / 'tables' / 'thesis_figure_manifest.tsv'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
