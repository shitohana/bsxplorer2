#!/usr/bin/env python
"""Summarize external DMR caller outputs and build simple diagnostic plots."""

from __future__ import annotations

import argparse
import html
from pathlib import Path

import pandas as pd

from dmr_validation_framework.core.io import read_table
from dmr_validation_framework.core.columns import normalize_region_table
from dmr_validation_framework.reports.palette import CONTEXT_COLORS, PALETTE, caller_color, color_for, legend_outside

DEFAULT_CONTEXTS = ("CG", "CHG", "CHH")
DEFAULT_CALLERS = ("DSS", "methylKit", "dmrseq", "BSmooth", "comb-p", "metilene")

CALLER_PATTERNS: dict[str, tuple[str, ...]] = {
    "DSS": ("dss/dss_dmrs_{context}.tsv",),
    "methylKit": ("methylkit/methylkit_dmrs_{context}.tsv",),
    "dmrseq": ("dmrseq/dmrseq_dmrs_{context}.tsv",),
    "BSmooth": ("bsmooth/bsmooth_dmrs_{context}.tsv", "bsseq/bsmooth_dmrs_{context}.tsv"),
    "comb-p": ("combp/combp_dmrs_{context}.tsv", "comb-p/combp_dmrs_{context}.tsv"),
    "metilene": ("metilene/metilene_dmrs_{context}.tsv",),
    "metilene_input": ("metilene_input/metilene_{context}.tsv",),
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-root", required=True, help="Directory with caller output folders.")
    parser.add_argument("--out-dir", required=True, help="Directory for summary tables and plots.")
    parser.add_argument("--contexts", default="CG,CHG,CHH")
    parser.add_argument("--callers", default="DSS,methylKit,dmrseq,BSmooth,comb-p,metilene")
    parser.add_argument("--status-table", help="Optional external_caller_run_status.tsv.")
    parser.add_argument("--max-plot-rows", type=int, default=200_000)
    parser.add_argument("--no-plots", action="store_true")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def first_existing(columns: pd.Index, candidates: tuple[str, ...]) -> str | None:
    lower = {str(col).lower(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return str(lower[candidate.lower()])
    return None


def discover_path(external_root: Path, caller: str, context: str) -> Path | None:
    patterns = CALLER_PATTERNS.get(caller, ())
    if caller == "metilene":
        patterns = (*patterns, *CALLER_PATTERNS["metilene_input"])
    for pattern in patterns:
        path = external_root / pattern.format(context=context)
        if path.exists():
            return path
    return None


def summarize_table(caller: str, context: str, path: Path | None) -> tuple[dict, pd.DataFrame]:
    base = {
        "caller": caller,
        "context": context,
        "path": str(path) if path else "",
        "exists": bool(path and path.exists()),
        "status": "missing",
        "n_rows": 0,
        "n_dmrs": 0,
        "n_q05": 0,
        "n_q10": 0,
        "n_hyper": 0,
        "n_hypo": 0,
        "median_length": pd.NA,
        "mean_abs_delta": pd.NA,
        "median_abs_delta": pd.NA,
        "notes": "",
    }
    if path is None or not path.exists():
        return base, pd.DataFrame()
    try:
        raw = read_table(path)
    except Exception as exc:  # noqa: BLE001
        base["status"] = "error"
        base["notes"] = str(exc)
        return base, pd.DataFrame()

    base["status"] = "ok"
    base["n_rows"] = len(raw)

    # Use the same canonical normalization layer as harmonization/validation.
    # This fixes caller-specific delta scales, e.g. BSmooth percent points -> fraction.
    canonical, norm_notes = normalize_region_table(raw, source=caller)
    df = canonical if not canonical.empty else raw.copy()

    if norm_notes:
        base["notes"] = "; ".join(norm_notes)

    is_dmr_table = {"chrom", "start", "end"}.issubset(set(df.columns))
    base["n_dmrs"] = len(df) if is_dmr_table and caller != "metilene_input" else 0

    q_col = first_existing(df.columns, ("q_value", "qvalue", "q", "fdr", "padj"))
    delta_col = first_existing(df.columns, ("delta", "meth.diff", "meth_diff", "diff"))

    if q_col:
        q = pd.to_numeric(df[q_col], errors="coerce")
        base["n_q05"] = int((q <= 0.05).sum())
        base["n_q10"] = int((q <= 0.10).sum())

    if "start" in df.columns and "end" in df.columns:
        length = pd.to_numeric(df["end"], errors="coerce") - pd.to_numeric(df["start"], errors="coerce")
        base["median_length"] = float(length.dropna().median()) if length.notna().any() else pd.NA

    if delta_col:
        delta = pd.to_numeric(df[delta_col], errors="coerce")
        base["n_hyper"] = int((delta > 0).sum())
        base["n_hypo"] = int((delta < 0).sum())
        abs_delta = delta.abs()
        base["mean_abs_delta"] = float(abs_delta.dropna().mean()) if abs_delta.notna().any() else pd.NA
        base["median_abs_delta"] = float(abs_delta.dropna().median()) if abs_delta.notna().any() else pd.NA
        df = df.assign(_caller=caller, _context=context, _abs_delta=abs_delta)
    else:
        df = df.assign(_caller=caller, _context=context)

    if "start" in df.columns and "end" in df.columns:
        df["_length"] = pd.to_numeric(df["end"], errors="coerce") - pd.to_numeric(df["start"], errors="coerce")

    return base, df


def load_status(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    try:
        return read_table(path)
    except Exception:
        return pd.DataFrame()


def merge_run_status(summary: pd.DataFrame, status: pd.DataFrame) -> pd.DataFrame:
    if status.empty or not {"caller", "context", "status"}.issubset(status.columns):
        summary["run_status"] = ""
        summary["run_message"] = ""
        return summary
    latest = status.groupby(["caller", "context"], dropna=False).tail(1)
    latest = latest.rename(columns={"status": "run_status", "message": "run_message"})
    return summary.merge(
        latest[["caller", "context", "run_status", "run_message"]],
        on=["caller", "context"],
        how="left",
    )


def write_plot_counts(summary: pd.DataFrame, out_dir: Path) -> list[str]:
    import matplotlib.pyplot as plt

    plot_files: list[str] = []
    if summary.empty:
        return plot_files
    pivot = summary.pivot_table(index="caller", columns="context", values="n_dmrs", aggfunc="sum", fill_value=0)
    colors = [color_for(ctx, CONTEXT_COLORS) for ctx in pivot.columns]
    ax = pivot.plot(kind="bar", figsize=(10, 5), color=colors)
    ax.set_ylabel("DMR count")
    ax.set_xlabel("Caller")
    ax.set_title("DMR count by caller and context")
    legend_outside(ax, title="Context")
    plt.tight_layout()
    path = out_dir / "dmr_count_by_caller_context.png"
    plt.savefig(path, dpi=160, bbox_inches="tight")
    plt.close()
    plot_files.append(path.name)
    return plot_files


def _faceted_box(
    plot_df: pd.DataFrame,
    value_col: str,
    ylabel: str,
    title: str,
    path: Path,
    *,
    log: bool = False,
    exclude_callers: set[str] | None = None,
) -> bool:
    """One subplot per context; one box per caller, colored and annotated with n."""
    import matplotlib.pyplot as plt

    df = plot_df.copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    df = df.dropna(subset=[value_col])
    if exclude_callers:
        df = df[~df["_caller"].astype(str).isin(exclude_callers)]
    if df.empty:
        return False
    contexts = sorted(df["_context"].astype(str).unique())
    fig, axes = plt.subplots(
        1,
        len(contexts),
        figsize=(max(6, len(contexts) * 4.2), 5),
        squeeze=False,
        sharey=True,
    )
    for ax, ctx in zip(axes[0], contexts):
        sub = df[df["_context"].astype(str) == ctx]
        callers = sorted(sub["_caller"].astype(str).unique())
        data = [sub.loc[sub["_caller"].astype(str) == c, value_col].to_numpy() for c in callers]
        bp = ax.boxplot(data, labels=callers, showfliers=False, patch_artist=True)
        for patch, caller in zip(bp["boxes"], callers):
            patch.set_facecolor(caller_color(caller))
            patch.set_alpha(0.75)
        for median in bp["medians"]:
            median.set_color(PALETTE["black"])
        if log:
            ax.set_yscale("log")
        ax.set_title(ctx)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=45)
        top = ax.get_ylim()[1]
        for i, values in enumerate(data, start=1):
            ax.text(i, top, f"n={len(values)}", ha="center", va="bottom", fontsize=7, color=PALETTE["gray"])
    fig.suptitle(title, y=1.03)
    fig.tight_layout()
    plt.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return True


def write_plot_boxplots(all_rows: pd.DataFrame, out_dir: Path, max_plot_rows: int) -> list[str]:
    plot_files: list[str] = []
    if all_rows.empty:
        return plot_files
    plot_df = all_rows.head(max_plot_rows).copy()

    if "_abs_delta" in plot_df.columns and plot_df["_abs_delta"].notna().any():
        if _faceted_box(
            plot_df,
            "_abs_delta",
            "|delta|",
            "Absolute effect size by context and caller",
            out_dir / "abs_delta_distribution_by_caller_context.png",
        ):
            plot_files.append("abs_delta_distribution_by_caller_context.png")
        positive = pd.to_numeric(plot_df["_abs_delta"], errors="coerce")
        positive = positive[positive > 0].dropna()
        if not positive.empty and positive.max() / max(float(positive.min()), 1e-12) > 100:
            if _faceted_box(
                plot_df,
                "_abs_delta",
                "|delta| (log scale)",
                "Absolute effect size by context and caller (log scale)",
                out_dir / "abs_delta_distribution_by_caller_context_logscale.png",
                log=True,
            ):
                plot_files.append("abs_delta_distribution_by_caller_context_logscale.png")

    if "_length" in plot_df.columns and plot_df["_length"].notna().any():
        # Exclude fixed-window callers (e.g. methylKit tiles): their length is
        # window-defined, not data-defined, so a near-zero-variance box would be
        # misleading next to data-defined callers.
        length_std = plot_df.groupby("_caller")["_length"].transform(lambda values: values.dropna().std(ddof=0))
        fixed_callers = set(
            plot_df.loc[length_std.notna() & (length_std <= 1.0), "_caller"].astype(str).unique()
        )
        title = "DMR length distribution by context (data-defined callers)"
        if fixed_callers:
            title += f"\nfixed-window callers excluded: {', '.join(sorted(fixed_callers))}"
        if _faceted_box(
            plot_df,
            "_length",
            "DMR length (bp)",
            title,
            out_dir / "dmr_length_distribution_by_caller_context.png",
            exclude_callers=fixed_callers,
        ):
            plot_files.append("dmr_length_distribution_by_caller_context.png")
    return plot_files


def write_html(summary: pd.DataFrame, out_dir: Path, plot_files: list[str]) -> None:
    image_tags = "\n".join(
        f'<h2>{html.escape(name)}</h2><img src="{html.escape(name)}" style="max-width:100%;height:auto">'
        for name in plot_files
    )
    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\"><title>DMR caller summary</title>",
            "<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:5px 7px;font-size:12px}th{background:#f4f4f4}</style>",
            "</head><body>",
            "<h1>DMR caller summary</h1>",
            "<p>This report summarizes caller-native DMR output tables. It is descriptive and does not rank caller quality by DMR count alone.</p>",
            summary.to_html(index=False, escape=True),
            image_tags,
            "</body></html>",
        ]
    )
    (out_dir / "caller_summary.html").write_text(body + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    external_root = Path(args.external_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    contexts = [ctx.upper() for ctx in split_csv(args.contexts)] or list(DEFAULT_CONTEXTS)
    callers = split_csv(args.callers) or list(DEFAULT_CALLERS)

    rows: list[dict] = []
    tables: list[pd.DataFrame] = []
    for caller in callers:
        for context in contexts:
            path = discover_path(external_root, caller, context)
            row, table = summarize_table(caller, context, path)
            rows.append(row)
            if not table.empty and (("_abs_delta" in table.columns) or ("_length" in table.columns)):
                tables.append(table)

    summary = pd.DataFrame(rows)
    status_path = Path(args.status_table) if args.status_table else external_root / "external_caller_run_status.tsv"
    summary = merge_run_status(summary, load_status(status_path))
    summary.to_csv(out_dir / "caller_context_summary.tsv", sep="\t", index=False)

    all_rows = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    if not all_rows.empty:
        all_rows.to_csv(out_dir / "caller_dmr_rows_for_plots.tsv", sep="\t", index=False)

    plot_files: list[str] = []
    if not args.no_plots:
        from dmr_validation_framework.reports.theme import apply_matplotlib_theme

        apply_matplotlib_theme()
        try:
            plot_files.extend(write_plot_counts(summary, out_dir))
            plot_files.extend(write_plot_boxplots(all_rows, out_dir, args.max_plot_rows))
        except Exception as exc:  # noqa: BLE001
            (out_dir / "plot_warnings.txt").write_text(str(exc) + "\n", encoding="utf-8")

    write_html(summary, out_dir, plot_files)
    print(f"wrote caller summary: {out_dir / 'caller_context_summary.tsv'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
