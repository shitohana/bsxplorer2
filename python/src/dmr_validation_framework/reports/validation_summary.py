#!/usr/bin/env python
"""Build a compact HTML/TSV summary from DMR validation outputs."""

from __future__ import annotations

import argparse
import html
from pathlib import Path

import pandas as pd

from dmr_validation_framework.core.io import read_table
from dmr_validation_framework.reports.palette import (
    CONFIDENCE_COLORS,
    PALETTE,
    STATUS_COLORS,
    color_for,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", required=True, help="Directory produced by run-all/critical-validation.")
    parser.add_argument("--out-dir", required=True, help="Directory for summary tables and plots.")
    parser.add_argument("--no-plots", action="store_true")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def safe_read(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return read_table(path)
    except Exception:
        return pd.DataFrame()


def count_column(df: pd.DataFrame, column: str) -> pd.DataFrame:
    if df.empty or column not in df.columns:
        return pd.DataFrame(columns=[column, "n"])
    out = df[column].astype(str).fillna("NA").value_counts(dropna=False).reset_index()
    out.columns = [column, "n"]
    return out


def summarize_file(path: Path, preferred_status_col: str = "status") -> dict:
    df = safe_read(path)
    row = {
        "file": path.name,
        "path": str(path),
        "exists": path.exists(),
        "n_rows": len(df),
        "status_counts": "",
        "notes": "",
    }
    if df.empty:
        row["notes"] = "missing_or_empty"
        return row
    status_col = preferred_status_col if preferred_status_col in df.columns else None
    if status_col is None:
        for candidate in ("final_confidence_class", "model_agreement_status", "coverage_qc_status"):
            if candidate in df.columns:
                status_col = candidate
                break
    if status_col:
        counts = df[status_col].astype(str).value_counts(dropna=False).to_dict()
        row["status_counts"] = "; ".join(f"{key}={value}" for key, value in counts.items())
    return row


def write_barplot(df: pd.DataFrame, x: str, y: str, title: str, path: Path, color_map: dict[str, str] | None = None) -> bool:
    if df.empty or x not in df.columns or y not in df.columns:
        return False
    import matplotlib.pyplot as plt

    categories = df[x].astype(str)
    colors = (
        [color_for(value, color_map) for value in categories]
        if color_map
        else PALETTE["blue"]
    )
    fig, ax = plt.subplots(figsize=(max(7, len(df) * 0.8), 4.5))
    ax.bar(categories, pd.to_numeric(df[y], errors="coerce").fillna(0), color=colors)
    ax.set_title(title)
    ax.set_ylabel(y)
    ax.tick_params(axis="x", rotation=35)
    plt.tight_layout()
    plt.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return True


def build_summary_tables(validation_dir: Path, out_dir: Path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    files = [
        "run_summary.tsv",
        "delta_weighting_summary.tsv",
        "delta_weighting_sensitivity.tsv",
        "delta_bootstrap_summary.tsv",
        "coverage_set_summary.tsv",
        "glm_glmm_status_summary.tsv",
        "model_status_audit.tsv",
        "direction_agreement.tsv",
        "overlap_threshold_sensitivity.tsv",
        "random_control_occupancy_summary.tsv",
    ]
    inventory = pd.DataFrame([summarize_file(validation_dir / name) for name in files])
    inventory.to_csv(out_dir / "validation_output_inventory.tsv", sep="\t", index=False)

    tables: dict[str, pd.DataFrame] = {}
    run_summary = safe_read(validation_dir / "run_summary.tsv")
    if not run_summary.empty:
        tables["check_status"] = count_column(run_summary, "status")
        tables["check_status"].to_csv(out_dir / "validation_check_status_counts.tsv", sep="\t", index=False)

    delta = safe_read(validation_dir / "delta_weighting_sensitivity.tsv")
    if not delta.empty:
        if "final_confidence_class" in delta.columns:
            tables["confidence"] = count_column(delta, "final_confidence_class")
            tables["confidence"].to_csv(out_dir / "final_confidence_class_counts.tsv", sep="\t", index=False)
        if "model_agreement_status" in delta.columns:
            tables["model_agreement"] = count_column(delta, "model_agreement_status")
            tables["model_agreement"].to_csv(out_dir / "model_agreement_status_counts.tsv", sep="\t", index=False)

    glm = safe_read(validation_dir / "glm_glmm_status_summary.tsv")
    if not glm.empty and "model_agreement_status" in glm.columns:
        tables["glm_glmm"] = count_column(glm, "model_agreement_status")
        tables["glm_glmm"].to_csv(out_dir / "glm_glmm_status_counts.tsv", sep="\t", index=False)

    coverage = safe_read(validation_dir / "coverage_set_summary.tsv")
    if not coverage.empty and "status" in coverage.columns:
        tables["coverage_status"] = count_column(coverage, "status")
        tables["coverage_status"].to_csv(out_dir / "coverage_status_counts.tsv", sep="\t", index=False)

    return inventory, tables


def write_html(out_dir: Path, inventory: pd.DataFrame, tables: dict[str, pd.DataFrame], plot_files: list[str]) -> None:
    sections = []
    for name, table in tables.items():
        sections.append(f"<h2>{html.escape(name)}</h2>")
        sections.append(table.to_html(index=False, escape=True))
    image_tags = "\n".join(
        f'<h2>{html.escape(name)}</h2><img src="{html.escape(name)}" style="max-width:100%;height:auto">'
        for name in plot_files
    )
    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\"><title>DMR validation summary</title>",
            "<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse;margin:12px 0}td,th{border:1px solid #ddd;padding:5px 7px;font-size:12px}th{background:#f4f4f4}</style>",
            "</head><body>",
            "<h1>DMR validation summary</h1>",
            "<p>This report summarizes framework validation outputs. Missing checks are shown explicitly in the inventory.</p>",
            "<h2>Output inventory</h2>",
            inventory.to_html(index=False, escape=True),
            *sections,
            image_tags,
            "</body></html>",
        ]
    )
    (out_dir / "validation_summary.html").write_text(body + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> int:
    validation_dir = Path(args.validation_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    inventory, tables = build_summary_tables(validation_dir, out_dir)
    plot_files: list[str] = []
    if not args.no_plots:
        from dmr_validation_framework.reports.theme import apply_matplotlib_theme

        apply_matplotlib_theme()
        try:
            plot_specs = [
                ("check_status", "status", "n", "Validation check status", "validation_check_status.png", None),
                ("confidence", "final_confidence_class", "n", "Final confidence class", "final_confidence_class.png", CONFIDENCE_COLORS),
                ("model_agreement", "model_agreement_status", "n", "Model agreement status", "model_agreement_status.png", STATUS_COLORS),
                ("coverage_status", "status", "n", "Coverage QC status", "coverage_status.png", None),
            ]
            for key, x, y, title, filename, color_map in plot_specs:
                if write_barplot(tables.get(key, pd.DataFrame()), x, y, title, out_dir / filename, color_map):
                    plot_files.append(filename)
        except Exception as exc:  # noqa: BLE001
            (out_dir / "plot_warnings.txt").write_text(str(exc) + "\n", encoding="utf-8")
    write_html(out_dir, inventory, tables, plot_files)
    print(f"wrote validation summary: {out_dir / 'validation_summary.html'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
