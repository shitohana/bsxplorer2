#!/usr/bin/env python3
"""Summarize agreement classes between aggregated GLM and CpG-level GLMM."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from dmr_validation_framework.core.io import ensure_out_dir, find_first, read_table, write_tsv


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def status_not_ok(value: object) -> bool:
    text = str(value).lower()
    if text in {"", "nan", "none"}:
        return False
    return any(
        token in text
        for token in [
            "fail",
            "error",
            "unavailable",
            "skipped",
            "insufficient",
            "not_tested",
            "no_lrt",
            "not_converged",
            "singular",
        ]
    )


def first_col(df: pd.DataFrame, names: list[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    path = find_first("glm_vs_glmm_comparison")
    if not path:
        skipped = [{"status": "SKIPPED", "notes": "GLM vs GLMM comparison table not found"}]
        write_tsv(out_dir / "glm_glmm_status_summary.tsv", skipped)
        write_tsv(out_dir / "glm_glmm_region_classes.tsv", skipped)
        return 0
    df = read_table(path)
    glm_q_col = first_col(df, ["glm_q_value", "glm_q", "aggregated_glm_q"])
    glmm_q_col = first_col(df, ["glmm_q_value", "glmm_q", "q_value"])
    glm_status_col = first_col(df, ["glm_status", "model_status", "aggregated_glm_status"])
    glmm_status_col = first_col(df, ["glmm_status", "model_status"])
    region_col = first_col(df, ["region_id", "dmr_id", "harmonized_region_id"]) or df.columns[0]
    if not glm_q_col or not glmm_q_col:
        skipped = [{"status": "SKIPPED", "notes": "glm/glmm q-value columns not found", "input_file": str(path)}]
        write_tsv(out_dir / "glm_glmm_status_summary.tsv", skipped)
        write_tsv(out_dir / "glm_glmm_region_classes.tsv", skipped)
        return 0
    out_rows: list[dict] = []
    for row in df.to_dict(orient="records"):
        glm_q = pd.to_numeric(row.get(glm_q_col), errors="coerce")
        glmm_q = pd.to_numeric(row.get(glmm_q_col), errors="coerce")
        glm_status = row.get(glm_status_col, "") if glm_status_col else ""
        glmm_status = row.get(glmm_status_col, "") if glmm_status_col else ""
        glm_sig = pd.notna(glm_q) and glm_q < args.alpha
        glmm_sig = pd.notna(glmm_q) and glmm_q < args.alpha and not status_not_ok(glmm_status)
        failed = status_not_ok(glm_status) or status_not_ok(glmm_status)
        if failed:
            klass = "failed/fallback/model_status_not_ok"
        elif glm_sig and glmm_sig:
            klass = "GLMM-confirmed"
        elif glm_sig and not glmm_sig:
            klass = "GLM-only"
        elif glmm_sig and not glm_sig:
            klass = "GLMM-only"
        else:
            klass = "neither"
        out_rows.append(
            {
                "region_id": row.get(region_col, ""),
                "chrom": row.get("chrom", ""),
                "start": row.get("start", ""),
                "end": row.get("end", ""),
                "context": row.get("context", ""),
                "glm_status": glm_status,
                "glmm_status": glmm_status,
                "glm_q": glm_q,
                "glmm_q": glmm_q,
                "glm_delta": row.get("region_delta_callus_minus_seedling", row.get("glm_delta", "")),
                "glmm_delta": row.get("delta_methylation", row.get("glmm_delta", "")),
                "class": klass,
                "notes": "CpG-level GLMM is confirmatory for selected regions, not genome-wide calling",
            }
        )
    classes = pd.Series([r["class"] for r in out_rows]).value_counts().to_dict()
    summary = [
        {
            "class": klass,
            "n_regions": count,
            "input_file": str(path),
            "status": "PASS",
            "notes": "q-values from GLM and GLMM are compared descriptively, not treated as equivalent",
        }
        for klass, count in classes.items()
    ]
    write_tsv(out_dir / "glm_glmm_region_classes.tsv", out_rows)
    write_tsv(out_dir / "glm_glmm_status_summary.tsv", summary)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
