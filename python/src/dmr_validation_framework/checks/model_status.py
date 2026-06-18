#!/usr/bin/env python3
"""Audit model status and convergence/fallback flags in GLM/GLMM tables."""

from __future__ import annotations

import argparse
from pathlib import Path

from dmr_validation_framework.core.io import classify_file_role, ensure_out_dir, find_candidate_files, read_table, write_tsv


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def classify_status(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"", "nan", "none"}:
        return "missing_status"
    if "fallback" in text or "quasi" in text:
        return "fallback"
    if any(token in text for token in ["fail", "error", "unavailable", "skipped", "insufficient"]):
        return "failed"
    if "converg" in text or "singular" in text:
        return "convergence_warning"
    if text.startswith("ok") or text in {"true", "1"}:
        return "ok"
    return "other"


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    rows: list[dict] = []
    for path in find_candidate_files():
        role = classify_file_role(path)
        if role not in {"glm_results", "glmm_results", "glm_vs_glmm_comparison"}:
            continue
        try:
            df = read_table(path)
        except Exception as exc:
            rows.append({"file": str(path), "role": role, "status_category": "read_error", "n_rows": 0, "notes": str(exc)})
            continue
        status_cols = [c for c in df.columns if "status" in c.lower() or c.lower() in {"converged", "singular", "warning"}]
        if not status_cols:
            rows.append(
                {
                    "file": str(path),
                    "role": role,
                    "status_column": "NA",
                    "status_value": "NA",
                    "status_category": "missing_status",
                    "n_rows": len(df),
                    "notes": "no model_status-like column found",
                }
            )
            continue
        for col in status_cols:
            counts = df[col].astype(str).fillna("").map(classify_status).value_counts()
            for category, count in counts.items():
                rows.append(
                    {
                        "file": str(path),
                        "role": role,
                        "status_column": col,
                        "status_value": category,
                        "status_category": category,
                        "n_rows": int(count),
                        "notes": "status categories are descriptive audit bins",
                    }
                )
    if not rows:
        rows = [{"status": "SKIPPED", "notes": "no GLM/GLMM result tables found"}]
    write_tsv(out_dir / "model_status_audit.tsv", rows)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
