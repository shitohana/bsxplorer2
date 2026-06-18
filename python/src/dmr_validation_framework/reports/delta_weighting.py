"""Delta-weighting report builders."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

def write_html_report(out_dir: Path, out: pd.DataFrame, summary: pd.DataFrame) -> None:
    display_cols = [
        "region_id",
        "delta_rep",
        "delta_pool",
        "abs_delta_shift",
        "relative_delta_shift",
        "loo_direction_changed",
        "ci_low",
        "ci_high",
        "coverage_imbalance_score",
        "model_agreement_status",
        "final_confidence_class",
    ]
    display_cols = [col for col in display_cols if col in out.columns]
    class_counts = out["final_confidence_class"].value_counts().reset_index()
    class_counts.columns = ["final_confidence_class", "n_regions"]
    order = {"LOW_CONFIDENCE": 0, "WEAK_CONFIDENCE": 1, "MODERATE_CONFIDENCE": 2, "HIGH_CONFIDENCE": 3}
    ranked = out.assign(_rank=out["final_confidence_class"].map(order).fillna(9))
    ranked = ranked.sort_values(["_rank", "abs_delta_shift"], ascending=[True, False]).head(50)
    body = "\n".join(
        [
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\"><title>Delta weighting sensitivity</title>",
            "<style>body{font-family:Arial,sans-serif;margin:24px;color:#222}table{border-collapse:collapse;margin:16px 0}th,td{border:1px solid #ddd;padding:6px 8px;font-size:13px}th{background:#f4f4f4}code{background:#f4f4f4;padding:2px 4px}</style>",
            "</head><body>",
            "<h1>Delta weighting sensitivity enhanced report</h1>",
            "<p>This report combines replicate-level vs pooled delta, leave-one-replicate-out, bootstrap CI, coverage imbalance, GLM/GLMM agreement, and confidence classes.</p>",
            "<h2>Summary</h2>",
            summary.to_html(index=False, escape=True),
            "<h2>Confidence classes</h2>",
            class_counts.to_html(index=False, escape=True),
            "<h2>Regions requiring most attention</h2>",
            ranked[display_cols].to_html(index=False, escape=True),
            "</body></html>",
        ]
    )
    (out_dir / "delta_weighting_summary.html").write_text(body + "\n", encoding="utf-8")
