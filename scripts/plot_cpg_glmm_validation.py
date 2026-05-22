#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.viz.glmm_validation_plots import (
    plot_example_cpg_consistency,
    plot_glm_vs_glmm_scatter,
    plot_glmm_confirmation_by_context,
    plot_glmm_status_summary,
    plot_top_cpg_contribution,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create diagnostic plots for confirmatory CpG-level GLMM validation."
    )
    parser.add_argument("--comparison", required=True)
    parser.add_argument("--glmm-results", required=True)
    parser.add_argument("--region-cpg-counts", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--design")
    parser.add_argument("--condition-column", default="condition")
    parser.add_argument("--case-label")
    parser.add_argument("--control-label")
    parser.add_argument("--example-region-id")
    parser.add_argument("--region-delta", type=float)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    comparison = pd.read_csv(args.comparison, sep="\t")
    glmm = pd.read_csv(args.glmm_results, sep="\t")
    cpg_counts = pd.read_csv(args.region_cpg_counts, sep="\t")
    design = pd.read_csv(args.design, sep=None, engine="python") if args.design else None

    plot_glm_vs_glmm_scatter(comparison, out_dir / "glm_vs_glmm_scatter.png")
    plot_glmm_status_summary(glmm, out_dir / "glmm_status_summary.png")
    plot_glmm_confirmation_by_context(comparison, out_dir / "glmm_confirmation_by_context.png")
    plot_top_cpg_contribution(cpg_counts, out_dir / "top_cpg_contribution.png")
    region_id = args.example_region_id
    if region_id is None and not cpg_counts.empty:
        region_id = str(cpg_counts["region_id"].iloc[0])
    if region_id:
        plot_example_cpg_consistency(
            cpg_counts,
            region_id=region_id,
            out_path=out_dir / f"example_cpg_consistency_region_{region_id}.png",
            design_df=design,
            condition_column=args.condition_column,
            case_label=args.case_label,
            control_label=args.control_label,
            region_delta=args.region_delta,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
