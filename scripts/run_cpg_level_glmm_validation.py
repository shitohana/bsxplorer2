#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis.cpg_level_glmm import (
    glmmTMB_available,
    read_design_table,
    read_dmr_evidence_for_glmm,
    read_region_cpg_counts,
    rscript_available,
    run_cpg_level_glmm_validation,
)


def sha256_file(path: Path) -> str | None:
    if not path.exists() or path.is_dir():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run optional confirmatory CpG-level GLMM validation for selected DMR candidates."
    )
    parser.add_argument("--region-cpg-counts", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--dmr-evidence")
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--condition-column", default="condition")
    parser.add_argument("--case-label")
    parser.add_argument("--control-label")
    parser.add_argument("--covariates", default="")
    parser.add_argument("--min-cpg", type=int, default=3)
    parser.add_argument("--min-replicates-per-group", type=int, default=2)
    parser.add_argument("--min-total", type=int, default=1)
    parser.add_argument("--max-zero-coverage-fraction", type=float, default=0.5)
    parser.add_argument("--qvalue-method", default="BH")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = out_dir / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    counts = read_region_cpg_counts(args.region_cpg_counts)
    design = read_design_table(args.design)
    evidence = read_dmr_evidence_for_glmm(args.dmr_evidence) if args.dmr_evidence else None
    covariates = [item.strip() for item in args.covariates.split(",") if item.strip()]

    results, warnings = run_cpg_level_glmm_validation(
        counts,
        design,
        evidence_df=evidence,
        top_n=args.top_n,
        condition_column=args.condition_column,
        case_label=args.case_label,
        control_label=args.control_label,
        covariates=covariates,
        min_cpg=args.min_cpg,
        min_replicates_per_group=args.min_replicates_per_group,
        min_total=args.min_total,
        max_zero_coverage_fraction=args.max_zero_coverage_fraction,
        qvalue_method=args.qvalue_method,
        temp_dir=temp_dir,
    )

    results_path = out_dir / "cpg_level_glmm_results.tsv"
    warnings_path = out_dir / "cpg_level_glmm_warnings.tsv"
    summary_path = out_dir / "cpg_level_glmm_summary.md"
    manifest_path = out_dir / "cpg_level_glmm_manifest.json"
    results.to_csv(results_path, sep="\t", index=False)
    warnings.to_csv(warnings_path, sep="\t", index=False)

    status_counts = results["model_status"].value_counts(dropna=False).to_dict() if not results.empty else {}
    summary_path.write_text(
        "\n".join([
            "# CpG-level GLMM Confirmatory Validation Summary",
            "",
            "This optional layer validates selected predefined DMR candidates; it is not a genome-wide DMR caller.",
            "",
            f"- top_n_requested: {args.top_n}",
            f"- n_regions_reported: {len(results)}",
            f"- Rscript_available: {rscript_available()}",
            f"- glmmTMB_available: {glmmTMB_available()}",
            f"- model_status_counts: {status_counts}",
            f"- warnings: {len(warnings)}",
            "",
        ]),
        encoding="utf-8",
    )
    output_files = [results_path, warnings_path, summary_path]
    manifest = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "region_cpg_counts": str(args.region_cpg_counts),
            "design": str(args.design),
            "dmr_evidence": str(args.dmr_evidence) if args.dmr_evidence else None,
        },
        "parameters": {
            "top_n": args.top_n,
            "condition_column": args.condition_column,
            "case_label": args.case_label,
            "control_label": args.control_label,
            "covariates": covariates,
            "min_cpg": args.min_cpg,
            "min_replicates_per_group": args.min_replicates_per_group,
            "min_total": args.min_total,
            "max_zero_coverage_fraction": args.max_zero_coverage_fraction,
            "qvalue_method": args.qvalue_method,
        },
        "outputs": {path.name: {"path": str(path), "sha256": sha256_file(path)} for path in output_files},
        "status_counts": status_counts,
        "warnings": warnings.to_dict(orient="records"),
        "note": "No raw/Bismark processing, external DMR caller, or DMR statistical model change was performed.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
