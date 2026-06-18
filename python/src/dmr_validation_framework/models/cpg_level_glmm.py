#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


from dmr_validation_framework.core.glm_glmm import (
    glmmTMB_available,
    read_design_table,
    read_dmr_evidence_for_glmm,
    read_region_cpg_counts,
    rscript_available,
    run_cpg_level_glmm_validation,
)
from dmr_validation_framework.core.coverage import aggregate_counts_common_cpg, build_cpg_coverage_qc


def sha256_file(path: Path) -> str | None:
    if not path.exists() or path.is_dir():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run optional confirmatory CpG-level GLMM validation for selected DMR candidates."
    )
    parser.add_argument("--region-cpg-counts", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--rscript", help="Optional explicit Rscript/Rscript.exe path for glmmTMB.")
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
    parser.add_argument("--coverage-set-mode", choices=("per_sample", "common"), default="per_sample")
    parser.add_argument("--min-coverage", type=int)
    parser.add_argument("--min-covered-per-group", type=int)
    parser.add_argument("--min-covered-A", type=int)
    parser.add_argument("--min-covered-B", type=int)
    parser.add_argument("--min-common-cpgs", type=int, default=3)
    parser.add_argument("--write-coverage-qc", action="store_true")
    parser.add_argument("--out-dir", required=True)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = out_dir / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    counts = read_region_cpg_counts(args.region_cpg_counts)
    design = read_design_table(args.design)
    evidence = read_dmr_evidence_for_glmm(args.dmr_evidence) if args.dmr_evidence else None
    covariates = [item.strip() for item in args.covariates.split(",") if item.strip()]
    min_covered_per_group = args.min_covered_per_group
    if args.min_covered_A is not None or args.min_covered_B is not None:
        min_covered_per_group = {
            "A": args.min_covered_A if args.min_covered_A is not None else args.min_covered_per_group or 1,
            "B": args.min_covered_B if args.min_covered_B is not None else args.min_covered_per_group or 1,
        }

    coverage_output_files = []
    if args.write_coverage_qc:
        common_cpg, region_qc = build_cpg_coverage_qc(
            counts,
            design,
            condition_col=args.condition_column,
            min_coverage=args.min_total if args.min_coverage is None else args.min_coverage,
            min_covered_per_group=min_covered_per_group,
            min_common_cpgs=args.min_common_cpgs,
        )
        common_region = aggregate_counts_common_cpg(
            counts,
            common_cpg,
            design,
            condition_col=args.condition_column,
            min_common_cpgs=args.min_common_cpgs,
        )
        common_cpg_path = out_dir / "coverage_set_common_cpg.tsv"
        region_qc_path = out_dir / "coverage_set_region_qc.tsv"
        common_region_path = out_dir / "coverage_set_common_region_counts.tsv"
        common_cpg.to_csv(common_cpg_path, sep="\t", index=False)
        region_qc.to_csv(region_qc_path, sep="\t", index=False)
        common_region.to_csv(common_region_path, sep="\t", index=False)
        coverage_output_files = [common_cpg_path, region_qc_path, common_region_path]

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
        rscript=args.rscript,
        coverage_set_mode=args.coverage_set_mode,
        min_coverage=args.min_coverage,
        min_covered_per_group=min_covered_per_group,
        min_common_cpgs=args.min_common_cpgs,
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
            f"- Rscript_available: {rscript_available(args.rscript)}",
            f"- glmmTMB_available: {glmmTMB_available(args.rscript)}",
            f"- model_status_counts: {status_counts}",
            f"- coverage_set_mode: {args.coverage_set_mode}",
            f"- min_coverage: {args.min_total if args.min_coverage is None else args.min_coverage}",
            f"- min_covered_per_group: {min_covered_per_group}",
            f"- min_common_cpgs: {args.min_common_cpgs}",
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
            "rscript": str(args.rscript) if args.rscript else None,
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
            "coverage_set_mode": args.coverage_set_mode,
            "min_coverage": args.min_total if args.min_coverage is None else args.min_coverage,
            "min_covered_per_group": min_covered_per_group,
            "min_common_cpgs": args.min_common_cpgs,
            "write_coverage_qc": args.write_coverage_qc,
        },
        "outputs": {path.name: {"path": str(path), "sha256": sha256_file(path)} for path in [*output_files, *coverage_output_files]},
        "status_counts": status_counts,
        "warnings": warnings.to_dict(orient="records"),
        "note": "No raw/Bismark processing, external DMR caller, or DMR statistical model change was performed.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
