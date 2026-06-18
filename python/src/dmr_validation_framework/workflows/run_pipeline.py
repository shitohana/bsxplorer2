#!/usr/bin/env python
"""End-to-end terminal UX for caller execution, validation, and reports."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from dmr_validation_framework.core.io import SEARCH_ROOTS_ENV
from dmr_validation_framework.reports import (
    caller_summary,
    interactive_report,
    upset_caller_support,
    validation_summary,
)
from dmr_validation_framework.workflows import (
    critical_validation,
    external_callers,
    glm_glmm_validation,
    run_external_callers,
    thesis_figures,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the DMR validation framework pipeline: external callers, optional harmonization, "
            "controlled GLM/GLMM confirmation, validation profile, and summary reports."
        )
    )
    parser.add_argument("--root", required=True, help="Dataset root with raw_cx/ and metadata/sample_manifest.tsv.")
    parser.add_argument("--rscript", default=run_external_callers.DEFAULT_RSCRIPT)
    parser.add_argument("--r-lib", action="append")
    parser.add_argument("--external-root", help="Defaults to <root>/external_callers.")
    parser.add_argument("--validation-out-dir", help="Defaults to <root>/validation_<profile>.")
    parser.add_argument("--report-dir", help="Defaults to <root>/dmr_validation_report.")
    parser.add_argument("--harmonized-out-dir", help="Defaults to <root>/external_callers_harmonized.")
    parser.add_argument("--glm-glmm-out-dir", help="Defaults to <root>/glm_glmm_validation.")
    parser.add_argument("--target-dmr-table", help="Optional target/internal DMR table for caller-support overlap.")
    parser.add_argument("--contexts", default="CG,CHG,CHH")
    parser.add_argument("--callers", default=run_external_callers.DEFAULT_CALLERS)
    parser.add_argument("--profile", default="core")
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--min-sites", type=int, default=3)
    parser.add_argument("--max-gap", type=int, default=1000)
    parser.add_argument("--q-threshold", type=float, default=0.10)
    parser.add_argument("--window-size", type=int, default=1000)
    parser.add_argument("--step-size", type=int, default=1000)
    parser.add_argument(
        "--methylkit-min-diff",
        type=float,
        default=10.0,
        help="methylKit effect-size threshold in percentage points (caller-specific; not comparable with DSS/metilene delta).",
    )
    parser.add_argument("--max-sites-per-context", type=int, default=0)
    parser.add_argument("--read-chunk-size", type=int, default=1_000_000)
    parser.add_argument(
        "--dmrseq-permutations",
        type=int,
        default=0,
        help="Forwarded to dmrseq::dmrseq(maxPerms=...). 0 keeps the package default.",
    )
    parser.add_argument("--metilene")
    parser.add_argument("--glm-glmm-caller", default="methylKit")
    parser.add_argument("--glm-glmm-top-n-per-context", type=int, default=10)
    parser.add_argument("--glm-glmm-condition-column", default="condition")
    parser.add_argument("--glm-glmm-case-label", default="treatment")
    parser.add_argument("--glm-glmm-control-label", default="control")
    parser.add_argument("--glm-glmm-covariates", default="")
    parser.add_argument("--glm-glmm-min-cpg", type=int, default=3)
    parser.add_argument("--glm-glmm-max-zero-coverage-fraction", type=float, default=0.8)
    parser.add_argument("--glm-glmm-coverage-set-mode", choices=("per_sample", "common"), default="per_sample")
    parser.add_argument("--glm-glmm-min-common-cpgs", type=int, default=3)
    parser.add_argument("--annotation-dir", help="Optional gene/TE/GO annotation package directory. Enables the annotation figure stage and the report Annotation tab.")
    parser.add_argument("--expression-table", help="Optional gene-level expression evidence table (gene_id-like column) for the annotation figures.")
    parser.add_argument("--thesis-figures-out-dir", help="Defaults to <report-dir>/thesis_figures.")
    # Metagene / annotation-figure tuning (forwarded to the thesis-figures stage).
    parser.add_argument("--promoter-window", type=int, default=2000, help="Promoter half-window (bp) for DMR-gene assignment.")
    parser.add_argument("--upstream-len", type=int, default=2000, help="Metagene upstream length (bp).")
    parser.add_argument("--downstream-len", type=int, default=2000, help="Metagene downstream length (bp).")
    parser.add_argument("--upstream-bins", type=int, default=20, help="Metagene upstream bin count.")
    parser.add_argument("--body-bins", type=int, default=100, help="Metagene gene-body bin count.")
    parser.add_argument("--downstream-bins", type=int, default=20, help="Metagene downstream bin count.")
    parser.add_argument("--random-control-iterations", type=int, default=100, help="Random-control occupancy iterations (0 disables).")
    parser.add_argument("--max-random-regions-per-iteration", type=int, default=5000)
    parser.add_argument("--max-metagene-dmrs", type=int, default=100_000, help="Cap on DMRs used for the lightweight metagene audit.")
    parser.add_argument("--metagene-seed", type=int, default=202715)
    parser.add_argument("--max-forest-rows", type=int, default=40, help="Max regions in the bootstrap CI forest figure.")
    parser.add_argument("--skip-extended-metagene", action="store_true", help="Skip the metagene occupancy / projection-sensitivity computation in the annotation stage.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip caller stage if external-root already has a completed manifest.")
    parser.add_argument("--skip-callers", action="store_true")
    parser.add_argument("--skip-harmonize", action="store_true")
    parser.add_argument("--skip-glm-glmm", action="store_true")
    parser.add_argument("--skip-validation", action="store_true")
    parser.add_argument("--skip-reports", action="store_true")
    parser.add_argument("--skip-thesis-figures", action="store_true", help="Do not build the annotation/downstream figure bundle even if --annotation-dir is given.")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--strict", action="store_true", help="Return non-zero if validation/harmonization/report stages fail.")
    parser.add_argument("--dry-run", action="store_true", help="Print/write manifests but do not execute R callers.")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def default_target(root: Path) -> Path:
    return root / "strict_all_windows" / "processed_geo_cx" / "run" / "dmr_evidence_scores.tsv"


def stage_row(stage: str, status: str, returncode: int = 0, notes: str = "") -> dict:
    return {
        "stage": stage,
        "status": status,
        "returncode": returncode,
        "notes": notes,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def write_manifest(report_dir: Path, args: argparse.Namespace, rows: list[dict], outputs: dict[str, str]) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "root": args.root,
        "contexts": split_csv(args.contexts),
        "callers": split_csv(args.callers),
        "profile": args.profile,
        "parameters": {
            "min_coverage": args.min_coverage,
            "min_sites": args.min_sites,
            "max_gap": args.max_gap,
            "q_threshold": args.q_threshold,
            "max_sites_per_context": args.max_sites_per_context,
            "read_chunk_size": args.read_chunk_size,
            "dmrseq_permutations": args.dmrseq_permutations,
            "glm_glmm_caller": args.glm_glmm_caller,
            "glm_glmm_top_n_per_context": args.glm_glmm_top_n_per_context,
        },
        "stages": rows,
        "outputs": outputs,
        "note": (
            "run-pipeline is an orchestration wrapper. Caller-native statistics remain caller-specific; "
            "validation outputs are robustness diagnostics, not a new genome-wide DMR caller."
        ),
    }
    (report_dir / "pipeline_manifest.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def caller_manifest_completed(external_root: Path) -> bool:
    path = external_root / "external_caller_runner_manifest.json"
    if not path.exists():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    return str(payload.get("status", "")).startswith("completed")


def run_callers(args: argparse.Namespace, root: Path, external_root: Path) -> int:
    caller_args = argparse.Namespace(
        root=str(root),
        rscript=args.rscript,
        r_lib=args.r_lib,
        external_root=str(external_root),
        contexts=args.contexts,
        callers=args.callers,
        min_coverage=args.min_coverage,
        min_sites=args.min_sites,
        max_gap=args.max_gap,
        q_threshold=args.q_threshold,
        window_size=args.window_size,
        step_size=args.step_size,
        methylkit_min_diff=args.methylkit_min_diff,
        max_sites_per_context=args.max_sites_per_context,
        read_chunk_size=args.read_chunk_size,
        dmrseq_permutations=args.dmrseq_permutations,
        metilene=args.metilene,
        dry_run=args.dry_run,
        skip_harmonize=True,
        harmonize_max_target_regions=200_000,
    )
    return run_external_callers.run(caller_args)


def run_harmonization(
    args: argparse.Namespace,
    root: Path,
    external_root: Path,
    harmonized_out_dir: Path,
    target_path: Path,
) -> int:
    harmonize_args = argparse.Namespace(
        root=str(root),
        target_dmr_table=str(target_path),
        external_root=str(external_root),
        out_dir=str(harmonized_out_dir),
        caller_output=[],
        caller_output_manifest=None,
        context=split_csv(args.contexts),
        contrast_id="control_vs_treatment",
        condition_a="control",
        condition_b="treatment",
        overlap_threshold=0.5,
        include_candidate_only=False,
        max_target_regions=200_000,
        no_auto_discover=False,
        write_run_plan=True,
    )
    return external_callers.run(harmonize_args)


def run_glm_glmm(args: argparse.Namespace, root: Path, external_root: Path, glm_glmm_out_dir: Path) -> int:
    glm_args = argparse.Namespace(
        root=str(root),
        external_root=str(external_root),
        out_dir=str(glm_glmm_out_dir),
        rscript=args.rscript,
        contexts=args.contexts,
        glm_caller=args.glm_glmm_caller,
        top_n_per_context=args.glm_glmm_top_n_per_context,
        q_threshold=args.q_threshold,
        candidate_q_threshold=None,
        condition_column=args.glm_glmm_condition_column,
        case_label=args.glm_glmm_case_label,
        control_label=args.glm_glmm_control_label,
        covariates=args.glm_glmm_covariates,
        min_cpg=args.glm_glmm_min_cpg,
        min_replicates_per_group=2,
        min_total=args.min_coverage,
        max_zero_coverage_fraction=args.glm_glmm_max_zero_coverage_fraction,
        coverage_set_mode=args.glm_glmm_coverage_set_mode,
        min_common_cpgs=args.glm_glmm_min_common_cpgs,
        no_coverage_qc=False,
        no_plots=args.no_plots,
    )
    return glm_glmm_validation.run(glm_args)


def run_validation(
    args: argparse.Namespace,
    root: Path,
    external_root: Path,
    validation_out_dir: Path,
    harmonized_out_dir: Path,
    glm_glmm_out_dir: Path | None = None,
) -> int:
    roots = []
    if glm_glmm_out_dir is not None:
        roots.extend([glm_glmm_out_dir, glm_glmm_out_dir / "glm_vs_glmm", glm_glmm_out_dir / "glmm"])
    roots.extend([external_root, harmonized_out_dir, validation_out_dir, root / "outputs", root])
    old_value = os.environ.get(SEARCH_ROOTS_ENV)
    os.environ[SEARCH_ROOTS_ENV] = os.pathsep.join(str(path) for path in roots if path.exists())
    try:
        validation_args = argparse.Namespace(out_dir=validation_out_dir, profile=args.profile, list_profiles=False)
        return critical_validation.run(validation_args)
    finally:
        if old_value is None:
            os.environ.pop(SEARCH_ROOTS_ENV, None)
        else:
            os.environ[SEARCH_ROOTS_ENV] = old_value


def run_thesis_figures(
    args: argparse.Namespace,
    root: Path,
    external_root: Path,
    validation_out_dir: Path,
    report_dir: Path,
    glm_glmm_out_dir: Path,
    thesis_out_dir: Path,
) -> int:
    """Build the annotation/downstream figure bundle (gene/TE/GO/metagene)."""
    thesis_args = argparse.Namespace(
        root=str(root),
        run_id="",
        external_root=str(external_root),
        report_dir=str(report_dir),
        validation_dir=str(validation_out_dir),
        glm_glmm_dir=str(glm_glmm_out_dir),
        annotation_dir=args.annotation_dir,
        out_dir=str(thesis_out_dir),
        figure_prefix="fig_dmr",
        start_index=1,
        max_forest_rows=args.max_forest_rows,
        max_table_copy_rows=200_000,
        expression_table=args.expression_table,
        promoter_window=args.promoter_window,
        skip_extended_metagene=args.skip_extended_metagene,
        random_control_iterations=args.random_control_iterations,
        max_random_regions_per_iteration=args.max_random_regions_per_iteration,
        max_metagene_dmrs=args.max_metagene_dmrs,
        metagene_seed=args.metagene_seed,
        upstream_len=args.upstream_len,
        downstream_len=args.downstream_len,
        upstream_bins=args.upstream_bins,
        body_bins=args.body_bins,
        downstream_bins=args.downstream_bins,
    )
    return thesis_figures.run(thesis_args)


def run_reports(
    args: argparse.Namespace,
    root: Path,
    external_root: Path,
    validation_out_dir: Path,
    report_dir: Path,
    glm_glmm_out_dir: Path | None = None,
    harmonized_out_dir: Path | None = None,
) -> int:
    caller_report_dir = report_dir / "caller_summary"
    validation_report_dir = report_dir / "validation_summary"
    caller_code = caller_summary.run(
        argparse.Namespace(
            external_root=str(external_root),
            out_dir=str(caller_report_dir),
            contexts=args.contexts,
            callers=args.callers,
            status_table=str(external_root / "external_caller_run_status.tsv"),
            max_plot_rows=200_000,
            no_plots=args.no_plots,
        )
    )
    validation_code = validation_summary.run(
        argparse.Namespace(
            validation_dir=str(validation_out_dir),
            out_dir=str(validation_report_dir),
            no_plots=args.no_plots,
        )
    )

    # UpSet caller-support intersections. Builds the consensus support matrix
    # from the caller rows when harmonization was skipped, then renders
    # ComplexUpset (R). Skips gracefully if R/ComplexUpset are unavailable.
    upset_code = 0
    upset_figures_dir: Path | None = None
    if not args.no_plots:
        upset_dir = report_dir / "upset"
        upset_code = upset_caller_support.run(
            argparse.Namespace(
                out_dir=str(upset_dir),
                support_matrix=None,
                caller_rows=str(caller_report_dir / "caller_dmr_rows_for_plots.tsv"),
                rscript=args.rscript,
                overlap_threshold=0.5,
                min_size=1,
            )
        )
        if upset_dir.exists() and any(upset_dir.glob("upset_caller_support_*.png")):
            upset_figures_dir = upset_dir

    # Optional annotation/downstream figure bundle. Runs after the static
    # report tables exist (it consumes caller_dmr_rows_for_plots.tsv) and before
    # the interactive report, which embeds its PNGs in the Annotation tab.
    thesis_code = 0
    annotation_figures_dir: Path | None = None
    if args.annotation_dir and not args.skip_thesis_figures and glm_glmm_out_dir is not None:
        thesis_out_dir = (
            Path(args.thesis_figures_out_dir)
            if args.thesis_figures_out_dir
            else report_dir / "thesis_figures"
        )
        thesis_code = run_thesis_figures(
            args, root, external_root, validation_out_dir, report_dir, glm_glmm_out_dir, thesis_out_dir
        )
        figures_png = thesis_out_dir / "figures_png"
        if figures_png.exists():
            annotation_figures_dir = figures_png

    report_code = 0
    if not args.no_plots:
        # Combined interactive report: read the tables the static reports just
        # wrote, plus validation/glm-glmm/consensus outputs discovered via the
        # search roots below.
        search_roots = [
            path
            for path in (
                caller_report_dir,
                validation_report_dir,
                validation_out_dir,
                report_dir / "upset",
                glm_glmm_out_dir,
                glm_glmm_out_dir / "glm_vs_glmm" if glm_glmm_out_dir else None,
                harmonized_out_dir,
                external_root,
            )
            if path is not None
        ]
        report_code = interactive_report.run(
            argparse.Namespace(
                out_dir=str(report_dir / "report"),
                caller_rows=str(caller_report_dir / "caller_dmr_rows_for_plots.tsv"),
                caller_summary=str(caller_report_dir / "caller_context_summary.tsv"),
                consensus=None,
                glm_glmm=None,
                bootstrap=None,
                overlap_sensitivity=None,
                confidence=None,
                annotation_figures_dir=str(annotation_figures_dir) if annotation_figures_dir else None,
                upset_figures_dir=str(upset_figures_dir) if upset_figures_dir else None,
                search_root=[str(path) for path in search_roots],
                title="DMR validation report",
                width=820,
                height=420,
            )
        )
    return max(caller_code, validation_code, upset_code, thesis_code, report_code)


def run(args: argparse.Namespace) -> int:
    root = Path(args.root)
    external_root = Path(args.external_root) if args.external_root else root / "external_callers"
    validation_out_dir = Path(args.validation_out_dir) if args.validation_out_dir else root / f"validation_{args.profile}"
    report_dir = Path(args.report_dir) if args.report_dir else root / "dmr_validation_report"
    harmonized_out_dir = Path(args.harmonized_out_dir) if args.harmonized_out_dir else root / "external_callers_harmonized"
    glm_glmm_out_dir = Path(args.glm_glmm_out_dir) if args.glm_glmm_out_dir else root / "glm_glmm_validation"
    target_path = Path(args.target_dmr_table) if args.target_dmr_table else default_target(root)

    rows: list[dict] = []
    outputs = {
        "external_root": str(external_root),
        "validation_out_dir": str(validation_out_dir),
        "report_dir": str(report_dir),
        "harmonized_out_dir": str(harmonized_out_dir),
        "glm_glmm_out_dir": str(glm_glmm_out_dir),
    }

    if args.skip_callers:
        rows.append(stage_row("run_external_callers", "SKIPPED", notes="--skip-callers"))
    elif args.skip_existing and caller_manifest_completed(external_root):
        rows.append(stage_row("run_external_callers", "SKIPPED", notes="existing completed caller manifest"))
    else:
        code = run_callers(args, root, external_root)
        rows.append(stage_row("run_external_callers", "PASS" if code == 0 else "FAIL", code))
        if code != 0 and args.strict:
            write_manifest(report_dir, args, rows, outputs)
            return code

    if args.skip_harmonize:
        rows.append(stage_row("harmonize_callers", "SKIPPED", notes="--skip-harmonize"))
    elif not target_path.exists():
        rows.append(stage_row("harmonize_callers", "SKIPPED", notes=f"target DMR table not found: {target_path}"))
    else:
        code = run_harmonization(args, root, external_root, harmonized_out_dir, target_path)
        rows.append(stage_row("harmonize_callers", "PASS" if code == 0 else "FAIL", code))
        if code != 0 and args.strict:
            write_manifest(report_dir, args, rows, outputs)
            return code

    if args.skip_glm_glmm:
        rows.append(stage_row("glm_glmm_validation", "SKIPPED", notes="--skip-glm-glmm"))
    elif args.dry_run:
        rows.append(stage_row("glm_glmm_validation", "SKIPPED", notes="--dry-run"))
    else:
        code = run_glm_glmm(args, root, external_root, glm_glmm_out_dir)
        rows.append(stage_row("glm_glmm_validation", "PASS" if code == 0 else "WARN", code, notes="controlled top-N confirmatory layer"))
        if code != 0 and args.strict:
            write_manifest(report_dir, args, rows, outputs)
            return code

    if args.skip_validation:
        rows.append(stage_row("validation_profile", "SKIPPED", notes="--skip-validation"))
    else:
        code = run_validation(args, root, external_root, validation_out_dir, harmonized_out_dir, glm_glmm_out_dir)
        status = "PASS" if code == 0 else "WARN"
        rows.append(stage_row("validation_profile", status, code, notes="non-zero usually means at least one check failed"))
        if code != 0 and args.strict:
            write_manifest(report_dir, args, rows, outputs)
            return code

    if args.skip_reports:
        rows.append(stage_row("reports", "SKIPPED", notes="--skip-reports"))
    else:
        code = run_reports(args, root, external_root, validation_out_dir, report_dir, glm_glmm_out_dir, harmonized_out_dir)
        rows.append(stage_row("reports", "PASS" if code == 0 else "FAIL", code))
        if code != 0 and args.strict:
            write_manifest(report_dir, args, rows, outputs)
            return code

    write_manifest(report_dir, args, rows, outputs)
    print(f"pipeline manifest: {report_dir / 'pipeline_manifest.json'}")
    if not args.skip_reports:
        print(f"combined interactive report: {report_dir / 'report' / 'index.html'}")
        print(f"caller summary: {report_dir / 'caller_summary' / 'caller_summary.html'}")
        print(f"validation summary: {report_dir / 'validation_summary' / 'validation_summary.html'}")
        print(f"upset figures: {report_dir / 'upset'}")
        if args.annotation_dir and not args.skip_thesis_figures:
            thesis_dir = Path(args.thesis_figures_out_dir) if args.thesis_figures_out_dir else report_dir / "thesis_figures"
            print(f"annotation figures: {thesis_dir / 'figures_pdf'}")
    if not args.skip_glm_glmm and not args.dry_run:
        print(f"glm/glmm plots: {glm_glmm_out_dir / 'glmm_plots'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except Exception as exc:  # noqa: BLE001
        print(f"run-pipeline failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
