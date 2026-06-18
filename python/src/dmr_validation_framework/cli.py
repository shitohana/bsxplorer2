"""Unified CLI for DMR validation framework checks."""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

from .checks import CHECKS, PROFILE_DESCRIPTIONS, checks_for_profile, get_check, profile_names

WORKFLOWS = {
    "external-callers": "dmr_validation_framework.workflows.external_callers",
    "run-external-callers": "dmr_validation_framework.workflows.run_external_callers",
    "run-pipeline": "dmr_validation_framework.workflows.run_pipeline",
    "glm-glmm-validation": "dmr_validation_framework.workflows.glm_glmm_validation",
    "consensus-tiers": "dmr_validation_framework.workflows.consensus_tiers",
    "critical-validation": "dmr_validation_framework.workflows.critical_validation",
    "region-cpg-counts": "dmr_validation_framework.workflows.region_cpg_counts",
    "region-signal-aggregation": "dmr_validation_framework.workflows.region_signal_aggregation",
    "expression-support": "dmr_validation_framework.workflows.expression_support",
    "plant-annotation": "dmr_validation_framework.workflows.plant_annotation",
    "thesis-figures": "dmr_validation_framework.workflows.thesis_figures",
}

REPORTS = {
    "caller-summary": "dmr_validation_framework.reports.caller_summary",
    "validation-summary": "dmr_validation_framework.reports.validation_summary",
    "interactive-report": "dmr_validation_framework.reports.interactive_report",
    "upset": "dmr_validation_framework.reports.upset_caller_support",
    "curve-bundle": "dmr_validation_framework.reports.dmr_curve_bundle",
    "glmm-plots": "dmr_validation_framework.reports.glmm_plots",
}


def add_common_run_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("check_args", nargs=argparse.REMAINDER, help="Arguments passed after '--' to the underlying check script.")


def normalize_passthrough(args: list[str]) -> list[str]:
    if args and args[0] == "--":
        return args[1:]
    return args


def extract_out_dir_from_remainder(out_dir: Path, check_args: list[str]) -> tuple[Path, list[str]]:
    cleaned: list[str] = []
    i = 0
    while i < len(check_args):
        if check_args[i] == "--out-dir":
            if i + 1 >= len(check_args):
                raise SystemExit("--out-dir requires a value")
            out_dir = Path(check_args[i + 1])
            i += 2
            continue
        cleaned.append(check_args[i])
        i += 1
    return out_dir, cleaned


def run_check(check_name: str, out_dir: Path, passthrough: list[str]) -> int:
    check = get_check(check_name)
    try:
        module = importlib.import_module(check.module)
    except Exception as exc:
        print(f"Could not import check module {check.module}: {exc}", file=sys.stderr)
        return 2
    argv = ["--out-dir", str(out_dir), *normalize_passthrough(passthrough)]
    try:
        return int(module.main(argv))
    except SystemExit as exc:
        return int(exc.code or 0)


def run_workflow(workflow_name: str, passthrough: list[str]) -> int:
    module_name = WORKFLOWS[workflow_name]
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:
        print(f"Could not import workflow module {module_name}: {exc}", file=sys.stderr)
        return 2
    try:
        return int(module.main(normalize_passthrough(passthrough)))
    except SystemExit as exc:
        return int(exc.code or 0)


def run_report(report_name: str, passthrough: list[str]) -> int:
    module_name = REPORTS[report_name]
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:
        print(f"Could not import report module {module_name}: {exc}", file=sys.stderr)
        return 2
    try:
        return int(module.main(normalize_passthrough(passthrough)))
    except SystemExit as exc:
        return int(exc.code or 0)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="DMR validation framework command runner.")
    sub = parser.add_subparsers(dest="command", required=True)

    list_parser = sub.add_parser("list", help="List available validation checks.")
    list_parser.add_argument("--profile", choices=profile_names())

    sub.add_parser("profiles", help="List validation profiles.")

    run_parser = sub.add_parser("run", help="Run one validation check.")
    run_parser.add_argument("check", choices=[check.name for check in CHECKS])
    add_common_run_args(run_parser)

    all_parser = sub.add_parser("run-all", help="Run checks from a validation profile.")
    all_parser.add_argument("--profile", choices=profile_names(), default="core")
    add_common_run_args(all_parser)

    workflow_parser = sub.add_parser("workflow", help="Run a validation framework workflow.")
    workflow_parser.add_argument("workflow", choices=sorted(WORKFLOWS))
    workflow_parser.add_argument(
        "workflow_args",
        nargs=argparse.REMAINDER,
        help="Arguments passed after '--' to the underlying workflow.",
    )

    report_parser = sub.add_parser("report", help="Build a DMR validation report or plot bundle.")
    report_parser.add_argument("report", choices=sorted(REPORTS))
    report_parser.add_argument(
        "report_args",
        nargs=argparse.REMAINDER,
        help="Arguments passed after '--' to the underlying report.",
    )

    args = parser.parse_args(argv)
    if args.command == "list":
        checks = checks_for_profile(args.profile) if args.profile else CHECKS
        for check in checks:
            profile_text = ",".join(check.profiles)
            print(f"{check.name}\t{profile_text}\t{check.layer}\t{check.entrypoint}\t{check.description}")
        return 0
    if args.command == "profiles":
        for name, description in PROFILE_DESCRIPTIONS.items():
            checks = ",".join(check.name for check in checks_for_profile(name))
            print(f"{name}\t{description}\tchecks={checks}")
        return 0
    if args.command == "run":
        args.out_dir, args.check_args = extract_out_dir_from_remainder(args.out_dir, args.check_args)
        return run_check(args.check, args.out_dir, args.check_args)
    if args.command == "run-all":
        args.out_dir, args.check_args = extract_out_dir_from_remainder(args.out_dir, args.check_args)
        failures = []
        for check in checks_for_profile(args.profile):
            code = run_check(check.name, args.out_dir, args.check_args)
            if code != 0:
                failures.append((check.name, code))
        if failures:
            for name, code in failures:
                print(f"{name}\tFAIL\treturncode={code}", file=sys.stderr)
            return 1
        return 0
    if args.command == "workflow":
        return run_workflow(args.workflow, args.workflow_args)
    if args.command == "report":
        return run_report(args.report, args.report_args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
