#!/usr/bin/env python
"""Run available external DMR callers and then harmonize their outputs."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from dmr_validation_framework.workflows import external_callers

DEFAULT_RSCRIPT = r"C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe"
DEFAULT_CALLERS = "DSS,methylKit,dmrseq,BSmooth,comb-p,metilene"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run third-party DMR callers from Bismark CX reports and harmonize outputs.",
    )
    parser.add_argument("--root", required=True, help="Dataset root with raw_cx/ and metadata/sample_manifest.tsv.")
    parser.add_argument("--rscript", default=DEFAULT_RSCRIPT, help="Path to Rscript.exe/Rscript.")
    parser.add_argument(
        "--r-lib",
        action="append",
        help="Extra R library path. Repeat if needed. Passed to .libPaths().",
    )
    parser.add_argument("--external-root", help="Defaults to <root>/external_callers.")
    parser.add_argument("--contexts", default="CG,CHG,CHH")
    parser.add_argument(
        "--callers",
        default=DEFAULT_CALLERS,
        help="Comma-separated subset: DSS,methylKit,dmrseq,BSmooth,comb-p,metilene,DMRcate.",
    )
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
    parser.add_argument(
        "--max-sites-per-context",
        type=int,
        default=0,
        help="Smoke-test limiter. 0 means no limit.",
    )
    parser.add_argument(
        "--read-chunk-size",
        type=int,
        default=1_000_000,
        help="Rows per CX-report chunk when reading in R.",
    )
    parser.add_argument(
        "--dmrseq-permutations",
        type=int,
        default=0,
        help="Forwarded to dmrseq::dmrseq(maxPerms=...). 0 keeps the package default.",
    )
    parser.add_argument("--metilene", help="Optional metilene executable path.")
    parser.add_argument("--dry-run", action="store_true", help="Print command and write manifest; do not run R.")
    parser.add_argument(
        "--skip-harmonize",
        action="store_true",
        help="Only run callers; do not run workflow external-callers afterwards.",
    )
    parser.add_argument(
        "--harmonize-max-target-regions",
        type=int,
        default=200_000,
        help="Forwarded to workflow external-callers after caller execution.",
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def split_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def resolve_rscript(path: str) -> str:
    candidate = Path(path)
    if candidate.exists():
        return str(candidate)
    found = shutil.which(path)
    if found:
        return found
    raise FileNotFoundError(f"Rscript not found: {path}")


def build_r_command(args: argparse.Namespace, r_script: Path, external_root: Path) -> list[str]:
    command = [
        resolve_rscript(args.rscript),
        str(r_script),
        "--root",
        str(Path(args.root)),
        "--out-root",
        str(external_root),
        "--contexts",
        args.contexts,
        "--callers",
        args.callers,
        "--min-coverage",
        str(args.min_coverage),
        "--min-sites",
        str(args.min_sites),
        "--max-gap",
        str(args.max_gap),
        "--q-threshold",
        str(args.q_threshold),
        "--window-size",
        str(args.window_size),
        "--step-size",
        str(args.step_size),
        "--methylkit-min-diff",
        str(args.methylkit_min_diff),
        "--max-sites-per-context",
        str(args.max_sites_per_context),
        "--read-chunk-size",
        str(args.read_chunk_size),
    ]
    if args.dmrseq_permutations > 0:
        command.extend(["--dmrseq-permutations", str(args.dmrseq_permutations)])
    for r_lib in args.r_lib or []:
        command.extend(["--r-lib", r_lib])
    if args.metilene:
        command.extend(["--metilene", args.metilene])
    return command


def write_run_manifest(out_dir: Path, command: list[str], args: argparse.Namespace, status: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "root": args.root,
        "contexts": split_csv(args.contexts),
        "callers": split_csv(args.callers),
        "command": command,
        "note": "This step invokes third-party R/Bioconductor callers; harmonization is handled separately.",
    }
    (out_dir / "external_caller_runner_manifest.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )


def run(args: argparse.Namespace) -> int:
    root = Path(args.root)
    external_root = Path(args.external_root) if args.external_root else root / "external_callers"
    r_script = Path(__file__).with_name("run_external_callers.R")
    command = build_r_command(args, r_script, external_root)
    write_run_manifest(external_root, command, args, "dry_run" if args.dry_run else "started")

    if args.dry_run:
        print(" ".join(command))
        return 0

    completed = subprocess.run(command, cwd=Path.cwd(), check=False)
    if completed.returncode != 0:
        write_run_manifest(external_root, command, args, f"failed:{completed.returncode}")
        return completed.returncode
    write_run_manifest(external_root, command, args, "completed")

    if args.skip_harmonize:
        return 0

    harmonize_args = argparse.Namespace(
        root=str(root),
        target_dmr_table=None,
        external_root=str(external_root),
        out_dir=None,
        caller_output=[],
        caller_output_manifest=None,
        context=split_csv(args.contexts),
        contrast_id="control_vs_treatment",
        condition_a="control",
        condition_b="treatment",
        overlap_threshold=0.5,
        include_candidate_only=False,
        max_target_regions=args.harmonize_max_target_regions,
        no_auto_discover=False,
        write_run_plan=True,
    )
    return external_callers.run(harmonize_args)


def main(argv: list[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except Exception as exc:  # noqa: BLE001
        print(f"run-external-callers failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
