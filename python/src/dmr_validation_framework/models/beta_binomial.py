#!/usr/bin/env python
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def _fallback(out: Path, flag: str, *, returncode: int | None = None, stderr: str = "") -> None:
    """Write a visible fallback marker and surface the failure.

    The R beta-binomial step can be unavailable (no Rscript) or fail; we keep the
    return code at 0 so the pipeline does not crash, but the failure must not be
    silent. The fallback row carries the flag and return code, the captured R
    stderr is written next to the output, and a warning is printed to stderr so
    the orchestrator/manifest layer can pick it up.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{"region_id": pd.NA, "qc_flag": flag, "converged": False, "returncode": returncode}]
    ).to_csv(out, sep="\t", index=False)
    if stderr:
        out.with_suffix(out.suffix + ".stderr.log").write_text(stderr, encoding="utf-8")
    print(
        f"WARNING: beta-binomial validation did not complete ({flag}"
        + (f", returncode={returncode}" if returncode is not None else "")
        + f"); wrote fallback marker to {out}",
        file=sys.stderr,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run R beta-binomial DMR validation with graceful fallback.")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--mode", default="aggregated")
    parser.add_argument("--full", default="condition + batch")
    parser.add_argument("--reduced", default="batch")
    parser.add_argument("--script", default="R/beta_binom_dmr.R")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    out = Path(args.out)
    rscript = shutil.which("Rscript")
    if rscript is None:
        _fallback(out, "beta_binom_rscript_unavailable")
        return 0
    cmd = [rscript, args.script, "--counts", args.counts, "--design", args.design, "--out", str(out), "--mode", args.mode, "--full", args.full, "--reduced", args.reduced]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        _fallback(out, "beta_binom_failed", returncode=result.returncode, stderr=result.stderr or "")
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
