#!/usr/bin/env python
from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd


def _fallback(out: Path, flag: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"region_id": pd.NA, "qc_flag": flag, "converged": False}]).to_csv(out, sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run R beta-binomial DMR validation with graceful fallback.")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--mode", default="aggregated")
    parser.add_argument("--full", default="condition + batch")
    parser.add_argument("--reduced", default="batch")
    parser.add_argument("--script", default="R/beta_binom_dmr.R")
    args = parser.parse_args()
    out = Path(args.out)
    rscript = shutil.which("Rscript")
    if rscript is None:
        _fallback(out, "beta_binom_failed")
        return 0
    cmd = [rscript, args.script, "--counts", args.counts, "--design", args.design, "--out", str(out), "--mode", args.mode, "--full", args.full, "--reduced", args.reduced]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        _fallback(out, "beta_binom_failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
