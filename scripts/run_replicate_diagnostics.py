#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis import compute_pairwise_replicate_consistency, compute_region_replicate_summary, write_replicate_diagnostics_outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Run region-level replicate diagnostics.")
    parser.add_argument("--region-counts", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    counts = pd.read_csv(args.region_counts, sep=None, engine="python")
    design = pd.read_csv(args.design, sep=None, engine="python")
    condition = compute_region_replicate_summary(counts, design)
    consistency = compute_pairwise_replicate_consistency(counts, design)
    write_replicate_diagnostics_outputs(condition, consistency, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
