#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from dmr_validation_framework.core.io import read_table

from bsx2.analysis import compute_pairwise_replicate_consistency, compute_region_replicate_summary, write_replicate_diagnostics_outputs


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run region-level replicate diagnostics.")
    parser.add_argument("--region-counts", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    counts = read_table(args.region_counts)
    design = read_table(args.design)
    condition = compute_region_replicate_summary(counts, design)
    consistency = compute_pairwise_replicate_consistency(counts, design)
    write_replicate_diagnostics_outputs(condition, consistency, args.out_dir)
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
