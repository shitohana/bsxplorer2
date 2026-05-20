#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis import RegionSignalConfig, aggregate_region_signal, read_seqname_aliases, write_region_signal_qc, write_region_signal_table


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate point methylation counts over interval regions.")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--context")
    parser.add_argument("--aliases")
    parser.add_argument("--strand-policy", default="ignore", choices=["ignore", "same", "opposite"])
    args = parser.parse_args()
    aliases = read_seqname_aliases(args.aliases) if args.aliases else None
    regions = pd.read_csv(args.regions, sep=None, engine="python")
    counts = pd.read_csv(args.counts, sep=None, engine="python")
    result = aggregate_region_signal(regions, counts, RegionSignalConfig(context=args.context, strand_policy=args.strand_policy, seqname_aliases=aliases))
    write_region_signal_table(result, args.out)
    write_region_signal_qc(result, args.qc_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
