#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis import check_coordinate_bounds, read_genome_sizes, read_seqname_aliases, write_assembly_compatibility_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Check interval coordinate compatibility against optional genome sizes.")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--genome-sizes")
    parser.add_argument("--aliases")
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary-out")
    args = parser.parse_args()
    aliases = read_seqname_aliases(args.aliases) if args.aliases else None
    regions = pd.read_csv(args.regions, sep=None, engine="python")
    genome = read_genome_sizes(args.genome_sizes) if args.genome_sizes else None
    report = check_coordinate_bounds(regions, genome, aliases)
    write_assembly_compatibility_report(report, args.out)
    if args.summary_out:
        Path(args.summary_out).write_text(f"# Assembly Compatibility Summary\n\nwarnings={(report['compatibility_status'] != 'ok').sum()}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
