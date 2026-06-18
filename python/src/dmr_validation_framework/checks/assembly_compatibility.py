#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from dmr_validation_framework.core.io import read_table

from bsx2.analysis import check_coordinate_bounds, read_genome_sizes, read_seqname_aliases, write_assembly_compatibility_report


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check interval coordinate compatibility against optional genome sizes.")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--genome-sizes")
    parser.add_argument("--aliases")
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary-out")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    aliases = read_seqname_aliases(args.aliases) if args.aliases else None
    regions = read_table(args.regions)
    genome = read_genome_sizes(args.genome_sizes) if args.genome_sizes else None
    report = check_coordinate_bounds(regions, genome, aliases)
    write_assembly_compatibility_report(report, args.out)
    if args.summary_out:
        Path(args.summary_out).write_text(f"# Assembly Compatibility Summary\n\nwarnings={(report['compatibility_status'] != 'ok').sum()}\n", encoding="utf-8")
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
