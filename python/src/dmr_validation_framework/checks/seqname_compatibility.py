#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from dmr_validation_framework.core.io import read_table

from bsx2.analysis import read_seqname_aliases, write_seqname_compatibility_report


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare seqname compatibility between two tables.")
    parser.add_argument("--left", required=True)
    parser.add_argument("--right", required=True)
    parser.add_argument("--left-column", required=True)
    parser.add_argument("--right-column", required=True)
    parser.add_argument("--aliases")
    parser.add_argument("--out", required=True)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    alias_map = read_seqname_aliases(args.aliases) if args.aliases else None
    left = read_table(args.left)
    right = read_table(args.right)
    write_seqname_compatibility_report(left[args.left_column], right[args.right_column], args.out, alias_map)
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
