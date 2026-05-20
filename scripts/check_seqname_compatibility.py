#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis import read_seqname_aliases, write_seqname_compatibility_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare seqname compatibility between two tables.")
    parser.add_argument("--left", required=True)
    parser.add_argument("--right", required=True)
    parser.add_argument("--left-column", required=True)
    parser.add_argument("--right-column", required=True)
    parser.add_argument("--aliases")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    alias_map = read_seqname_aliases(args.aliases) if args.aliases else None
    left = pd.read_csv(args.left, sep=None, engine="python")
    right = pd.read_csv(args.right, sep=None, engine="python")
    write_seqname_compatibility_report(left[args.left_column], right[args.right_column], args.out, alias_map)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
