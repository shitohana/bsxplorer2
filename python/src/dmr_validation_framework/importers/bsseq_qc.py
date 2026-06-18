#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


from bsx2.analysis import classify_bsseq_qc, compute_counts_qc, parse_bismark_alignment_report, parse_bismark_dedup_report, parse_bismark_mbias_report, write_bsseq_qc_outputs


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Import BS-seq QC reports without running raw pipeline tools.")
    parser.add_argument("--counts")
    parser.add_argument("--alignment-report")
    parser.add_argument("--dedup-report")
    parser.add_argument("--mbias-report")
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-warnings", required=True)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    rows = []
    warnings = []
    if args.counts:
        qc = compute_counts_qc(pd.read_csv(args.counts, sep=None, engine="python"))
        rows.append({"source": args.counts, **qc, **classify_bsseq_qc(qc)})
    if args.alignment_report:
        row = parse_bismark_alignment_report(args.alignment_report); rows.append(row)
        if row.get("warning"): warnings.append({"source": args.alignment_report, "warning": row["warning"]})
    if args.dedup_report:
        row = parse_bismark_dedup_report(args.dedup_report); rows.append(row)
        if row.get("warning"): warnings.append({"source": args.dedup_report, "warning": row["warning"]})
    if args.mbias_report:
        mbias = parse_bismark_mbias_report(args.mbias_report)
        rows.append({"source": args.mbias_report, "mbias_rows": len(mbias)})
    write_bsseq_qc_outputs(rows, warnings, args.out_summary, args.out_warnings)
    return 0



def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
