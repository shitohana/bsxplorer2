#!/usr/bin/env python
"""UpSet plots of DMR caller-support intersections (ComplexUpset, R).

Builds the caller-support membership table from a consensus support matrix (or,
when none exists yet, directly from the per-DMR caller rows the caller summary
writes) and renders ComplexUpset plots via a small R helper. The R step is
optional: if Rscript or the R packages are unavailable, the membership table is
still written and the step exits cleanly with a note.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from dmr_validation_framework.core.io import find_preferred_file, read_table

SUPPORT_MATRIX_NAMES = [
    "dmr_caller_support_matrix.tsv",
    "caller_support_matrix.tsv",
    "dmr_method_support_matrix.tsv",
    "external_caller_support_matrix.tsv",
]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, help="Directory for the membership table and UpSet figures.")
    parser.add_argument("--support-matrix", help="Consensus support matrix TSV (auto-discovered if omitted).")
    parser.add_argument("--caller-rows", help="caller_dmr_rows_for_plots.tsv, used to build a support matrix if none exists.")
    parser.add_argument("--rscript", help="Path to Rscript/Rscript.exe for the ComplexUpset step.")
    parser.add_argument("--overlap-threshold", type=float, default=0.5)
    parser.add_argument("--min-size", type=int, default=1, help="Minimum intersection size shown in the UpSet plot.")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def _read(path: Path | None) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    try:
        return read_table(path)
    except Exception:
        return pd.DataFrame()


def support_matrix_from_caller_rows(caller_rows: pd.DataFrame, overlap_threshold: float) -> pd.DataFrame:
    """Run the canonical consensus clustering on the flat per-DMR caller rows."""
    from bsx2.analysis.dmr_harmonization import build_caller_support_matrix

    df = caller_rows.copy()
    # _caller / _context are the authoritative per-row caller and context written
    # by the caller summary; prefer them over any canonical columns.
    if "_caller" in df.columns:
        df["source_caller"] = df["_caller"]
    elif "source_caller" not in df.columns and "caller" in df.columns:
        df["source_caller"] = df["caller"]
    if "_context" in df.columns:
        df["context"] = df["_context"]
    required = {"chrom", "start", "end", "context", "source_caller"}
    if not required.issubset(df.columns):
        return pd.DataFrame()
    return build_caller_support_matrix([df], overlap_threshold=overlap_threshold)


def membership_from_support(support: pd.DataFrame) -> pd.DataFrame:
    """One-hot caller membership per consensus region, preserving real names."""
    if support.empty or "supporting_callers" not in support.columns:
        return pd.DataFrame()
    rows: list[dict] = []
    callers: list[str] = []
    for value in support["supporting_callers"].dropna().astype(str):
        for caller in (item.strip() for item in value.split(",")):
            if caller and caller not in callers:
                callers.append(caller)
    if not callers:
        return pd.DataFrame()
    for idx, row in support.reset_index(drop=True).iterrows():
        present = {item.strip() for item in str(row.get("supporting_callers", "")).split(",") if item.strip()}
        record: dict = {"region_id": str(row.get("region_id", row.get("harmonized_region_id", idx)))}
        if "context" in support.columns:
            record["context"] = str(row.get("context", "NA"))
        for caller in callers:
            record[caller] = int(caller in present)
        rows.append(record)
    return pd.DataFrame(rows)


def resolve_rscript(path: str | None) -> str | None:
    if path:
        candidate = Path(path)
        if candidate.exists():
            return str(candidate)
        found = shutil.which(path)
        if found:
            return found
    return shutil.which("Rscript")


def run(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Source precedence (deterministic): an explicit support matrix, then the
    # current caller rows (built fresh so the UpSet reflects exactly the callers
    # of this run), and only as a last resort an auto-discovered matrix. Auto-
    # discovery is for standalone CLI use; in a pipeline it can otherwise pick up
    # a stale/target matrix (e.g. one with an "internal" set), which silently
    # mislabels the plot.
    support = pd.DataFrame()
    if args.support_matrix:
        support = _read(Path(args.support_matrix))
    elif args.caller_rows and Path(args.caller_rows).exists():
        caller_rows = _read(Path(args.caller_rows))
        if not caller_rows.empty:
            support = support_matrix_from_caller_rows(caller_rows, args.overlap_threshold)
            if not support.empty:
                support.to_csv(out_dir / "dmr_caller_support_matrix.tsv", sep="\t", index=False)
    if support.empty:
        support = _read(find_preferred_file(SUPPORT_MATRIX_NAMES))

    membership = membership_from_support(support)
    if membership.empty:
        (out_dir / "upset_status.txt").write_text(
            "no caller-support membership could be built (no support matrix or caller rows).\n",
            encoding="utf-8",
        )
        print("upset: no membership table to plot", file=sys.stderr)
        return 0

    membership_path = out_dir / "upset_membership.tsv"
    membership.to_csv(membership_path, sep="\t", index=False)

    # Shared palette: one stable color per caller, matching the rest of the reports.
    from dmr_validation_framework.reports.palette import caller_color

    set_cols = [c for c in membership.columns if c not in {"region_id", "context"}]
    colors_path = out_dir / "upset_colors.tsv"
    pd.DataFrame(
        {"caller": set_cols, "color": [caller_color(c) for c in set_cols]}
    ).to_csv(colors_path, sep="\t", index=False)

    rscript = resolve_rscript(args.rscript)
    if rscript is None:
        (out_dir / "upset_status.txt").write_text(
            "Rscript not found; wrote membership table only. Install R + ComplexUpset to render UpSet plots.\n",
            encoding="utf-8",
        )
        print("upset: Rscript unavailable; wrote membership table only", file=sys.stderr)
        return 0

    r_script = Path(__file__).with_name("upset_caller_support.R")
    command = [
        rscript,
        str(r_script),
        "--membership",
        str(membership_path),
        "--out-dir",
        str(out_dir),
        "--colors",
        str(colors_path),
        "--min-size",
        str(args.min_size),
    ]
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        print(f"upset: ComplexUpset step exited with {completed.returncode}; see {out_dir / 'upset_status.txt'}", file=sys.stderr)
    else:
        print(f"upset figures: {out_dir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
