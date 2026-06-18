#!/usr/bin/env python3
"""Check delta-methylation direction agreement for matched DMR candidates."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np

from dmr_validation_framework.core.intervals import matched_pairs
from dmr_validation_framework.core.io import ensure_out_dir, load_caller_tables, write_tsv


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--thresholds", default="0.3,0.5,0.8")
    parser.add_argument("--delta-thresholds", default="0,0.05,0.1")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    thresholds = [float(x) for x in args.thresholds.split(",") if x.strip()]
    delta_thresholds = [float(x) for x in args.delta_thresholds.split(",") if x.strip()]
    tables = load_caller_tables()
    rows: list[dict] = []
    if len(tables) < 2:
        rows.append(
            {
                "threshold": "NA",
                "delta_threshold": "NA",
                "caller_a": "NA",
                "caller_b": "NA",
                "n_matched": 0,
                "n_tested_after_delta_threshold": 0,
                "n_same_direction": 0,
                "n_opposite_direction": 0,
                "n_zero_or_below_threshold": 0,
                "fraction_same_direction": "NA",
                "status": "SKIPPED",
                "notes": "fewer than two caller DMR tables found",
            }
        )
    else:
        for threshold in thresholds:
            for a, b in combinations(tables, 2):
                if "delta" not in a.data.columns or "delta" not in b.data.columns:
                    for delta0 in delta_thresholds:
                        rows.append(
                            {
                                "threshold": threshold,
                                "delta_threshold": delta0,
                                "caller_a": a.name,
                                "caller_b": b.name,
                                "n_matched": 0,
                                "n_tested_after_delta_threshold": 0,
                                "n_same_direction": 0,
                                "n_opposite_direction": 0,
                                "n_zero_or_below_threshold": 0,
                                "fraction_same_direction": "NA",
                                "status": "SKIPPED",
                                "notes": "delta column missing in at least one caller table",
                            }
                        )
                    continue
                pairs = matched_pairs(a.data, b.data, threshold)
                for delta0 in delta_thresholds:
                    if pairs.empty:
                        n_tested = n_same = n_opp = n_excluded = 0
                        frac = "NA"
                    else:
                        valid = pairs.dropna(subset=["delta_a", "delta_b"]).copy()
                        delta_a = valid["delta_a"].astype(float)
                        delta_b = valid["delta_b"].astype(float)
                        if delta0 == 0:
                            tested_mask = (delta_a != 0) & (delta_b != 0)
                        else:
                            tested_mask = (delta_a.abs() >= delta0) & (delta_b.abs() >= delta0)
                        tested = valid[tested_mask].copy()
                        n_tested = len(tested)
                        n_excluded = len(pairs) - n_tested
                        products = tested["delta_a"].astype(float) * tested["delta_b"].astype(float)
                        n_same = int((products > 0).sum())
                        n_opp = int((products < 0).sum())
                        frac = n_same / n_tested if n_tested else "NA"
                    excluded_fraction = n_excluded / len(pairs) if len(pairs) else 0
                    rows.append(
                        {
                            "threshold": threshold,
                            "delta_threshold": delta0,
                            "caller_a": a.name,
                            "caller_b": b.name,
                            "n_matched": len(pairs),
                            "n_tested_after_delta_threshold": n_tested,
                            "n_same_direction": n_same,
                            "n_opposite_direction": n_opp,
                            "n_zero_or_below_threshold": n_excluded,
                            "fraction_same_direction": frac,
                            "status": "WARN" if excluded_fraction > 0.5 and len(pairs) else ("PASS" if len(pairs) else "WARN"),
                            "notes": "direction agreement is evaluated only for effects above |Delta| threshold; q-values are not treated as equivalent",
                        }
                    )
    write_tsv(out_dir / "direction_agreement.tsv", rows)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
