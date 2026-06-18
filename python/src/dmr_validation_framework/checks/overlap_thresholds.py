#!/usr/bin/env python3
"""Check caller overlap sensitivity across reciprocal-overlap thresholds."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

from dmr_validation_framework.core.intervals import matched_pairs
from dmr_validation_framework.core.io import ensure_out_dir, load_caller_tables, write_tsv


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--thresholds", default="0.3,0.5,0.8")
    parser.add_argument("--matching-policy", default="many_to_many,best_reciprocal,one_to_one_greedy")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def apply_matching_policy(pairs, policy: str):
    if pairs.empty or policy == "many_to_many":
        return pairs.copy()
    ordered = pairs.sort_values("reciprocal_overlap", ascending=False).copy()
    if policy == "best_reciprocal":
        best_b_for_a = ordered.drop_duplicates("region_id_a", keep="first").set_index("region_id_a")["region_id_b"].to_dict()
        best_a_for_b = ordered.drop_duplicates("region_id_b", keep="first").set_index("region_id_b")["region_id_a"].to_dict()
        mask = [
            best_b_for_a.get(row.region_id_a) == row.region_id_b
            and best_a_for_b.get(row.region_id_b) == row.region_id_a
            for row in ordered.itertuples()
        ]
        return ordered.loc[mask].copy()
    if policy == "one_to_one_greedy":
        used_a: set[str] = set()
        used_b: set[str] = set()
        rows = []
        for row in ordered.to_dict("records"):
            a_id = str(row["region_id_a"])
            b_id = str(row["region_id_b"])
            if a_id in used_a or b_id in used_b:
                continue
            used_a.add(a_id)
            used_b.add(b_id)
            rows.append(row)
        return pairs.__class__(rows)
    raise ValueError(f"unsupported matching policy: {policy}")


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    thresholds = [float(x) for x in args.thresholds.split(",") if x.strip()]
    policies = [x.strip() for x in args.matching_policy.split(",") if x.strip()]
    tables = load_caller_tables(include_support_matrix=False)
    if len(tables) < 2:
        row = {
            "threshold": "NA",
            "matching_policy": "NA",
            "caller_a": "NA",
            "caller_b": "NA",
            "n_a": 0,
            "n_b": 0,
            "n_pairs": 0,
            "n_matched_a": 0,
            "n_matched_b": 0,
            "jaccard_like": "NA",
            "reciprocal_overlap_definition": "strict_min_overlap_over_each_interval",
            "notes": (
                "SKIPPED: fewer than two caller-native DMR tables found; harmonized support matrices "
                "are excluded because they make overlap-threshold sensitivity degenerate"
            ),
        }
        write_tsv(out_dir / "overlap_threshold_sensitivity.tsv", [row])
        write_tsv(out_dir / "overlap_threshold_pairwise.tsv", [row])
        return 0

    summary: list[dict] = []
    pair_rows: list[dict] = []
    min_threshold = min(thresholds) if thresholds else 0.0
    for a, b in combinations(tables, 2):
        candidate_pairs = matched_pairs(a.data, b.data, min_threshold)
        for threshold in thresholds:
            if candidate_pairs.empty:
                all_pairs = candidate_pairs
            else:
                all_pairs = candidate_pairs[candidate_pairs["reciprocal_overlap"] >= threshold].copy()
            for policy in policies:
                pairs = apply_matching_policy(all_pairs, policy)
                matched_a = pairs["region_id_a"].nunique() if not pairs.empty else 0
                matched_b = pairs["region_id_b"].nunique() if not pairs.empty else 0
                denom = max(1, len(a.data) + len(b.data) - matched_a - matched_b + len(pairs))
                note = f"{a.note}; {b.note}"
                if policy == "many_to_many":
                    note += "; many-to-many can inflate pair counts when one interval overlaps multiple intervals"
                note += "; filtering uses strict reciprocal overlap min(O/L_a, O/L_b); min_length_overlap retained only for diagnostic comparison"
                summary.append(
                    {
                        "threshold": threshold,
                        "matching_policy": policy,
                        "caller_a": a.name,
                        "caller_b": b.name,
                        "n_a": len(a.data),
                        "n_b": len(b.data),
                        "n_pairs": len(pairs),
                        "n_matched_a": matched_a,
                        "n_matched_b": matched_b,
                        "jaccard_like": len(pairs) / denom,
                        "caller_a_source": str(a.path),
                        "caller_b_source": str(b.path),
                        "caller_a_source_kind": a.source_kind,
                        "caller_b_source_kind": b.source_kind,
                        "reciprocal_overlap_definition": "strict_min_overlap_over_each_interval",
                        "notes": note,
                    }
                )
                for row in pairs.to_dict("records"):
                    row.update({"threshold": threshold, "matching_policy": policy, "caller_a": a.name, "caller_b": b.name})
                    pair_rows.append(row)
    write_tsv(out_dir / "overlap_threshold_sensitivity.tsv", summary)
    write_tsv(out_dir / "overlap_threshold_pairwise.tsv", pair_rows or [{"status": "SKIPPED", "notes": "no overlapping pairs"}])
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
