#!/usr/bin/env python3
"""Check DMR occupancy density sensitivity to denominator choice."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from dmr_validation_framework.core.columns import bin_cols, section_for_bin
from dmr_validation_framework.core.io import (
    ensure_out_dir,
    find_first,
    read_table,
    write_tsv,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--upstream-bins", type=int, default=20)
    parser.add_argument("--body-bins", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=20)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    occ_path = find_first("occupancy_matrix")
    mapped_path = find_first("mapped_dmr_centers")
    if not occ_path:
        skipped = [{"status": "SKIPPED", "notes": "occupancy_matrix_* input not found"}]
        write_tsv(out_dir / "density_denominator_sensitivity.tsv", skipped)
        write_tsv(out_dir / "density_by_section.tsv", skipped)
        write_tsv(out_dir / "top_density_bins.tsv", skipped)
        return 0

    occ = read_table(occ_path)
    bins = bin_cols(occ)
    values = occ[bins].apply(pd.to_numeric, errors="coerce").fillna(0).clip(lower=0)
    values = (values > 0).astype(int)
    n_all = len(values)
    linked_mask = values.sum(axis=1) > 0
    n_linked = int(linked_mask.sum())
    denom_rows: list[dict] = []
    density_tables = {
        "all_genes": values.sum(axis=0) / max(1, n_all),
        "dmr_linked_genes": values[linked_mask].sum(axis=0) / max(1, n_linked),
    }
    for denom, density in density_tables.items():
        denom_rows.append(
            {
                "denominator": denom,
                "input_file": str(occ_path),
                "N_all_genes": n_all,
                "N_dmr_linked_genes": n_linked,
                "mean_density": float(density.mean()),
                "max_density": float(density.max()),
                "occupied_bins": int((density > 0).sum()),
                "status": "PASS",
                "notes": "density_b = sum_g O[g,b] / denominator",
            }
        )

    section_rows: list[dict] = []
    top_rows: list[dict] = []
    for denom, density in density_tables.items():
        for idx, (bin_name, value) in enumerate(density.items()):
            section = section_for_bin(idx, args.upstream_bins, args.body_bins)
            top_rows.append(
                {
                    "denominator": denom,
                    "bin": bin_name,
                    "bin_index": idx,
                    "section": section,
                    "density": float(value),
                }
            )
        top_df = pd.DataFrame(top_rows)
        sub = top_df[top_df["denominator"] == denom]
        for section, sec_df in sub.groupby("section"):
            section_rows.append(
                {
                    "denominator": denom,
                    "section": section,
                    "N_all_genes": n_all,
                    "N_dmr_linked_genes": n_linked,
                    "mean_density": float(sec_df["density"].mean()),
                    "max_density": float(sec_df["density"].max()),
                    "n_bins": len(sec_df),
                    "status": "PASS",
                    "notes": "section inferred from bin index",
                }
            )
    top_rows = sorted(top_rows, key=lambda r: r["density"], reverse=True)[: args.top_n]
    if mapped_path:
        denom_rows.append(
            {
                "denominator": "mapped_dmr_centers",
                "input_file": str(mapped_path),
                "N_all_genes": n_all,
                "N_dmr_linked_genes": n_linked,
                "mean_density": "NA",
                "max_density": "NA",
                "occupied_bins": "NA",
                "status": "PASS",
                "notes": "mapped centers available for section-level cross-check",
            }
        )
    write_tsv(out_dir / "density_denominator_sensitivity.tsv", denom_rows)
    write_tsv(out_dir / "density_by_section.tsv", section_rows)
    write_tsv(out_dir / "top_density_bins.tsv", top_rows)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
