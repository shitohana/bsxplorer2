#!/usr/bin/env python3
"""Compare center and interval-overlap metagene DMR projection strategies."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from dmr_validation_framework.core.columns import bin_cols, section_for_bin
from dmr_validation_framework.core.io import (
    ensure_out_dir,
    find_first,
    find_preferred_file,
    read_table,
    write_tsv,
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--upstream-len", type=int, default=2000)
    parser.add_argument("--downstream-len", type=int, default=2000)
    parser.add_argument("--upstream-bins", type=int, default=20)
    parser.add_argument("--body-bins", type=int, default=100)
    parser.add_argument("--downstream-bins", type=int, default=20)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def read_genes() -> pd.DataFrame:
    gene_path = find_preferred_file(
        ["metagene_gene_regions.bed", "gene_regions.bed", "genes.tsv", "gene_annotation.tsv", "genes.gff3", "annotation.gff3"]
    ) or find_first("gene_annotation")
    if not gene_path:
        return pd.DataFrame()
    if gene_path.suffix == ".bed":
        df = pd.read_csv(gene_path, sep="\t", header=None)
        df = df.iloc[:, :6]
        df.columns = ["chrom", "start", "end", "gene_id", "score", "strand"]
        df["start"] = pd.to_numeric(df["start"], errors="coerce") + 1
        df["end"] = pd.to_numeric(df["end"], errors="coerce")
        return df[["gene_id", "chrom", "start", "end", "strand"]].dropna()
    df = read_table(gene_path)
    if {"gene_id", "chrom", "start", "end", "strand"}.issubset(df.columns):
        return df[["gene_id", "chrom", "start", "end", "strand"]].dropna()
    return pd.DataFrame()


def map_pos_to_bin(pos: float, gene: pd.Series, args: argparse.Namespace) -> int | None:
    start = float(gene.start)
    end = float(gene.end)
    strand = str(gene.strand)
    if strand == "+":
        if start - args.upstream_len <= pos < start:
            rel = (pos - (start - args.upstream_len)) / max(1, args.upstream_len)
            return min(args.upstream_bins - 1, max(0, int(rel * args.upstream_bins)))
        if start <= pos <= end:
            rel = (pos - start) / max(1, end - start + 1)
            return args.upstream_bins + min(args.body_bins - 1, max(0, int(rel * args.body_bins)))
        if end < pos <= end + args.downstream_len:
            rel = (pos - end) / max(1, args.downstream_len)
            return args.upstream_bins + args.body_bins + min(args.downstream_bins - 1, max(0, int(rel * args.downstream_bins)))
    elif strand == "-":
        if end < pos <= end + args.upstream_len:
            rel = ((end + args.upstream_len) - pos) / max(1, args.upstream_len)
            return min(args.upstream_bins - 1, max(0, int(rel * args.upstream_bins)))
        if start <= pos <= end:
            rel = (end - pos) / max(1, end - start + 1)
            return args.upstream_bins + min(args.body_bins - 1, max(0, int(rel * args.body_bins)))
        if start - args.downstream_len <= pos < start:
            rel = (start - pos) / max(1, args.downstream_len)
            return args.upstream_bins + args.body_bins + min(args.downstream_bins - 1, max(0, int(rel * args.downstream_bins)))
    return None


def bin_interval(gene: pd.Series, bin_idx: int, args: argparse.Namespace) -> tuple[float, float]:
    start = float(gene.start)
    end = float(gene.end)
    strand = str(gene.strand)
    if bin_idx < args.upstream_bins:
        frac0 = bin_idx / args.upstream_bins
        frac1 = (bin_idx + 1) / args.upstream_bins
        if strand == "+":
            return start - args.upstream_len + frac0 * args.upstream_len, start - args.upstream_len + frac1 * args.upstream_len
        return end + (1 - frac1) * args.upstream_len, end + (1 - frac0) * args.upstream_len
    if bin_idx < args.upstream_bins + args.body_bins:
        i = bin_idx - args.upstream_bins
        frac0 = i / args.body_bins
        frac1 = (i + 1) / args.body_bins
        if strand == "+":
            return start + frac0 * (end - start + 1), start + frac1 * (end - start + 1)
        return end - frac1 * (end - start + 1), end - frac0 * (end - start + 1)
    i = bin_idx - args.upstream_bins - args.body_bins
    frac0 = i / args.downstream_bins
    frac1 = (i + 1) / args.downstream_bins
    if strand == "+":
        return end + frac0 * args.downstream_len, end + frac1 * args.downstream_len
    return start - frac1 * args.downstream_len, start - frac0 * args.downstream_len


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    mapped_path = find_first("mapped_dmr_centers")
    occ_path = find_first("occupancy_matrix")
    genes = read_genes()
    if not mapped_path or genes.empty:
        skipped = [{"status": "SKIPPED", "notes": "mapped_dmr_centers or gene annotation not found"}]
        write_tsv(out_dir / "projection_sensitivity_summary.tsv", skipped)
        write_tsv(out_dir / "projection_sensitivity_density.tsv", skipped)
        write_tsv(out_dir / "projection_changed_regions.tsv", skipped)
        return 0
    dmr = read_table(mapped_path)
    required = {"chrom", "start", "end", "gene_id"}
    if not required.issubset(dmr.columns):
        skipped = [{"status": "SKIPPED", "notes": "mapped DMR center table lacks required columns"}]
        write_tsv(out_dir / "projection_sensitivity_summary.tsv", skipped)
        write_tsv(out_dir / "projection_sensitivity_density.tsv", skipped)
        write_tsv(out_dir / "projection_changed_regions.tsv", skipped)
        return 0
    genes = genes.drop_duplicates("gene_id").set_index("gene_id")
    n_bins = args.upstream_bins + args.body_bins + args.downstream_bins
    if occ_path:
        occ = read_table(occ_path)
        n_bins = len(bin_cols(occ)) or n_bins
    center_density = np.zeros(n_bins)
    interval_density = np.zeros(n_bins)
    changed: list[dict] = []
    mapped = 0
    for row in dmr.itertuples(index=False):
        gene_id = str(row.gene_id)
        if gene_id not in genes.index:
            continue
        gene = genes.loc[gene_id]
        center = (float(row.start) + float(row.end)) / 2.0
        center_bin = map_pos_to_bin(center, gene, args)
        interval_bins: list[int] = []
        overlap_by_bin: list[tuple[int, float]] = []
        for b in range(n_bins):
            b0, b1 = bin_interval(gene, b, args)
            lo = max(min(b0, b1), float(row.start))
            hi = min(max(b0, b1), float(row.end))
            overlap = max(0.0, hi - lo + 1)
            if overlap > 0:
                interval_bins.append(b)
                overlap_by_bin.append((b, overlap))
        if center_bin is not None:
            center_density[center_bin] += 1
        for b in set(interval_bins):
            interval_density[b] += 1
        if center_bin is not None or interval_bins:
            mapped += 1
        dominant_interval_bin = max(overlap_by_bin, key=lambda x: x[1])[0] if overlap_by_bin else None
        center_section = section_for_bin(center_bin, args.upstream_bins, args.body_bins) if center_bin is not None else "unmapped"
        interval_section = section_for_bin(dominant_interval_bin, args.upstream_bins, args.body_bins) if dominant_interval_bin is not None else "unmapped"
        if center_bin != dominant_interval_bin or center_section != interval_section:
            changed.append(
                {
                    "gene_id": gene_id,
                    "chrom": row.chrom,
                    "start": row.start,
                    "end": row.end,
                    "center_bin": center_bin,
                    "dominant_interval_bin": dominant_interval_bin,
                    "center_section": center_section,
                    "dominant_interval_section": interval_section,
                    "n_interval_bins": len(set(interval_bins)),
                }
            )
    denom = max(1, dmr["gene_id"].nunique())
    center_density = center_density / denom
    interval_density = interval_density / denom
    corr = float(np.corrcoef(center_density, interval_density)[0, 1]) if n_bins > 1 and np.isfinite(center_density).all() else np.nan
    summary = [
        {
            "n_dmr_input": len(dmr),
            "n_mapped_any_method": mapped,
            "n_changed_regions": len(changed),
            "density_profile_correlation": corr,
            "status": "PASS",
            "notes": "interval-overlap projection is a sensitivity check; center projection remains the primary companion figure rule",
        }
    ]
    density_rows = []
    for b in range(n_bins):
        density_rows.append(
            {
                "bin": f"bin_{b}",
                "bin_index": b,
                "section": section_for_bin(b, args.upstream_bins, args.body_bins),
                "density_center": center_density[b],
                "density_interval_overlap": interval_density[b],
                "density_difference_interval_minus_center": interval_density[b] - center_density[b],
            }
        )
    write_tsv(out_dir / "projection_sensitivity_summary.tsv", summary)
    write_tsv(out_dir / "projection_sensitivity_density.tsv", density_rows)
    write_tsv(out_dir / "projection_changed_regions.tsv", changed or [{"status": "PASS", "notes": "no changed regions"}])
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
