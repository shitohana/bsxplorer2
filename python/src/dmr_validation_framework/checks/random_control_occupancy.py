#!/usr/bin/env python3
"""Random-control occupancy audit for DMR metagene density profiles."""

from __future__ import annotations

import argparse
import os
import random
from pathlib import Path

import numpy as np
import pandas as pd

from dmr_validation_framework.core.columns import bin_cols, section_for_bin
from dmr_validation_framework.core.io import (
    SKIP_DIRS,
    default_search_roots,
    ensure_out_dir,
    find_first,
    find_preferred_file,
    read_table,
    write_tsv,
)
from dmr_validation_framework.core.stats import bh_qvalues


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/validation_audit"))
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=202715)
    parser.add_argument("--max-regions-per-iteration", type=int, default=1000)
    parser.add_argument("--upstream-len", type=int, default=2000)
    parser.add_argument("--downstream-len", type=int, default=2000)
    parser.add_argument("--upstream-bins", type=int, default=20)
    parser.add_argument("--body-bins", type=int, default=100)
    parser.add_argument("--downstream-bins", type=int, default=20)
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def find_fai() -> Path | None:
    for root in default_search_roots():
        if not root.exists():
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            current = Path(dirpath)
            try:
                rel_depth = len(current.relative_to(root).parts)
            except ValueError:
                rel_depth = 0
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".")]
            if rel_depth >= 6:
                dirnames[:] = []
            for filename in filenames:
                if filename.endswith(".fai") or filename.endswith(".fa.fai"):
                    return current / filename
    return None


def read_fai(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return sizes


def read_genes() -> pd.DataFrame:
    gene_path = find_preferred_file(
        ["metagene_gene_regions.bed", "gene_regions.bed", "genes.tsv", "gene_annotation.tsv", "genes.gff3", "annotation.gff3"]
    ) or find_first("gene_annotation")
    if not gene_path:
        return pd.DataFrame()
    if gene_path.suffix == ".bed":
        df = pd.read_csv(gene_path, sep="\t", header=None).iloc[:, :6]
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
    if strand == "-":
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


def build_gene_index(genes: pd.DataFrame, args: argparse.Namespace):
    index = {}
    max_window = 0
    for chrom, sub in genes.groupby("chrom"):
        sub = sub.copy()
        sub["window_start"] = sub["start"] - args.downstream_len
        sub["window_end"] = sub["end"] + args.upstream_len
        sub = sub.sort_values("window_start").reset_index(drop=True)
        max_window = max(max_window, int((sub["window_end"] - sub["window_start"]).max()))
        index[str(chrom)] = (sub, sub["window_start"].to_numpy())
    return index, max_window


def project_center(chrom: str, center: float, gene_index, max_window: int, args: argparse.Namespace):
    if chrom not in gene_index:
        return None
    sub, starts = gene_index[chrom]
    lo = np.searchsorted(starts, center - max_window - 1, side="left")
    hi = np.searchsorted(starts, center, side="right")
    best = None
    best_dist = None
    for gene in sub.iloc[lo:hi].itertuples(index=False):
        if not (gene.window_start <= center <= gene.window_end):
            continue
        bin_index = map_pos_to_bin(center, gene, args)
        if bin_index is None:
            continue
        tss = float(gene.start if gene.strand == "+" else gene.end)
        dist = abs(center - tss)
        if best is None or dist < best_dist:
            best = (str(gene.gene_id), bin_index)
            best_dist = dist
    return best


def run(args: argparse.Namespace) -> int:
    out_dir = ensure_out_dir(args.out_dir)
    mapped_path = find_first("mapped_dmr_centers")
    occ_path = find_first("occupancy_matrix")
    fai_path = find_fai()
    genes = read_genes()
    if not mapped_path or not occ_path or not fai_path or genes.empty:
        skipped = [
            {
                "status": "SKIPPED",
                "notes": "random control requires mapped DMR centers, occupancy matrix, gene annotation and genome sizes/faidx; at least one input was not found",
                "mapped_dmr_centers": str(mapped_path) if mapped_path else "missing",
                "occupancy_matrix": str(occ_path) if occ_path else "missing",
                "genome_fai": str(fai_path) if fai_path else "missing",
                "genes": "found" if not genes.empty else "missing",
            }
        ]
        write_tsv(out_dir / "random_control_occupancy_summary.tsv", skipped)
        write_tsv(out_dir / "random_control_density_by_section.tsv", skipped)
        write_tsv(out_dir / "random_control_bin_empirical_p.tsv", skipped)
        return 0

    occ = read_table(occ_path)
    bins = bin_cols(occ)
    n_bins = len(bins)
    n_genes = len(occ)
    observed_values = (occ[bins].apply(pd.to_numeric, errors="coerce").fillna(0) > 0).astype(int)
    observed_density = observed_values.sum(axis=0).to_numpy(dtype=float) / max(1, n_genes)
    sizes = read_fai(fai_path)
    gene_index, max_window = build_gene_index(genes, args)
    observed = read_table(mapped_path)
    observed["length"] = pd.to_numeric(observed["end"], errors="coerce") - pd.to_numeric(observed["start"], errors="coerce") + 1
    observed = observed[observed["chrom"].astype(str).isin(sizes)].dropna(subset=["length"]).copy()
    if len(observed) > args.max_regions_per_iteration:
        observed = observed.sample(args.max_regions_per_iteration, random_state=args.seed)
        sampling_note = f"observed DMR lengths downsampled to {args.max_regions_per_iteration} per iteration for lightweight audit"
    else:
        sampling_note = "all observed DMR lengths used"

    rng = random.Random(args.seed)
    ge_counts = np.zeros(n_bins, dtype=int)
    random_global_max: list[float] = []
    random_section_rows: list[dict] = []
    section_indices = {
        section: [i for i in range(n_bins) if section_for_bin(i, args.upstream_bins, args.body_bins) == section]
        for section in ["upstream", "gene_body", "downstream"]
    }
    for iteration in range(args.iterations):
        occupied: set[tuple[str, int]] = set()
        for row in observed.itertuples(index=False):
            chrom = str(row.chrom)
            chrom_size = sizes.get(chrom)
            length = int(max(1, row.length))
            if not chrom_size or chrom_size <= length + 2:
                continue
            start = rng.randint(1, chrom_size - length)
            center = start + (length - 1) / 2
            mapped = project_center(chrom, center, gene_index, max_window, args)
            if mapped is not None:
                occupied.add(mapped)
        random_density = np.zeros(n_bins, dtype=float)
        for _, bin_index in occupied:
            random_density[bin_index] += 1
        random_density = random_density / max(1, n_genes)
        ge_counts += random_density >= observed_density
        random_global_max.append(float(random_density.max()) if len(random_density) else 0.0)
        for section, idx in section_indices.items():
            random_section_rows.append(
                {
                    "iteration": iteration,
                    "section": section,
                    "random_mean_density": float(random_density[idx].mean()) if idx else 0.0,
                    "n_projected_gene_bins": len(occupied),
                }
            )

    empirical_p = (1 + ge_counts) / (args.iterations + 1)
    empirical_q = bh_qvalues(empirical_p)
    observed_global_max_density = float(observed_density.max()) if len(observed_density) else 0.0
    global_profile_p = (1 + sum(x >= observed_global_max_density for x in random_global_max)) / (args.iterations + 1)
    summary = [
        {
            "status": "PASS" if sampling_note.startswith("all") else "WARN",
            "n_observed_dmr_lengths_used": len(observed),
            "iterations_requested": args.iterations,
            "n_genes_denominator": n_genes,
            "genome_fai": str(fai_path),
            "observed_global_max_density": observed_global_max_density,
            "global_profile_p": global_profile_p,
            "multiple_testing_note": "bin-level empirical p-values are BH-corrected; global p uses max-density statistic across bins",
            "notes": sampling_note,
        }
    ]
    bin_rows = [
        {
            "bin": f"bin_{i}",
            "bin_index": i,
            "section": section_for_bin(i, args.upstream_bins, args.body_bins),
            "observed_density": float(observed_density[i]),
            "empirical_p_random_ge_observed": float(empirical_p[i]),
            "empirical_q_BH": float(empirical_q[i]) if not pd.isna(empirical_q[i]) else np.nan,
            "global_profile_p": global_profile_p,
            "global_statistic": "max_density",
            "observed_global_max_density": observed_global_max_density,
            "iterations": args.iterations,
        }
        for i in range(n_bins)
    ]
    write_tsv(out_dir / "random_control_occupancy_summary.tsv", summary)
    write_tsv(out_dir / "random_control_density_by_section.tsv", random_section_rows)
    write_tsv(out_dir / "random_control_bin_empirical_p.tsv", bin_rows)
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
