#!/usr/bin/env python
from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis import (
    RegionSignalConfig,
    aggregate_region_signal,
    available_region_signal_backends,
    read_seqname_aliases,
    write_region_signal_qc,
    write_region_signal_table,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate point methylation counts over interval regions.")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--summary-out")
    parser.add_argument("--context")
    parser.add_argument("--aliases")
    parser.add_argument(
        "--strand-policy",
        default="both",
        choices=["both", "ignore", "same", "region_strand", "opposite", "plus", "minus", "+", "-"],
    )
    parser.add_argument("--backend", default="auto", choices=["auto", "rust", "pandas"])
    parser.add_argument("--chunk-size", type=int, default=10000)
    parser.add_argument("--sample-id")
    empty = parser.add_mutually_exclusive_group()
    empty.add_argument("--include-empty-regions", dest="include_empty_regions", action="store_true", default=True)
    empty.add_argument("--drop-empty-regions", dest="include_empty_regions", action="store_false")
    args = parser.parse_args()
    aliases = read_seqname_aliases(args.aliases) if args.aliases else None
    regions = pd.read_csv(args.regions, sep=None, engine="python")
    counts_path = Path(args.counts)
    backend_info = available_region_signal_backends()
    backend_used = args.backend
    warnings: list[str] = []
    counts: pd.DataFrame | str
    if args.backend == "rust" or (args.backend == "auto" and counts_path.suffix == ".bsx" and backend_info["rust"]):
        counts = str(counts_path)
        backend_used = "rust"
    else:
        if args.backend == "auto" and counts_path.suffix == ".bsx" and not backend_info["rust"]:
            raise SystemExit("backend='auto' received a .bsx file, but the Rust binding is unavailable")
        if args.backend == "rust":
            warnings.append("Rust backend requested but unavailable; command will fail in API layer.")
        counts = pd.read_csv(args.counts, sep=None, engine="python")
        backend_used = "pandas"

    started = time.perf_counter()
    result = aggregate_region_signal(
        regions,
        counts,
        RegionSignalConfig(
            context=args.context,
            strand_policy=args.strand_policy,
            seqname_aliases=aliases,
            backend=args.backend,
            chunk_size=args.chunk_size,
            include_empty_regions=args.include_empty_regions,
        ),
        sample_id=args.sample_id,
    )
    runtime_seconds = time.perf_counter() - started
    write_region_signal_table(result, args.out)
    write_region_signal_qc(result, args.qc_out)
    n_zero = int((result["coverage_qc"].isin(["zero_coverage", "no_records"])).sum()) if not result.empty else 0
    if args.summary_out:
        summary_path = Path(args.summary_out)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            "\n".join([
                "# Region Signal Aggregation Summary",
                "",
                f"- backend_used: {backend_used}",
                f"- n_regions: {len(regions)}",
                f"- n_output_rows: {len(result)}",
                f"- n_zero_coverage_regions: {n_zero}",
                f"- runtime_seconds: {runtime_seconds:.6f}",
                f"- warnings: {'; '.join(warnings) if warnings else 'none'}",
                "",
            ]),
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
