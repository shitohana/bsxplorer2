#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python" / "src"))

from bsx2.analysis.region_cpg_counts import (
    extract_region_cpg_counts,
    region_cpg_missing_region_warnings,
    rust_region_cpg_extractor_available,
    summarize_region_cpg_counts,
)


def sha256_file(path: Path) -> str | None:
    if not path.exists() or path.is_dir():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract per-CpG methylation counts for predefined regions."
    )
    parser.add_argument("--methylation-path", required=True, help=".bsx file for Rust backend or TSV/CSV counts table for pandas backend")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--sample-id")
    parser.add_argument("--context", default=None)
    parser.add_argument("--strand-policy", default="both", choices=["both", "ignore", "same", "region_strand", "plus", "minus", "+", "-"])
    parser.add_argument("--min-total", type=int, default=0)
    parser.add_argument("--chunk-size", type=int, default=10_000)
    parser.add_argument("--top-n", type=int)
    parser.add_argument("--backend", default="auto", choices=["auto", "rust", "pandas"])
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    regions = pd.read_csv(args.regions, sep=None, engine="python")
    methylation: str | Path | pd.DataFrame = Path(args.methylation_path)
    backend_used = args.backend
    if args.backend == "pandas":
        methylation = pd.read_csv(args.methylation_path, sep=None, engine="python")
        backend_used = "pandas"
    elif args.backend == "auto":
        if Path(args.methylation_path).suffix == ".bsx":
            backend_used = "rust" if rust_region_cpg_extractor_available() else "unavailable"
        else:
            methylation = pd.read_csv(args.methylation_path, sep=None, engine="python")
            backend_used = "pandas"
    else:
        backend_used = "rust"

    rows = extract_region_cpg_counts(
        methylation,
        regions,
        sample_id=args.sample_id,
        context=args.context,
        strand_policy=args.strand_policy,
        min_total=args.min_total,
        chunk_size=args.chunk_size,
        top_n=args.top_n,
        backend=args.backend,
    )
    warnings = region_cpg_missing_region_warnings(
        regions.head(args.top_n) if args.top_n is not None else regions,
        rows,
    )

    counts_path = out_dir / "region_cpg_counts.tsv"
    summary_path = out_dir / "region_cpg_counts_summary.md"
    manifest_path = out_dir / "region_cpg_counts_manifest.json"
    warnings_path = out_dir / "warnings.tsv"
    rows.to_csv(counts_path, sep="\t", index=False)
    warnings.to_csv(warnings_path, sep="\t", index=False)

    summary = summarize_region_cpg_counts(rows)
    summary_path.write_text(
        "\n".join([
            "# Region CpG Counts Extraction Summary",
            "",
            f"- backend_used: {backend_used}",
            f"- methylation_path: {args.methylation_path}",
            f"- regions_path: {args.regions}",
            f"- n_requested_regions: {len(regions.head(args.top_n)) if args.top_n is not None else len(regions)}",
            f"- n_output_rows: {summary['n_rows']}",
            f"- n_regions_with_records: {summary['n_regions']}",
            f"- n_unique_cpg: {summary['n_cpg']}",
            f"- n_samples: {summary['n_samples']}",
            f"- warnings: {len(warnings)}",
            "",
            "This extraction is a predefined-region confirmatory data layer, not a genome-wide DMR caller.",
            "",
        ]),
        encoding="utf-8",
    )
    output_files = [counts_path, summary_path, warnings_path]
    manifest = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "methylation_path": str(args.methylation_path),
            "regions": str(args.regions),
        },
        "parameters": {
            "sample_id": args.sample_id,
            "context": args.context,
            "strand_policy": args.strand_policy,
            "min_total": args.min_total,
            "chunk_size": args.chunk_size,
            "top_n": args.top_n,
            "backend": args.backend,
        },
        "backend_used": backend_used,
        "outputs": {path.name: {"path": str(path), "sha256": sha256_file(path)} for path in output_files},
        "warnings": warnings.to_dict(orient="records"),
        "note": "No raw/Bismark processing and no DMR statistical model changes were performed.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
