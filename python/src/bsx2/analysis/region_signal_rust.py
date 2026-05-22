"""Rust-backed region signal aggregation wrapper.

The Rust binding works on BSX methylation files and returns pandas-friendly
rows. The public RegionSignal API keeps pandas as a fallback.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .seqname_harmonization import normalize_seqname

try:  # pragma: no cover - availability depends on the compiled extension
    from bsx2._bsx2 import aggregate_region_counts_rust as _aggregate_region_counts_rust
except Exception:  # pragma: no cover
    _aggregate_region_counts_rust = None


REGION_SIGNAL_SCHEMA = [
    "region_id",
    "seqname",
    "chrom",
    "start",
    "end",
    "strand",
    "context",
    "sample_id",
    "mC",
    "uC",
    "total",
    "n_cytosines",
    "mean_methylation",
    "coverage_qc",
]


def rust_region_aggregator_available() -> bool:
    return _aggregate_region_counts_rust is not None


def aggregate_region_signal_rust(
    methylation_path: str | Path,
    regions_df: pd.DataFrame,
    *,
    sample_id: Optional[str] = None,
    context: Optional[str] = None,
    strand_policy: str = "both",
    min_total: int = 0,
    chunk_size: int = 10_000,
    include_empty_regions: bool = True,
    seqname_aliases: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    if _aggregate_region_counts_rust is None:
        raise RuntimeError("Rust region aggregation binding is not available")
    if str(strand_policy).lower() == "opposite":
        raise ValueError("strand_policy='opposite' is not supported by Rust backend; use backend='pandas' or backend='auto'.")

    regions = regions_df.copy()
    if "seqname" not in regions.columns and "chrom" in regions.columns:
        regions["seqname"] = regions["chrom"]
    if "chrom" not in regions.columns and "seqname" in regions.columns:
        regions["chrom"] = regions["seqname"]
    if seqname_aliases:
        seq_col = "seqname" if "seqname" in regions.columns else "chrom"
        regions[seq_col] = regions[seq_col].map(
            lambda value: normalize_seqname(str(value), seqname_aliases)
        )
        regions["chrom"] = regions[seq_col]
        regions["seqname"] = regions[seq_col]

    region_rows = regions.to_dict(orient="records")
    rows = _aggregate_region_counts_rust(
        str(methylation_path),
        region_rows,
        sample_id,
        context,
        strand_policy,
        int(min_total),
        int(chunk_size),
        bool(include_empty_regions),
    )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=REGION_SIGNAL_SCHEMA)
    for column in REGION_SIGNAL_SCHEMA:
        if column not in out.columns:
            out[column] = pd.NA
    return out[REGION_SIGNAL_SCHEMA]
