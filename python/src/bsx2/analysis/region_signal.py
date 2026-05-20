"""Unified region-level methylation signal aggregation.

Purpose:
    Aggregate methylated/unmethylated counts over arbitrary interval regions,
    including genes, promoters, DMRs, BED intervals, and external caller
    candidates.

Input assumptions:
    Regions are interval-like tables with chromosome/seqname, start, and end
    aliases. Counts are point methylation calls with sample_id, position, and
    methylated plus unmethylated or total counts.

Limitations:
    The MVP implementation is pandas-based. Production-scale WGBS can later use
    a chunked or Rust backend with the same public schema.

Stability:
    Stable additive API. It does not replace the existing Regional Evidence
    Model; it makes the interval aggregation contract reusable.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .seqname_harmonization import normalize_seqname


@dataclass(frozen=True)
class RegionSignalConfig:
    context: Optional[str] = None
    strand_policy: str = "ignore"
    min_total: int = 0
    seqname_aliases: Optional[dict[str, str]] = None
    region_id_column: str = "region_id"


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def normalize_region_table(regions_df: pd.DataFrame) -> pd.DataFrame:
    chrom_col = _find_column(regions_df, ("chrom", "chr", "chromosome", "seqname"))
    start_col = _find_column(regions_df, ("start", "start_bp", "begin"))
    end_col = _find_column(regions_df, ("end", "end_bp", "stop"))
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError("regions_df must contain chrom/seqname, start, and end columns")
    out = pd.DataFrame({
        "chrom": regions_df[chrom_col].astype(str),
        "start": pd.to_numeric(regions_df[start_col], errors="coerce"),
        "end": pd.to_numeric(regions_df[end_col], errors="coerce"),
    })
    if "region_id" in regions_df.columns:
        out["region_id"] = regions_df["region_id"].astype(str)
    elif "name" in regions_df.columns:
        out["region_id"] = regions_df["name"].astype(str)
    else:
        out["region_id"] = [f"region_{i + 1}" for i in range(len(out))]
    if "strand" in regions_df.columns:
        out["strand"] = regions_df["strand"].astype(str)
    if "context" in regions_df.columns:
        out["context"] = regions_df["context"].astype(str).str.upper()
    return out[["region_id", "chrom", "start", "end"] + [c for c in ("strand", "context") if c in out.columns]]


def normalize_counts_table(counts_df: pd.DataFrame) -> pd.DataFrame:
    chrom_col = _find_column(counts_df, ("chrom", "chr", "chromosome", "seqname"))
    pos_col = _find_column(counts_df, ("position", "pos"))
    sample_col = _find_column(counts_df, ("sample_id", "sample"))
    mc_col = _find_column(counts_df, ("mC", "methylated_count", "count_m"))
    uc_col = _find_column(counts_df, ("uC", "unmethylated_count"))
    total_col = _find_column(counts_df, ("total", "coverage", "count_total"))
    if chrom_col is None or pos_col is None or sample_col is None or mc_col is None:
        raise ValueError("counts_df must contain chrom/seqname, position, sample_id, and methylated count")
    out = pd.DataFrame({
        "chrom": counts_df[chrom_col].astype(str),
        "position": pd.to_numeric(counts_df[pos_col], errors="coerce"),
        "sample_id": counts_df[sample_col].astype(str),
        "mC": pd.to_numeric(counts_df[mc_col], errors="coerce").fillna(0),
    })
    if uc_col is not None:
        out["uC"] = pd.to_numeric(counts_df[uc_col], errors="coerce").fillna(0)
        out["total"] = out["mC"] + out["uC"]
    elif total_col is not None:
        out["total"] = pd.to_numeric(counts_df[total_col], errors="coerce").fillna(0)
        out["uC"] = out["total"] - out["mC"]
    else:
        raise ValueError("counts_df must contain uC/unmethylated_count or total/coverage")
    if (out[["mC", "uC", "total"]] < 0).any().any():
        raise ValueError("counts_df contains negative counts")
    if "context" in counts_df.columns:
        out["context"] = counts_df["context"].astype(str).str.upper()
    if "strand" in counts_df.columns:
        out["strand"] = counts_df["strand"].astype(str)
    return out


def aggregate_region_signal(regions_df: pd.DataFrame, counts_df: pd.DataFrame, config: RegionSignalConfig | None = None) -> pd.DataFrame:
    config = config or RegionSignalConfig()
    regions = normalize_region_table(regions_df)
    counts = normalize_counts_table(counts_df)
    if config.seqname_aliases:
        regions["chrom"] = regions["chrom"].map(lambda v: normalize_seqname(v, config.seqname_aliases))
        counts["chrom"] = counts["chrom"].map(lambda v: normalize_seqname(v, config.seqname_aliases))
    rows: list[dict[str, object]] = []
    samples = sorted(counts["sample_id"].dropna().unique())
    for _, region in regions.iterrows():
        base = counts[
            (counts["chrom"] == region["chrom"])
            & (counts["position"] >= region["start"])
            & (counts["position"] <= region["end"])
        ]
        context = config.context or region.get("context", None)
        if context is not None and "context" in base.columns:
            base = base[base["context"].str.upper() == str(context).upper()]
        if config.strand_policy == "same" and "strand" in region and "strand" in base.columns:
            base = base[base["strand"] == region["strand"]]
        elif config.strand_policy == "opposite" and "strand" in region and "strand" in base.columns:
            base = base[base["strand"] != region["strand"]]
        elif config.strand_policy != "ignore":
            raise ValueError("strand_policy must be ignore, same, or opposite")
        for sample in samples:
            sample_rows = base[base["sample_id"] == sample]
            mc = float(sample_rows["mC"].sum())
            uc = float(sample_rows["uC"].sum())
            total = mc + uc
            if len(sample_rows) == 0:
                qc = "no_cytosines"
            elif total <= 0:
                qc = "zero_coverage"
            elif total < config.min_total:
                qc = "below_min_total"
            else:
                qc = "ok"
            rows.append({
                "region_id": region["region_id"],
                "sample_id": sample,
                "chrom": region["chrom"],
                "start": region["start"],
                "end": region["end"],
                "mC": mc,
                "uC": uc,
                "total": total,
                "n_cytosines": int(len(sample_rows)),
                "mean_methylation": np.nan if total <= 0 else mc / total,
                "coverage_qc": qc,
            })
    return pd.DataFrame(rows)


def write_region_signal_table(df: pd.DataFrame, out_path: str | Path) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)


def write_region_signal_qc(df: pd.DataFrame, out_path: str | Path) -> pd.DataFrame:
    qc = df.groupby("coverage_qc", dropna=False).size().reset_index(name="n_rows") if not df.empty else pd.DataFrame(columns=["coverage_qc", "n_rows"])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(out_path, sep="\t", index=False)
    return qc
