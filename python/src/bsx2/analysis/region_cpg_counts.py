"""Per-CpG count extraction for confirmatory DMR validation.

This module extracts methylation observations inside predefined regions.  It is
not a DMR caller and does not perform statistical testing.  The Rust path uses
indexed BSX region queries when the compiled binding is available; the pandas
path is kept for tests and small synthetic fixtures.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

import pandas as pd

from .region_signal import normalize_counts_table, normalize_region_table
from .seqname_harmonization import normalize_seqname

try:  # pragma: no cover - depends on compiled extension availability
    from bsx2._bsx2 import extract_region_cpg_counts_rust as _extract_region_cpg_counts_rust
except Exception:  # pragma: no cover
    _extract_region_cpg_counts_rust = None


REGION_CPG_COUNTS_COLUMNS = [
    "region_id",
    "cpg_id",
    "seqname",
    "chrom",
    "position",
    "strand",
    "context",
    "sample_id",
    "mC",
    "uC",
    "total",
    "coverage_qc",
]


def rust_region_cpg_extractor_available() -> bool:
    return _extract_region_cpg_counts_rust is not None


def _normalize_context(context: Optional[str]) -> Optional[str]:
    if context is None or str(context).lower() == "all":
        return None
    normalized = str(context).upper()
    if normalized not in {"CG", "CHG", "CHH"}:
        raise ValueError("context must be CG, CHG, CHH, all, or None")
    return normalized


def _normalize_strand_policy(policy: str) -> str:
    aliases = {
        "both": "both",
        "ignore": "both",
        "region_strand": "region_strand",
        "same": "region_strand",
        "plus": "plus",
        "+": "plus",
        "minus": "minus",
        "-": "minus",
    }
    key = str(policy).lower()
    if key not in aliases:
        raise ValueError("strand_policy must be both/ignore, region_strand/same, plus, or minus")
    return aliases[key]


def normalize_cpg_regions(
    regions_df: pd.DataFrame,
    *,
    top_n: Optional[int] = None,
    seqname_aliases: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    regions = normalize_region_table(regions_df)
    if seqname_aliases:
        regions["seqname"] = regions["seqname"].map(lambda v: normalize_seqname(str(v), seqname_aliases))
        regions["chrom"] = regions["seqname"]
    if top_n is not None:
        regions = regions.head(int(top_n)).copy()
    return regions


def _coverage_qc(total: float, min_total: int) -> str:
    if total <= 0:
        return "zero_coverage"
    if total < min_total:
        return "low_coverage"
    return "ok"


def extract_region_cpg_counts_pandas(
    counts_df: pd.DataFrame,
    regions_df: pd.DataFrame,
    *,
    sample_id: Optional[str] = None,
    context: Optional[str] = None,
    strand_policy: str = "both",
    min_total: int = 0,
    top_n: Optional[int] = None,
    seqname_aliases: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    regions = normalize_cpg_regions(regions_df, top_n=top_n, seqname_aliases=seqname_aliases)
    counts = normalize_counts_table(counts_df)
    if seqname_aliases:
        counts["seqname"] = counts["seqname"].map(lambda v: normalize_seqname(str(v), seqname_aliases))
        counts["chrom"] = counts["seqname"]
    context = _normalize_context(context)
    strand_policy = _normalize_strand_policy(strand_policy)
    if sample_id is not None:
        counts = counts[counts["sample_id"].astype(str) == str(sample_id)]
    if context is not None and "context" in counts.columns:
        counts = counts[counts["context"].astype(str).str.upper() == context]
    if min_total > 0:
        counts = counts[counts["total"] >= min_total]

    rows: list[dict[str, object]] = []
    for _, region in regions.iterrows():
        base = counts[
            (counts["seqname"] == region["seqname"])
            & (counts["position"] >= region["start"])
            & (counts["position"] <= region["end"])
        ]
        region_context = context or region.get("context", None)
        if region_context is not None and "context" in base.columns:
            base = base[base["context"].astype(str).str.upper() == str(region_context).upper()]
        if strand_policy == "region_strand" and "strand" in region and "strand" in base.columns:
            base = base[base["strand"].astype(str) == str(region["strand"])]
        elif strand_policy == "plus" and "strand" in base.columns:
            base = base[base["strand"].astype(str).isin(["+", "plus", "Forward"])]
        elif strand_policy == "minus" and "strand" in base.columns:
            base = base[base["strand"].astype(str).isin(["-", "minus", "Reverse"])]

        for _, record in base.iterrows():
            total = float(record["total"])
            strand = str(record.get("strand", "."))
            record_context = str(record.get("context", region_context or "NA")).upper()
            seqname = str(record["seqname"])
            position = int(record["position"])
            rows.append({
                "region_id": region["region_id"],
                "cpg_id": f"{seqname}:{position}:{strand}:{record_context}",
                "seqname": seqname,
                "chrom": seqname,
                "position": position,
                "strand": strand,
                "context": record_context,
                "sample_id": str(record.get("sample_id", sample_id or "")),
                "mC": float(record["mC"]),
                "uC": float(record["uC"]),
                "total": total,
                "coverage_qc": _coverage_qc(total, int(min_total)),
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=REGION_CPG_COUNTS_COLUMNS)
    for column in REGION_CPG_COUNTS_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    return out[REGION_CPG_COUNTS_COLUMNS]


def extract_region_cpg_counts_rust(
    methylation_path: str | Path,
    regions_df: pd.DataFrame,
    *,
    sample_id: Optional[str] = None,
    context: Optional[str] = None,
    strand_policy: str = "both",
    min_total: int = 0,
    chunk_size: int = 10_000,
    top_n: Optional[int] = None,
    seqname_aliases: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    if _extract_region_cpg_counts_rust is None:
        raise RuntimeError("Rust per-CpG extraction binding is not available")
    context = _normalize_context(context)
    strand_policy = _normalize_strand_policy(strand_policy)
    regions = normalize_cpg_regions(regions_df, top_n=top_n, seqname_aliases=seqname_aliases)
    rows = _extract_region_cpg_counts_rust(
        str(methylation_path),
        regions.to_dict(orient="records"),
        sample_id,
        context,
        strand_policy,
        int(min_total),
        int(chunk_size),
    )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=REGION_CPG_COUNTS_COLUMNS)
    for column in REGION_CPG_COUNTS_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    return out[REGION_CPG_COUNTS_COLUMNS]


def extract_region_cpg_counts(
    methylation: str | Path | pd.DataFrame,
    regions_df: pd.DataFrame,
    *,
    sample_id: Optional[str] = None,
    context: Optional[str] = None,
    strand_policy: str = "both",
    min_total: int = 0,
    chunk_size: int = 10_000,
    top_n: Optional[int] = None,
    seqname_aliases: Optional[dict[str, str]] = None,
    backend: Literal["auto", "rust", "pandas"] = "auto",
) -> pd.DataFrame:
    if backend not in {"auto", "rust", "pandas"}:
        raise ValueError("backend must be auto, rust, or pandas")
    if backend == "pandas" or isinstance(methylation, pd.DataFrame):
        if not isinstance(methylation, pd.DataFrame):
            methylation = pd.read_csv(methylation, sep=None, engine="python")
        return extract_region_cpg_counts_pandas(
            methylation,
            regions_df,
            sample_id=sample_id,
            context=context,
            strand_policy=strand_policy,
            min_total=min_total,
            top_n=top_n,
            seqname_aliases=seqname_aliases,
        )
    if backend == "rust" or rust_region_cpg_extractor_available():
        return extract_region_cpg_counts_rust(
            methylation,
            regions_df,
            sample_id=sample_id,
            context=context,
            strand_policy=strand_policy,
            min_total=min_total,
            chunk_size=chunk_size,
            top_n=top_n,
            seqname_aliases=seqname_aliases,
        )
    raise RuntimeError("backend='auto' requires Rust binding for .bsx input or a pandas DataFrame/table input")


def region_cpg_missing_region_warnings(
    regions_df: pd.DataFrame,
    cpg_counts_df: pd.DataFrame,
) -> pd.DataFrame:
    regions = normalize_cpg_regions(regions_df)
    observed = set(cpg_counts_df["region_id"].astype(str)) if "region_id" in cpg_counts_df.columns else set()
    rows = [
        {
            "warning_type": "region_no_cpg_records",
            "region_id": region_id,
            "severity": "warning",
            "message": "No CpG-level records were extracted for this predefined region.",
        }
        for region_id in regions["region_id"].astype(str)
        if region_id not in observed
    ]
    return pd.DataFrame(rows, columns=["warning_type", "region_id", "severity", "message"])


def summarize_region_cpg_counts(cpg_counts_df: pd.DataFrame) -> dict[str, object]:
    return {
        "n_rows": int(len(cpg_counts_df)),
        "n_regions": int(cpg_counts_df["region_id"].nunique()) if "region_id" in cpg_counts_df.columns else 0,
        "n_cpg": int(cpg_counts_df["cpg_id"].nunique()) if "cpg_id" in cpg_counts_df.columns else 0,
        "n_samples": int(cpg_counts_df["sample_id"].nunique()) if "sample_id" in cpg_counts_df.columns else 0,
        "n_zero_coverage": int((cpg_counts_df.get("coverage_qc", pd.Series(dtype=str)) == "zero_coverage").sum()),
    }
