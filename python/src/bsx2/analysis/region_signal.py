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
from typing import Literal, Optional

import numpy as np
import pandas as pd

from .seqname_harmonization import normalize_seqname


@dataclass(frozen=True)
class RegionSignalConfig:
    context: Optional[str] = None
    strand_policy: str = "both"
    min_total: int = 0
    seqname_aliases: Optional[dict[str, str]] = None
    region_id_column: str = "region_id"
    backend: Literal["auto", "rust", "pandas"] = "auto"
    methylation_path: Optional[str | Path] = None
    chunk_size: int = 10_000
    include_empty_regions: bool = True


REGION_SIGNAL_OUTPUT_COLUMNS = [
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
    if out[["start", "end"]].isna().any().any():
        raise ValueError("regions_df contains invalid start/end coordinates")
    if (out["start"] > out["end"]).any():
        raise ValueError("regions_df contains regions with start > end")
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
    out["seqname"] = out["chrom"]
    return out[["region_id", "seqname", "chrom", "start", "end"] + [c for c in ("strand", "context") if c in out.columns]]


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
    out["seqname"] = out["chrom"]
    return out


def available_region_signal_backends() -> dict[str, object]:
    try:
        from .region_signal_rust import rust_region_aggregator_available

        rust_available = rust_region_aggregator_available()
    except Exception:
        rust_available = False
    return {
        "rust": rust_available,
        "pandas": True,
        "default": "rust" if rust_available else "pandas",
    }


def _resolve_backend(
    backend: Literal["auto", "rust", "pandas"],
    methylation_path: str | Path | None,
) -> Literal["rust", "pandas"]:
    if backend not in {"auto", "rust", "pandas"}:
        raise ValueError("backend must be auto, rust, or pandas")
    if backend == "pandas":
        return "pandas"

    rust_available = bool(available_region_signal_backends()["rust"])
    if backend == "rust":
        if not rust_available:
            raise RuntimeError("backend='rust' requested, but the Rust binding is unavailable")
        if methylation_path is None:
            raise ValueError("backend='rust' requires methylation_path or a BSX path counts_df")
        return "rust"

    if rust_available and methylation_path is not None:
        return "rust"
    return "pandas"


def _normalize_strand_policy(policy: str) -> str:
    normalized = str(policy).lower()
    aliases = {
        "ignore": "both",
        "both": "both",
        "same": "region_strand",
        "region_strand": "region_strand",
        "opposite": "opposite",
        "plus": "plus",
        "+": "plus",
        "minus": "minus",
        "-": "minus",
    }
    if normalized not in aliases:
        raise ValueError("strand_policy must be both/ignore, region_strand/same, opposite, plus, or minus")
    return aliases[normalized]


def _aggregate_region_signal_pandas(
    regions_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    config: RegionSignalConfig,
) -> pd.DataFrame:
    config = config or RegionSignalConfig()
    regions = normalize_region_table(regions_df)
    counts = normalize_counts_table(counts_df)
    if config.seqname_aliases:
        regions["chrom"] = regions["chrom"].map(lambda v: normalize_seqname(v, config.seqname_aliases))
        regions["seqname"] = regions["chrom"]
        counts["chrom"] = counts["chrom"].map(lambda v: normalize_seqname(v, config.seqname_aliases))
        counts["seqname"] = counts["chrom"]
    if config.min_total > 0:
        counts = counts[counts["total"] >= config.min_total]
    strand_policy = _normalize_strand_policy(config.strand_policy)
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
        if strand_policy == "region_strand" and "strand" in region and "strand" in base.columns:
            base = base[base["strand"] == region["strand"]]
        elif strand_policy == "opposite" and "strand" in region and "strand" in base.columns:
            base = base[base["strand"] != region["strand"]]
        elif strand_policy == "plus" and "strand" in base.columns:
            base = base[base["strand"] == "+"]
        elif strand_policy == "minus" and "strand" in base.columns:
            base = base[base["strand"] == "-"]
        for sample in samples:
            sample_rows = base[base["sample_id"] == sample]
            if len(sample_rows) == 0 and not config.include_empty_regions:
                continue
            mc = float(sample_rows["mC"].sum())
            uc = float(sample_rows["uC"].sum())
            total = mc + uc
            if len(sample_rows) == 0:
                qc = "no_records"
            elif total <= 0:
                qc = "zero_coverage"
            elif total < config.min_total:
                qc = "low_coverage"
            else:
                qc = "ok"
            rows.append({
                "region_id": region["region_id"],
                "seqname": region["seqname"],
                "chrom": region["chrom"],
                "start": region["start"],
                "end": region["end"],
                "strand": region.get("strand", "."),
                "context": str(context).upper() if context is not None else "all",
                "sample_id": sample,
                "mC": mc,
                "uC": uc,
                "total": total,
                "n_cytosines": int(len(sample_rows)),
                "mean_methylation": np.nan if total <= 0 else mc / total,
                "coverage_qc": qc,
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=REGION_SIGNAL_OUTPUT_COLUMNS)
    for column in REGION_SIGNAL_OUTPUT_COLUMNS:
        if column not in out.columns:
            out[column] = pd.NA
    return out[REGION_SIGNAL_OUTPUT_COLUMNS]


def aggregate_region_signal(
    regions_df: pd.DataFrame,
    counts_df: pd.DataFrame | str | Path | None = None,
    config: RegionSignalConfig | None = None,
    *,
    backend: Literal["auto", "rust", "pandas"] | None = None,
    methylation_path: str | Path | None = None,
    sample_id: str | None = None,
) -> pd.DataFrame:
    config = config or RegionSignalConfig()
    selected_backend = backend or config.backend
    selected_methylation_path = methylation_path or config.methylation_path
    if selected_methylation_path is None and isinstance(counts_df, (str, Path)):
        selected_methylation_path = counts_df

    backend_used = _resolve_backend(selected_backend, selected_methylation_path)
    if backend_used == "rust":
        from .region_signal_rust import aggregate_region_signal_rust

        return aggregate_region_signal_rust(
            selected_methylation_path,
            normalize_region_table(regions_df),
            sample_id=sample_id,
            context=config.context,
            strand_policy=_normalize_strand_policy(config.strand_policy),
            min_total=config.min_total,
            chunk_size=config.chunk_size,
            include_empty_regions=config.include_empty_regions,
            seqname_aliases=config.seqname_aliases,
        )

    if counts_df is None or isinstance(counts_df, (str, Path)):
        raise ValueError("pandas backend requires counts_df as a pandas DataFrame")
    return _aggregate_region_signal_pandas(regions_df, counts_df, config)


def write_region_signal_table(df: pd.DataFrame, out_path: str | Path) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)


def write_region_signal_qc(df: pd.DataFrame, out_path: str | Path) -> pd.DataFrame:
    qc = df.groupby("coverage_qc", dropna=False).size().reset_index(name="n_rows") if not df.empty else pd.DataFrame(columns=["coverage_qc", "n_rows"])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(out_path, sep="\t", index=False)
    return qc
