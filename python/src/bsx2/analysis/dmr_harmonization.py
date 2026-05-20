"""Canonical schema and overlap helpers for external DMR harmonization.

Purpose:
    Normalize already produced internal/external DMR candidate tables into a
    shared interval schema and compare caller support by reciprocal overlap.

Limitations:
    The layer standardizes schema and interval matching only. It does not make
    caller-specific p-values or q-values statistically equivalent and does not
    execute DSS, methylKit, dmrseq, metilene, or any other external caller.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


CANONICAL_DMR_COLUMNS = [
    "dmr_id", "source_caller", "caller_version", "contrast_id", "condition_a", "condition_b",
    "context", "chrom", "start", "end", "strand", "direction", "delta", "p_value", "q_value",
    "n_sites", "n_cytosines", "mean_methylation_a", "mean_methylation_b", "source_file",
    "source_status", "method_notes",
]


def canonical_dmr_columns() -> list[str]:
    return list(CANONICAL_DMR_COLUMNS)


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def validate_canonical_dmr_schema(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for column in CANONICAL_DMR_COLUMNS:
        rows.append({"column": column, "present": column in df.columns, "status": "ok" if column in df.columns else "missing"})
    required = ["chrom", "start", "end"]
    if any(c not in df.columns for c in required):
        rows.append({"column": "required_coordinates", "present": False, "status": "missing_required_coordinate"})
    return pd.DataFrame(rows)


def normalize_dmr_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    chrom_col = _find_column(df, ("chrom", "chr", "chromosome", "seqname"))
    start_col = _find_column(df, ("start", "start_bp", "begin"))
    end_col = _find_column(df, ("end", "end_bp", "stop"))
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError("DMR table requires chrom/chr, start, and end columns")
    out = df.copy()
    out["chrom"] = out[chrom_col].astype(str)
    out["start"] = pd.to_numeric(out[start_col], errors="coerce")
    out["end"] = pd.to_numeric(out[end_col], errors="coerce")
    return out


def assign_missing_dmr_ids(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    out = df.copy()
    if "dmr_id" not in out.columns:
        out["dmr_id"] = [f"{prefix}_{i + 1}" for i in range(len(out))]
    out["dmr_id"] = out["dmr_id"].fillna(pd.Series([f"{prefix}_{i + 1}" for i in range(len(out))], index=out.index)).astype(str)
    return out


def write_schema_validation_report(df: pd.DataFrame, out_path: str | Path) -> pd.DataFrame:
    report = validate_canonical_dmr_schema(df)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out_path, sep="\t", index=False)
    return report


def _reciprocal_overlap(a_start: float, a_end: float, b_start: float, b_end: float) -> float:
    overlap = max(0.0, min(a_end, b_end) - max(a_start, b_start))
    a_len = max(float(a_end - a_start), 1.0)
    b_len = max(float(b_end - b_start), 1.0)
    return min(overlap / a_len, overlap / b_len)


def build_caller_support_matrix(canonical_tables: list[pd.DataFrame], overlap_threshold: float = 0.5) -> pd.DataFrame:
    all_rows = []
    for table in canonical_tables:
        if table.empty:
            continue
        all_rows.extend(table.to_dict("records"))
    clusters: list[dict[str, Any]] = []
    for row in all_rows:
        matched = None
        for cluster in clusters:
            same = row.get("chrom") == cluster["chrom"] and str(row.get("context", "")) == str(cluster.get("context", ""))
            if same and _reciprocal_overlap(float(row.get("start", 0)), float(row.get("end", 0)), float(cluster["start"]), float(cluster["end"])) >= overlap_threshold:
                matched = cluster
                break
        if matched is None:
            matched = {"chrom": row.get("chrom"), "start": row.get("start"), "end": row.get("end"), "context": row.get("context"), "callers": set(), "q_values": [], "deltas": [], "evidence_class": row.get("evidence_class", pd.NA)}
            clusters.append(matched)
        matched["callers"].add(str(row.get("source_caller", "external")))
        if pd.notna(row.get("q_value")):
            matched["q_values"].append(float(row["q_value"]))
        if pd.notna(row.get("delta")):
            matched["deltas"].append(float(row["delta"]))
    rows = []
    for i, cluster in enumerate(clusters, start=1):
        callers = cluster["callers"]
        rows.append({
            "harmonized_region_id": f"harmonized_{i}",
            "chrom": cluster["chrom"],
            "start": cluster["start"],
            "end": cluster["end"],
            "context": cluster["context"],
            "internal_support": "internal" in callers,
            "DSS_support": "DSS" in callers,
            "methylKit_support": "methylKit" in callers,
            "dmrseq_support": "dmrseq" in callers,
            "metilene_support": "metilene" in callers,
            "n_callers_supporting": len(callers),
            "best_q_value": min(cluster["q_values"]) if cluster["q_values"] else pd.NA,
            "max_abs_delta": max(abs(v) for v in cluster["deltas"]) if cluster["deltas"] else pd.NA,
            "evidence_class": cluster["evidence_class"],
        })
    return pd.DataFrame(rows)
