"""Safe DMR visualization table loaders and column normalization.

This module is intentionally small and dependency-light. It accepts existing
BSX2 or harmonized external DMR TSV outputs and adds stable normalized columns
where possible while preserving the original input columns.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "chrom": ("chrom", "chr", "chromosome", "seqname"),
    "start": ("start", "begin", "start_bp"),
    "end": ("end", "stop", "end_bp"),
    "region_id": ("region_id", "dmr_id", "id", "harmonized_region_id"),
    "q_value": ("q_value", "q", "fdr", "padj", "qvalue", "qval", "region_q_value"),
    "p_value": ("p_value", "p", "pval", "pvalue", "region_p_value"),
    "delta": ("delta", "mean_delta", "methylation_delta", "diff", "region_delta", "beta_binom_delta"),
    "context": ("context", "methylation_context"),
    "evidence_class": ("evidence_class", "class", "support_class"),
    "Y": ("Y", "mC", "methylated", "methylated_count", "count_m"),
    "m": ("m", "total", "coverage", "count_total"),
    "sample_id": ("sample_id", "sample", "sample_name"),
    "condition": ("condition", "group", "treatment"),
}


def _key(column: object) -> str:
    return str(column).strip().lower().replace("-", "_").replace(".", "_").replace(" ", "_")


def normalize_dmr_curve_columns(df: pd.DataFrame, require_coordinates: bool = True) -> tuple[pd.DataFrame, list[str]]:
    """Return a copy with stable DMR visualization aliases added where possible."""

    out = df.copy()
    warnings: list[str] = []
    keyed = {_key(col): col for col in out.columns}
    for canonical, aliases in COLUMN_ALIASES.items():
        if canonical in out.columns:
            continue
        source = None
        for alias in aliases:
            source = keyed.get(_key(alias))
            if source is not None:
                break
        if source is not None:
            out[canonical] = out[source]

    for numeric in ("start", "end", "q_value", "p_value", "delta", "Y", "m"):
        if numeric in out.columns:
            out[numeric] = pd.to_numeric(out[numeric], errors="coerce")

    if "region_id" not in out.columns and {"chrom", "start", "end"}.issubset(out.columns):
        out["region_id"] = (
            out["chrom"].astype(str) + ":" + out["start"].astype("Int64").astype(str) + "-" + out["end"].astype("Int64").astype(str)
        )
        warnings.append("region_id was synthesized from chrom/start/end")

    missing_core = [col for col in ("chrom", "start", "end") if col not in out.columns]
    if missing_core and require_coordinates:
        warnings.append("missing coordinate columns: " + ",".join(missing_core))
    return out, warnings


def read_table(
    path_or_file: Any,
    *,
    max_rows: int | None = None,
    require_coordinates: bool = True,
) -> tuple[pd.DataFrame, list[str]]:
    """Read a TSV-like table and normalize common DMR visualization columns."""

    warnings: list[str] = []
    try:
        df = pd.read_csv(path_or_file, sep="\t", nrows=max_rows)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(), ["empty table"]
    except Exception as exc:
        return pd.DataFrame(), [f"read_error={exc}"]

    normalized, norm_warnings = normalize_dmr_curve_columns(df, require_coordinates=require_coordinates)
    warnings.extend(norm_warnings)
    return normalized, warnings


def read_dmr_table(path_or_file: Any, *, max_rows: int | None = None) -> tuple[pd.DataFrame, list[str]]:
    return read_table(path_or_file, max_rows=max_rows, require_coordinates=True)


def read_region_counts(path_or_file: Any, *, max_rows: int | None = None) -> tuple[pd.DataFrame, list[str]]:
    df, warnings = read_table(path_or_file, max_rows=max_rows, require_coordinates=False)
    required = [col for col in ("region_id", "sample_id", "Y", "m") if col not in df.columns]
    if required:
        warnings.append("region_counts missing columns: " + ",".join(required))
    return df, warnings


def read_design_table(path_or_file: Any, *, max_rows: int | None = None) -> tuple[pd.DataFrame, list[str]]:
    df, warnings = read_table(path_or_file, max_rows=max_rows, require_coordinates=False)
    required = [col for col in ("sample_id", "condition") if col not in df.columns]
    if required:
        warnings.append("design missing columns: " + ",".join(required))
    return df, warnings


def read_beta_binom_table(path_or_file: Any, *, max_rows: int | None = None) -> tuple[pd.DataFrame, list[str]]:
    return read_table(path_or_file, max_rows=max_rows, require_coordinates=False)


def read_caller_support_table(path_or_file: Any, *, max_rows: int | None = None) -> tuple[pd.DataFrame, list[str]]:
    return read_table(path_or_file, max_rows=max_rows, require_coordinates=False)


def read_annotation_table(path_or_file: Any, *, max_rows: int | None = None) -> tuple[pd.DataFrame, list[str]]:
    return read_table(path_or_file, max_rows=max_rows, require_coordinates=False)


def coerce_table(value: Any, reader=read_table, *, max_rows: int | None = None) -> tuple[pd.DataFrame | None, list[str]]:
    """Coerce a DataFrame/path/file-like value into a normalized DataFrame."""

    if value is None:
        return None, []
    if isinstance(value, pd.DataFrame):
        return normalize_dmr_curve_columns(value, require_coordinates=reader in {read_table, read_dmr_table})
    if isinstance(value, (str, Path)) or hasattr(value, "read"):
        return reader(value, max_rows=max_rows)
    return None, [f"unsupported table input type: {type(value).__name__}"]


__all__ = [
    "COLUMN_ALIASES",
    "coerce_table",
    "normalize_dmr_curve_columns",
    "read_annotation_table",
    "read_beta_binom_table",
    "read_caller_support_table",
    "read_design_table",
    "read_dmr_table",
    "read_region_counts",
    "read_table",
]
