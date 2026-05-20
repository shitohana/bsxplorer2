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

ANNOTATION_COLUMNS = ("annotation", "feature_type", "region_class", "genomic_feature", "feature", "region_type", "class")
COUNT_COLUMNS = (
    "count",
    "n_regions",
    "n_dmrs",
    "n",
    "frequency",
    "overlap_count",
    "observed_dmr_count",
    "dmr_count",
)


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


def classify_dmr_table(df: pd.DataFrame) -> tuple[str, list[str]]:
    """Classify whether a DMR-like table is region-level or only a summary."""

    warnings: list[str] = []
    if df is None or df.empty:
        return "unsupported", ["dmr_table_empty"]
    normalized, norm_warnings = normalize_dmr_curve_columns(df, require_coordinates=False)
    warnings.extend(norm_warnings)
    has_region = "region_id" in normalized.columns or {"chrom", "start", "end"}.issubset(normalized.columns)
    if not has_region:
        return "unsupported", warnings + ["dmr_table_missing_region_id_or_coordinates"]
    if len(normalized) < 20:
        return "summary_table", warnings + ["dmr_table_has_few_rows_check_summary_input"]
    if "chrom" in normalized.columns and normalized["chrom"].nunique(dropna=True) == 1:
        warnings.append("single_chromosome_input")
    return ("full_region_level_dmr" if len(normalized) >= 1000 else "top_region_level_dmr"), warnings


def classify_annotation_table(df: pd.DataFrame) -> tuple[str, list[str]]:
    """Classify annotation input as region-level, enrichment/feature summary, or unsupported."""

    warnings: list[str] = []
    if df is None or df.empty:
        return "unsupported", ["annotation_table_empty"]
    normalized, norm_warnings = normalize_dmr_curve_columns(df, require_coordinates=False)
    warnings.extend(norm_warnings)
    lower = {_key(col): col for col in normalized.columns}
    has_region = "region_id" in normalized.columns or {"chrom", "start", "end"}.issubset(normalized.columns)
    ann_cols = [lower[_key(col)] for col in ANNOTATION_COLUMNS if _key(col) in lower]
    count_cols = [lower[_key(col)] for col in COUNT_COLUMNS if _key(col) in lower]
    enrichment_cols = [col for col in normalized.columns if _key(col) in {"p_value", "pvalue", "pval", "q_value", "qvalue", "fdr", "enrichment", "odds_ratio", "fraction"}]
    if has_region and ann_cols:
        return "region_level_annotation", warnings
    if ann_cols and (count_cols or enrichment_cols):
        return "enrichment_summary", warnings + ["annotation_input_is_summary_not_region_level"]
    if ann_cols:
        return "feature_summary", warnings + ["annotation_summary_missing_region_ids_and_count_column"]
    return "unsupported", warnings + ["annotation_table_missing_annotation_like_column"]


def classify_region_counts_table(df: pd.DataFrame) -> tuple[str, list[str]]:
    """Classify whether a table can support DMR sample-level methylation curves."""

    warnings: list[str] = []
    if df is None or df.empty:
        return "unsupported", ["region_counts_empty"]
    normalized, norm_warnings = normalize_dmr_curve_columns(df, require_coordinates=False)
    warnings.extend(norm_warnings)
    required = [col for col in ("region_id", "sample_id", "Y", "m") if col not in normalized.columns]
    if required:
        return "unsupported", warnings + ["region_counts_missing_columns:" + ",".join(required)]
    n_samples = normalized["sample_id"].nunique(dropna=True)
    if n_samples < 2:
        warnings.append("region_counts_have_fewer_than_two_samples")
    return "region_sample_counts", warnings


def classify_design_table(df: pd.DataFrame, region_counts_df: pd.DataFrame | None = None) -> tuple[str, list[str]]:
    """Classify design table and optionally check sample overlap with region counts."""

    warnings: list[str] = []
    if df is None or df.empty:
        return "unsupported", ["design_table_empty"]
    normalized, norm_warnings = normalize_dmr_curve_columns(df, require_coordinates=False)
    warnings.extend(norm_warnings)
    required = [col for col in ("sample_id", "condition") if col not in normalized.columns]
    if required:
        return "unsupported", warnings + ["design_missing_columns:" + ",".join(required)]
    if normalized["sample_id"].duplicated().any():
        warnings.append("design_has_duplicate_sample_id")
    if region_counts_df is not None and not region_counts_df.empty:
        counts_norm, _ = normalize_dmr_curve_columns(region_counts_df, require_coordinates=False)
        if "sample_id" in counts_norm.columns:
            counts_samples = set(counts_norm["sample_id"].dropna().astype(str))
            design_samples = set(normalized["sample_id"].dropna().astype(str))
            missing_design = sorted(counts_samples - design_samples)
            missing_counts = sorted(design_samples - counts_samples)
            if missing_design:
                warnings.append("samples_in_counts_missing_from_design:" + ",".join(missing_design[:10]))
            if missing_counts:
                warnings.append("samples_in_design_missing_from_counts:" + ",".join(missing_counts[:10]))
    return "sample_design", warnings


def find_best_dmr_table(start_path: str | Path) -> tuple[Path | None, str]:
    """Find a fuller DMR table near ``start_path`` for distribution-style curves."""

    start = Path(start_path)
    search_roots = [start.parent]
    for parent in start.parents:
        if parent.name in {"run", "processed_geo_cx"}:
            search_roots.append(parent)
    candidates: list[tuple[tuple[int, int, int, int], Path]] = []
    for root in dict.fromkeys(search_roots):
        if not root.exists():
            continue
        for name in ("dmr_evidence_scores.tsv", "dmr_regions.tsv", "dmr_region_count_tests.tsv"):
            path = root / name
            if not path.exists():
                continue
            try:
                with path.open("r", encoding="utf-8", errors="replace") as fh:
                    header = fh.readline().rstrip("\n").split("\t")
                    n_rows = sum(1 for _ in fh)
            except Exception:
                continue
            low = {_key(col) for col in header}
            has_coord = int(bool({"chrom", "chr", "seqname"} & low) and bool({"start", "start_bp"} & low) and bool({"end", "end_bp"} & low))
            has_evidence = int("evidence_class" in low)
            is_subset = int(any(token in path.name.lower() for token in ("top", "head", "subset", "smoke")))
            score = (has_coord, n_rows, has_evidence, -is_subset)
            candidates.append((score, path))
    if not candidates:
        return None, "no_candidate_dmr_table_found"
    candidates.sort(reverse=True, key=lambda item: item[0])
    selected = candidates[0][1]
    return selected, f"selected_by_has_coordinates_row_count_evidence_columns:{selected}"


__all__ = [
    "COLUMN_ALIASES",
    "classify_annotation_table",
    "classify_design_table",
    "classify_dmr_table",
    "classify_region_counts_table",
    "coerce_table",
    "find_best_dmr_table",
    "normalize_dmr_curve_columns",
    "read_annotation_table",
    "read_beta_binom_table",
    "read_caller_support_table",
    "read_design_table",
    "read_dmr_table",
    "read_region_counts",
    "read_table",
]
