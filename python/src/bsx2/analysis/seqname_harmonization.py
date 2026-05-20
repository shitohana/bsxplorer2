"""Sequence-name harmonization helpers for non-model organism workflows.

Purpose:
    Compare and harmonize arbitrary sequence names such as chromosomes,
    scaffolds, contigs, or RefSeq accessions.

Input assumptions:
    Alias tables contain one of ``input_name/canonical_name``,
    ``source/target``, or ``raw_name/normalized_name``. Names are matched as
    strings and arbitrary scaffold names are preserved by default.

Limitations:
    This module harmonizes labels only. It does not perform liftover,
    coordinate conversion, assembly conversion, or orthology mapping.

Stability:
    Stable lightweight utility layer used by CLI checks and interval
    aggregation code.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd


_ALIAS_COLUMN_PAIRS = (
    ("input_name", "canonical_name"),
    ("source", "target"),
    ("raw_name", "normalized_name"),
)


def _read_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    suffix = table_path.suffix.lower()
    sep = "\t" if suffix in {".tsv", ".tab"} else "," if suffix == ".csv" else None
    return pd.read_csv(table_path, sep=sep, engine="python")


def read_seqname_aliases(path: str | Path) -> dict[str, str]:
    df = _read_table(path)
    lower = {str(c).lower(): c for c in df.columns}
    for source, target in _ALIAS_COLUMN_PAIRS:
        if source in lower and target in lower:
            return {
                str(raw): str(norm)
                for raw, norm in zip(df[lower[source]], df[lower[target]])
                if pd.notna(raw) and pd.notna(norm)
            }
    raise ValueError("Alias table must contain input_name/canonical_name, source/target, or raw_name/normalized_name")


def normalize_seqname(name: Any, alias_map: Mapping[str, str] | None = None, strip_chr: bool = False) -> str:
    value = "" if pd.isna(name) else str(name)
    if strip_chr and value.lower().startswith("chr"):
        value = value[3:]
    if alias_map is None:
        return value
    return str(alias_map.get(value, value))


def apply_seqname_aliases(
    df: pd.DataFrame,
    column: str,
    alias_map: Mapping[str, str] | None,
    output_column: str | None = None,
    *,
    copy: bool = True,
) -> pd.DataFrame:
    if column not in df.columns:
        raise ValueError(f"Column not found: {column}")
    out = df.copy() if copy else df
    target = output_column or column
    out[target] = out[column].map(lambda value: normalize_seqname(value, alias_map))
    return out


def compare_seqname_sets(
    left_names: Iterable[Any],
    right_names: Iterable[Any],
    alias_map: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    left = {normalize_seqname(v, alias_map) for v in left_names if pd.notna(v)}
    right = {normalize_seqname(v, alias_map) for v in right_names if pd.notna(v)}
    common = left & right
    denominator = max(len(left), len(right), 1)
    fraction = len(common) / denominator
    return {
        "n_left": len(left),
        "n_right": len(right),
        "n_common": len(common),
        "left_only": ";".join(sorted(left - right)),
        "right_only": ";".join(sorted(right - left)),
        "compatibility_status": "pass" if fraction >= 0.95 else "warning",
        "matching_fraction": fraction,
    }


def write_seqname_compatibility_report(
    left_names: Iterable[Any],
    right_names: Iterable[Any],
    out_path: str | Path,
    alias_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    report = pd.DataFrame([compare_seqname_sets(left_names, right_names, alias_map)])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out_path, sep="\t", index=False)
    return report


def validate_seqname_compatibility(
    left_names: Iterable[Any],
    right_names: Iterable[Any],
    min_fraction: float = 0.95,
) -> dict[str, Any]:
    result = compare_seqname_sets(left_names, right_names)
    result["compatibility_status"] = "pass" if result["matching_fraction"] >= min_fraction else "warning"
    result["min_fraction"] = min_fraction
    return result
