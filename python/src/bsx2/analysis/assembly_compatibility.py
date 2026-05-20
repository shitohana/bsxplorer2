"""Assembly and coordinate compatibility diagnostics.

Purpose:
    Check whether interval tables use compatible sequence names and coordinate
    bounds relative to optional genome sizes.

Limitations:
    These checks do not perform liftover, assembly conversion, or infer
    assembly equivalence. Missing genome sizes produce warnings, not crashes.

Stability:
    Diagnostic hardening layer for thesis workflows.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd

from .seqname_harmonization import normalize_seqname


def _find_column(df: pd.DataFrame, names: tuple[str, ...]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name in lower:
            return lower[name]
    return None


def read_genome_sizes(path: str | Path) -> pd.DataFrame:
    rows = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                rows.append({"seqname": fields[0], "length": int(float(fields[1]))})
    return pd.DataFrame(rows)


def normalize_coordinate_table(df: pd.DataFrame, alias_map: Mapping[str, str] | None = None) -> pd.DataFrame:
    chrom_col = _find_column(df, ("seqname", "chrom", "chr", "chromosome"))
    start_col = _find_column(df, ("start", "start_bp", "begin"))
    end_col = _find_column(df, ("end", "end_bp", "stop"))
    if chrom_col is None or start_col is None or end_col is None:
        raise ValueError("coordinate table requires seqname/chrom, start, and end")
    return pd.DataFrame({
        "seqname": df[chrom_col].map(lambda v: normalize_seqname(v, alias_map)),
        "start": pd.to_numeric(df[start_col], errors="coerce"),
        "end": pd.to_numeric(df[end_col], errors="coerce"),
    })


def check_coordinate_bounds(
    regions_df: pd.DataFrame,
    genome_sizes_df: pd.DataFrame | None = None,
    alias_map: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    regions = normalize_coordinate_table(regions_df, alias_map)
    length_map = None
    if genome_sizes_df is not None and not genome_sizes_df.empty:
        genome = genome_sizes_df.copy()
        if "seqname" not in genome.columns:
            genome = normalize_coordinate_table(genome, alias_map)
        genome["seqname"] = genome["seqname"].map(lambda v: normalize_seqname(v, alias_map))
        length_col = _find_column(genome, ("length", "size"))
        if length_col is not None:
            length_map = dict(zip(genome["seqname"].astype(str), pd.to_numeric(genome[length_col], errors="coerce")))
    rows = []
    for idx, row in regions.iterrows():
        warnings = []
        start, end, seqname = row["start"], row["end"], row["seqname"]
        seq_length = pd.NA
        if pd.isna(start) or pd.isna(end):
            warnings.append("non_numeric_coordinate")
        else:
            if start < 0:
                warnings.append("negative_start")
            if end < start:
                warnings.append("end_before_start")
            if length_map is None:
                warnings.append("genome_sizes_unavailable")
            elif seqname not in length_map:
                warnings.append("missing_in_genome_sizes")
            else:
                seq_length = int(length_map[seqname])
                if end > seq_length:
                    warnings.append("end_exceeds_seq_length")
        rows.append({
            "row_index": idx,
            "seqname": seqname,
            "start": start,
            "end": end,
            "seq_length": seq_length,
            "out_of_bounds": any(w in warnings for w in ("negative_start", "end_before_start", "end_exceeds_seq_length")),
            "warnings": ";".join(warnings),
            "compatibility_status": "ok" if not warnings else "warning",
        })
    return pd.DataFrame(rows)


def compare_coordinate_sources(*tables: Any, genome_sizes_df: pd.DataFrame | None = None, alias_map: Mapping[str, str] | None = None) -> pd.DataFrame:
    rows = []
    genome_names = set(genome_sizes_df["seqname"].astype(str)) if genome_sizes_df is not None and "seqname" in genome_sizes_df.columns else None
    for index, table in enumerate(tables):
        source_name, df = table if isinstance(table, tuple) else (f"source_{index + 1}", table)
        coords = normalize_coordinate_table(df, alias_map)
        seqnames = set(coords["seqname"].astype(str))
        bounds = check_coordinate_bounds(coords, genome_sizes_df, alias_map)
        missing = sorted(seqnames - genome_names) if genome_names is not None else []
        out_count = int(bounds["out_of_bounds"].sum()) if not bounds.empty else 0
        rows.append({
            "source_name": source_name,
            "n_regions": len(coords),
            "n_seqnames": len(seqnames),
            "min_start": coords["start"].min() if not coords.empty else pd.NA,
            "max_end": coords["end"].max() if not coords.empty else pd.NA,
            "missing_in_genome_sizes": ";".join(missing),
            "out_of_bounds_count": out_count,
            "compatibility_status": "warning_no_genome_sizes" if genome_sizes_df is None else "warning" if missing or out_count else "ok",
        })
    return pd.DataFrame(rows)


def write_assembly_compatibility_report(report: pd.DataFrame, out_path: str | Path) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out_path, sep="\t", index=False)
