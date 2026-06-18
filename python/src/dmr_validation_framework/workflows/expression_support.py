#!/usr/bin/env python
"""Normalize gene-level expression evidence for DMR downstream figures."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from dmr_validation_framework.core.io import read_table


DEFAULT_GENE_REGEX = r"(Solyc\d{2}g\d{6}(?:\.\d+)?)"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build expression_support.tsv from a gene-level RNA-seq/expression table, "
            "or from a GEO series matrix plus optional platform/probe annotation."
        )
    )
    parser.add_argument("--root", type=Path, help="Dataset root. Used only for default output path.")
    parser.add_argument("--out-dir", type=Path, help="Defaults to <root>/expression_support.")
    parser.add_argument("--expression-table", type=Path, help="Generic gene-level expression table, e.g. TPM/count matrix.")
    parser.add_argument("--series-matrix", type=Path, help="GEO series_matrix.txt or .gz expression matrix.")
    parser.add_argument("--platform-annot", type=Path, help="Optional GEO platform annot file for probe-to-gene mapping.")
    parser.add_argument("--probe-gene-map", type=Path, help="Optional two-column probe_id/gene_id mapping table.")
    parser.add_argument("--gene-id-column", help="Gene id column in --expression-table.")
    parser.add_argument("--probe-id-column", help="Probe id column in --probe-gene-map or platform table.")
    parser.add_argument("--gene-map-column", help="Gene id column in --probe-gene-map.")
    parser.add_argument("--gene-id-regex", default=DEFAULT_GENE_REGEX)
    parser.add_argument("--min-expression", type=float, default=0.0, help="Detected if max expression is greater than this value.")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _clean(value: object) -> str:
    return str(value).strip().strip('"')


def _read_any_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return read_table(path)


def _canonical_gene(value: object) -> str:
    text = _clean(value)
    match = re.search(r"(Solyc\d{2}g\d{6})(?:\.\d+)?", text, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return text


def _extract_gene_ids(text: object, pattern: str) -> list[str]:
    value = _clean(text)
    return [_canonical_gene(match) for match in re.findall(pattern, value, flags=re.IGNORECASE)]


def read_geo_series_matrix(path: Path) -> pd.DataFrame:
    rows: list[str] = []
    in_table = False
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                rows.append(line)
    if not rows:
        return pd.DataFrame()
    from io import StringIO

    df = pd.read_csv(StringIO("\n".join(rows)), sep="\t")
    df.columns = [_clean(col) for col in df.columns]
    if "ID_REF" in df.columns:
        df = df.rename(columns={"ID_REF": "feature_id"})
    elif df.columns.size:
        df = df.rename(columns={df.columns[0]: "feature_id"})
    df["feature_id"] = df["feature_id"].map(_clean)
    return df


def read_geo_platform_table(path: Path) -> pd.DataFrame:
    rows: list[str] = []
    in_table = False
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                rows.append(line)
    if not rows:
        return pd.DataFrame()
    from io import StringIO

    df = pd.read_csv(StringIO("\n".join(rows)), sep="\t")
    df.columns = [_clean(col) for col in df.columns]
    return df


def detect_gene_column(df: pd.DataFrame, explicit: str | None = None) -> str | None:
    if explicit and explicit in df.columns:
        return explicit
    lowered = {str(col).lower(): col for col in df.columns}
    for key in ("gene_id", "gene", "id", "locus", "locus_id"):
        if key in lowered:
            return lowered[key]
    for col in df.columns:
        sample = " ".join(df[col].dropna().astype(str).head(200).tolist())
        if re.search(DEFAULT_GENE_REGEX, sample, flags=re.IGNORECASE):
            return str(col)
    return None


def numeric_sample_columns(df: pd.DataFrame, exclude: set[str]) -> list[str]:
    columns: list[str] = []
    for col in df.columns:
        if col in exclude:
            continue
        converted = pd.to_numeric(df[col], errors="coerce")
        if converted.notna().any():
            columns.append(str(col))
    return columns


def read_probe_gene_map(args: argparse.Namespace) -> pd.DataFrame:
    if args.probe_gene_map:
        mapping = _read_any_table(args.probe_gene_map)
        probe_col = args.probe_id_column or next((c for c in mapping.columns if str(c).lower() in {"probe_id", "feature_id", "id", "id_ref"}), None)
        gene_col = args.gene_map_column or detect_gene_column(mapping)
        if probe_col and gene_col:
            return mapping[[probe_col, gene_col]].rename(columns={probe_col: "feature_id", gene_col: "gene_id"}).dropna()
        return pd.DataFrame()
    if not args.platform_annot:
        return pd.DataFrame()
    platform = read_geo_platform_table(args.platform_annot)
    if platform.empty:
        return pd.DataFrame()
    probe_col = args.probe_id_column or ("ID" if "ID" in platform.columns else platform.columns[0])
    rows: list[dict[str, str]] = []
    for values in platform.fillna("").to_dict("records"):
        feature_id = _clean(values.get(probe_col, ""))
        if not feature_id:
            continue
        found: set[str] = set()
        for value in values.values():
            found.update(_extract_gene_ids(value, args.gene_id_regex))
        for gene_id in sorted(found):
            rows.append({"feature_id": feature_id, "gene_id": gene_id})
    return pd.DataFrame(rows)


def summarize_expression(matrix: pd.DataFrame, gene_col: str, feature_col: str | None, min_expression: float) -> pd.DataFrame:
    exclude = {gene_col}
    if feature_col:
        exclude.add(feature_col)
    sample_cols = numeric_sample_columns(matrix, exclude)
    if not sample_cols:
        return pd.DataFrame(columns=["gene_id", "n_expression_features", "mean_expression", "max_expression", "n_samples", "detected", "source_feature_ids"])
    work = matrix.copy()
    for col in sample_cols:
        work[col] = pd.to_numeric(work[col], errors="coerce")
    work["gene_id"] = work[gene_col].map(_canonical_gene)
    work = work[work["gene_id"].astype(str).ne("")]
    if work.empty:
        return pd.DataFrame()
    row_mean = work[sample_cols].mean(axis=1, skipna=True)
    row_max = work[sample_cols].max(axis=1, skipna=True)
    work["_row_mean_expression"] = row_mean
    work["_row_max_expression"] = row_max
    if feature_col:
        work["_source_feature_id"] = work[feature_col].astype(str)
    else:
        work["_source_feature_id"] = work["gene_id"].astype(str)
    summary = (
        work.groupby("gene_id", dropna=False)
        .agg(
            n_expression_features=("_source_feature_id", "nunique"),
            mean_expression=("_row_mean_expression", "mean"),
            max_expression=("_row_max_expression", "max"),
            source_feature_ids=("_source_feature_id", lambda x: ",".join(sorted(set(map(str, x)))[:20])),
        )
        .reset_index()
    )
    summary["n_samples"] = len(sample_cols)
    summary["detected"] = pd.to_numeric(summary["max_expression"], errors="coerce").fillna(0) > min_expression
    return summary.sort_values(["detected", "mean_expression", "gene_id"], ascending=[False, False, True])


def run(args: argparse.Namespace) -> int:
    out_dir = args.out_dir or ((args.root / "expression_support") if args.root else Path("outputs/expression_support"))
    out_dir.mkdir(parents=True, exist_ok=True)
    sources: list[str] = []
    matrix = pd.DataFrame()
    feature_col: str | None = None
    status = "PASS"
    notes = ""

    if args.expression_table:
        matrix = _read_any_table(args.expression_table)
        sources.append(str(args.expression_table))
        gene_col = detect_gene_column(matrix, args.gene_id_column)
        feature_col = None
        if gene_col is None:
            status = "SKIPPED"
            notes = "expression table has no gene_id-like column"
            summary = pd.DataFrame()
        else:
            summary = summarize_expression(matrix, gene_col, feature_col, args.min_expression)
    elif args.series_matrix:
        matrix = read_geo_series_matrix(args.series_matrix)
        sources.append(str(args.series_matrix))
        if matrix.empty:
            status = "SKIPPED"
            notes = "GEO series matrix table is empty or not found"
            summary = pd.DataFrame()
        else:
            mapping = read_probe_gene_map(args)
            if args.platform_annot:
                sources.append(str(args.platform_annot))
            if args.probe_gene_map:
                sources.append(str(args.probe_gene_map))
            if mapping.empty:
                status = "WARN"
                notes = "no probe-to-gene mapping found; provide --probe-gene-map or a platform annotation containing target gene IDs"
                summary = pd.DataFrame()
            else:
                mapping["feature_id"] = mapping["feature_id"].map(_clean)
                mapping["gene_id"] = mapping["gene_id"].map(_canonical_gene)
                merged = matrix.merge(mapping, on="feature_id", how="inner")
                feature_col = "feature_id"
                summary = summarize_expression(merged, "gene_id", feature_col, args.min_expression)
                if summary.empty:
                    status = "WARN"
                    notes = "probe mapping was found, but no expression rows matched mapped probes"
    else:
        status = "SKIPPED"
        notes = "provide --expression-table or --series-matrix"
        summary = pd.DataFrame()

    if not matrix.empty:
        matrix.head(200000).to_csv(out_dir / "expression_matrix_preview.tsv", sep="\t", index=False)
    if summary.empty:
        summary = pd.DataFrame(
            columns=["gene_id", "n_expression_features", "mean_expression", "max_expression", "n_samples", "detected", "source_feature_ids"]
        )
    summary.to_csv(out_dir / "expression_support.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "status": status,
                "n_expression_support_genes": len(summary),
                "n_detected_genes": int(summary["detected"].sum()) if "detected" in summary.columns else 0,
                "sources": ";".join(sources),
                "notes": notes,
            }
        ]
    ).to_csv(out_dir / "expression_support_summary.tsv", sep="\t", index=False)
    print(f"expression support: {out_dir / 'expression_support.tsv'}")
    print(f"summary: {out_dir / 'expression_support_summary.tsv'}")
    if status != "PASS":
        print(f"status: {status}; {notes}")
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
