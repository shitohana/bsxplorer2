"""BS-seq QC report importers and lightweight counts QC.

Purpose:
    Import selected Bismark alignment, deduplication, and M-bias report fields,
    and compute simple QC summaries from methylation count tables.

Limitations:
    Missing reports are warnings, not failures. Bisulfite conversion failure is
    not estimated unless a dedicated spike-in or conversion report is provided;
    this module reports it as unavailable.

Stability:
    Diagnostic importer layer. It does not run Bismark or any raw pipeline.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _read_text(path: str | Path) -> tuple[str, str]:
    path = Path(path)
    if not path.exists():
        return "", f"missing_file:{path}"
    return path.read_text(encoding="utf-8", errors="replace"), ""


def _number_after(pattern: str, text: str) -> float | None:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    return float(match.group(1).replace(",", "")) if match else None


def parse_bismark_alignment_report(path: str | Path) -> dict[str, Any]:
    text, warning = _read_text(path)
    if warning:
        return {"source": str(path), "warning": warning}
    total = _number_after(r"Sequences analysed in total:\s*([0-9,]+)", text) or _number_after(r"total reads[:\s]+([0-9,]+)", text)
    aligned = _number_after(r"Number of alignments with a unique best hit.*?:\s*([0-9,]+)", text) or _number_after(r"aligned reads[:\s]+([0-9,]+)", text)
    rate = _number_after(r"Mapping efficiency:\s*([0-9.]+)%", text) or ((aligned / total * 100) if total and aligned else None)
    return {
        "source": str(path),
        "total_reads": total,
        "aligned_reads": aligned,
        "alignment_rate": rate,
        "methylation_CG": _number_after(r"C methylated in CpG context:\s*([0-9.]+)%", text),
        "methylation_CHG": _number_after(r"C methylated in CHG context:\s*([0-9.]+)%", text),
        "methylation_CHH": _number_after(r"C methylated in CHH context:\s*([0-9.]+)%", text),
        "warning": "",
    }


def parse_bismark_dedup_report(path: str | Path) -> dict[str, Any]:
    text, warning = _read_text(path)
    if warning:
        return {"source": str(path), "warning": warning}
    dedup = _number_after(r"Total number of alignments analysed in total:\s*([0-9,]+)", text) or _number_after(r"deduplicated reads[:\s]+([0-9,]+)", text)
    duplicate_rate = _number_after(r"Duplicated alignments removed:\s*[0-9,]+\s*\(([0-9.]+)%\)", text) or _number_after(r"duplicate rate[:\s]+([0-9.]+)", text)
    return {"source": str(path), "deduplicated_reads": dedup, "duplicate_rate": duplicate_rate, "retained_rate": None if duplicate_rate is None else 100 - duplicate_rate, "warning": ""}


def parse_bismark_mbias_report(path: str | Path) -> pd.DataFrame:
    text, warning = _read_text(path)
    if warning:
        return pd.DataFrame([{"source": str(path), "warning": warning}])
    rows = []
    for line in text.splitlines():
        fields = re.split(r"\s+|\t", line.strip())
        if len(fields) >= 5 and fields[0].upper() in {"CG", "CHG", "CHH"}:
            rows.append({"context": fields[0].upper(), "position": int(float(fields[1])), "methylation_percentage": float(fields[2]), "count_methylated": float(fields[3]), "count_unmethylated": float(fields[4]), "warning": ""})
    return pd.DataFrame(rows)


def _find_column(df: pd.DataFrame, names: tuple[str, ...]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def _balance(series: pd.Series) -> str:
    counts = series.astype(str).value_counts(dropna=False)
    total = counts.sum()
    return ";".join(f"{k}:{v / total:.3f}" for k, v in counts.items()) if total else ""


def compute_counts_qc(counts_df: pd.DataFrame) -> dict[str, Any]:
    if counts_df.empty:
        return {"n_rows": 0, "total_mC": 0, "total_uC": 0, "total_coverage": 0, "mean_coverage": np.nan, "median_coverage": np.nan, "context_balance": "", "strand_balance": "", "retained_cytosines": 0, "high_coverage_fraction": np.nan, "zero_coverage_fraction": np.nan}
    mc_col = _find_column(counts_df, ("mC", "methylated_count", "count_m"))
    uc_col = _find_column(counts_df, ("uC", "unmethylated_count"))
    total_col = _find_column(counts_df, ("total", "coverage", "count_total"))
    if mc_col is None:
        raise ValueError("counts_df missing methylated count column")
    mc = pd.to_numeric(counts_df[mc_col], errors="coerce").fillna(0).clip(lower=0)
    if uc_col is not None:
        uc = pd.to_numeric(counts_df[uc_col], errors="coerce").fillna(0).clip(lower=0)
    elif total_col is not None:
        uc = (pd.to_numeric(counts_df[total_col], errors="coerce").fillna(0).clip(lower=0) - mc).clip(lower=0)
    else:
        raise ValueError("counts_df must contain uC/unmethylated_count or total/coverage")
    total = pd.to_numeric(counts_df[total_col], errors="coerce").fillna(0).clip(lower=0) if total_col else mc + uc
    return {
        "n_rows": int(len(counts_df)),
        "total_mC": float(mc.sum()),
        "total_uC": float(uc.sum()),
        "total_coverage": float(total.sum()),
        "mean_coverage": float(total.mean()),
        "median_coverage": float(total.median()),
        "context_balance": _balance(counts_df["context"]) if "context" in counts_df.columns else "",
        "strand_balance": _balance(counts_df["strand"]) if "strand" in counts_df.columns else "",
        "retained_cytosines": int((total > 0).sum()),
        "high_coverage_fraction": float((total >= 100).mean()),
        "zero_coverage_fraction": float((total == 0).mean()),
    }


def classify_bsseq_qc(qc_row: dict[str, Any] | pd.Series) -> dict[str, str]:
    row = dict(qc_row)
    reasons = []
    status = "usable"
    if float(row.get("zero_coverage_fraction", 0) or 0) > 0.5:
        status = "fail"; reasons.append("zero_coverage_fraction_above_0.5")
    if status != "fail" and float(row.get("high_coverage_fraction", 0) or 0) > 0.2:
        status = "caution"; reasons.append("high_coverage_fraction_above_0.2")
    if not reasons:
        reasons.append("conversion_failure_unavailable_without_spikein_or_report")
    return {"qc_status": status, "warning_reasons": ";".join(reasons)}


def write_bsseq_qc_outputs(summary_rows: list[dict[str, Any]], warnings: list[dict[str, Any]], out_summary: str | Path, out_warnings: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.DataFrame(summary_rows)
    warnings_df = pd.DataFrame(warnings) if warnings else pd.DataFrame(columns=["source", "warning"])
    Path(out_summary).parent.mkdir(parents=True, exist_ok=True)
    Path(out_warnings).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_summary, sep="\t", index=False)
    warnings_df.to_csv(out_warnings, sep="\t", index=False)
    return summary, warnings_df
