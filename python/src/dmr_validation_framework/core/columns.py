"""Column-name normalization and region-table helpers."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence

import pandas as pd

ALIASES = {
    "region_id": ["region_id", "dmr_id", "harmonized_region_id", "id"],
    "caller": ["source_caller", "caller", "method", "source"],
    "chrom": ["chrom", "seqname", "seqnames", "chr", "chromosome"],
    "start": ["start", "start0", "begin"],
    "end": ["end", "stop"],
    "context": ["context", "ctx"],
    "delta": [
        "delta",
        "delta_methylation",
        "meth.diff",
        "meth_diff",
        "diff.methy",
        "diff_methy",
        "mean_diff",
        "meandiff",
        "maxdiff",
        "beta",
        "effect",
        "region_delta",
        "delta_callus_minus_seedling",
        "region_delta_callus_minus_seedling",
    ],
    "q_value": ["q_value", "qvalue", "q.value", "qval", "q", "padj", "fdr", "minfdr", "region_q_value"],
}

def first_existing(columns: Sequence[str], aliases: Sequence[str]) -> str | None:
    lower_map = {str(col).lower(): str(col) for col in columns}
    for alias in aliases:
        if alias.lower() in lower_map:
            return lower_map[alias.lower()]
    return None

def has_columns(df: pd.DataFrame, required: set[str]) -> tuple[bool, str]:
    missing = sorted(required - set(df.columns))
    return len(missing) == 0, ",".join(missing)


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(column).strip() for column in out.columns]
    return out


def find_col(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    return first_existing(df.columns, list(candidates))


def require_columns(df: pd.DataFrame, required: Iterable[str], *, table_name: str = "table") -> None:
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"{table_name} missing required columns: {', '.join(missing)}")


def normalize_delta_scale(df: pd.DataFrame, source: str = "", notes: list[str] | None = None) -> pd.DataFrame:
    """Normalize methylation delta to biological fraction scale.

    Canonical convention:
    - signed delta is in [-1, 1]
    - abs_delta is in [0, 1]

    Some callers, especially BSmooth/bsseq-derived outputs, may report
    methylation differences in percent points, e.g. 12.5 instead of 0.125.
    """
    if "delta" not in df.columns:
        return df

    out = df.copy()
    delta = pd.to_numeric(out["delta"], errors="coerce")
    abs_delta = delta.abs().dropna()

    out["delta_scale_original"] = "fraction_or_unknown"
    out["delta_scale_normalized"] = "unchanged"

    if abs_delta.empty:
        out["abs_delta"] = delta.abs()
        return out

    caller_hint = str(source)
    if "caller" in out.columns and out["caller"].notna().any():
        caller_hint += " " + " ".join(out["caller"].dropna().astype(str).unique()[:5])
    caller_hint = caller_hint.lower()

    q95_abs = float(abs_delta.quantile(0.95))
    max_abs = float(abs_delta.max())

    is_known_percent_caller = ("bsmooth" in caller_hint) or ("bsseq" in caller_hint)

    # Known BSmooth/bsseq-like outputs are usually percent-point deltas.
    # Do not require max_abs <= 100 here, because raw outputs may contain outliers.
    looks_like_percent_known = is_known_percent_caller and q95_abs > 1.0

    # For unknown callers, keep auto-detection conservative.
    looks_like_percent_auto = (not is_known_percent_caller) and q95_abs > 1.0 and max_abs <= 100.0

    if looks_like_percent_known:
        delta = delta / 100.0
        out["delta_scale_original"] = "percent"
        out["delta_scale_normalized"] = "percent_to_fraction"
        if notes is not None:
            notes.append("delta normalized from percent to fraction for BSmooth/bsseq-like caller")
    elif looks_like_percent_auto:
        delta = delta / 100.0
        out["delta_scale_original"] = "percent_auto_detected"
        out["delta_scale_normalized"] = "percent_to_fraction"
        if notes is not None:
            notes.append("delta auto-normalized from percent-like scale to fraction scale")
    else:
        out["delta_scale_original"] = "fraction"
        out["delta_scale_normalized"] = "unchanged"

    # Biological sanity: methylation fraction difference cannot exceed [-1, 1].
    too_large_after = delta.abs() > 1.0
    if bool(too_large_after.fillna(False).any()):
        if notes is not None:
            notes.append("delta values outside [-1,1] clipped after scale normalization")
        delta = delta.clip(-1.0, 1.0)

    out["delta"] = delta
    out["abs_delta"] = delta.abs()
    return out


def normalize_region_table(df: pd.DataFrame, source: str = "") -> tuple[pd.DataFrame, list[str]]:
    notes: list[str] = []
    out = pd.DataFrame(index=df.index)
    for canonical, aliases in ALIASES.items():
        col = first_existing(df.columns, aliases)
        if col is not None:
            out[canonical] = df[col]
    if "region_id" not in out.columns:
        out["region_id"] = [f"{source}:{i}" for i in range(len(df))]
        notes.append("region_id synthesized from row index")
    if "chrom" not in out.columns or "start" not in out.columns or "end" not in out.columns:
        rid_col = first_existing(df.columns, ["region_id", "dmr_id", "harmonized_region_id"])
        if rid_col:
            parsed = df[rid_col].astype(str).str.extract(r"^([^:]+):(\d+)-(\d+)(?::([^:]+))?")
            if "chrom" not in out.columns:
                out["chrom"] = parsed[0]
            if "start" not in out.columns:
                out["start"] = parsed[1]
            if "end" not in out.columns:
                out["end"] = parsed[2]
            if "context" not in out.columns:
                out["context"] = parsed[3]
            notes.append("coordinates parsed from region_id-like column")
    if "start" not in out.columns or "end" not in out.columns:
        notes.append("coordinate columns missing; table skipped for interval analysis")
        return pd.DataFrame(), notes
    for col in ("start", "end", "delta", "q_value"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    if "context" not in out.columns:
        out["context"] = "NA"
        notes.append("context missing; set to NA")
    if "caller" not in out.columns:
        out["caller"] = source
    out = normalize_delta_scale(out, source=source, notes=notes)
    out["chrom"] = out.get("chrom", pd.Series(["NA"] * len(out), index=out.index)).astype(str)
    out["context"] = out["context"].astype(str)
    out = out.dropna(subset=["start", "end"]).copy()
    out = out[out["start"] < out["end"]].copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    if "q_value" in out.columns:
        out["significant_by_q"] = out["q_value"] < 0.05
    return out.reset_index(drop=True), notes

def significant_or_all(df: pd.DataFrame, max_rows: int = 25_000) -> tuple[pd.DataFrame, str]:
    note = "all rows used"
    if "q_value" in df.columns and df["q_value"].notna().any():
        sig = df[df["q_value"] < 0.05].copy()
        if not sig.empty:
            if len(sig) > max_rows:
                sig = sig.sort_values("q_value", na_position="last").head(max_rows)
                return sig.reset_index(drop=True), f"filtered to q_value < 0.05 and capped to {max_rows} rows"
            return sig.reset_index(drop=True), "filtered to q_value < 0.05"
    bool_cols = [c for c in df.columns if "significant" in str(c).lower()]
    if bool_cols:
        mask = pd.Series(False, index=df.index)
        for col in bool_cols:
            mask |= df[col].astype(str).str.lower().isin(["true", "1", "yes"])
        sig = df[mask].copy()
        if not sig.empty:
            if len(sig) > max_rows:
                sig = sig.head(max_rows)
                return sig.reset_index(drop=True), f"filtered by significant columns and capped to {max_rows} rows: {','.join(bool_cols)}"
            return sig.reset_index(drop=True), f"filtered by significant columns: {','.join(bool_cols)}"
    if len(df) > max_rows:
        if "q_value" in df.columns:
            df = df.sort_values("q_value", na_position="last").head(max_rows)
            note = f"capped to first {max_rows} rows after q_value sort"
        else:
            df = df.head(max_rows)
            note = f"capped to first {max_rows} rows"
    return df.reset_index(drop=True), note

def canonical_caller_name(value: str) -> str:
    lower = value.lower()
    if "methylkit" in lower:
        return "methylKit"
    if "dss" in lower:
        return "DSS"
    if "regional" in lower or "bsx2" in lower or "evidence" in lower:
        return "Regional Evidence"
    if "dmrseq" in lower:
        return "dmrseq"
    if "radmeth" in lower:
        return "RADMeth"
    if "moabs" in lower:
        return "MOABS"
    if "methylsig" in lower:
        return "methylSig"
    if "bsmooth" in lower or "bsseq" in lower:
        return "BSmooth"
    if "biseq" in lower:
        return "BiSeq"
    if "dmrcate" in lower:
        return "DMRcate"
    if "methylasso" in lower:
        return "MethyLasso"
    if "hmm-dm" in lower or "hmmdm" in lower:
        return "HMM-DM"
    if "comb-p" in lower or "combp" in lower:
        return "comb-p"
    if "metilene" in lower:
        return "metilene"
    return value

def section_for_bin(bin_index: int, upstream_bins: int = 20, body_bins: int = 100) -> str:
    if bin_index < upstream_bins:
        return "upstream"
    if bin_index < upstream_bins + body_bins:
        return "gene_body"
    return "downstream"

def bin_cols(df: pd.DataFrame) -> list[str]:
    cols = [str(c) for c in df.columns if re.fullmatch(r"bin_\d+", str(c))]
    return sorted(cols, key=lambda c: int(c.split("_", 1)[1]))
