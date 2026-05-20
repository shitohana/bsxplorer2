"""Adapters for importing already produced external DMR caller outputs.

Purpose:
    Read generic BED, DSS-like, methylKit-like, dmrseq-like, and metilene-like
    tabular outputs into the BSX2 canonical DMR schema.

Limitations:
    Adapters do not run external callers and do not validate caller-specific
    statistical assumptions. Harmonized q-values remain caller-specific.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .dmr_harmonization import assign_missing_dmr_ids, canonical_dmr_columns, normalize_dmr_coordinates


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".bed":
        return pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand", "context", "delta", "p_value", "q_value"], engine="python")
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else "," if path.suffix.lower() == ".csv" else None
    return pd.read_csv(path, sep=sep, engine="python")


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


class ExternalDmrAdapter:
    source_caller = "external"
    caller_version = "unknown"

    def __init__(self, path: str | Path, *, contrast_id: str = "", condition_a: str = "", condition_b: str = "", context: str | None = None) -> None:
        self.path = Path(path)
        self.contrast_id = contrast_id
        self.condition_a = condition_a
        self.condition_b = condition_b
        self.context = context
        self.warnings: list[str] = []

    def read_raw(self) -> pd.DataFrame:
        return _read_table(self.path)

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        return raw.copy()

    def read(self) -> pd.DataFrame:
        mapped = normalize_dmr_coordinates(self.map_columns(self.read_raw()))
        mapped = assign_missing_dmr_ids(mapped, self.source_caller.lower())
        for column in canonical_dmr_columns():
            if column not in mapped.columns:
                mapped[column] = pd.NA
        mapped["source_caller"] = self.source_caller
        mapped["caller_version"] = self.caller_version
        mapped["contrast_id"] = self.contrast_id
        mapped["condition_a"] = self.condition_a
        mapped["condition_b"] = self.condition_b
        if self.context is not None:
            mapped["context"] = self.context
        mapped["source_file"] = str(self.path)
        mapped["source_status"] = "ok" if not self.warnings else "warning"
        mapped["method_notes"] = ";".join(self.warnings)
        return mapped[canonical_dmr_columns()]


class GenericBedAdapter(ExternalDmrAdapter):
    source_caller = "generic_bed"

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        out = raw.copy()
        if "name" in out.columns:
            out["dmr_id"] = out["name"]
        return out


class DSSAdapter(ExternalDmrAdapter):
    source_caller = "DSS"

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        out = raw.copy()
        for target, aliases in {
            "chrom": ("chr", "chrom"),
            "p_value": ("pvalue", "p.value", "pval"),
            "q_value": ("fdr", "qvalue", "q_value"),
            "delta": ("diff.Methy", "diff.methy", "delta"),
            "n_cytosines": ("nCG", "ncg", "n_cytosines"),
        }.items():
            col = _find_column(out, aliases)
            if col is not None:
                out[target] = out[col]
        return out


class MethylKitAdapter(ExternalDmrAdapter):
    source_caller = "methylKit"

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        out = raw.copy()
        for target, aliases in {
            "chrom": ("chr", "chrom"),
            "p_value": ("pvalue", "p.value", "pval"),
            "q_value": ("qvalue", "q.value", "qval"),
            "delta": ("meth.diff", "meth_diff", "delta"),
        }.items():
            col = _find_column(out, aliases)
            if col is not None:
                out[target] = out[col]
        return out


class DMRseqAdapter(ExternalDmrAdapter):
    source_caller = "dmrseq"


class MetileneAdapter(ExternalDmrAdapter):
    source_caller = "metilene"


ADAPTERS = {
    "generic_bed": GenericBedAdapter,
    "DSS": DSSAdapter,
    "dss": DSSAdapter,
    "methylKit": MethylKitAdapter,
    "methylkit": MethylKitAdapter,
    "dmrseq": DMRseqAdapter,
    "metilene": MetileneAdapter,
}


def adapter_for_caller(caller: str) -> type[ExternalDmrAdapter]:
    try:
        return ADAPTERS[caller]
    except KeyError as exc:
        raise ValueError(f"Unsupported caller: {caller}") from exc
