"""Adapters for importing already produced external DMR caller outputs.

Purpose:
    Read generic BED, DSS-like, methylKit-like, dmrseq-like, metilene-like,
    and other common DMR caller tabular outputs into the BSX2 canonical DMR
    schema.

Limitations:
    Adapters do not run external callers and do not validate caller-specific
    statistical assumptions. Harmonized q-values remain caller-specific.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .dmr_harmonization import (
    CALLER_MODEL_FAMILY,
    COORD_ONE_BASED_INCLUSIVE,
    COORD_ZERO_BASED_HALF_OPEN,
    assign_missing_dmr_ids,
    caller_delta_sign,
    canonical_dmr_columns,
    harmonize_delta_sign,
    normalize_dmr_coordinates,
    normalize_interval_convention,
)


COMMON_DMR_ALIASES = {
    "dmr_id": ("dmr_id", "region_id", "id", "name", "cluster", "cluster_id", "block_id", "segment_id", "state_id"),
    "chrom": ("chrom", "chr", "chromosome", "seqname", "seqnames", "seqid", "contig"),
    "start": ("start", "begin", "start_bp", "pos_start", "idxstart", "idx_start", "locus_start"),
    "end": ("end", "stop", "end_bp", "pos_end", "idxend", "idx_end", "locus_end"),
    "strand": ("strand", "str"),
    "context": ("context", "ctx", "methylation_context", "sequence_context"),
    "direction": ("direction", "meth_status", "change", "state"),
    "delta": (
        "delta",
        "diff",
        "difference",
        "meth.diff",
        "meth_diff",
        "diff.methy",
        "diff_methy",
        "mean_diff",
        "meandiff",
        "meanDiff",
        "maxdiff",
        "max_diff",
        "beta",
        "estimate",
        "effect",
        "methylation_difference",
        "mean_methylation_difference",
        "areaStat",
        "value",
    ),
    "p_value": (
        "p_value",
        "pvalue",
        "p.value",
        "pval",
        "p",
        "p_val",
        "rawp",
        "raw_p",
        "region_p",
        "region_p_value",
        "comb_p",
        "slk_p",
        "fisher",
        "stouffer",
    ),
    "q_value": (
        "q_value",
        "qvalue",
        "q.value",
        "qval",
        "q",
        "fdr",
        "padj",
        "adj_p",
        "adj.p",
        "adjusted_p",
        "minfdr",
        "min_fdr",
        "region_q",
        "region_q_value",
        "region_slk_sidak_p",
        "slk_sidak_p",
        "sidak_p",
        "fwer",
    ),
    "n_sites": ("n_sites", "nsites", "n_cpg", "n_cpgs", "n_cytosines", "ncpg", "no.cpgs", "l", "L"),
    "n_cytosines": ("n_cytosines", "n_sites", "nsites", "n_cpg", "n_cpgs", "ncpg", "no.cpgs", "l", "L"),
    "mean_methylation_a": ("mean_methylation_a", "mean_a", "meth_a", "meth1", "mean1", "mu1", "group1", "control_methylation"),
    "mean_methylation_b": ("mean_methylation_b", "mean_b", "meth_b", "meth2", "mean2", "mu2", "group2", "case_methylation"),
}


FAMILY_NOTES = {
    "generic_bed": "Generic interval import; no caller-specific statistical model is inferred.",
    "fisher_logistic": "Simple count/proportion model family; useful as a baseline for simple designs.",
    "beta_binomial": "Count-based model family for methylated/unmethylated counts with overdispersion diagnostics.",
    "smoothing": "Smoothing family; regional signal is inferred by borrowing information across neighboring CpGs.",
    "smoothing_beta_regression": "Smoothing plus beta-regression family for regional methylation signal.",
    "limma_kernel_smoothing": "Limma-style statistics with regional kernel smoothing.",
    "segmentation": "Segmentation family; adjacent CpGs are grouped into candidate regions or methylation states.",
    "hmm": "Hidden Markov model family; spatial dependence between neighboring CpGs is part of the caller model.",
    "pvalue_combination": "Regional enrichment family; CpG-level p-values are combined into candidate regions.",
    "regional_fdr": "Strict region-level family; significance is interpreted at the region/bump level.",
    "unknown": "Schema harmonization only; no caller-specific statistical model is inferred.",
}


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".bed":
        return pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand", "context", "delta", "p_value", "q_value"], engine="python")
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else "," if path.suffix.lower() == ".csv" else None
    return pd.read_csv(path, sep=sep, engine="python")


def _find_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    normalized = {str(c).lower().replace("-", "_").replace(".", "_").replace(" ", "_"): c for c in df.columns}
    for candidate in candidates:
        key = candidate.lower()
        if key in lower:
            return lower[key]
        normalized_key = key.replace("-", "_").replace(".", "_").replace(" ", "_")
        if normalized_key in normalized:
            return normalized[normalized_key]
    return None


def _apply_aliases(out: pd.DataFrame, aliases: dict[str, tuple[str, ...]] | None = None) -> pd.DataFrame:
    aliases = aliases or COMMON_DMR_ALIASES
    for target, candidates in aliases.items():
        col = _find_column(out, candidates)
        if col is not None:
            out[target] = out[col]
    return out


class ExternalDmrAdapter:
    source_caller = "external"
    caller_version = "unknown"
    model_family = "generic_bed"
    # Native coordinate convention of this caller's output. Normalized to
    # BED-style 0-based half-open at import (see normalize_interval_convention).
    coordinate_system = COORD_ZERO_BASED_HALF_OPEN

    def __init__(
        self,
        path: str | Path,
        *,
        contrast_id: str = "",
        condition_a: str = "",
        condition_b: str = "",
        context: str | None = None,
        coordinate_system: str | None = None,
        delta_sign: int | None = None,
    ) -> None:
        self.path = Path(path)
        self.contrast_id = contrast_id
        self.condition_a = condition_a
        self.condition_b = condition_b
        self.context = context
        # Per-dataset overrides for callers whose group coding or coordinate
        # convention differs from the documented defaults.
        if coordinate_system is not None:
            self.coordinate_system = coordinate_system
        self.delta_sign = delta_sign if delta_sign is not None else caller_delta_sign(self.source_caller)
        self.warnings: list[str] = []

    def read_raw(self) -> pd.DataFrame:
        return _read_table(self.path)

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        return _apply_aliases(raw.copy())

    def read(self) -> pd.DataFrame:
        mapped = normalize_dmr_coordinates(self.map_columns(self.read_raw()))
        mapped = normalize_interval_convention(mapped, self.coordinate_system)
        mapped = harmonize_delta_sign(mapped, self.delta_sign)
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
        notes = [f"model_family={self.model_family}", FAMILY_NOTES.get(self.model_family, "Schema harmonization only.")]
        notes.extend(self.warnings)
        mapped["method_notes"] = ";".join(notes)
        return mapped[canonical_dmr_columns()]


class GenericBedAdapter(ExternalDmrAdapter):
    source_caller = "generic_bed"
    model_family = "generic_bed"

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        out = _apply_aliases(raw.copy())
        if "name" in out.columns:
            out["dmr_id"] = out["name"]
        return out


class DSSAdapter(ExternalDmrAdapter):
    source_caller = "DSS"
    model_family = CALLER_MODEL_FAMILY["DSS"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        return _apply_aliases(raw.copy())


class MethylKitAdapter(ExternalDmrAdapter):
    source_caller = "methylKit"
    model_family = CALLER_MODEL_FAMILY["methylKit"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE

    def map_columns(self, raw: pd.DataFrame) -> pd.DataFrame:
        return _apply_aliases(raw.copy())


class DMRseqAdapter(ExternalDmrAdapter):
    source_caller = "dmrseq"
    model_family = CALLER_MODEL_FAMILY["dmrseq"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE


class RADMethAdapter(ExternalDmrAdapter):
    source_caller = "RADMeth"
    model_family = CALLER_MODEL_FAMILY["RADMeth"]


class MOABSAdapter(ExternalDmrAdapter):
    source_caller = "MOABS"
    model_family = CALLER_MODEL_FAMILY["MOABS"]


class MethylSigAdapter(ExternalDmrAdapter):
    source_caller = "methylSig"
    model_family = CALLER_MODEL_FAMILY["methylSig"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE


class BSmoothAdapter(ExternalDmrAdapter):
    source_caller = "BSmooth"
    model_family = CALLER_MODEL_FAMILY["BSmooth"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE


class BiSeqAdapter(ExternalDmrAdapter):
    source_caller = "BiSeq"
    model_family = CALLER_MODEL_FAMILY["BiSeq"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE


class DMRcateAdapter(ExternalDmrAdapter):
    source_caller = "DMRcate"
    model_family = CALLER_MODEL_FAMILY["DMRcate"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE


class MethyLassoAdapter(ExternalDmrAdapter):
    source_caller = "MethyLasso"
    model_family = CALLER_MODEL_FAMILY["MethyLasso"]
    coordinate_system = COORD_ONE_BASED_INCLUSIVE


class HmmdmAdapter(ExternalDmrAdapter):
    source_caller = "HMM-DM"
    model_family = CALLER_MODEL_FAMILY["HMM-DM"]


class CombPAdapter(ExternalDmrAdapter):
    source_caller = "comb-p"
    model_family = CALLER_MODEL_FAMILY["comb-p"]


class MetileneAdapter(ExternalDmrAdapter):
    source_caller = "metilene"
    model_family = CALLER_MODEL_FAMILY["metilene"]


ADAPTERS = {
    "generic_bed": GenericBedAdapter,
    "DSS": DSSAdapter,
    "dss": DSSAdapter,
    "methylKit": MethylKitAdapter,
    "methylkit": MethylKitAdapter,
    "dmrseq": DMRseqAdapter,
    "RADMeth": RADMethAdapter,
    "radmeth": RADMethAdapter,
    "MOABS": MOABSAdapter,
    "moabs": MOABSAdapter,
    "methylSig": MethylSigAdapter,
    "methylsig": MethylSigAdapter,
    "methylSig2": MethylSigAdapter,
    "methylsig2": MethylSigAdapter,
    "BSmooth": BSmoothAdapter,
    "bsmooth": BSmoothAdapter,
    "bsseq": BSmoothAdapter,
    "BiSeq": BiSeqAdapter,
    "biseq": BiSeqAdapter,
    "DMRcate": DMRcateAdapter,
    "dmrcate": DMRcateAdapter,
    "MethyLasso": MethyLassoAdapter,
    "methylasso": MethyLassoAdapter,
    "methyLasso": MethyLassoAdapter,
    "HMM-DM": HmmdmAdapter,
    "hmm-dm": HmmdmAdapter,
    "hmmdm": HmmdmAdapter,
    "comb-p": CombPAdapter,
    "combp": CombPAdapter,
    "CombP": CombPAdapter,
    "metilene": MetileneAdapter,
}


def adapter_for_caller(caller: str) -> type[ExternalDmrAdapter]:
    try:
        return ADAPTERS[caller]
    except KeyError as exc:
        raise ValueError(f"Unsupported caller: {caller}") from exc
