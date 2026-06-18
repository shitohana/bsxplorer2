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

from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from bsx2.analysis.intervals import reciprocal_overlap


CANONICAL_DMR_COLUMNS = [
    "dmr_id", "source_caller", "caller_version", "contrast_id", "condition_a", "condition_b",
    "context", "chrom", "start", "end", "strand", "direction", "delta", "p_value", "q_value",
    "n_sites", "n_cytosines", "mean_methylation_a", "mean_methylation_b", "source_file",
    "source_status", "method_notes",
]


CALLER_MODEL_FAMILY = {
    "dmrseq": "regional_fdr",
    "RADMeth": "beta_binomial",
    "MOABS": "beta_binomial",
    "methylSig": "beta_binomial",
    "methylSig2": "beta_binomial",
    "BSmooth": "smoothing",
    "BiSeq": "smoothing_beta_regression",
    "DMRcate": "limma_kernel_smoothing",
    "MethyLasso": "segmentation",
    "HMM-DM": "hmm",
    "comb-p": "pvalue_combination",
    "DSS": "beta_binomial",
    "methylKit": "fisher_logistic",
    "metilene": "segmentation",
    "generic_bed": "unknown",
    "internal": "internal",
}

# Coordinate conventions. Every adapter declares the convention of its native
# output; coordinates are normalized to BED-style 0-based half-open at import so
# reciprocal-overlap math (bsx2.analysis.intervals) is consistent everywhere.
COORD_ZERO_BASED_HALF_OPEN = "zero_based_half_open"
COORD_ONE_BASED_INCLUSIVE = "one_based_inclusive"

# Per-caller sign correction mapping each caller's native methylation
# difference onto the canonical ``condition_b - condition_a`` orientation
# (positive == higher methylation in condition_b == "hyper").
#
# Assumption: condition_a is the first/reference group passed to the caller.
# Callers that natively report ``group1 - group2`` (reference - other) therefore
# need a sign flip. DSS (``diff.Methy = mu1 - mu2``) and metilene
# (``mean_difference = mean_g1 - mean_g2``) are the well-documented reversed
# cases; all other callers default to +1. Override per dataset via the adapter's
# ``delta_sign`` argument if your group coding differs. When both
# ``mean_methylation_a`` and ``mean_methylation_b`` are present, direction is
# always derived from ``mean_b - mean_a`` and this map is not consulted.
CALLER_DELTA_SIGN = {
    "DSS": -1,
    "metilene": -1,
}

TIER_DESCRIPTIONS = {
    "TIER_1_CONSENSUS_VALIDATED": "Region is supported by multiple callers from multiple model families, has consistent direction, and passes validation audit.",
    "TIER_2_SINGLE_CALLER_VALIDATED": "Region is supported by one caller only, but has strong effect size and passes validation audit.",
    "TIER_3_MULTI_CALLER_QC_LIMITED": "Region is supported by multiple callers, but has validation or QC warnings.",
    "TIER_4_CALLER_DISCORDANT": "Callers disagree in effect direction or region interpretation.",
    "TIER_5_REJECTED_BY_AUDIT": "Region is unstable according to validation audit: direction changes, LOO changes direction, CI includes zero, or confidence class is LOW_CONFIDENCE.",
    "CANDIDATE_ONLY": "Region has insufficient support for strong interpretation.",
}

SUPPORT_CALLERS = [
    "internal",
    "DSS",
    "methylKit",
    "dmrseq",
    "metilene",
    "RADMeth",
    "MOABS",
    "methylSig",
    "BSmooth",
    "BiSeq",
    "DMRcate",
    "MethyLasso",
    "HMM-DM",
    "comb-p",
]

VALIDATION_DEFAULTS = {
    "final_confidence_class": "NOT_EVALUATED",
    "model_agreement_status": "NOT_EVALUATED",
    "direction_changed": False,
    "loo_direction_changed": False,
    "ci_includes_zero": False,
    "coverage_imbalance_score": pd.NA,
    "fraction_common_cpg": pd.NA,
    "relative_delta_shift": pd.NA,
}


def canonical_dmr_columns() -> list[str]:
    return list(CANONICAL_DMR_COLUMNS)


def caller_model_family(caller: Any) -> str:
    if caller is None or pd.isna(caller):
        return "unknown"
    text = str(caller)
    if text in CALLER_MODEL_FAMILY:
        return CALLER_MODEL_FAMILY[text]
    lower = {key.lower(): value for key, value in CALLER_MODEL_FAMILY.items()}
    return lower.get(text.lower(), "unknown")


def caller_delta_sign(caller: Any) -> int:
    """Sign correction that orients a caller's delta to ``condition_b - condition_a``."""
    if caller is None:
        return 1
    try:
        if pd.isna(caller):
            return 1
    except TypeError:
        pass
    text = str(caller)
    if text in CALLER_DELTA_SIGN:
        return CALLER_DELTA_SIGN[text]
    lower = {key.lower(): value for key, value in CALLER_DELTA_SIGN.items()}
    return lower.get(text.lower(), 1)


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


def normalize_interval_convention(df: pd.DataFrame, coordinate_system: str) -> pd.DataFrame:
    """Normalize ``start``/``end`` to BED-style 0-based half-open coordinates.

    A 1-based inclusive interval ``[start, end]`` maps to the 0-based half-open
    interval ``[start - 1, end)``; 0-based half-open input is returned unchanged.
    """
    if coordinate_system == COORD_ZERO_BASED_HALF_OPEN:
        return df
    if coordinate_system == COORD_ONE_BASED_INCLUSIVE:
        out = df.copy()
        out["start"] = pd.to_numeric(out["start"], errors="coerce") - 1
        return out
    raise ValueError(f"Unknown coordinate system: {coordinate_system}")


def harmonize_delta_sign(df: pd.DataFrame, sign: int) -> pd.DataFrame:
    """Flip the ``delta`` column onto the canonical ``condition_b - condition_a`` orientation."""
    if sign == 1 or "delta" not in df.columns:
        return df
    out = df.copy()
    out["delta"] = pd.to_numeric(out["delta"], errors="coerce") * sign
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


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value: Any) -> int:
    number = _safe_float(value)
    if number is None:
        return 0
    return int(number)


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except TypeError:
        pass
    text = str(value).strip().lower()
    return text in {"1", "true", "t", "yes", "y", "warn", "warning"}


def _direction_from_delta(delta: Any) -> str:
    number = _safe_float(delta)
    if number is None or number == 0:
        return "unknown"
    return "hyper" if number > 0 else "hypo"


def _normalize_direction(direction: Any, delta: Any = None) -> str:
    try:
        if direction is not None and not pd.isna(direction):
            text = str(direction).strip().lower()
            if text in {"hyper", "hypermethylated", "positive", "pos", "up", "+", "1", "a_gt_b", "b_lt_a"}:
                return "hyper"
            if text in {"hypo", "hypomethylated", "negative", "neg", "down", "-", "-1", "a_lt_b", "b_gt_a"}:
                return "hypo"
    except TypeError:
        pass
    return _direction_from_delta(delta)


def _effect_delta(row: dict[str, Any]) -> float | None:
    mean_a = _safe_float(row.get("mean_methylation_a"))
    mean_b = _safe_float(row.get("mean_methylation_b"))
    if mean_a is not None and mean_b is not None:
        return mean_b - mean_a
    return _safe_float(row.get("delta"))


def _support_column(caller: str) -> str:
    clean = "".join(ch if ch.isalnum() else "_" for ch in caller)
    clean = "_".join(part for part in clean.split("_") if part)
    return f"{clean}_support"


# Public aliases. These canonical helpers are the single implementation reused by
# the validation-framework workflows; do not duplicate the logic there.
safe_float = _safe_float
support_column = _support_column
direction_from_delta = _direction_from_delta


def _unique_sorted(values: list[Any]) -> list[str]:
    out = []
    for value in values:
        if value is None:
            continue
        try:
            if pd.isna(value):
                continue
        except TypeError:
            pass
        text = str(value)
        if text and text not in out:
            out.append(text)
    return sorted(out)


def compute_multi_caller_evidence_score(row: pd.Series | dict[str, Any]) -> int:
    score = 0
    n_callers = _safe_int(row.get("n_callers_supporting"))
    n_families = _safe_int(row.get("n_model_families_supporting"))
    fraction_same = _safe_float(row.get("fraction_same_direction"))
    n_q05 = _safe_int(row.get("n_callers_q05"))
    max_abs_delta = _safe_float(row.get("max_abs_delta"))

    if n_callers >= 3:
        score += 2
    elif n_callers == 2:
        score += 1

    if n_families >= 3:
        score += 2
    elif n_families == 2:
        score += 1

    if fraction_same == 1.0:
        score += 2
    elif fraction_same is not None and fraction_same >= 0.7:
        score += 1
    else:
        score -= 2

    if n_q05 >= 2:
        score += 2
    elif n_q05 == 1:
        score += 1

    if max_abs_delta is not None and max_abs_delta >= 0.2:
        score += 2
    elif max_abs_delta is not None and max_abs_delta >= 0.1:
        score += 1

    if _truthy(row.get("caller_conflict_flag")):
        score -= 3
    return score


def multi_caller_evidence_class(score: int, row: pd.Series | dict[str, Any]) -> str:
    if _truthy(row.get("caller_conflict_flag")):
        return "CONFLICTING_DIRECTIONS"
    if score >= 8:
        return "STRONG_MULTI_CALLER_EVIDENCE"
    if score >= 5:
        return "MODERATE_MULTI_CALLER_EVIDENCE"
    if score >= 2:
        return "LIMITED_MULTI_CALLER_EVIDENCE"
    return "WEAK_MULTI_CALLER_EVIDENCE"


def compute_validation_robustness_score(row: pd.Series | dict[str, Any]) -> int:
    score = 0
    confidence = str(row.get("final_confidence_class", "")).strip()
    model_status = str(row.get("model_agreement_status", "")).strip()
    coverage_imbalance = _safe_float(row.get("coverage_imbalance_score"))
    relative_shift = _safe_float(row.get("relative_delta_shift"))
    fraction_common = _safe_float(row.get("fraction_common_cpg"))

    if confidence == "HIGH_CONFIDENCE":
        score += 3
    elif confidence == "MODERATE_CONFIDENCE":
        score += 2
    elif confidence == "WEAK_CONFIDENCE":
        score += 1
    elif confidence == "LOW_CONFIDENCE":
        score -= 3

    if model_status == "GLMM_CONFIRMED":
        score += 2
    elif model_status == "GLM_ONLY":
        score += 1
    elif model_status in {"DIRECTION_DISCORDANT", "NOT_SIGNIFICANT"}:
        score -= 2

    if _truthy(row.get("direction_changed")):
        score -= 4
    if _truthy(row.get("loo_direction_changed")):
        score -= 4
    if _truthy(row.get("ci_includes_zero")):
        score -= 3

    if coverage_imbalance is not None and coverage_imbalance > 3:
        score -= 1
    if coverage_imbalance is not None and coverage_imbalance > 10:
        score -= 2
    if relative_shift is not None and relative_shift > 0.5:
        score -= 2
    if fraction_common is not None and fraction_common < 0.5:
        score -= 2
    return score


def assign_final_dmr_tier(row: pd.Series | dict[str, Any]) -> str:
    confidence = str(row.get("final_confidence_class", "NOT_EVALUATED")).strip()
    model_status = str(row.get("model_agreement_status", "NOT_EVALUATED")).strip()
    direction_consensus = str(row.get("direction_consensus", "unknown")).strip()
    n_callers = _safe_int(row.get("n_callers_supporting"))
    n_families = _safe_int(row.get("n_model_families_supporting"))
    max_abs_delta = _safe_float(row.get("max_abs_delta"))

    if _truthy(row.get("caller_conflict_flag")) or direction_consensus == "mixed":
        return "TIER_4_CALLER_DISCORDANT"

    if (
        confidence == "LOW_CONFIDENCE"
        or _truthy(row.get("direction_changed"))
        or _truthy(row.get("loo_direction_changed"))
        or _truthy(row.get("ci_includes_zero"))
    ):
        return "TIER_5_REJECTED_BY_AUDIT"

    if (
        n_callers >= 2
        and n_families >= 2
        and direction_consensus == "same"
        and confidence in {"HIGH_CONFIDENCE", "MODERATE_CONFIDENCE"}
        and model_status in {"GLMM_CONFIRMED", "GLM_ONLY"}
    ):
        return "TIER_1_CONSENSUS_VALIDATED"

    if (
        n_callers == 1
        and confidence in {"HIGH_CONFIDENCE", "MODERATE_CONFIDENCE"}
        and max_abs_delta is not None
        and max_abs_delta >= 0.2
    ):
        return "TIER_2_SINGLE_CALLER_VALIDATED"

    if n_callers >= 2:
        return "TIER_3_MULTI_CALLER_QC_LIMITED"

    return "CANDIDATE_ONLY"


def build_tier_reasons(row: pd.Series | dict[str, Any]) -> str:
    reasons = [
        f"supported_by_{_safe_int(row.get('n_callers_supporting'))}_callers",
        f"{_safe_int(row.get('n_model_families_supporting'))}_model_families",
        str(row.get("direction_consensus", "unknown")) + "_direction",
    ]
    confidence = str(row.get("final_confidence_class", "NOT_EVALUATED"))
    model_status = str(row.get("model_agreement_status", "NOT_EVALUATED"))
    if model_status and model_status != "NOT_EVALUATED":
        reasons.append(model_status)
    if confidence and confidence != "NOT_EVALUATED":
        reasons.append(confidence)
    if _truthy(row.get("caller_conflict_flag")):
        reasons.append("caller_direction_conflict")
    if _truthy(row.get("direction_changed")):
        reasons.append("direction_changed")
    if _truthy(row.get("loo_direction_changed")):
        reasons.append("loo_direction_changed")
    if _truthy(row.get("ci_includes_zero")):
        reasons.append("ci_includes_zero")
    coverage_imbalance = _safe_float(row.get("coverage_imbalance_score"))
    if coverage_imbalance is not None and coverage_imbalance > 10:
        reasons.append("very_high_coverage_imbalance")
    elif coverage_imbalance is not None and coverage_imbalance > 3:
        reasons.append("high_coverage_imbalance")
    relative_shift = _safe_float(row.get("relative_delta_shift"))
    if relative_shift is not None and relative_shift > 0.5:
        reasons.append("large_relative_delta_shift")
    fraction_common = _safe_float(row.get("fraction_common_cpg"))
    if fraction_common is not None and fraction_common < 0.5:
        reasons.append("low_common_cpg_fraction")
    return ";".join(reasons)


def _validation_lookup(validation_df: pd.DataFrame) -> dict[str, dict[str, Any]]:
    if validation_df.empty:
        return {}
    key_columns = [c for c in ("region_id", "harmonized_region_id", "dmr_id", "source_region_id", "representative_region_id") if c in validation_df.columns]
    lookup: dict[str, dict[str, Any]] = {}
    for row in validation_df.to_dict("records"):
        for column in key_columns:
            value = row.get(column)
            if value is None:
                continue
            try:
                if pd.isna(value):
                    continue
            except TypeError:
                pass
            text = str(value)
            if text:
                lookup.setdefault(text, row)
    return lookup


def _match_validation_row(row: pd.Series, lookup: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    for column in ("region_id", "harmonized_region_id", "representative_region_id"):
        value = row.get(column)
        if value is not None:
            match = lookup.get(str(value))
            if match is not None:
                return match
    source_ids = str(row.get("source_region_ids", "")).split(",")
    for source_id in source_ids:
        match = lookup.get(source_id.strip())
        if match is not None:
            return match
    return None


def add_dmr_tiers(support_df: pd.DataFrame, validation_df: pd.DataFrame | None = None) -> pd.DataFrame:
    out = support_df.copy()
    if out.empty:
        for column, value in VALIDATION_DEFAULTS.items():
            out[column] = value
        out["multi_caller_evidence_score"] = pd.Series(dtype="int64")
        out["validation_robustness_score"] = pd.Series(dtype="int64")
        out["final_dmr_tier"] = pd.Series(dtype="object")
        out["tier_reasons"] = pd.Series(dtype="object")
        return out

    if "region_id" not in out.columns:
        out["region_id"] = out.get("harmonized_region_id", pd.Series([f"harmonized_{i + 1}" for i in range(len(out))], index=out.index))

    for column, value in VALIDATION_DEFAULTS.items():
        if column not in out.columns:
            out[column] = value

    if validation_df is not None and not validation_df.empty:
        lookup = _validation_lookup(validation_df)
        validation_columns = [c for c in validation_df.columns if c not in {"chrom", "start", "end", "context"}]
        for column in validation_columns:
            if column not in out.columns:
                out[column] = pd.Series([pd.NA] * len(out), index=out.index, dtype="object")
            elif column in VALIDATION_DEFAULTS:
                out[column] = out[column].astype("object")
        for idx, row in out.iterrows():
            match = _match_validation_row(row, lookup)
            if match is None:
                continue
            for column in validation_columns:
                if column in {"region_id", "harmonized_region_id"}:
                    continue
                out.at[idx, column] = match.get(column)
            for column, value in VALIDATION_DEFAULTS.items():
                if column in match:
                    out.at[idx, column] = match.get(column)

    out["multi_caller_evidence_score"] = [compute_multi_caller_evidence_score(row) for _, row in out.iterrows()]
    out["multi_caller_evidence_class"] = [multi_caller_evidence_class(int(row["multi_caller_evidence_score"]), row) for _, row in out.iterrows()]
    out["validation_robustness_score"] = [compute_validation_robustness_score(row) for _, row in out.iterrows()]
    out["final_dmr_tier"] = [assign_final_dmr_tier(row) for _, row in out.iterrows()]
    out["tier_reasons"] = [build_tier_reasons(row) for _, row in out.iterrows()]
    return out


class _DisjointSet:
    """Union-find with path compression and union by rank."""

    def __init__(self, n: int) -> None:
        self._parent = list(range(n))
        self._rank = [0] * n

    def find(self, x: int) -> int:
        root = x
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[x] != root:
            self._parent[x], x = root, self._parent[x]
        return root

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self._rank[ra] < self._rank[rb]:
            ra, rb = rb, ra
        self._parent[rb] = ra
        if self._rank[ra] == self._rank[rb]:
            self._rank[ra] += 1


def cluster_records_by_overlap(
    records: list[dict[str, Any]],
    overlap_threshold: float,
) -> list[dict[str, Any]]:
    """Partition DMR records into connected components of the reciprocal-overlap graph.

    The graph ``G_tau`` connects two *original* intervals (same chrom+context)
    iff ``reciprocal_overlap >= overlap_threshold``; clusters are its connected
    components, computed with a per-group sweep line plus union-find. The result
    is deterministic and independent of input order, and the cluster envelope is
    derived only *after* the partition, so it never influences membership (no
    chaining via a growing envelope). See
    ``docs/dmr_validation_math_problems.md`` task 3 for the formal argument.
    """
    n = len(records)
    dsu = _DisjointSet(n)

    groups: dict[tuple[str, str], list[int]] = defaultdict(list)
    for idx, record in enumerate(records):
        groups[(str(record.get("chrom")), str(record.get("context", "")))].append(idx)

    edge_overlaps: dict[int, list[float]] = defaultdict(list)
    edges: list[tuple[int, int, float]] = []
    for member_indices in groups.values():
        member_indices.sort(key=lambda i: (records[i]["_start"], records[i]["_end"]))
        active: list[int] = []
        for i in member_indices:
            s_i = records[i]["_start"]
            e_i = records[i]["_end"]
            active = [j for j in active if records[j]["_end"] > s_i]
            for j in active:
                overlap = reciprocal_overlap(records[j]["_start"], records[j]["_end"], s_i, e_i)
                if overlap >= overlap_threshold and overlap > 0:
                    dsu.union(i, j)
                    edges.append((i, j, overlap))
            active.append(i)

    for i, _j, overlap in edges:
        edge_overlaps[dsu.find(i)].append(overlap)

    # Stable ordering: clusters are numbered by the first record (in flatten
    # order) that belongs to each component, keeping output deterministic.
    component_order: list[int] = []
    component_position: dict[int, int] = {}
    for idx in range(n):
        root = dsu.find(idx)
        if root not in component_position:
            component_position[root] = len(component_order)
            component_order.append(root)
    members: list[list[int]] = [[] for _ in component_order]
    for idx in range(n):
        members[component_position[dsu.find(idx)]].append(idx)

    clusters: list[dict[str, Any]] = []
    for position, root in enumerate(component_order):
        member_indices = members[position]
        first = records[member_indices[0]]
        clusters.append(
            {
                "chrom": first.get("chrom"),
                "context": first.get("context"),
                "start": min(records[i]["_start"] for i in member_indices),
                "end": max(records[i]["_end"] for i in member_indices),
                "records": [records[i]["_payload"] for i in member_indices],
                "overlaps": edge_overlaps.get(root, []),
                "evidence_class": first.get("evidence_class", pd.NA),
            }
        )
    return clusters


def build_caller_support_matrix(canonical_tables: list[pd.DataFrame], overlap_threshold: float = 0.5) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for table in canonical_tables:
        if table.empty:
            continue
        for row in table.to_dict("records"):
            row_start = _safe_float(row.get("start"))
            row_start = 0.0 if row_start is None else row_start
            row_end = _safe_float(row.get("end"))
            row_end = row_start if row_end is None else row_end
            caller = str(row.get("source_caller", "external"))
            delta = _effect_delta(row)
            records.append(
                {
                    "chrom": row.get("chrom"),
                    "context": row.get("context"),
                    "_start": row_start,
                    "_end": row_end,
                    "evidence_class": row.get("evidence_class", pd.NA),
                    "_payload": {
                        "caller": caller,
                        "family": caller_model_family(caller),
                        "dmr_id": row.get("dmr_id", pd.NA),
                        "q_value": _safe_float(row.get("q_value")),
                        "delta": delta,
                        "direction": _normalize_direction(row.get("direction", pd.NA), delta),
                    },
                }
            )
    clusters = cluster_records_by_overlap(records, overlap_threshold)
    rows = []
    for i, cluster in enumerate(clusters, start=1):
        records = cluster["records"]
        callers = _unique_sorted([record["caller"] for record in records])
        families = _unique_sorted([record["family"] for record in records])
        q_by_caller: dict[str, float] = {}
        for record in records:
            q_value = record["q_value"]
            if q_value is None:
                continue
            caller = record["caller"]
            q_by_caller[caller] = min(q_by_caller.get(caller, q_value), q_value)
        q_values = list(q_by_caller.values())
        deltas = [record["delta"] for record in records if record["delta"] is not None]
        abs_deltas = [abs(delta) for delta in deltas]
        direction_by_caller = {}
        for caller in callers:
            caller_directions = {record["direction"] for record in records if record["caller"] == caller and record["direction"] != "unknown"}
            if len(caller_directions) > 1:
                direction_by_caller[caller] = "mixed"
            elif len(caller_directions) == 1:
                direction_by_caller[caller] = next(iter(caller_directions))
            else:
                direction_by_caller[caller] = "unknown"
        known_directions = [direction for direction in direction_by_caller.values() if direction in {"hyper", "hypo"}]
        direction_counts = {direction: known_directions.count(direction) for direction in set(known_directions)}
        has_unknown_direction = any(direction == "unknown" for direction in direction_by_caller.values())
        has_within_caller_mixed_direction = any(direction == "mixed" for direction in direction_by_caller.values())
        if not direction_counts:
            direction_consensus = "unknown"
            fraction_same = pd.NA
            n_opposite = 0
        elif has_within_caller_mixed_direction or len(direction_counts) > 1:
            direction_consensus = "mixed"
            fraction_same = max(direction_counts.values()) / max(len(direction_by_caller), 1)
            n_opposite = min(direction_counts.values())
        elif has_unknown_direction:
            direction_consensus = "unknown"
            fraction_same = max(direction_counts.values()) / max(len(direction_by_caller), 1)
            n_opposite = 0
        else:
            direction_consensus = "same"
            fraction_same = 1.0
            n_opposite = 0

        source_ids = _unique_sorted([record["dmr_id"] for record in records])
        representative_region_id = source_ids[0] if source_ids else f"harmonized_{i}"
        row_out = {
            "harmonized_region_id": f"harmonized_{i}",
            "region_id": f"harmonized_{i}",
            "representative_region_id": representative_region_id,
            "source_region_ids": ",".join(source_ids),
            "chrom": cluster["chrom"],
            "start": cluster["start"],
            "end": cluster["end"],
            "context": cluster["context"],
            "supporting_callers": ",".join(callers),
            "supporting_model_families": ",".join(families),
            "n_callers_supporting": len(callers),
            "n_model_families_supporting": len(families),
            "direction_consensus": direction_consensus,
            "fraction_same_direction": fraction_same,
            "n_opposite_direction": n_opposite,
            "best_q_value": min(q_values) if q_values else pd.NA,
            "n_callers_q05": sum(1 for value in q_values if value <= 0.05),
            "n_callers_q10": sum(1 for value in q_values if value <= 0.10),
            "max_abs_delta": max(abs_deltas) if abs_deltas else pd.NA,
            "median_abs_delta": pd.Series(abs_deltas).median() if abs_deltas else pd.NA,
            "mean_abs_delta": pd.Series(abs_deltas).mean() if abs_deltas else pd.NA,
            "caller_conflict_flag": direction_consensus == "mixed",
            "mean_reciprocal_overlap": pd.Series(cluster["overlaps"]).mean() if cluster["overlaps"] else pd.NA,
            "min_reciprocal_overlap": min(cluster["overlaps"]) if cluster["overlaps"] else pd.NA,
            "evidence_class": cluster["evidence_class"],
        }
        for caller in SUPPORT_CALLERS:
            row_out[_support_column(caller)] = caller in callers
        rows.append(row_out)
    out = pd.DataFrame(rows)
    if out.empty:
        return add_dmr_tiers(out)
    support_columns = [column for column in out.columns if column.endswith("_support")]
    out[support_columns] = out[support_columns].fillna(False).astype(bool)
    return add_dmr_tiers(out)
