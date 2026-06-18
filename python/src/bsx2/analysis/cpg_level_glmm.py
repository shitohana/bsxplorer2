"""CpG-level GLMM confirmatory validation for selected DMR candidates.

The layer is intentionally confirmatory.  It fits region-wise models only for a
selected top-N set of predefined DMR candidates and does not replace the main
DMR statistical model or external genome-wide DMR callers.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from .coverage_set_qc import build_cpg_coverage_qc


CPG_LEVEL_GLMM_COLUMNS = [
    "region_id",
    "context",
    "n_cpg",
    "n_samples",
    "n_rows",
    "effect_logit",
    "delta_methylation",
    "logLik_full",
    "logLik_null",
    "lrt_stat",
    "p_value",
    "q_value",
    "converged",
    "singular",
    "model_status",
    "warning",
    "coverage_set_mode",
    "n_cpg_union",
    "n_cpg_common",
    "fraction_common",
    "coverage_qc_status",
]


def _find_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def read_region_cpg_counts(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"region_id", "cpg_id", "sample_id", "mC", "uC", "total"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"region_cpg_counts table is missing required columns: {', '.join(missing)}")
    if "context" not in df.columns:
        df["context"] = "NA"
    for column in ("mC", "uC", "total"):
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0)
    return df


def read_design_table(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    sample_col = _find_column(df, ("sample_id", "sample"))
    if sample_col is None:
        raise ValueError("design table must contain sample_id")
    if sample_col != "sample_id":
        df = df.rename(columns={sample_col: "sample_id"})
    df["sample_id"] = df["sample_id"].astype(str)
    return df


def read_dmr_evidence_for_glmm(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    region_col = _find_column(df, ("region_id", "dmr_id"))
    if region_col is None:
        raise ValueError("DMR evidence table must contain region_id or dmr_id")
    if region_col != "region_id":
        df = df.rename(columns={region_col: "region_id"})
    return df


def resolve_rscript(rscript: str | Path | None = None) -> str | None:
    if rscript:
        candidate = Path(rscript)
        if candidate.exists():
            return str(candidate)
        found = shutil.which(str(rscript))
        if found:
            return found
        return None
    return shutil.which("Rscript")


def rscript_available(rscript: str | Path | None = None) -> bool:
    return resolve_rscript(rscript) is not None


def glmmTMB_available(rscript: str | Path | None = None) -> bool:
    resolved = resolve_rscript(rscript)
    if resolved is None:
        return False
    command = [
        resolved,
        "-e",
        "suppressPackageStartupMessages(library(glmmTMB)); cat('ok\\n')",
    ]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    except Exception:
        return False
    return result.returncode == 0 and "ok" in result.stdout


def bh_qvalues(p_values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(p_values, errors="coerce")
    q_values = pd.Series(pd.NA, index=p_values.index, dtype="Float64")
    valid = numeric.dropna()
    m = len(valid)
    if m == 0:
        return q_values
    ordered = valid.sort_values(ascending=False)
    running = 1.0
    for rank_from_high, (idx, p_value) in enumerate(ordered.items(), start=1):
        rank = m - rank_from_high + 1
        running = min(running, float(p_value) * m / rank)
        q_values.loc[idx] = min(running, 1.0)
    return q_values


def _evidence_sort_key(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    q_col = _find_column(out, ("region_q_value", "q_value", "qvalue", "padj", "fdr"))
    score_col = _find_column(out, ("final_evidence_score", "evidence_score", "score"))
    class_col = _find_column(out, ("evidence_class",))
    priority = {
        "strong": 0,
        "moderate": 1,
        "weak": 2,
        "candidate_only": 3,
        "qc_limited": 4,
    }
    out["_class_priority"] = out[class_col].map(priority).fillna(99) if class_col else 99
    out["_q"] = pd.to_numeric(out[q_col], errors="coerce").fillna(1.0) if q_col else 1.0
    out["_score"] = -pd.to_numeric(out[score_col], errors="coerce").fillna(0.0) if score_col else 0.0
    return out.sort_values(["_class_priority", "_score", "_q"])


def select_top_regions(
    cpg_counts_df: pd.DataFrame,
    evidence_df: Optional[pd.DataFrame] = None,
    *,
    top_n: int = 100,
) -> list[str]:
    if evidence_df is not None and not evidence_df.empty:
        evidence = evidence_df.copy()
        if "region_id" not in evidence.columns and "dmr_id" in evidence.columns:
            evidence = evidence.rename(columns={"dmr_id": "region_id"})
        evidence = _evidence_sort_key(evidence)
        ordered = [str(v) for v in evidence["region_id"].dropna().astype(str).drop_duplicates()]
        return ordered[: int(top_n)]
    return list(cpg_counts_df["region_id"].astype(str).drop_duplicates().head(int(top_n)))


def _delta_methylation(
    region_df: pd.DataFrame,
    design_df: pd.DataFrame,
    condition_column: str,
    case_label: Optional[str],
    control_label: Optional[str],
) -> float | None:
    merged = region_df.merge(design_df[["sample_id", condition_column]], on="sample_id", how="left")
    if merged[condition_column].isna().any():
        return None
    groups = list(merged[condition_column].astype(str).drop_duplicates())
    if case_label is None or control_label is None:
        if len(groups) != 2:
            return None
        control_label, case_label = sorted(groups)[:2]
    case = merged[merged[condition_column].astype(str) == str(case_label)]
    control = merged[merged[condition_column].astype(str) == str(control_label)]
    if case.empty or control.empty:
        return None
    case_total = case["total"].sum()
    control_total = control["total"].sum()
    if case_total <= 0 or control_total <= 0:
        return None
    return float(case["mC"].sum() / case_total - control["mC"].sum() / control_total)


def _validate_region(
    region_df: pd.DataFrame,
    design_df: pd.DataFrame,
    condition_column: str,
    *,
    min_cpg: int,
    min_replicates_per_group: int,
    min_total: int,
    max_zero_coverage_fraction: float,
) -> tuple[str, str]:
    if condition_column not in design_df.columns:
        return "missing_condition", f"design table is missing {condition_column}"
    filtered = region_df[region_df["total"] >= min_total].copy()
    if filtered["cpg_id"].nunique() < min_cpg:
        return "insufficient_cpg", f"fewer than {min_cpg} CpG records passed filters"
    zero_fraction = float((region_df["total"] <= 0).mean()) if len(region_df) else 1.0
    if zero_fraction > max_zero_coverage_fraction:
        return "too_many_zero_coverage", "zero-coverage fraction exceeds threshold"
    merged = filtered.merge(design_df[["sample_id", condition_column]], on="sample_id", how="left")
    if merged[condition_column].isna().any():
        return "missing_condition", "some CpG rows do not match design samples"
    replicate_counts = merged[["sample_id", condition_column]].drop_duplicates().groupby(condition_column).size()
    if len(replicate_counts) < 2 or (replicate_counts < min_replicates_per_group).any():
        return "insufficient_replicates", f"requires at least {min_replicates_per_group} samples per group"
    return "ready", ""


def _unavailable_rows(
    counts_df: pd.DataFrame,
    selected_regions: list[str],
    status: str,
    warning: str,
    design_df: pd.DataFrame,
    condition_column: str,
    case_label: Optional[str],
    control_label: Optional[str],
    coverage_set_mode: str = "per_sample",
    coverage_region_qc: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    qc = coverage_region_qc.set_index("region_id") if coverage_region_qc is not None and not coverage_region_qc.empty else None
    rows = []
    for region_id in selected_regions:
        region = counts_df[counts_df["region_id"].astype(str) == region_id]
        delta = _delta_methylation(region, design_df, condition_column, case_label, control_label)
        qc_row = qc.loc[region_id] if qc is not None and region_id in qc.index else {}
        rows.append({
            "region_id": region_id,
            "context": ";".join(sorted(region["context"].astype(str).dropna().unique())) if not region.empty else "NA",
            "n_cpg": int(region["cpg_id"].nunique()) if not region.empty else 0,
            "n_samples": int(region["sample_id"].nunique()) if not region.empty else 0,
            "n_rows": int(len(region)),
            "effect_logit": pd.NA,
            "delta_methylation": delta,
            "logLik_full": pd.NA,
            "logLik_null": pd.NA,
            "lrt_stat": pd.NA,
            "p_value": pd.NA,
            "q_value": pd.NA,
            "converged": False,
            "singular": pd.NA,
            "model_status": status,
            "warning": warning,
            "coverage_set_mode": coverage_set_mode,
            "n_cpg_union": qc_row.get("n_cpg_union", int(region["cpg_id"].nunique()) if not region.empty else 0),
            "n_cpg_common": qc_row.get("n_cpg_common", pd.NA),
            "fraction_common": qc_row.get("fraction_common", pd.NA),
            "coverage_qc_status": qc_row.get("coverage_qc_status", "not_evaluated"),
        })
    return pd.DataFrame(rows, columns=CPG_LEVEL_GLMM_COLUMNS)


def _run_region_glmm(
    region_df: pd.DataFrame,
    design_df: pd.DataFrame,
    *,
    region_id: str,
    condition_column: str,
    covariates: list[str],
    case_label: Optional[str],
    control_label: Optional[str],
    temp_dir: Path,
    rscript: str | Path | None = None,
) -> dict[str, object]:
    resolved_rscript = resolve_rscript(rscript)
    if resolved_rscript is None:
        raise RuntimeError("Rscript is unavailable")
    temp_dir.mkdir(parents=True, exist_ok=True)
    model_df = region_df.merge(design_df, on="sample_id", how="left").copy()
    model_df["condition_factor"] = model_df[condition_column].astype(str)
    for covariate in covariates:
        if covariate not in model_df.columns:
            raise ValueError(f"covariate is missing from design table: {covariate}")
    input_path = temp_dir / f"{region_id}.glmm_input.tsv"
    output_path = temp_dir / f"{region_id}.glmm_output.tsv"
    model_df.to_csv(input_path, sep="\t", index=False)
    covariate_terms = " + ".join(covariates)
    rhs_null = covariate_terms if covariate_terms else "1"
    rhs_full = "condition_factor" + (f" + {covariate_terms}" if covariate_terms else "")
    r_code = "; ".join(
        [
            "suppressPackageStartupMessages(library(glmmTMB))",
            f"d <- read.delim('{input_path.as_posix()}', check.names=FALSE)",
            "d$condition_factor <- factor(d$condition_factor)",
            (
                f"full <- glmmTMB(cbind(mC, uC) ~ {rhs_full} + "
                "(1 | cpg_id) + (1 | sample_id), "
                "family=betabinomial(link='logit'), data=d)"
            ),
            (
                f"null <- glmmTMB(cbind(mC, uC) ~ {rhs_null} + "
                "(1 | cpg_id) + (1 | sample_id), "
                "family=betabinomial(link='logit'), data=d)"
            ),
            "ll_full <- as.numeric(logLik(full))",
            "ll_null <- as.numeric(logLik(null))",
            "df_full <- attr(logLik(full), 'df')",
            "df_null <- attr(logLik(null), 'df')",
            "lrt <- 2 * (ll_full - ll_null)",
            "p <- pchisq(lrt, df=max(1, df_full - df_null), lower.tail=FALSE)",
            "coefs <- summary(full)$coefficients$cond",
            "effect <- NA",
            "cond_rows <- grep('^condition_factor', rownames(coefs))",
            "if (length(cond_rows) > 0) effect <- coefs[cond_rows[1], 'Estimate']",
            "conv <- isTRUE(full$fit$convergence == 0)",
            "singular <- !isTRUE(full$sdr$pdHess)",
            (
                "out <- data.frame(effect_logit=effect, logLik_full=ll_full, "
                "logLik_null=ll_null, lrt_stat=lrt, p_value=p, converged=conv, "
                "singular=singular)"
            ),
            (
                f"write.table(out, file='{output_path.as_posix()}', "
                "sep='\\t', quote=FALSE, row.names=FALSE)"
            ),
        ]
    )
    command = [resolved_rscript, "-e", r_code]
    result = subprocess.run(command, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip() or "glmmTMB failed")
    out = pd.read_csv(output_path, sep="\t").iloc[0].to_dict()
    out["delta_methylation"] = _delta_methylation(region_df, design_df, condition_column, case_label, control_label)
    return out


def run_cpg_level_glmm_validation(
    cpg_counts_df: pd.DataFrame,
    design_df: pd.DataFrame,
    *,
    evidence_df: Optional[pd.DataFrame] = None,
    top_n: int = 100,
    condition_column: str = "condition",
    case_label: Optional[str] = None,
    control_label: Optional[str] = None,
    covariates: Optional[list[str]] = None,
    min_cpg: int = 3,
    min_replicates_per_group: int = 2,
    min_total: int = 1,
    max_zero_coverage_fraction: float = 0.5,
    qvalue_method: str = "BH",
    temp_dir: str | Path | None = None,
    force_glmm_unavailable: bool = False,
    rscript: str | Path | None = None,
    coverage_set_mode: str = "per_sample",
    min_coverage: int | None = None,
    min_covered_per_group: int | dict[str, int] | None = None,
    min_common_cpgs: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if qvalue_method.upper() != "BH":
        raise ValueError("only BH q-value correction is currently supported")
    if coverage_set_mode not in {"per_sample", "common"}:
        raise ValueError("coverage_set_mode must be per_sample or common")
    covariates = covariates or []
    counts = cpg_counts_df.copy()
    counts["region_id"] = counts["region_id"].astype(str)
    counts["sample_id"] = counts["sample_id"].astype(str)
    selected = select_top_regions(counts, evidence_df, top_n=top_n)
    warnings: list[dict[str, object]] = []
    coverage_region_qc: pd.DataFrame | None = None

    if coverage_set_mode == "common":
        common_cpg, coverage_region_qc = build_cpg_coverage_qc(
            counts,
            design_df,
            condition_col=condition_column,
            min_coverage=min_total if min_coverage is None else min_coverage,
            min_covered_per_group=min_covered_per_group,
            min_common_cpgs=min_common_cpgs,
        )
        common_keys = common_cpg[common_cpg["is_common"].astype(bool)][["region_id", "cpg_id"]]
        counts = counts.merge(common_keys, on=["region_id", "cpg_id"], how="inner")
        warnings.append(
            {
                "warning_type": "coverage_set_mode_common",
                "severity": "info",
                "message": "CpG-level GLMM input was filtered to common CpG set before validation.",
            }
        )

    if force_glmm_unavailable or not rscript_available(rscript) or not glmmTMB_available(rscript):
        warning = "Rscript/glmmTMB unavailable; CpG-level GLMM confirmatory model was not fitted."
        warnings.append({"warning_type": "glmmTMB_unavailable", "severity": "warning", "message": warning})
        return (
            _unavailable_rows(
                counts,
                selected,
                "glmmTMB_unavailable",
                warning,
                design_df,
                condition_column,
                case_label,
                control_label,
                coverage_set_mode=coverage_set_mode,
                coverage_region_qc=coverage_region_qc,
            ),
            pd.DataFrame(warnings),
        )

    if temp_dir is None:
        warning = "temp_dir is required for GLMM execution to keep runtime files under the caller-controlled output directory."
        warnings.append({"warning_type": "glmm_temp_dir_missing", "severity": "error", "message": warning})
        return (
            _unavailable_rows(
                counts,
                selected,
                "glmm_temp_dir_missing",
                warning,
                design_df,
                condition_column,
                case_label,
                control_label,
                coverage_set_mode=coverage_set_mode,
                coverage_region_qc=coverage_region_qc,
            ),
            pd.DataFrame(warnings),
        )

    rows: list[dict[str, object]] = []
    qc_index = coverage_region_qc.set_index("region_id") if coverage_region_qc is not None and not coverage_region_qc.empty else None
    for region_id in selected:
        region = counts[counts["region_id"].astype(str) == region_id].copy()
        context = ";".join(sorted(region["context"].astype(str).dropna().unique())) if not region.empty else "NA"
        qc_row = qc_index.loc[region_id] if qc_index is not None and region_id in qc_index.index else {}
        status, warning = _validate_region(
            region,
            design_df,
            condition_column,
            min_cpg=min_cpg,
            min_replicates_per_group=min_replicates_per_group,
            min_total=min_total,
            max_zero_coverage_fraction=max_zero_coverage_fraction,
        )
        base = {
            "region_id": region_id,
            "context": context,
            "n_cpg": int(region["cpg_id"].nunique()) if not region.empty else 0,
            "n_samples": int(region["sample_id"].nunique()) if not region.empty else 0,
            "n_rows": int(len(region)),
            "effect_logit": pd.NA,
            "delta_methylation": _delta_methylation(region, design_df, condition_column, case_label, control_label),
            "logLik_full": pd.NA,
            "logLik_null": pd.NA,
            "lrt_stat": pd.NA,
            "p_value": pd.NA,
            "q_value": pd.NA,
            "converged": False,
            "singular": pd.NA,
            "model_status": status,
            "warning": warning,
            "coverage_set_mode": coverage_set_mode,
            "n_cpg_union": qc_row.get("n_cpg_union", int(region["cpg_id"].nunique()) if not region.empty else 0),
            "n_cpg_common": qc_row.get("n_cpg_common", pd.NA),
            "fraction_common": qc_row.get("fraction_common", pd.NA),
            "coverage_qc_status": qc_row.get(
                "coverage_qc_status",
                "per_sample_not_evaluated" if coverage_set_mode == "per_sample" else "not_evaluable",
            ),
        }
        if coverage_set_mode == "common" and base["coverage_qc_status"] != "pass_common_cpg":
            base["model_status"] = "insufficient_common_cpgs"
            base["warning"] = "region does not pass common-CpG coverage QC"
            rows.append(base)
            warnings.append(
                {
                    "warning_type": "insufficient_common_cpgs",
                    "region_id": region_id,
                    "severity": "warning",
                    "message": base["warning"],
                }
            )
            continue
        if status != "ready":
            rows.append(base)
            warnings.append({"warning_type": status, "region_id": region_id, "severity": "warning", "message": warning})
            continue
        try:
            result = _run_region_glmm(
                region[region["total"] >= min_total],
                design_df,
                region_id=region_id,
                condition_column=condition_column,
                covariates=covariates,
                case_label=case_label,
                control_label=control_label,
                temp_dir=Path(temp_dir),
                rscript=rscript,
            )
            base.update(result)
            p_value = pd.to_numeric(pd.Series([base.get("p_value")]), errors="coerce").iloc[0]
            lrt_stat = pd.to_numeric(pd.Series([base.get("lrt_stat")]), errors="coerce").iloc[0]
            converged = bool(base.get("converged"))
            singular = bool(base.get("singular"))
            if pd.isna(p_value) or pd.isna(lrt_stat):
                base["model_status"] = "model_no_lrt"
                base["warning"] = "GLMM fitted but likelihood-ratio statistic or p-value is unavailable"
            elif not converged:
                base["model_status"] = "model_not_converged"
                base["warning"] = "GLMM fit did not converge"
            elif singular:
                base["model_status"] = "model_singular"
                base["warning"] = "GLMM fit has a singular/non-positive-definite Hessian"
            else:
                base["model_status"] = "ok"
                base["warning"] = ""
            rows.append(base)
        except Exception as exc:
            base["model_status"] = "model_error"
            base["warning"] = str(exc)
            rows.append(base)
            warnings.append({"warning_type": "model_error", "region_id": region_id, "severity": "warning", "message": str(exc)})

    out = pd.DataFrame(rows, columns=CPG_LEVEL_GLMM_COLUMNS)
    out["q_value"] = pd.NA
    for _context, idx in out.groupby("context", dropna=False).groups.items():
        out.loc[idx, "q_value"] = bh_qvalues(out.loc[idx, "p_value"])
    return out, pd.DataFrame(warnings)
