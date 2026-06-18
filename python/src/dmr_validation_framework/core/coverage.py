"""Coverage-set and common-CpG QC helpers for DMR validation."""

from pathlib import Path

import pandas as pd

from bsx2.analysis.coverage_set_qc import *  # noqa: F401,F403
from bsx2.analysis.coverage_set_qc import aggregate_counts_per_sample, normalize_cpg_counts_table
from dmr_validation_framework.core.columns import find_col
from dmr_validation_framework.core.io import read_table


def aggregate_counts_per_sample_chunked(
    cpg_counts: str | Path | pd.DataFrame,
    *,
    design: pd.DataFrame | None = None,
    condition_column: str = "condition",
    min_coverage: int = 1,
    chunk_size: int = 1_000_000,
) -> pd.DataFrame:
    """Aggregate region/sample CpG counts without loading large text files at once.

    Parquet/Feather inputs are read as a whole because pandas does not expose a
    simple chunked reader for them. TSV/CSV/text inputs are streamed by
    ``chunk_size`` rows.
    """
    if isinstance(cpg_counts, pd.DataFrame):
        counts = normalize_cpg_counts_table(cpg_counts)
        if design is not None:
            sample_col = find_col(design, ["sample_id", "sample"])
            cond_col = find_col(design, [condition_column, "condition", "group", "treatment"])
            if not sample_col or not cond_col:
                raise ValueError("design table lacks sample/condition columns")
            design_norm = design.rename(columns={sample_col: "sample_id", cond_col: "condition"})[["sample_id", "condition"]]
            counts = counts.drop(columns=["condition"]).merge(design_norm, on="sample_id", how="left")
        return aggregate_counts_per_sample(counts, min_coverage=min_coverage)

    path = Path(cpg_counts)
    suffix = path.suffix.lower()
    if suffix in {".parquet", ".feather", ".arrow"}:
        reader_df = read_table(path)
        return aggregate_counts_per_sample_chunked(
            reader_df,
            design=design,
            condition_column=condition_column,
            min_coverage=min_coverage,
            chunk_size=chunk_size,
        )

    sep = "," if suffix == ".csv" or ".csv" in [s.lower() for s in path.suffixes] else "\t"
    design_norm = None
    if design is not None:
        sample_col = find_col(design, ["sample_id", "sample"])
        cond_col = find_col(design, [condition_column, "condition", "group", "treatment"])
        if not sample_col or not cond_col:
            raise ValueError("design table lacks sample/condition columns")
        design_norm = design.rename(columns={sample_col: "sample_id", cond_col: "condition"})[["sample_id", "condition"]]

    count_parts: list[pd.DataFrame] = []
    cpg_key_parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep=sep, chunksize=int(chunk_size), low_memory=False):
        counts = normalize_cpg_counts_table(chunk)
        if design_norm is not None:
            counts = counts.drop(columns=["condition"]).merge(design_norm, on="sample_id", how="left")
        covered = counts[pd.to_numeric(counts["total"], errors="coerce").fillna(0) >= int(min_coverage)]
        if covered.empty:
            continue
        part = (
            covered.groupby(["region_id", "sample_id", "condition"], dropna=False)
            .agg(
                M_per_sample=("mC", "sum"),
                U_per_sample=("uC", "sum"),
                N_per_sample=("total", "sum"),
            )
            .reset_index()
        )
        count_parts.append(part)
        cpg_key_parts.append(
            covered[["region_id", "sample_id", "condition", "cpg_id"]].drop_duplicates()
        )
    if not count_parts:
        return pd.DataFrame(
            columns=["region_id", "sample_id", "condition", "M_per_sample", "U_per_sample", "N_per_sample", "n_cpg_per_sample"]
        )
    combined = pd.concat(count_parts, ignore_index=True)
    out = (
        combined.groupby(["region_id", "sample_id", "condition"], dropna=False)
        .agg(
            M_per_sample=("M_per_sample", "sum"),
            U_per_sample=("U_per_sample", "sum"),
            N_per_sample=("N_per_sample", "sum"),
        )
        .reset_index()
    )
    cpg_counts = (
        pd.concat(cpg_key_parts, ignore_index=True)
        .drop_duplicates()
        .groupby(["region_id", "sample_id", "condition"], dropna=False)
        .size()
        .rename("n_cpg_per_sample")
        .reset_index()
    )
    out = out.merge(cpg_counts, on=["region_id", "sample_id", "condition"], how="left")
    out["mu_hat_per_sample"] = out["M_per_sample"] / out["N_per_sample"].where(out["N_per_sample"] > 0)
    return out
