"""Streamlit viewer for BSX2 DMR evidence tables.

This app is a lightweight thesis/demo frontend. It reads existing TSV outputs,
summarizes and filters them, and never runs DMR calling, raw read processing,
or external caller execution.
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:  # Streamlit is optional for static tests and source-only delivery.
    import streamlit as st
except ModuleNotFoundError:  # pragma: no cover - exercised by import smoke only
    st = None


APP_TITLE = "BSX2 DMR Evidence Viewer"
APP_DESCRIPTION = (
    "Demo frontend for viewing BSX2 DMR evidence outputs. "
    "It does not run DMR calling."
)
MAX_FULL_READ_BYTES = 50 * 1024 * 1024
PREVIEW_ROWS = 20_000

CANONICAL_ALIASES = {
    "chrom": {"chrom", "chr", "chromosome", "seqname"},
    "start": {"start", "start_bp", "begin"},
    "end": {"end", "end_bp", "stop"},
    "delta": {"delta", "region_delta", "mean_delta", "beta_binom_delta"},
    "p_value": {"p_value", "pvalue", "pval", "region_p_value"},
    "q_value": {"q_value", "qvalue", "qval", "region_q_value"},
    "region_id": {"region_id", "dmr_id", "harmonized_region_id"},
}


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with common DMR/evidence column aliases normalized."""

    out = df.copy()
    lower_to_column = {str(col).lower(): col for col in out.columns}
    renames: dict[str, str] = {}
    for canonical, aliases in CANONICAL_ALIASES.items():
        if canonical in out.columns:
            continue
        for alias in aliases:
            source = lower_to_column.get(alias)
            if source is not None:
                renames[source] = canonical
                break
    return out.rename(columns=renames)


def read_tsv_source(source: Any, *, source_name: str = "uploaded file") -> tuple[pd.DataFrame, str]:
    """Read a TSV from a Streamlit upload object or a path-like value."""

    warning = ""
    if isinstance(source, (str, Path)):
        path = Path(source).expanduser()
        size = path.stat().st_size
        nrows = PREVIEW_ROWS if size > MAX_FULL_READ_BYTES else None
        if nrows:
            warning = f"{source_name} is large; loaded first {PREVIEW_ROWS} rows only."
        df = pd.read_csv(path, sep="\t", nrows=nrows)
        return normalize_columns(df), warning

    if hasattr(source, "seek"):
        source.seek(0)
    size = getattr(source, "size", None)
    nrows = PREVIEW_ROWS if size and size > MAX_FULL_READ_BYTES else None
    if nrows:
        warning = f"{source_name} is large; loaded first {PREVIEW_ROWS} rows only."
    df = pd.read_csv(source, sep="\t", nrows=nrows)
    return normalize_columns(df), warning


def numeric_series(df: pd.DataFrame, column: str | None) -> pd.Series:
    if column is None or column not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def first_present(df: pd.DataFrame | None, columns: list[str]) -> str | None:
    if df is None:
        return None
    for column in columns:
        if column in df.columns:
            return column
    return None


def negative_log10(values: pd.Series) -> pd.Series:
    finite = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan)
    positive = finite[finite > 0]
    floor = positive.min() if not positive.empty else 1e-300
    clipped = finite.clip(lower=max(float(floor), 1e-300))
    return -np.log10(clipped)


def apply_filters(
    df: pd.DataFrame,
    *,
    contexts: list[str] | None = None,
    evidence_classes: list[str] | None = None,
    q_threshold: float | None = None,
    abs_delta_threshold: float | None = None,
    caller_support_min: int | None = None,
    top_n: int | None = None,
) -> pd.DataFrame:
    """Apply viewer filters without mutating the source table."""

    out = df.copy()
    if contexts and "context" in out.columns:
        out = out[out["context"].astype(str).isin(contexts)]
    if evidence_classes and "evidence_class" in out.columns:
        out = out[out["evidence_class"].astype(str).isin(evidence_classes)]
    if q_threshold is not None:
        q_col = first_present(out, ["q_value", "region_q_value", "beta_binom_q_value"])
        if q_col:
            out = out[numeric_series(out, q_col) <= q_threshold]
    if abs_delta_threshold is not None:
        delta_col = first_present(out, ["delta", "region_delta", "mean_delta", "beta_binom_delta"])
        if delta_col:
            out = out[numeric_series(out, delta_col).abs() >= abs_delta_threshold]
    if caller_support_min is not None and "n_callers_supporting" in out.columns:
        out = out[numeric_series(out, "n_callers_supporting") >= caller_support_min]
    sort_col = first_present(out, ["q_value", "region_q_value", "beta_binom_q_value"])
    if sort_col:
        out = out.assign(_sort_q=numeric_series(out, sort_col)).sort_values("_sort_q").drop(columns="_sort_q")
    if top_n:
        out = out.head(int(top_n))
    return out


def table_to_tsv(df: pd.DataFrame) -> bytes:
    buffer = StringIO()
    df.to_csv(buffer, sep="\t", index=False)
    return buffer.getvalue().encode("utf-8")


def summary_metrics(
    dmr_df: pd.DataFrame | None,
    beta_df: pd.DataFrame | None = None,
    support_df: pd.DataFrame | None = None,
) -> dict[str, int | str]:
    if dmr_df is None:
        return {
            "total_regions": "not available",
            "significant_q05": "not available",
            "strong_regions": "not available",
            "moderate_regions": "not available",
            "weak_regions": "not available",
            "beta_binomial_confirmed": "not available",
            "external_supported": "not available",
        }
    q_col = first_present(dmr_df, ["q_value", "region_q_value"])
    evidence_col = first_present(dmr_df, ["evidence_class"])
    metrics: dict[str, int | str] = {
        "total_regions": int(len(dmr_df)),
        "significant_q05": int((numeric_series(dmr_df, q_col) < 0.05).sum()) if q_col else "not available",
        "strong_regions": "not available",
        "moderate_regions": "not available",
        "weak_regions": "not available",
        "beta_binomial_confirmed": "not available",
        "external_supported": "not available",
    }
    if evidence_col:
        classes = dmr_df[evidence_col].astype(str).str.lower()
        metrics["strong_regions"] = int((classes == "strong").sum())
        metrics["moderate_regions"] = int((classes == "moderate").sum())
        metrics["weak_regions"] = int((classes == "weak").sum())
    if beta_df is not None:
        if "significant_beta_binom" in beta_df.columns:
            metrics["beta_binomial_confirmed"] = int(beta_df["significant_beta_binom"].astype(str).str.lower().isin({"true", "1"}).sum())
        else:
            beta_q = first_present(beta_df, ["q_value", "beta_binom_q_value"])
            beta_delta = first_present(beta_df, ["delta", "beta_binom_delta"])
            if beta_q and beta_delta:
                metrics["beta_binomial_confirmed"] = int(
                    ((numeric_series(beta_df, beta_q) < 0.05) & (numeric_series(beta_df, beta_delta).abs() >= 0.2)).sum()
                )
    if support_df is not None and "n_callers_supporting" in support_df.columns:
        metrics["external_supported"] = int((numeric_series(support_df, "n_callers_supporting") > 0).sum())
    elif dmr_df is not None and "n_callers_supporting" in dmr_df.columns:
        metrics["external_supported"] = int((numeric_series(dmr_df, "n_callers_supporting") > 0).sum())
    return metrics


def render_metric_cards(metrics: dict[str, int | str]) -> None:
    cols = st.columns(4)
    cols[0].metric("Total regions", metrics["total_regions"])
    cols[1].metric("q < 0.05", metrics["significant_q05"])
    cols[2].metric("Strong", metrics["strong_regions"])
    cols[3].metric("Beta-binomial confirmed", metrics["beta_binomial_confirmed"])
    cols2 = st.columns(3)
    cols2[0].metric("Moderate", metrics["moderate_regions"])
    cols2[1].metric("Weak", metrics["weak_regions"])
    cols2[2].metric("External supported", metrics["external_supported"])


def render_bar_counts(df: pd.DataFrame, column: str, label: str) -> None:
    counts = df[column].astype(str).value_counts().rename_axis(label).reset_index(name="count")
    try:
        import plotly.express as px

        fig = px.bar(counts, x=label, y="count", template="plotly_dark")
        st.plotly_chart(fig, use_container_width=True)
    except Exception:
        st.bar_chart(counts.set_index(label))


def render_histogram(values: pd.Series, label: str) -> None:
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if clean.empty:
        st.info(f"{label} is not available.")
        return
    try:
        import plotly.express as px

        fig = px.histogram(pd.DataFrame({label: clean}), x=label, nbins=40, template="plotly_dark")
        st.plotly_chart(fig, use_container_width=True)
    except Exception:
        hist, edges = np.histogram(clean, bins=min(40, max(5, len(clean))))
        plot_df = pd.DataFrame({"bin": edges[:-1], "count": hist}).set_index("bin")
        st.bar_chart(plot_df)


def load_optional_table(label: str, help_text: str) -> tuple[pd.DataFrame | None, str]:
    upload = st.file_uploader(label, type=["tsv", "txt", "csv"], help=help_text)
    path_text = st.text_input(f"Optional path for {label}", value="", placeholder="Leave empty unless running locally")
    if upload is None and not path_text.strip():
        return None, ""
    try:
        source = upload if upload is not None else path_text.strip()
        return read_tsv_source(source, source_name=label)
    except Exception as exc:
        st.error(f"Could not read {label}: {exc}")
        return None, str(exc)


def inject_dark_theme() -> None:
    st.markdown(
        """
        <style>
        :root { color-scheme: dark; }
        .stApp { background: #101418; color: #eef2f5; }
        [data-testid="stSidebar"] { background: #151b22; }
        h1, h2, h3 { color: #f5f7fa; letter-spacing: 0; }
        .stMetric {
            background: #1b232c;
            border: 1px solid #2c3844;
            border-radius: 8px;
            padding: 12px;
        }
        div[data-testid="stDataFrame"] {
            border: 1px solid #2c3844;
            border-radius: 8px;
        }
        .muted-note { color: #aeb8c2; }
        </style>
        """,
        unsafe_allow_html=True,
    )


def main() -> None:
    if st is None:
        raise RuntimeError("Streamlit is required to run this viewer.")
    st.set_page_config(page_title=APP_TITLE, layout="wide", initial_sidebar_state="expanded")
    inject_dark_theme()

    st.sidebar.header("Inputs")
    dmr_df, dmr_warning = load_optional_table(
        "DMR/evidence TSV",
        "Examples: dmr_evidence_scores.tsv or dmr_region_count_tests.tsv",
    )
    beta_df, beta_warning = load_optional_table("Optional beta-binomial TSV", "Example: dmr_beta_binomial_tests.tsv")
    support_df, support_warning = load_optional_table("Optional caller support matrix TSV", "Example: dmr_caller_support_matrix.tsv")
    annotation_df, annotation_warning = load_optional_table("Optional annotation/enrichment TSV", "Optional contextual annotation table")

    for warning in [dmr_warning, beta_warning, support_warning, annotation_warning]:
        if warning:
            st.sidebar.warning(warning)

    st.title(APP_TITLE)
    st.caption(APP_DESCRIPTION)
    st.markdown(
        '<p class="muted-note">Upload existing BSX2 evidence outputs to inspect regions, filters, '
        "summary metrics, and caller support. The viewer does not write uploaded data to the repository.</p>",
        unsafe_allow_html=True,
    )

    if dmr_df is None:
        st.info("Upload a DMR/evidence TSV to start.")
        st.stop()

    context_values = sorted(dmr_df["context"].dropna().astype(str).unique().tolist()) if "context" in dmr_df.columns else []
    evidence_values = (
        sorted(dmr_df["evidence_class"].dropna().astype(str).unique().tolist()) if "evidence_class" in dmr_df.columns else []
    )
    st.sidebar.header("Filters")
    selected_contexts = st.sidebar.multiselect("Context", context_values, default=context_values)
    selected_classes = st.sidebar.multiselect("Evidence class", evidence_values, default=evidence_values)
    q_threshold = st.sidebar.slider("q-value threshold", min_value=0.0, max_value=1.0, value=0.05, step=0.01)
    abs_delta_threshold = st.sidebar.slider("abs(delta) threshold", min_value=0.0, max_value=1.0, value=0.2, step=0.05)
    caller_support_min = st.sidebar.number_input("Minimum caller support", min_value=0, value=0, step=1)
    top_n = st.sidebar.number_input("Top N", min_value=10, max_value=100_000, value=500, step=10)

    filtered = apply_filters(
        dmr_df,
        contexts=selected_contexts if context_values else None,
        evidence_classes=selected_classes if evidence_values else None,
        q_threshold=q_threshold,
        abs_delta_threshold=abs_delta_threshold,
        caller_support_min=int(caller_support_min),
        top_n=int(top_n),
    )

    overview, table_tab, distributions, caller_support, method_notes = st.tabs(
        ["Overview", "DMR Table", "Evidence Distributions", "Caller Support", "Method Notes"]
    )

    with overview:
        render_metric_cards(summary_metrics(filtered, beta_df, support_df))
        st.subheader("Filtered Preview")
        st.dataframe(filtered.head(50), use_container_width=True)

    with table_tab:
        st.dataframe(filtered, use_container_width=True, height=520)
        st.download_button(
            "Download filtered TSV",
            data=table_to_tsv(filtered),
            file_name="bsx2_filtered_dmr_evidence.tsv",
            mime="text/tab-separated-values",
        )

    with distributions:
        if "evidence_class" in filtered.columns:
            st.subheader("Evidence Class Counts")
            render_bar_counts(filtered, "evidence_class", "evidence_class")
        delta_col = first_present(filtered, ["delta", "region_delta", "mean_delta", "beta_binom_delta"])
        st.subheader("Delta Distribution")
        render_histogram(numeric_series(filtered, delta_col), "delta")
        q_col = first_present(filtered, ["q_value", "region_q_value", "beta_binom_q_value"])
        st.subheader("-log10(q-value) Distribution")
        render_histogram(negative_log10(numeric_series(filtered, q_col)), "-log10(q_value)")

    with caller_support:
        support_source = support_df if support_df is not None else filtered
        if support_source is not None and "n_callers_supporting" in support_source.columns:
            st.subheader("Caller Support Counts")
            render_bar_counts(support_source, "n_callers_supporting", "n_callers_supporting")
            st.dataframe(support_source.head(500), use_container_width=True)
        else:
            st.info("Caller support matrix is not available.")
        if annotation_df is not None:
            st.subheader("Annotation/Enrichment Preview")
            st.dataframe(annotation_df.head(500), use_container_width=True)

    with method_notes:
        st.markdown(
            """
            - This frontend is a viewer only; it does not run DMR calling.
            - It does not run raw read processing, Bismark, DSS, methylKit, dmrseq, or metilene.
            - q-values may come from different methods and should be interpreted with their source method.
            - External caller harmonization standardizes schema, not caller-specific statistical assumptions.
            - Beta-binomial aggregated mode is a GLM validation layer, not a random-effect GLMM.
            - Uploaded files are read in memory for display and are not stored in the repository.
            """
        )


if __name__ == "__main__":
    main()
