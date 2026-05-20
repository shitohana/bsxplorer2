"""Streamlit viewer for BSX2 DMR evidence tables.

This app is a lightweight thesis/demo frontend. It reads existing TSV outputs,
summarizes and filters them, and never runs DMR calling, raw read processing,
or external caller execution.
"""

from __future__ import annotations

from html import escape
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
    "region_id": {"region_id", "dmr_id", "harmonized_region_id", "id", "name"},
    "chrom": {"chrom", "chr", "chromosome", "seqname"},
    "start": {"start", "start_bp", "begin"},
    "end": {"end", "end_bp", "stop"},
    "context": {"context", "methylation_context"},
    "delta": {
        "delta",
        "region_delta",
        "mean_delta",
        "beta_binom_delta",
        "meth_diff",
        "meth_diff_percent",
        "mean_methylation_difference",
    },
    "p_value": {"p_value", "p", "pvalue", "pval", "region_p_value", "beta_binom_p_value"},
    "q_value": {"q_value", "q", "fdr", "qvalue", "qval", "region_q_value", "beta_binom_q_value", "padj"},
    "evidence_class": {"evidence_class", "class", "evidence", "evidence_level"},
    "n_callers_supporting": {"n_callers_supporting", "caller_support_count", "external_caller_support_count"},
    "source_caller": {"source_caller", "caller", "dmr_caller"},
}

INPUT_GUIDANCE = {
    "main": {
        "label": "DMR/evidence TSV",
        "required": True,
        "purpose": "Primary table for region browsing, filtering, summary cards, and plots.",
        "recommended": ["dmr_evidence_scores.tsv", "dmr_regions.tsv", "dmr_region_count_tests.tsv"],
        "best": "dmr_evidence_scores.tsv",
    },
    "beta": {
        "label": "Beta-binomial validation TSV",
        "required": False,
        "purpose": "Adds complementary beta-binomial confirmation for regions.",
        "recommended": [
            "dmr_beta_binomial_tests.tsv",
            "dmr_beta_binomial_tests_top1000_adjusted.tsv",
            "dmr_beta_binomial_tests_top100_real.tsv",
        ],
        "best": "dmr_beta_binomial_tests.tsv",
    },
    "support": {
        "label": "Caller support matrix TSV",
        "required": False,
        "purpose": "Shows which callers support each region, such as BSX2, DSS, methylKit, dmrseq, or metilene.",
        "recommended": ["dmr_caller_support_matrix.tsv", "external_caller_support_matrix.tsv"],
        "best": "dmr_caller_support_matrix.tsv",
    },
    "annotation": {
        "label": "Annotation / enrichment TSV",
        "required": False,
        "purpose": "Adds biological annotation such as promoter, gene_body, intergenic, or enrichment summaries.",
        "recommended": [
            "differential_methylated_features.tsv",
            "dmr_annotation_enrichment.tsv",
            "plant_region_methylation_summary.tsv",
        ],
        "best": "differential_methylated_features.tsv",
    },
}


def _column_key(column: object) -> str:
    return str(column).strip().lower().replace("-", "_").replace(".", "_").replace(" ", "_")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with common DMR/evidence column aliases normalized."""

    out = df.copy()
    keyed_columns = {_column_key(col): col for col in out.columns}
    renames: dict[str, str] = {}
    existing_keys = {_column_key(col) for col in out.columns}
    for canonical, aliases in CANONICAL_ALIASES.items():
        if canonical in out.columns or canonical in existing_keys:
            continue
        for alias in aliases:
            source = keyed_columns.get(_column_key(alias))
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
    keyed = {_column_key(col): col for col in df.columns}
    for column in columns:
        direct = keyed.get(_column_key(column))
        if direct is not None:
            return direct
    return None


def negative_log10(values: pd.Series) -> pd.Series:
    finite = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan)
    positive = finite[finite > 0]
    floor = positive.min() if not positive.empty else 1e-300
    clipped = finite.clip(lower=max(float(floor), 1e-300))
    return -np.log10(clipped)


def display_metric_value(value: int | str) -> str:
    if isinstance(value, str) and value == "not available":
        return "N/A"
    return str(value)


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
    context_col = first_present(out, ["context", "methylation_context"])
    evidence_col = first_present(out, ["evidence_class", "class"])
    if contexts and context_col:
        out = out[out[context_col].astype(str).isin(contexts)]
    if evidence_classes and evidence_col:
        out = out[out[evidence_col].astype(str).isin(evidence_classes)]
    if q_threshold is not None:
        q_col = first_present(out, ["q_value", "q", "fdr", "region_q_value", "beta_binom_q_value"])
        if q_col:
            out = out[numeric_series(out, q_col) <= q_threshold]
    if abs_delta_threshold is not None:
        delta_col = first_present(out, ["delta", "mean_delta", "region_delta", "beta_binom_delta"])
        if delta_col:
            out = out[numeric_series(out, delta_col).abs() >= abs_delta_threshold]
    support_col = first_present(out, ["n_callers_supporting", "caller_support_count"])
    if caller_support_min is not None and support_col:
        out = out[numeric_series(out, support_col) >= caller_support_min]
    sort_col = first_present(out, ["q_value", "q", "fdr", "region_q_value", "beta_binom_q_value"])
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
    *,
    q_threshold: float = 0.05,
) -> dict[str, int | str]:
    if dmr_df is None:
        return {
            "total_regions": "not available",
            "significant_q": "not available",
            "strong_regions": "not available",
            "moderate_regions": "not available",
            "weak_regions": "not available",
            "beta_binomial_confirmed": "not available",
            "external_supported": "not available",
        }
    q_col = first_present(dmr_df, ["q_value", "q", "fdr", "region_q_value"])
    evidence_col = first_present(dmr_df, ["evidence_class", "class"])
    metrics: dict[str, int | str] = {
        "total_regions": int(len(dmr_df)),
        "significant_q": int((numeric_series(dmr_df, q_col) < q_threshold).sum()) if q_col else "not available",
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
        sig_col = first_present(beta_df, ["significant_beta_binom"])
        if sig_col:
            metrics["beta_binomial_confirmed"] = int(beta_df[sig_col].astype(str).str.lower().isin({"true", "1", "yes"}).sum())
        else:
            beta_q = first_present(beta_df, ["q_value", "q", "fdr", "beta_binom_q_value"])
            beta_delta = first_present(beta_df, ["delta", "mean_delta", "beta_binom_delta"])
            if beta_q and beta_delta:
                metrics["beta_binomial_confirmed"] = int(
                    ((numeric_series(beta_df, beta_q) < q_threshold) & (numeric_series(beta_df, beta_delta).abs() >= 0.2)).sum()
                )
    support_col = first_present(support_df, ["n_callers_supporting", "caller_support_count"])
    dmr_support_col = first_present(dmr_df, ["n_callers_supporting", "caller_support_count"])
    if support_df is not None and support_col:
        metrics["external_supported"] = int((numeric_series(support_df, support_col) > 0).sum())
    elif dmr_support_col:
        metrics["external_supported"] = int((numeric_series(dmr_df, dmr_support_col) > 0).sum())
    return metrics


def render_metric_cards(metrics: dict[str, int | str], *, q_threshold: float) -> None:
    cols = st.columns(4)
    cols[0].metric("Total regions", display_metric_value(metrics["total_regions"]))
    cols[1].metric(f"q < {q_threshold:g}", display_metric_value(metrics["significant_q"]))
    cols[2].metric("Strong", display_metric_value(metrics["strong_regions"]))
    cols[3].metric("Moderate", display_metric_value(metrics["moderate_regions"]))
    cols2 = st.columns(3)
    cols2[0].metric("Weak", display_metric_value(metrics["weak_regions"]))
    cols2[1].metric("Beta-binomial confirmed", display_metric_value(metrics["beta_binomial_confirmed"]))
    cols2[2].metric("External caller supported", display_metric_value(metrics["external_supported"]))


def render_bar_counts(df: pd.DataFrame, column: str, label: str) -> None:
    counts = df[column].astype(str).value_counts().rename_axis(label).reset_index(name="count")
    try:
        import plotly.express as px

        fig = px.bar(counts, x=label, y="count", template="plotly_dark", color_discrete_sequence=["#5cc8ff"])
        fig.update_layout(paper_bgcolor="#111820", plot_bgcolor="#111820", font_color="#edf5ff")
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

        fig = px.histogram(
            pd.DataFrame({label: clean}),
            x=label,
            nbins=40,
            template="plotly_dark",
            color_discrete_sequence=["#7b8cff"],
        )
        fig.update_layout(paper_bgcolor="#111820", plot_bgcolor="#111820", font_color="#edf5ff")
        st.plotly_chart(fig, use_container_width=True)
    except Exception:
        hist, edges = np.histogram(clean, bins=min(40, max(5, len(clean))))
        plot_df = pd.DataFrame({"bin": edges[:-1], "count": hist}).set_index("bin")
        st.bar_chart(plot_df)


def file_pills_html(files: list[str], best: str | None = None) -> str:
    pills = []
    for filename in files:
        label = escape(filename)
        marker = '<span class="file-pill-tag">recommended</span>' if filename == best else ""
        pills.append(f'<span class="file-pill">{label}{marker}</span>')
    return "\n".join(pills)


def render_input_card(config: dict[str, Any]) -> None:
    st.markdown(
        f"""
        <div class="input-card">
          <div class="input-card-title">{escape(config["label"])}</div>
          <div class="input-card-purpose"><b>Purpose:</b> {escape(config["purpose"])}</div>
          <div class="input-card-subtitle">Recommended</div>
          <div class="file-list">
            {file_pills_html(config["recommended"], config.get("best"))}
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def load_table_input(key: str, *, use_local_paths: bool) -> tuple[pd.DataFrame | None, str]:
    config = INPUT_GUIDANCE[key]
    label = config["label"]
    render_input_card(config)
    upload = None
    path_text = ""
    if use_local_paths:
        path_text = st.text_input(
            f"Optional local path for {label}",
            value="",
            key=f"{key}_path",
            placeholder="Local run only: paste a TSV path",
            help="Local path mode is intended only for local runs.",
        )
        st.caption("Local path mode is intended only for local runs.")
    else:
        upload = st.file_uploader(
            label,
            type=["tsv", "txt", "csv"],
            key=f"{key}_upload",
            label_visibility="collapsed",
        )
    if upload is None and not path_text.strip():
        return None, ""
    try:
        source = upload if upload is not None else path_text.strip()
        return read_tsv_source(source, source_name=label)
    except Exception as exc:
        st.error(f"Could not read {label}: {exc}")
        return None, str(exc)


def render_guide_block() -> None:
    st.markdown(
        """
        <div class="guide-card">
          <div class="guide-title">How to use this viewer</div>
          <ol>
            <li>Upload the main DMR/evidence TSV.</li>
            <li>Optionally upload beta-binomial results.</li>
            <li>Optionally upload caller support matrix.</li>
            <li>Optionally upload annotation/enrichment results.</li>
            <li>Use filters to inspect regions.</li>
          </ol>
          <div class="guide-grid">
            <div><b>Main DMR/evidence</b><br><code>dmr_evidence_scores.tsv</code> recommended<br><code>dmr_regions.tsv</code><br><code>dmr_region_count_tests.tsv</code></div>
            <div><b>Beta-binomial</b><br><code>dmr_beta_binomial_tests.tsv</code><br><code>dmr_beta_binomial_tests_top1000_adjusted.tsv</code></div>
            <div><b>Caller support</b><br><code>dmr_caller_support_matrix.tsv</code><br><code>external_caller_support_matrix.tsv</code></div>
            <div><b>Annotation</b><br><code>differential_methylated_features.tsv</code><br><code>dmr_annotation_enrichment.tsv</code><br><code>plant_region_methylation_summary.tsv</code></div>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def inject_dark_theme() -> None:
    st.markdown(
        """
        <style>
        :root {
            color-scheme: dark;
            --bsx-bg: #070B12;
            --bsx-bg-2: #0B1220;
            --bsx-panel: #111827;
            --bsx-panel-2: #0F172A;
            --bsx-border: #263449;
            --bsx-border-2: #36506a;
            --bsx-text: #E5EDF7;
            --bsx-muted: #9FB0C5;
            --bsx-accent: #38BDF8;
            --bsx-accent-2: #8b7bff;
            --bsx-accent-muted: #1E3A5F;
            --bsx-warn: #f0ba55;
            --bsx-good: #57d49b;
        }

        * {
            box-sizing: border-box;
        }

        html, body, #root, .stApp,
        [data-testid="stAppViewContainer"] {
            background: radial-gradient(circle at 20% 0%, rgba(92, 200, 255, 0.10), transparent 28%),
                        linear-gradient(180deg, #070B12 0%, #0B1220 48%, #070B12 100%) !important;
            color: var(--bsx-text) !important;
        }

        [data-testid="stHeader"],
        [data-testid="stToolbar"] {
            background: rgba(7, 11, 18, 0.96) !important;
            color: var(--bsx-text) !important;
        }

        #MainMenu, footer, [data-testid="stDecoration"] {
            visibility: hidden !important;
            height: 0 !important;
        }

        .block-container {
            padding-top: 2.2rem !important;
            padding-bottom: 3rem !important;
            max-width: 1280px !important;
        }

        [data-testid="stSidebar"] {
            background: linear-gradient(180deg, #0B1220 0%, #0F172A 100%) !important;
            border-right: 1px solid var(--bsx-border) !important;
        }

        [data-testid="stSidebar"] * {
            max-width: 100%;
        }

        [data-testid="stSidebarContent"] {
            background: transparent !important;
            color: var(--bsx-text) !important;
        }

        h1, h2, h3, h4, h5, h6,
        [data-testid="stMarkdownContainer"],
        label, p, li, span {
            color: var(--bsx-text);
            letter-spacing: 0;
            overflow-wrap: anywhere;
        }

        p, li, span, div {
            overflow-wrap: anywhere;
        }

        a { color: var(--bsx-accent); }
        pre,
        code {
            color: #d8f2ff !important;
            background: #111827 !important;
            border: 1px solid #31465c;
            border-radius: 5px;
            padding: 0.05rem 0.3rem;
            white-space: pre-wrap;
            overflow-wrap: anywhere;
        }

        .muted-note {
            color: var(--bsx-muted);
            margin-top: -0.4rem;
            margin-bottom: 1rem;
        }

        .guide-card {
            background: linear-gradient(135deg, rgba(24, 34, 49, 0.96), rgba(17, 24, 34, 0.96));
            border: 1px solid var(--bsx-border-2);
            border-radius: 8px;
            padding: 1rem 1.1rem;
            margin: 1rem 0 1.2rem 0;
            box-shadow: 0 12px 28px rgba(0, 0, 0, 0.24);
        }

        .guide-title {
            color: var(--bsx-text);
            font-size: 1.1rem;
            font-weight: 700;
            margin-bottom: 0.5rem;
        }

        .guide-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 0.75rem;
            margin-top: 0.8rem;
        }

        .guide-grid > div {
            background: rgba(8, 12, 18, 0.46);
            border: 1px solid var(--bsx-border);
            border-radius: 8px;
            padding: 0.75rem;
            color: var(--bsx-muted);
        }

        .sidebar-section-title {
            margin-top: 1rem;
            margin-bottom: 0.45rem;
            color: #d8f2ff;
            font-weight: 700;
            font-size: 0.95rem;
            text-transform: uppercase;
            letter-spacing: 0.04em;
        }

        .input-card {
            width: 100%;
            max-width: 100%;
            background: linear-gradient(180deg, #111827 0%, #0F172A 100%);
            border: 1px solid var(--bsx-border);
            border-radius: 8px;
            padding: 0.72rem 0.78rem;
            margin: 0.55rem 0 0.42rem 0;
            box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.025);
        }

        .input-card-title {
            color: var(--bsx-text);
            font-size: 0.96rem;
            font-weight: 750;
            line-height: 1.25;
            margin-bottom: 0.35rem;
        }

        .input-card-purpose {
            color: var(--bsx-muted);
            font-size: 0.82rem;
            line-height: 1.35;
            margin-bottom: 0.55rem;
        }

        .input-card-purpose b {
            color: #c8d8ea;
        }

        .input-card-subtitle {
            color: #c8d8ea;
            font-size: 0.74rem;
            font-weight: 700;
            line-height: 1.2;
            margin-bottom: 0.34rem;
            text-transform: uppercase;
            letter-spacing: 0.04em;
        }

        .file-list {
            display: flex;
            flex-direction: column;
            flex-wrap: wrap;
            gap: 0.34rem;
            width: 100%;
            max-width: 100%;
        }

        .file-pill {
            display: block;
            width: 100%;
            max-width: 100%;
            color: #d8f2ff;
            background: rgba(30, 58, 95, 0.60);
            border: 1px solid rgba(56, 189, 248, 0.24);
            border-radius: 7px;
            padding: 0.34rem 0.45rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
            font-size: 0.72rem;
            line-height: 1.35;
            white-space: normal;
            overflow-wrap: anywhere;
            word-break: break-word;
        }

        .file-pill-tag {
            display: inline-block;
            color: #08111c;
            background: var(--bsx-accent);
            border-radius: 999px;
            padding: 0.05rem 0.32rem;
            margin-left: 0.35rem;
            font-family: inherit;
            font-size: 0.62rem;
            line-height: 1.2;
            font-weight: 800;
            vertical-align: middle;
        }

        .upload-start-card {
            background: rgba(92, 200, 255, 0.08);
            border: 1px solid rgba(92, 200, 255, 0.38);
            border-radius: 8px;
            padding: 1rem 1.1rem;
            color: #d8f2ff;
        }

        [data-testid="stMetric"] {
            background: linear-gradient(180deg, #162130 0%, #121a24 100%);
            border: 1px solid var(--bsx-border);
            border-radius: 8px;
            padding: 0.9rem 1rem;
            box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.03);
        }

        [data-testid="stMetric"] label,
        [data-testid="stMetric"] [data-testid="stMetricLabel"] {
            color: var(--bsx-muted) !important;
        }

        [data-testid="stMetric"] [data-testid="stMetricValue"] {
            color: var(--bsx-text) !important;
            font-weight: 700;
        }

        [data-testid="stFileUploader"] {
            background: #0F172A !important;
            border: 1px dashed var(--bsx-border) !important;
            border-radius: 8px !important;
            padding: 0.32rem !important;
            margin-bottom: 0.85rem !important;
        }

        [data-testid="stFileUploaderDropzone"],
        [data-testid="stFileUploader"] section {
            background: #0B1220 !important;
            border: 1px solid var(--bsx-border) !important;
            border-radius: 8px !important;
            color: var(--bsx-text) !important;
            min-height: 4.3rem !important;
            padding: 0.52rem !important;
        }

        [data-testid="stFileUploaderDropzone"] * {
            color: var(--bsx-text) !important;
            overflow-wrap: anywhere;
        }

        [data-testid="stFileUploader"] small,
        [data-testid="stFileUploader"] span,
        [data-testid="stFileUploader"] p {
            color: var(--bsx-muted) !important;
            font-size: 0.78rem !important;
        }

        input, textarea,
        [data-baseweb="input"] > div,
        [data-baseweb="select"] > div,
        [data-baseweb="base-input"],
        [data-testid="stTextInput"] input,
        [data-testid="stNumberInput"] input {
            background: #101923 !important;
            color: var(--bsx-text) !important;
            border-color: var(--bsx-border) !important;
            caret-color: var(--bsx-accent) !important;
        }

        [data-testid="stTextInput"] {
            margin-bottom: 0.75rem !important;
        }

        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] span,
        [data-testid="stCheckbox"] p {
            color: var(--bsx-text) !important;
        }

        [data-testid="stCheckbox"] [data-testid="stMarkdownContainer"] p {
            font-size: 0.86rem !important;
        }

        [data-baseweb="tag"] {
            background: rgba(92, 200, 255, 0.18) !important;
            color: var(--bsx-text) !important;
            border: 1px solid rgba(92, 200, 255, 0.35) !important;
        }

        [data-testid="stSlider"] [role="slider"] {
            background: var(--bsx-accent) !important;
            border-color: #d8f2ff !important;
        }

        [data-testid="stTabs"] {
            background: transparent !important;
        }

        [data-testid="stTabs"] button {
            color: var(--bsx-muted) !important;
            background: #111a25 !important;
            border-radius: 8px 8px 0 0 !important;
            border: 1px solid var(--bsx-border) !important;
            margin-right: 0.18rem !important;
        }

        [data-testid="stTabs"] button[aria-selected="true"] {
            color: var(--bsx-text) !important;
            background: linear-gradient(180deg, #1d2b3b, #142031) !important;
            border-bottom-color: var(--bsx-accent) !important;
        }

        [data-testid="stDataFrame"],
        div[data-testid="stDataFrame"] {
            background: #0e151f !important;
            border: 1px solid var(--bsx-border) !important;
            border-radius: 8px !important;
            overflow: hidden;
        }

        button,
        [data-testid="stDownloadButton"] button,
        [data-testid="stBaseButton-secondary"],
        [data-testid="stBaseButton-primary"] {
            background: linear-gradient(180deg, #22364b, #18283a) !important;
            color: var(--bsx-text) !important;
            border: 1px solid var(--bsx-border-2) !important;
            border-radius: 8px !important;
        }

        button:hover,
        [data-testid="stDownloadButton"] button:hover {
            border-color: var(--bsx-accent) !important;
            color: #ffffff !important;
        }

        [data-testid="stAlert"] {
            background: #121d29 !important;
            color: var(--bsx-text) !important;
            border: 1px solid var(--bsx-border-2) !important;
            border-radius: 8px !important;
        }

        [data-testid="stAlert"] * {
            color: var(--bsx-text) !important;
        }

        hr {
            border-color: var(--bsx-border) !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_caller_support_summary(support_df: pd.DataFrame | None, filtered_df: pd.DataFrame) -> None:
    support_source = support_df if support_df is not None else filtered_df
    support_col = first_present(support_source, ["n_callers_supporting", "caller_support_count"])
    caller_col = first_present(support_source, ["source_caller", "caller"])
    support_flag_cols = [col for col in support_source.columns if _column_key(col).endswith("_support")]

    if support_col:
        st.subheader("Caller Support Counts")
        render_bar_counts(support_source, support_col, support_col)
        values = numeric_series(support_source, support_col)
        c1, c2, c3 = st.columns(3)
        c1.metric("Regions in support table", len(support_source))
        c2.metric("Any caller support", int((values > 0).sum()))
        c3.metric("Max callers", display_metric_value(int(values.max()) if not values.dropna().empty else "not available"))
    elif support_flag_cols:
        st.subheader("Support Flag Summary")
        summary = pd.DataFrame(
            {"caller": support_flag_cols, "supported_regions": [int(support_source[col].astype(bool).sum()) for col in support_flag_cols]}
        )
        st.dataframe(summary, use_container_width=True)
    else:
        st.info("Caller support matrix is not available or has no recognized support columns.")

    if caller_col:
        st.subheader("Source Caller Summary")
        render_bar_counts(support_source, caller_col, caller_col)

    st.subheader("Support Table Preview")
    st.dataframe(support_source.head(500), use_container_width=True)


def main() -> None:
    if st is None:
        raise RuntimeError("Streamlit is required to run this viewer.")
    st.set_page_config(page_title=APP_TITLE, layout="wide", initial_sidebar_state="expanded")
    inject_dark_theme()

    with st.sidebar:
        st.title("BSX2 Viewer")
        st.markdown('<div class="sidebar-section-title">Required input</div>', unsafe_allow_html=True)
        use_local_paths = st.checkbox("Use local file paths instead of uploads", value=False)
        dmr_df, dmr_warning = load_table_input("main", use_local_paths=use_local_paths)

        st.markdown('<div class="sidebar-section-title">Optional supporting inputs</div>', unsafe_allow_html=True)
        beta_df, beta_warning = load_table_input("beta", use_local_paths=use_local_paths)
        support_df, support_warning = load_table_input("support", use_local_paths=use_local_paths)
        annotation_df, annotation_warning = load_table_input("annotation", use_local_paths=use_local_paths)

        for warning in [dmr_warning, beta_warning, support_warning, annotation_warning]:
            if warning:
                st.warning(warning)

    st.title(APP_TITLE)
    st.caption(APP_DESCRIPTION)
    st.markdown(
        '<p class="muted-note">Upload existing BSX2 evidence outputs to inspect regions, filters, '
        "summary metrics, and caller support. The viewer does not write uploaded data to the repository.</p>",
        unsafe_allow_html=True,
    )
    render_guide_block()

    if dmr_df is None:
        st.markdown(
            '<div class="upload-start-card"><b>Upload a DMR/evidence TSV to start.</b><br>'
            "Best first file: <code>dmr_evidence_scores.tsv</code>. "
            "You can also use <code>dmr_regions.tsv</code> or <code>dmr_region_count_tests.tsv</code>.</div>",
            unsafe_allow_html=True,
        )
        st.stop()

    context_col = first_present(dmr_df, ["context", "methylation_context"])
    evidence_col = first_present(dmr_df, ["evidence_class", "class"])
    context_values = sorted(dmr_df[context_col].dropna().astype(str).unique().tolist()) if context_col else []
    evidence_values = sorted(dmr_df[evidence_col].dropna().astype(str).unique().tolist()) if evidence_col else []

    with st.sidebar:
        st.markdown('<div class="sidebar-section-title">Filters</div>', unsafe_allow_html=True)
        selected_contexts = st.multiselect("Context", context_values, default=context_values, disabled=not bool(context_values))
        selected_classes = st.multiselect(
            "Evidence class",
            evidence_values,
            default=evidence_values,
            disabled=not bool(evidence_values),
        )
        q_threshold = st.slider("q-value threshold", min_value=0.0, max_value=1.0, value=0.05, step=0.01)
        abs_delta_threshold = st.slider("abs(delta) threshold", min_value=0.0, max_value=1.0, value=0.2, step=0.05)
        caller_support_min = st.number_input("Minimum caller support", min_value=0, value=0, step=1)
        top_n = st.number_input("Top N", min_value=10, max_value=100_000, value=500, step=10)

    recognized_notes = []
    for required_label, column_group in [
        ("q-value", ["q_value", "q", "fdr", "region_q_value"]),
        ("delta", ["delta", "mean_delta", "region_delta"]),
        ("evidence class", ["evidence_class", "class"]),
    ]:
        if not first_present(dmr_df, column_group):
            recognized_notes.append(required_label)
    if recognized_notes:
        st.warning("Some optional columns are not available: " + ", ".join(recognized_notes) + ". Related filters or plots are disabled.")

    filtered = apply_filters(
        dmr_df,
        contexts=selected_contexts if context_values else None,
        evidence_classes=selected_classes if evidence_values else None,
        q_threshold=q_threshold,
        abs_delta_threshold=abs_delta_threshold,
        caller_support_min=int(caller_support_min),
        top_n=int(top_n),
    )

    st.subheader("Summary")
    render_metric_cards(summary_metrics(filtered, beta_df, support_df, q_threshold=q_threshold), q_threshold=q_threshold)

    overview, table_tab, distributions, caller_support, method_notes = st.tabs(
        ["Overview", "DMR Table", "Evidence Distributions", "Caller Support", "Method Notes"]
    )

    with overview:
        st.subheader("Filtered Preview")
        if filtered.empty:
            st.info("No regions pass the current filters.")
        else:
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
        if evidence_col and evidence_col in filtered.columns:
            st.subheader("Evidence Class Counts")
            render_bar_counts(filtered, evidence_col, "evidence_class")
        else:
            st.info("Evidence class counts are not available because no evidence_class/class column was found.")

        delta_col = first_present(filtered, ["delta", "region_delta", "mean_delta", "beta_binom_delta"])
        st.subheader("Delta Distribution")
        render_histogram(numeric_series(filtered, delta_col), "delta")

        q_col = first_present(filtered, ["q_value", "q", "fdr", "region_q_value", "beta_binom_q_value"])
        st.subheader("-log10(q-value) Distribution")
        render_histogram(negative_log10(numeric_series(filtered, q_col)), "-log10(q_value)")

    with caller_support:
        render_caller_support_summary(support_df, filtered)
        if annotation_df is not None:
            st.subheader("Annotation/Enrichment Preview")
            st.dataframe(annotation_df.head(500), use_container_width=True)

    with method_notes:
        st.markdown(
            """
            - Viewer only; no DMR calling.
            - It does not run raw read processing, Bismark, DSS, methylKit, dmrseq, or metilene.
            - q-values may come from different methods and should be interpreted with their source method.
            - External caller harmonization standardizes schema, not statistical assumptions.
            - Beta-binomial validation is a complementary evidence layer.
            - Beta-binomial aggregated mode is a GLM validation layer, not a random-effect GLMM.
            - Uploaded files are read in memory for display and are not written into the repository.
            """
        )


if __name__ == "__main__":
    main()
