#!/usr/bin/env python
"""Single interactive HTML report for the DMR validation framework.

Builds one self-contained ``index.html`` (Panel + HoloViews/Plotly) that
gathers every plot the framework produces into a tabbed dashboard:

* Callers         - DMR counts, |delta| distribution, length distribution
                    (fixed-window callers are detected and annotated, never
                    silently shown as a zero-spread box).
* Consensus       - caller-combination support, final DMR tiers,
                    reciprocal-overlap threshold sensitivity.
* Confirmatory    - aggregated GLM vs CpG-level GLMM (hover = region_id),
                    Brown-corrected comb-p calibration QQ.
* Robustness      - replicate-delta bootstrap caterpillar with CIs.
* Validation      - final confidence-class distribution.

The report consumes tables already written by ``caller-summary``,
``critical-validation`` and ``glm-glmm-validation``. Each section degrades
gracefully: a missing input becomes an explanatory note instead of an error,
so the report is always produced.

This is a presentation layer only; it does not recompute statistics or call
any external DMR caller.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from dmr_validation_framework.core.io import SEARCH_ROOTS_ENV, find_preferred_file, read_table
from dmr_validation_framework.reports.concordance import compute_caller_concordance
from dmr_validation_framework.reports.funnel import compute_funnel_stages


# ---------------------------------------------------------------------------
# Input discovery
# ---------------------------------------------------------------------------

INPUT_NAMES = {
    "caller_rows": ["caller_dmr_rows_for_plots.tsv"],
    "caller_summary": ["caller_context_summary.tsv"],
    "consensus": [
        "dmr_caller_support_matrix.tsv",
        "dmr_method_support_matrix.tsv",
        "caller_support_matrix.tsv",
        "consensus_tiers.tsv",
        "dmr_regions_annotated.tsv",
    ],
    "membership": ["upset_membership.tsv"],
    "glm_glmm": ["glm_vs_glmm_comparison.tsv"],
    "bootstrap": ["delta_bootstrap_ci.tsv"],
    "overlap_sensitivity": ["overlap_threshold_sensitivity.tsv"],
    "confidence": ["delta_weighting_sensitivity.tsv"],
}


# Filename tokens that mark a thesis-figure PNG as annotation/downstream
# (as opposed to the core caller/validation/GLMM figures, which the interactive
# tabs already render natively).
ANNOTATION_TOKENS = (
    "annotation",
    "te_",
    "te_like",
    "te_enrichment",
    "go_",
    "occupancy",
    "projection",
    "metagene",
    "feature_class",
    "context_distribution",
    "seqname",
    "expression",
    "loci",
)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True, help="Directory for the combined report.")
    parser.add_argument("--caller-rows", help="caller_dmr_rows_for_plots.tsv (per-DMR rows).")
    parser.add_argument("--caller-summary", help="caller_context_summary.tsv.")
    parser.add_argument("--consensus", help="Caller support / consensus tier table.")
    parser.add_argument("--glm-glmm", help="glm_vs_glmm_comparison.tsv.")
    parser.add_argument("--bootstrap", help="delta_bootstrap_ci.tsv.")
    parser.add_argument("--overlap-sensitivity", help="overlap_threshold_sensitivity.tsv.")
    parser.add_argument("--confidence", help="delta_weighting_sensitivity.tsv (confidence classes).")
    parser.add_argument(
        "--annotation-figures-dir",
        help="figures_png directory of a thesis-figures bundle; annotation/downstream PNGs are embedded in an Annotation tab.",
    )
    parser.add_argument(
        "--upset-figures-dir",
        help="Directory with upset_caller_support_*.png; the overall UpSet plot is embedded in the Consensus tab.",
    )
    parser.add_argument(
        "--search-root",
        action="append",
        default=[],
        help="Extra directory to auto-discover input tables. Repeat as needed.",
    )
    parser.add_argument("--title", default="DMR validation report")
    parser.add_argument(
        "--width", type=int, default=820, help="Default plot width in pixels."
    )
    parser.add_argument(
        "--height", type=int, default=420, help="Default plot height in pixels."
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def resolve_input(explicit: str | None, key: str) -> Path | None:
    if explicit:
        path = Path(explicit)
        return path if path.exists() else None
    return find_preferred_file(INPUT_NAMES[key])


def safe_read(path: Path | None) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    try:
        return read_table(path)
    except Exception:
        return pd.DataFrame()


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


# ---------------------------------------------------------------------------
# Plot builders. Each returns a HoloViews element (or None when not possible).
# ---------------------------------------------------------------------------


def _opts(element, *, width: int, height: int, **kwargs):
    """Apply only backend-agnostic options so plotly rendering never errors."""
    safe = {k: v for k, v in kwargs.items() if v is not None}
    return element.opts(width=width, height=height, **safe)


def plot_funnel(stages: list[dict], *, width: int, height: int):
    import holoviews as hv

    if not stages:
        return None
    df = pd.DataFrame(stages)
    df["stage"] = [f"{i + 1}. {label}" for i, label in enumerate(df["label"])]
    bars = hv.Bars(df, kdims=["stage"], vdims=["n", "note", "nested"])
    return _opts(
        bars,
        width=width,
        height=height,
        title="Candidate selection funnel",
        xlabel="stage",
        ylabel="regions",
    ).opts(invert_axes=True)


def plot_caller_counts(summary: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    if summary.empty or not {"caller", "context", "n_dmrs"}.issubset(summary.columns):
        return None
    df = summary[["caller", "context", "n_dmrs"]].copy()
    df["n_dmrs"] = _num(df["n_dmrs"]).fillna(0)
    bars = hv.Bars(df, kdims=["caller", "context"], vdims=["n_dmrs"])
    return _opts(
        bars,
        width=width,
        height=height,
        title="Reported DMRs by caller and context",
        ylabel="DMR count",
        xlabel="caller",
    )


def plot_volcano(rows: pd.DataFrame, *, width: int, height: int, q_threshold: float = 0.1):
    import holoviews as hv
    import numpy as np

    from dmr_validation_framework.reports.palette import caller_cmap

    if rows.empty or "delta" not in rows.columns or "q_value" not in rows.columns:
        return None
    df = rows.copy()
    df["_x"] = _num(df["delta"])
    df["_y"] = -np.log10(_num(df["q_value"]).clip(lower=1e-300))
    df = df.dropna(subset=["_x", "_y"])
    if df.empty:
        return None
    ctx_col = "_context" if "_context" in df.columns else ("context" if "context" in df.columns else None)
    caller_col = "_caller" if "_caller" in df.columns else ("caller" if "caller" in df.columns else None)
    if caller_col is None:
        return None
    cmap = caller_cmap(df[caller_col].astype(str).tolist())
    sig_line = float(-np.log10(max(q_threshold, 1e-300)))
    sub_w = max(280, width // 3)

    def panel(sub: pd.DataFrame, label: str):
        points = hv.Points(sub, kdims=["_x", "_y"], vdims=[caller_col]).opts(
            color=caller_col, cmap=cmap, size=5, width=sub_w, height=height,
            xlabel="delta methylation", ylabel="-log10 q", title=label, show_legend=True,
        )
        return points * hv.HLine(sig_line).opts(color="gray", line_dash="dot")

    if ctx_col is None:
        return panel(df, "all contexts")
    panels = [panel(sub, str(ctx)) for ctx, sub in df.groupby(ctx_col)]
    if not panels:
        return None
    layout = panels[0]
    for item in panels[1:]:
        layout = layout + item
    return layout.cols(3)


def plot_abs_delta(rows: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    if rows.empty or "_abs_delta" not in rows.columns:
        return None
    df = rows[["_caller", "_context", "_abs_delta"]].copy()
    df["_abs_delta"] = _num(df["_abs_delta"])
    df = df.dropna(subset=["_abs_delta"])
    if df.empty:
        return None
    box = hv.BoxWhisker(df, kdims=["_context", "_caller"], vdims=["_abs_delta"])
    return _opts(
        box,
        width=width,
        height=height,
        title="Absolute effect-size |delta| by context and caller",
        ylabel="|delta|",
        xlabel="context / caller",
    )


def detect_fixed_window_callers(rows: pd.DataFrame) -> dict[str, float]:
    """Callers whose DMR length is effectively constant (e.g. methylKit tiles).

    These are window-defined, not data-defined, so their length distribution is
    not comparable with the other callers and must not be shown as a normal box.
    """
    fixed: dict[str, float] = {}
    if rows.empty or "_length" not in rows.columns or "_caller" not in rows.columns:
        return fixed
    length = _num(rows["_length"])
    for caller, sub in rows.assign(_length=length).groupby("_caller"):
        values = sub["_length"].dropna()
        if len(values) >= 3 and float(values.std(ddof=0)) <= 1.0:
            fixed[str(caller)] = float(values.median())
    return fixed


def plot_length(rows: pd.DataFrame, fixed: dict[str, float], *, width: int, height: int):
    import holoviews as hv

    if rows.empty or "_length" not in rows.columns:
        return None
    df = rows[["_caller", "_context", "_length"]].copy()
    df["_length"] = _num(df["_length"])
    df = df.dropna(subset=["_length"])
    df = df[~df["_caller"].astype(str).isin(fixed.keys())]
    if df.empty:
        return None
    box = hv.BoxWhisker(df, kdims=["_context", "_caller"], vdims=["_length"])
    return _opts(
        box,
        width=width,
        height=height,
        title="DMR length distribution (data-defined callers only)",
        ylabel="length (bp)",
        xlabel="context / caller",
    )


def plot_concordance(membership: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    long_df = compute_caller_concordance(membership)
    if long_df.empty:
        return None
    heatmap = hv.HeatMap(long_df, kdims=["caller_b", "caller_a"], vdims=["jaccard", "n_intersect"])
    return _opts(
        heatmap,
        width=width,
        height=height,
        title="Caller concordance (Jaccard overlap of supported regions)",
        xlabel="caller",
        ylabel="caller",
        colorbar=True,
        cmap="viridis",
    )


def plot_caller_support(consensus: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    if consensus.empty or "supporting_callers" not in consensus.columns:
        return None
    combos = (
        consensus["supporting_callers"].astype(str).str.strip().replace("", np.nan).dropna()
    )
    if combos.empty:
        return None
    counts = combos.value_counts().reset_index()
    counts.columns = ["caller_combination", "n_regions"]
    counts = counts.head(25)
    bars = hv.Bars(counts, kdims=["caller_combination"], vdims=["n_regions"])
    return _opts(
        bars,
        width=width,
        height=height,
        title="Consensus regions by supporting-caller combination",
        ylabel="consensus regions",
        xlabel="supporting callers",
    )


def plot_tiers(consensus: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    if consensus.empty or "final_dmr_tier" not in consensus.columns:
        return None
    df = consensus.copy()
    df["context"] = df["context"].astype(str) if "context" in df.columns else "NA"
    table = (
        df.groupby(["context", "final_dmr_tier"]).size().reset_index(name="n_regions")
    )
    bars = hv.Bars(table, kdims=["context", "final_dmr_tier"], vdims=["n_regions"])
    return _opts(
        bars,
        width=width,
        height=height,
        title="Final DMR tiers by context",
        ylabel="consensus regions",
        xlabel="context / tier",
    )


def plot_overlap_sensitivity(table: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    if table.empty:
        return None
    x_col = next((c for c in table.columns if "threshold" in str(c).lower()), None)
    y_candidates = [
        c
        for c in table.columns
        if c != x_col
        and any(token in str(c).lower() for token in ("n_", "count", "regions", "consensus", "matched"))
    ]
    if x_col is None or not y_candidates:
        return None
    df = table[[x_col, *y_candidates]].copy()
    df[x_col] = _num(df[x_col])
    df = df.dropna(subset=[x_col]).sort_values(x_col)
    curves = []
    for y_col in y_candidates:
        sub = df[[x_col, y_col]].copy()
        sub[y_col] = _num(sub[y_col])
        sub = sub.dropna()
        if sub.empty:
            continue
        sub = sub.rename(columns={x_col: "threshold", y_col: "value"})
        sub["metric"] = str(y_col)
        curves.append(hv.Curve(sub, kdims=["threshold"], vdims=["value", "metric"], label=str(y_col)))
    if not curves:
        return None
    overlay = curves[0]
    for curve in curves[1:]:
        overlay = overlay * curve
    return _opts(
        overlay,
        width=width,
        height=height,
        title="Reciprocal-overlap threshold sensitivity",
        ylabel="count",
        xlabel="reciprocal-overlap threshold",
    )


def plot_glm_vs_glmm(comparison: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    if comparison.empty or not {"glm_q_value", "glmm_q_value"}.issubset(comparison.columns):
        return None
    df = comparison.copy()
    df["x"] = -np.log10(_num(df["glm_q_value"]).clip(lower=1e-300))
    df["y"] = -np.log10(_num(df["glmm_q_value"]).clip(lower=1e-300))
    df = df.dropna(subset=["x", "y"])
    if df.empty:
        return None

    def categorize(row: pd.Series) -> str:
        if bool(row.get("confirmed_by_glmm", False)):
            return "GLMM_CONFIRMED"
        if bool(row.get("glm_only_candidate", False)):
            return "GLM_ONLY"
        return "OTHER"

    df["category"] = df.apply(categorize, axis=1)
    vdims = ["category"]
    for extra in ("region_id", "context", "glm_q_value", "glmm_q_value"):
        if extra in df.columns:
            vdims.append(extra)
    points = hv.Points(df, kdims=["x", "y"], vdims=vdims)
    points = _opts(
        points,
        width=width,
        height=height,
        title="Aggregated GLM vs CpG-level GLMM (hover = region_id)",
        xlabel="-log10 GLM q",
        ylabel="-log10 GLMM q",
        color="category",
        cmap="Category10",
    )
    lim = float(max(df["x"].max(), df["y"].max()))
    diagonal = hv.Curve([(0, 0), (lim, lim)]).opts(color="gray", dash="dot")
    return points * diagonal


def plot_brown_qq(combp_pvalues: list[float], *, width: int, height: int):
    import holoviews as hv

    p = np.array([v for v in combp_pvalues if np.isfinite(v) and 0 < v <= 1], dtype=float)
    if p.size < 5:
        return None
    p.sort()
    n = p.size
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    df = pd.DataFrame({"expected": expected, "observed": observed})
    points = hv.Scatter(df, kdims=["expected"], vdims=["observed"])
    lim = float(max(expected.max(), observed.max()))
    diagonal = hv.Curve([(0, 0), (lim, lim)]).opts(color="gray", dash="dot")
    overlay = points * diagonal
    return _opts(
        overlay,
        width=width,
        height=height,
        title="comb-p regional p-value calibration (Brown-corrected Fisher)",
        xlabel="expected -log10 p (uniform)",
        ylabel="observed -log10 p",
    )


def plot_bootstrap_caterpillar(bootstrap: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    needed = {"delta_rep", "ci_low", "ci_high"}
    if bootstrap.empty or not needed.issubset(bootstrap.columns):
        return None
    df = bootstrap.copy()
    for col in ("delta_rep", "ci_low", "ci_high"):
        df[col] = _num(df[col])
    df = df.dropna(subset=["delta_rep", "ci_low", "ci_high"]).reset_index(drop=True)
    if df.empty:
        return None
    df = df.sort_values("delta_rep").reset_index(drop=True)
    df["rank"] = np.arange(len(df))
    if "ci_includes_zero" in df.columns:
        df["stable"] = ~df["ci_includes_zero"].map(
            lambda v: str(v).strip().lower() in {"true", "1", "yes"}
        )
    else:
        df["stable"] = (df["ci_low"] > 0) | (df["ci_high"] < 0)
    segs = hv.Segments(df, kdims=["rank", "ci_low", "rank", "ci_high"]).opts(color="lightgray")
    pts = hv.Points(df, kdims=["rank", "delta_rep"], vdims=["stable"]).opts(
        color="stable", cmap={True: "#1f77b4", False: "#d62728"}
    )
    zero = hv.Curve([(0, 0), (len(df), 0)]).opts(color="black", dash="dot")
    overlay = segs * zero * pts
    return _opts(
        overlay,
        width=width,
        height=height,
        title="Replicate-delta bootstrap CIs (red = CI includes 0)",
        xlabel="region (sorted by delta)",
        ylabel="delta_rep (b - a)",
    )


def plot_manhattan(rows: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    from dmr_validation_framework.reports.genome import compute_manhattan_layout
    from dmr_validation_framework.reports.palette import CONTEXT_COLORS

    layout = compute_manhattan_layout(rows)
    if not layout:
        return None
    df = layout["data"]
    ctx_col = "_context" if "_context" in df.columns else ("context" if "context" in df.columns else None)
    points = hv.Points(df, kdims=["_gpos", "_y"], vdims=[ctx_col] if ctx_col else [])
    opts: dict = dict(
        width=max(width, 1000),
        height=height,
        title="Genome-wide DMR distribution",
        xlabel="chromosome",
        ylabel=layout["y_label"],
        xticks=[(pos, label) for pos, label in layout["ticks"]],
    )
    if ctx_col:
        opts["color"] = ctx_col
        opts["cmap"] = dict(CONTEXT_COLORS)
    return _opts(points, **opts)


def plot_confidence(confidence: pd.DataFrame, *, width: int, height: int):
    import holoviews as hv

    if confidence.empty or "final_confidence_class" not in confidence.columns:
        return None
    counts = (
        confidence["final_confidence_class"].astype(str).value_counts().reset_index()
    )
    counts.columns = ["final_confidence_class", "n_regions"]
    bars = hv.Bars(counts, kdims=["final_confidence_class"], vdims=["n_regions"])
    return _opts(
        bars,
        width=width,
        height=height,
        title="Final confidence-class distribution",
        ylabel="regions",
        xlabel="confidence class",
    )


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def collect_annotation_images(figures_dir: Path | None) -> list[Path]:
    if figures_dir is None or not Path(figures_dir).exists():
        return []
    images = [
        path
        for path in sorted(Path(figures_dir).glob("*.png"))
        if any(token in path.name.lower() for token in ANNOTATION_TOKENS)
    ]
    return images


def _figure_titles(figures_dir: Path | None) -> dict[str, str]:
    """Map figure_id -> human title from the thesis-figure manifest, if present."""
    if figures_dir is None:
        return {}
    manifest = Path(figures_dir).parent / "tables" / "thesis_figure_manifest.tsv"
    df = safe_read(manifest)
    if df.empty or "figure_id" not in df.columns or "title" not in df.columns:
        return {}
    return {str(row["figure_id"]): str(row["title"]) for _, row in df.iterrows()}


def build_annotation_tab(panel_module, figures_dir: Path | None):
    import panel as pn

    images = collect_annotation_images(figures_dir)
    if not images:
        return pn.Column(
            pn.pane.Markdown("## Annotation"),
            pn.pane.HTML(
                "<p style='color:#666;font-size:13px'>No annotation figures available. "
                "Pass <code>--annotation-dir</code> (and optionally <code>--expression-table</code>) "
                "to the pipeline to plug in a gene/TE/GO annotation and generate these.</p>"
            ),
            sizing_mode="stretch_width",
        )
    titles = _figure_titles(figures_dir)
    blocks = [
        pn.pane.Markdown("## Annotation"),
        pn.pane.HTML(
            "<p style='color:#666;font-size:13px'>Gene / TE / GO / metagene-occupancy figures "
            "built from the supplied annotation. Rendered from the thesis-figure bundle.</p>"
        ),
    ]
    for path in images:
        caption = titles.get(path.stem, path.stem.replace("_", " "))
        blocks.append(pn.pane.Markdown(f"**{caption}**"))
        blocks.append(pn.pane.PNG(str(path), sizing_mode="scale_width", max_width=900))
    return pn.Column(*blocks, sizing_mode="stretch_width")


def _discover_combp_pvalues() -> list[float]:
    values: list[float] = []
    for context in ("CG", "CHG", "CHH"):
        path = find_preferred_file([f"combp_dmrs_{context}.tsv"])
        df = safe_read(path)
        if df.empty:
            continue
        col = next((c for c in df.columns if str(c).lower() in {"p_value", "pvalue", "p"}), None)
        if col:
            values.extend(_num(df[col]).dropna().tolist())
    return values


def _section(panel_module, title: str, items: list):
    import panel as pn

    note_style = "<p style='color:#666;font-size:13px'>{}</p>"
    blocks = [pn.pane.Markdown(f"## {title}")]
    rendered = 0
    for label, element in items:
        if element is None:
            blocks.append(pn.pane.HTML(note_style.format(f"{label}: input not available.")))
            continue
        blocks.append(pn.pane.HoloViews(element, backend="plotly"))
        rendered += 1
    return pn.Column(*blocks, sizing_mode="stretch_width"), rendered


def build_report(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.search_root:
        existing = os.environ.get(SEARCH_ROOTS_ENV, "")
        roots = [r for r in args.search_root] + ([existing] if existing else [])
        os.environ[SEARCH_ROOTS_ENV] = os.pathsep.join(roots)

    try:
        import holoviews as hv
        import panel as pn
    except ModuleNotFoundError as exc:
        print(
            f"interactive report requires holoviews and panel ({exc.name} missing); "
            "install them or use the static reports instead.",
            file=sys.stderr,
        )
        return 2

    hv.extension("plotly")
    pn.extension("plotly")
    from dmr_validation_framework.reports.theme import apply_holoviews_theme

    apply_holoviews_theme()

    width, height = args.width, args.height

    caller_rows = safe_read(resolve_input(args.caller_rows, "caller_rows"))
    caller_summary = safe_read(resolve_input(args.caller_summary, "caller_summary"))
    consensus = safe_read(resolve_input(args.consensus, "consensus"))
    glm_glmm = safe_read(resolve_input(args.glm_glmm, "glm_glmm"))
    bootstrap = safe_read(resolve_input(args.bootstrap, "bootstrap"))
    overlap = safe_read(resolve_input(args.overlap_sensitivity, "overlap_sensitivity"))
    confidence = safe_read(resolve_input(args.confidence, "confidence"))
    membership = safe_read(resolve_input(None, "membership"))
    fixed = detect_fixed_window_callers(caller_rows)

    def guarded(builder, *builder_args, **builder_kwargs):
        try:
            return builder(*builder_args, **builder_kwargs)
        except Exception as exc:  # noqa: BLE001
            print(f"plot '{builder.__name__}' skipped: {exc}", file=sys.stderr)
            return None

    funnel_stages = compute_funnel_stages(
        consensus=consensus if not consensus.empty else None,
        glm_glmm=glm_glmm if not glm_glmm.empty else None,
        caller_rows=caller_rows if not caller_rows.empty else None,
    )
    overview_tab, _ = _section(
        pn,
        "Overview",
        [("Candidate selection funnel", guarded(plot_funnel, funnel_stages, width=width, height=height))],
    )
    if funnel_stages:
        notes = "; ".join(f"{s['label']}={s['n']} ({s['note']})" for s in funnel_stages)
        overview_tab.append(
            pn.pane.HTML(f"<p style='color:#666;font-size:12px'>{notes}</p>")
        )

    callers_tab, n_callers = _section(
        pn,
        "Callers",
        [
            ("DMR counts", guarded(plot_caller_counts, caller_summary, width=width, height=height)),
            ("Volcano (delta vs -log10 q)", guarded(plot_volcano, caller_rows, width=width, height=height)),
            ("|delta| distribution", guarded(plot_abs_delta, caller_rows, width=width, height=height)),
            ("DMR length", guarded(plot_length, caller_rows, fixed, width=width, height=height)),
        ],
    )
    if fixed:
        note = "; ".join(f"{caller} (~{int(length)} bp constant)" for caller, length in fixed.items())
        callers_tab.append(
            pn.pane.HTML(
                "<p style='color:#a33;font-size:13px'>Fixed-window callers excluded from the "
                f"length plot (window-defined, not data-defined): {note}.</p>"
            )
        )

    upset_dir = Path(args.upset_figures_dir) if args.upset_figures_dir else None
    upset_images = sorted((upset_dir.glob("upset_caller_support_*.png")) if upset_dir and upset_dir.exists() else [])
    # The UpSet plot is the preferred view of caller-support intersections; the
    # bar of caller combinations is only a fallback when UpSet is unavailable.
    consensus_items = [
        ("Caller concordance", guarded(plot_concordance, membership, width=width, height=height)),
        ("DMR tiers", guarded(plot_tiers, consensus, width=width, height=height)),
        ("Overlap sensitivity", guarded(plot_overlap_sensitivity, overlap, width=width, height=height)),
    ]
    if not upset_images:
        consensus_items.insert(
            0, ("Caller support", guarded(plot_caller_support, consensus, width=width, height=height))
        )
    consensus_tab, _ = _section(pn, "Consensus", consensus_items)
    if upset_images:
        upset_images.sort(key=lambda p: (0 if "overall" in p.name else 1, p.name))
        consensus_tab.insert(1, pn.pane.Markdown("### Caller-support intersections (UpSet)"))
        for offset, path in enumerate(upset_images, start=2):
            consensus_tab.insert(offset, pn.pane.PNG(str(path), sizing_mode="scale_width", max_width=900))

    confirmatory_tab, _ = _section(
        pn,
        "Confirmatory (GLM/GLMM)",
        [
            ("GLM vs GLMM", guarded(plot_glm_vs_glmm, glm_glmm, width=width, height=height)),
            ("comb-p Brown QQ", guarded(plot_brown_qq, _discover_combp_pvalues(), width=width, height=height)),
        ],
    )

    robustness_tab, _ = _section(
        pn,
        "Robustness",
        [("Bootstrap CIs", guarded(plot_bootstrap_caterpillar, bootstrap, width=width, height=height))],
    )

    validation_tab, _ = _section(
        pn,
        "Validation status",
        [("Confidence classes", guarded(plot_confidence, confidence, width=width, height=height))],
    )

    genome_tab, _ = _section(
        pn,
        "Genome",
        [("Genome-wide DMR distribution", guarded(plot_manhattan, caller_rows, width=width, height=height))],
    )

    annotation_dir = Path(args.annotation_figures_dir) if args.annotation_figures_dir else None
    annotation_tab = build_annotation_tab(pn, annotation_dir)

    header = pn.pane.Markdown(
        f"# {args.title}\n"
        "Interactive summary of the DMR validation framework. This is a presentation "
        "layer over caller-native outputs and robustness diagnostics; it is not a "
        "genome-wide DMR caller."
    )
    tabs = pn.Tabs(
        ("Overview", overview_tab),
        ("Callers", callers_tab),
        ("Consensus", consensus_tab),
        ("Confirmatory", confirmatory_tab),
        ("Robustness", robustness_tab),
        ("Validation", validation_tab),
        ("Genome", genome_tab),
        ("Annotation", annotation_tab),
    )
    report = pn.Column(header, tabs, sizing_mode="stretch_width")

    out_path = out_dir / "index.html"
    try:
        report.save(str(out_path), embed=True, resources="inline")
    except TypeError:
        report.save(str(out_path))
    print(f"wrote interactive report: {out_path}")
    return 0


def run(args: argparse.Namespace) -> int:
    return build_report(args)


def main(argv: list[str] | None = None) -> int:
    return run(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
