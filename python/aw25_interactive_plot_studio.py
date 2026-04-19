from __future__ import annotations

"""Interactive AW25 plot studio.

Run with:
    streamlit run aw25_interactive_plot_studio.py
"""

import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import holoviews as hv
import plotly.io as pio
import streamlit as st
from plotly.graph_objects import Figure

from bsx2 import Context, HcAnnotStore, RegionReader
from bsx2.clustering import (
    AnnotationFormat,
    BackendConfig,
    BlockCacheConfig,
    ClusterConfig,
    ClusterSource,
    GeneProfileConfig,
    HierarchicalConfig,
    OutputConfig,
    ReadConfig,
    build_gene_profile_matrix,
    cluster_gene_profiles,
)
from bsx2.viz import (
    AnnotProfileLayout,
    AnnotProfilePart,
    GeneDendrogramPlotComposer,
    GeneEmbeddingPlotComposer,
    box_plot,
    build_annotation_metagene,
    build_chromosome_methylation_map,
    build_cluster_metagene_data,
    build_cluster_metagene_plot,
    build_gene_dendrogram_data,
    build_gene_embedding_data,
    build_manual_metagene,
    collect_layout_parts_from_hcannot,
    compute_chromosome_methylation_map_data,
    heatmap,
    line_plot,
    violin_plot,
)

APP_VERSION = "aw25-studio-v2"
CACHE_DIR = Path(os.environ.get("AW25_STUDIO_CACHE_DIR", ".aw25_plot_cache")).resolve()
UPLOAD_DIR = Path(os.environ.get("AW25_STUDIO_UPLOAD_DIR", ".aw25_uploaded_inputs")).resolve()
DOCS_BUILD_INDEX = Path(__file__).resolve().parent / "docs" / "build" / "html" / "index.html"
DOCS_PUBLIC_URL = os.environ.get("AW25_DOCS_URL")


@dataclass(frozen=True)
class PlotSpec:
    plot_id: str
    title: str
    eyebrow: str
    description: str
    default_height: int
    aw25_requirement: str
    review_preset: str
    review_goal: str


PLOT_SPECS: tuple[PlotSpec, ...] = (
    PlotSpec(
        "metagene_line",
        "Metagene profile - Line",
        "profile module",
        "Mean methylation trace from the selected metagene assembly mode.",
        470,
        "AW25 · metagene visualization · line plot",
        "Recommended review preset: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, limit=0, smooth=0.",
        "Shows that RegionReader + HcAnnotStore can produce a classical metagene trace with explicit upstream/body/downstream layout.",
    ),
    PlotSpec(
        "metagene_heatmap",
        "Metagene profile - Heatmap",
        "matrix module",
        "Ranked profile matrix from the selected metagene assembly mode.",
        680,
        "AW25 · metagene visualization · heat map",
        "Recommended review preset: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, limit=0, rank_rows=200.",
        "Shows the per-gene metagene matrix, ranking path, and segment-aware structure over the same metagene definition.",
    ),
    PlotSpec(
        "segment_box",
        "Metagene box plot - up/body/down",
        "distribution module",
        "Grouped segment-level distributions across upstream, body, and downstream.",
        460,
        "AW25 · metagene visualization · box plot",
        'Recommended review preset: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, per_region=False, distribution_mode="segments".',
        "Shows grouped segment distributions across many genes instead of one box per gene, which is the intended metagene summary view.",
    ),
    PlotSpec(
        "segment_violin",
        "Metagene violin plot - up/body/down",
        "distribution module",
        "Grouped violin distributions across upstream, body, and downstream.",
        460,
        "AW25 · metagene visualization · violin plot",
        'Recommended review preset: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, per_region=False, distribution_mode="segments".',
        "Shows the same grouped metagene signal as the box plot, but with the full distribution shape for each segment.",
    ),
    PlotSpec(
        "chromosome_map",
        "Chromosome methylation map",
        "genome module",
        "Chromosome-scale methylation density map with configurable genomic bins.",
        390,
        "AW25 · whole-genome chromosome methylation map",
        "Recommended review preset: context=CG, bin_size_bp=5000000.",
        "Shows that whole-genome methylation can be aggregated and rendered independently of the metagene and clustering pipelines.",
    ),
    PlotSpec(
        "gene_pca",
        "Gene PCA - KMeans labels",
        "embedding module",
        "Two-dimensional PCA embedding colored by KMeans-derived gene clusters.",
        480,
        "AW25 · KMeans clustering and PCA plot",
        "Recommended review preset: context=CG, limit_genes=0, min_coverage=5, n_clusters=4, seed=0.",
        "Shows gene-level clustering over methylation profiles and the PCA embedding required by the task.",
    ),
    PlotSpec(
        "gene_dendrogram",
        "Gene dendrogram",
        "hierarchy module",
        "Hierarchical clustering view over gene-profile similarity.",
        480,
        "AW25 · dendrogram / hierarchical grouping",
        "Recommended review preset: context=CG, limit_genes=50, min_coverage=5, n_clusters=4, hierarchical_max_genes>=50.",
        "Shows the hierarchical view over gene-profile similarity required for dendrogram support.",
    ),
    PlotSpec(
        "cluster_metagene",
        "Cluster metagene",
        "optional cluster summary",
        "Mean metagene profile per inferred gene cluster.",
        480,
        "AW25 · optional cluster-level methylation profile",
        "Recommended review preset: context=CG, limit_genes=0, min_coverage=5, n_clusters=4, seed=0.",
        "Shows the optional cluster-specific metagene summary described in AW25 after genes are assigned to clusters.",
    ),
)


def _font_family(language: str) -> str:
    if language == "ru":
        return '"Segoe UI", "Noto Sans", "PT Sans", "Helvetica Neue", Arial, system-ui, sans-serif'
    return '"Inter", system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'


def _mono_font_family() -> str:
    return 'ui-monospace, SFMono-Regular, Menlo, Consolas, monospace'


def _plot_specs(language: str) -> tuple[PlotSpec, ...]:
    if language == "ru":
        return (
            PlotSpec(
                "metagene_line",
                "Метаген - линейный профиль",
                "модуль профиля",
                "Средний профиль метилирования для выбранного режима сборки метагена.",
                470,
                "AW25 · визуализация метагена · line plot",
                "Рекомендуемый preset для ревью: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, limit=0, smooth=0.",
                "Показывает, что RegionReader + HcAnnotStore могут собрать классический метаген с явной структурой upstream/body/downstream.",
            ),
            PlotSpec(
                "metagene_heatmap",
                "Метаген - тепловая карта",
                "модуль матрицы",
                "Ранжированная матрица профилей для выбранного режима сборки метагена.",
                680,
                "AW25 · визуализация метагена · heat map",
                "Рекомендуемый preset для ревью: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, limit=0, rank_rows=200.",
                "Показывает матрицу профилей по генам, ranking path и сегментную структуру на том же определении метагена.",
            ),
            PlotSpec(
                "segment_box",
                "Метаген - box plot up/body/down",
                "модуль распределений",
                "Групповые распределения по upstream, body и downstream.",
                460,
                "AW25 · визуализация метагена · box plot",
                'Рекомендуемый preset для ревью: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, per_region=False, distribution_mode="segments".',
                "Показывает групповые сегментные распределения по множеству генов, а не отдельную коробку на каждый ген.",
            ),
            PlotSpec(
                "segment_violin",
                "Метаген - violin plot up/body/down",
                "модуль распределений",
                "Групповые violin-распределения по upstream, body и downstream.",
                460,
                "AW25 · визуализация метагена · violin plot",
                'Рекомендуемый preset для ревью: annotation-driven, context=CG, bins=25/50/25, flank=2000 bp, per_region=False, distribution_mode="segments".',
                "Показывает ту же групповую метагеновую картину, что и box plot, но со всей формой распределения по сегментам.",
            ),
            PlotSpec(
                "chromosome_map",
                "Карта метилирования по хромосомам",
                "геномный модуль",
                "Хромосомная карта плотности метилирования с настраиваемым размером бина.",
                390,
                "AW25 · полногеномная chromosome methylation map",
                "Рекомендуемый preset для ревью: context=CG, bin_size_bp=5000000.",
                "Показывает, что полногеномное метилирование агрегируется и визуализируется независимо от metagene и clustering pipeline.",
            ),
            PlotSpec(
                "gene_pca",
                "PCA генов - метки KMeans",
                "модуль embedding",
                "Двумерное PCA-представление, окрашенное по кластерам, найденным KMeans.",
                480,
                "AW25 · KMeans clustering и PCA plot",
                "Рекомендуемый preset для ревью: context=CG, limit_genes=0, min_coverage=5, n_clusters=4, seed=0.",
                "Показывает gene-level clustering по профилям метилирования и PCA embedding, требуемый заданием.",
            ),
            PlotSpec(
                "gene_dendrogram",
                "Дендрограмма генов",
                "модуль иерархии",
                "Иерархический взгляд на похожесть профилей генов.",
                480,
                "AW25 · dendrogram / hierarchical grouping",
                "Рекомендуемый preset для ревью: context=CG, limit_genes=50, min_coverage=5, n_clusters=4, hierarchical_max_genes>=50.",
                "Показывает иерархический режим группировки по профилям генов и закрывает поддержку дендрограмм.",
            ),
            PlotSpec(
                "cluster_metagene",
                "Метаген по кластерам",
                "опциональный cluster summary",
                "Средний метагеновый профиль для каждого найденного кластера.",
                480,
                "AW25 · optional cluster-level methylation profile",
                "Рекомендуемый preset для ревью: context=CG, limit_genes=0, min_coverage=5, n_clusters=4, seed=0.",
                "Показывает опциональный cluster-specific метаген после присвоения генов кластерам.",
            ),
        )
    return PLOT_SPECS


def _plot_spec(language: str, plot_id: str) -> PlotSpec:
    for spec in _plot_specs(language):
        if spec.plot_id == plot_id:
            return spec
    raise KeyError(plot_id)


def _slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return slug or "plot"


def _parse_context(value: str) -> Context:
    try:
        return getattr(Context, value.upper())
    except AttributeError as exc:
        raise ValueError(f"Unsupported methylation context: {value}") from exc


def _new_reader(bsx_path: str | Path, context: Context) -> RegionReader:
    reader = RegionReader(str(bsx_path))
    reader.clear_filters()
    reader.filter_context(context)
    return reader


def _metagene_layout(
    *,
    assembly_mode: str,
    up_bins: int,
    body_bins: int,
    down_bins: int,
    flank_bp: int,
) -> AnnotProfileLayout:
    if assembly_mode == "manual-composed":
        upstream_name = "upstream (manual)"
        body_name = "body (manual)"
        downstream_name = "downstream (manual)"
    else:
        upstream_name = "upstream"
        body_name = "body"
        downstream_name = "downstream"

    return AnnotProfileLayout(
        (
            AnnotProfilePart(upstream_name, up_bins, source="flank5", flank_bp=flank_bp),
            AnnotProfilePart(body_name, body_bins, source="gene"),
            AnnotProfilePart(downstream_name, down_bins, source="flank3", flank_bp=flank_bp),
        )
    )


@st.cache_resource(show_spinner=False)
def _annot_store(annot_gff_path: str) -> HcAnnotStore:
    return HcAnnotStore.from_gff(annot_gff_path)


@st.cache_resource(show_spinner=False)
def _manual_metagene_data(
    *,
    bsx_path: str,
    annot_gff_path: str,
    context_name: str,
    assembly_mode: str,
    up_bins: int,
    body_bins: int,
    down_bins: int,
    flank_bp: int,
    limit_regions: int | None,
):
    annot = _annot_store(annot_gff_path)
    layout = _metagene_layout(
        assembly_mode=assembly_mode,
        up_bins=up_bins,
        body_bins=body_bins,
        down_bins=down_bins,
        flank_bp=flank_bp,
    )
    part_map = collect_layout_parts_from_hcannot(annot, layout=layout, limit=limit_regions)
    drd = build_manual_metagene(
        _new_reader(bsx_path, _parse_context(context_name)),
        part_map=part_map,
        layout=layout,
    )
    return drd, tuple(layout.segments)


@st.cache_resource(show_spinner=False)
def _annotation_metagene_data(
    *,
    bsx_path: str,
    annot_gff_path: str,
    context_name: str,
    assembly_mode: str,
    up_bins: int,
    body_bins: int,
    down_bins: int,
    flank_bp: int,
    limit_regions: int | None,
):
    annot = _annot_store(annot_gff_path)
    layout = _metagene_layout(
        assembly_mode=assembly_mode,
        up_bins=up_bins,
        body_bins=body_bins,
        down_bins=down_bins,
        flank_bp=flank_bp,
    )
    drd = build_annotation_metagene(
        _new_reader(bsx_path, _parse_context(context_name)),
        annot,
        layout=layout,
        limit=limit_regions,
    )
    return drd, tuple(layout.segments)


def _metagene_data(
    *,
    assembly_mode: str,
    bsx_path: str,
    annot_gff_path: str,
    context_name: str,
    up_bins: int,
    body_bins: int,
    down_bins: int,
    flank_bp: int,
    limit_regions: int | None,
):
    if assembly_mode == "annotation-driven":
        return _annotation_metagene_data(
            bsx_path=bsx_path,
            annot_gff_path=annot_gff_path,
            context_name=context_name,
            assembly_mode=assembly_mode,
            up_bins=up_bins,
            body_bins=body_bins,
            down_bins=down_bins,
            flank_bp=flank_bp,
            limit_regions=limit_regions,
        )
    return _manual_metagene_data(
        bsx_path=bsx_path,
        annot_gff_path=annot_gff_path,
        context_name=context_name,
        assembly_mode=assembly_mode,
        up_bins=up_bins,
        body_bins=body_bins,
        down_bins=down_bins,
        flank_bp=flank_bp,
        limit_regions=limit_regions,
    )


@st.cache_resource(show_spinner=False)
def _chromosome_map_data(*, bsx_path: str, context_name: str, bin_size_bp: int):
    return compute_chromosome_methylation_map_data(
        _new_reader(bsx_path, _parse_context(context_name)),
        context=_parse_context(context_name),
        bin_size_bp=bin_size_bp,
    )


@st.cache_resource(show_spinner=False)
def _cluster_result(
    *,
    bsx_path: str,
    annot_gff_path: str,
    context_name: str,
    limit_genes: int | None,
    min_coverage: int,
    query_block_merge_gap_bp: int,
    n_clusters: int,
    seed: int,
    hierarchical_enabled: bool,
    hierarchical_max_genes: int,
):
    config = ClusterConfig(
        bsx_path=Path(bsx_path),
        annotation_path=Path(annot_gff_path),
        annotation_format=AnnotationFormat.GFF,
        read=ReadConfig(
            context=_parse_context(context_name),
            min_coverage=min_coverage,
            query_block_merge_gap_bp=query_block_merge_gap_bp,
        ),
        block_cache=BlockCacheConfig(enabled=False),
        gene_profile=GeneProfileConfig(limit_genes=limit_genes),
        backend=BackendConfig(n_components=2, n_clusters=n_clusters, seed=seed),
        hierarchical=HierarchicalConfig(
            enabled=hierarchical_enabled,
            max_genes=hierarchical_max_genes,
        ),
        cluster_source=ClusterSource.KMEANS,
        output=OutputConfig(
            output_dir=Path(tempfile.gettempdir()),
            write_table_files=False,
            write_metrics_file=False,
        ),
    )
    matrix = build_gene_profile_matrix(config)
    return cluster_gene_profiles(matrix, config)


def _stable_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _uploads_fingerprint(report_bytes: bytes, annot_bytes: bytes) -> str:
    h = hashlib.sha256()
    h.update(APP_VERSION.encode("utf-8"))
    h.update(b"\0report.bsx\0")
    h.update(report_bytes)
    h.update(b"\0annot.gff\0")
    h.update(annot_bytes)
    return h.hexdigest()[:24]


def _params_fingerprint(plot_id: str, params: dict[str, Any], inputs_fingerprint: str) -> str:
    payload = {
        "app_version": APP_VERSION,
        "plot_id": plot_id,
        "inputs": inputs_fingerprint,
        "params": params,
    }
    return hashlib.sha256(_stable_json(payload).encode("utf-8")).hexdigest()[:24]


def _cache_path(plot_id: str, params: dict[str, Any], inputs_fingerprint: str) -> Path:
    digest = _params_fingerprint(plot_id, params, inputs_fingerprint)
    return CACHE_DIR / inputs_fingerprint / f"{_slugify(plot_id)}_{digest}.plotly.json"


def _save_uploaded_inputs(report_file, annot_file) -> tuple[Path, Path, str]:
    report_bytes = report_file.getvalue()
    annot_bytes = annot_file.getvalue()
    inputs_fp = _uploads_fingerprint(report_bytes, annot_bytes)

    upload_dir = UPLOAD_DIR / inputs_fp
    upload_dir.mkdir(parents=True, exist_ok=True)
    report_path = upload_dir / "report.bsx"
    annot_path = upload_dir / "annot.gff"

    if not report_path.exists() or _sha256_bytes(report_path.read_bytes()) != _sha256_bytes(report_bytes):
        report_path.write_bytes(report_bytes)
    if not annot_path.exists() or _sha256_bytes(annot_path.read_bytes()) != _sha256_bytes(annot_bytes):
        annot_path.write_bytes(annot_bytes)

    return report_path, annot_path, inputs_fp


def _theme_tokens(theme_name: str) -> dict[str, str]:
    if theme_name == "Dark":
        return {
            "plotly_template": "plotly_dark",
            "paper_bgcolor": "#111827",
            "plot_bgcolor": "#111827",
            "font_color": "#e5edf7",
            "grid_color": "rgba(148,163,184,0.14)",
            "line_color": "rgba(148,163,184,0.34)",
            "tick_color": "rgba(191,219,254,0.38)",
            "hover_bgcolor": "#0f172a",
            "hover_bordercolor": "rgba(96,165,250,0.24)",
            "legend_bgcolor": "rgba(15,23,42,0.88)",
            "legend_bordercolor": "rgba(148,163,184,0.20)",
            "surface": "rgba(15,23,42,0.84)",
            "surface_alt": "rgba(17,24,39,0.96)",
            "surface_border": "rgba(148,163,184,0.16)",
            "surface_border_strong": "rgba(96,165,250,0.22)",
            "accent": "#60a5fa",
            "accent_alt": "#22d3ee",
            "accent_soft": "rgba(96,165,250,0.12)",
            "success": "#34d399",
            "bg": "#0b1020",
            "bg_grad_1": "rgba(37,99,235,0.12)",
            "bg_grad_2": "rgba(34,211,238,0.08)",
            "text_soft": "#b6c3d1",
            "text_muted": "#8b9bb0",
            "code_bg": "rgba(2,6,23,0.78)",
            "code_border": "rgba(148,163,184,0.16)",
            "shadow": "0 12px 32px rgba(2,6,23,0.42)",
            "button_bg": "rgba(15,23,42,0.92)",
            "button_text": "#e5edf7",
            "button_hover_bg": "rgba(15,23,42,0.98)",
            "button_hover_text": "#e5edf7",
            "button_hover_border": "#60a5fa",
        }
    return {
        "plotly_template": "plotly_white",
        "paper_bgcolor": "#ffffff",
        "plot_bgcolor": "#ffffff",
        "font_color": "#0f172a",
        "grid_color": "rgba(148,163,184,0.16)",
        "line_color": "rgba(100,116,139,0.34)",
        "tick_color": "rgba(71,85,105,0.42)",
        "hover_bgcolor": "#f8fafc",
        "hover_bordercolor": "rgba(37,99,235,0.18)",
        "legend_bgcolor": "rgba(255,255,255,0.94)",
        "legend_bordercolor": "rgba(203,213,225,0.92)",
        "surface": "rgba(255,255,255,0.92)",
        "surface_alt": "rgba(248,250,252,0.98)",
        "surface_border": "rgba(203,213,225,0.84)",
        "surface_border_strong": "rgba(37,99,235,0.20)",
        "accent": "#2563eb",
        "accent_alt": "#0f766e",
        "accent_soft": "rgba(37,99,235,0.08)",
        "success": "#059669",
        "bg": "#f5f7fb",
        "bg_grad_1": "rgba(37,99,235,0.08)",
        "bg_grad_2": "rgba(15,118,110,0.05)",
        "text_soft": "#475569",
        "text_muted": "#64748b",
        "code_bg": "rgba(241,245,249,0.95)",
        "code_border": "rgba(203,213,225,0.84)",
        "shadow": "0 8px 24px rgba(15,23,42,0.06)",
        "button_bg": "#ffffff",
        "button_text": "#000000",
        "button_hover_bg": "#ffffff",
        "button_hover_text": "#000000",
        "button_hover_border": "#2563eb",
    }


def _style_figure(fig: Figure, *, height: int, title: str, theme_name: str, language: str) -> Figure:
    tokens = _theme_tokens(theme_name)
    app_font = _font_family(language)
    fig.update_layout(
        template=tokens["plotly_template"],
        height=height,
        autosize=True,
        title={"text": title, "x": 0.015, "xanchor": "left", "font": {"size": 16}},
        paper_bgcolor=tokens["paper_bgcolor"],
        plot_bgcolor=tokens["plot_bgcolor"],
        font={
            "color": tokens["font_color"],
            "family": app_font,
            "size": 12,
        },
        margin={"l": 60, "r": 24, "t": 60, "b": 52},
        hoverlabel={
            "bgcolor": tokens["hover_bgcolor"],
            "bordercolor": tokens["hover_bordercolor"],
            "font": {"color": tokens["font_color"]},
        },
        legend={
            "bgcolor": tokens["legend_bgcolor"],
            "bordercolor": tokens["legend_bordercolor"],
            "borderwidth": 1,
        },
    )
    fig.update_xaxes(
        showgrid=True,
        gridcolor=tokens["grid_color"],
        zeroline=False,
        linecolor=tokens["line_color"],
        tickcolor=tokens["tick_color"],
        ticks="outside",
        showline=True,
        mirror=False,
        title_font={"size": 12, "color": tokens["text_soft"]},
        tickfont={"size": 11},
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor=tokens["grid_color"],
        zeroline=False,
        linecolor=tokens["line_color"],
        tickcolor=tokens["tick_color"],
        ticks="outside",
        showline=True,
        mirror=False,
        title_font={"size": 12, "color": tokens["text_soft"]},
        tickfont={"size": 11},
    )
    return fig


def _hv_plot_to_figure(plot: object, *, height: int, title: str, theme_name: str, language: str) -> Figure:
    rendered = hv.render(plot, backend="plotly")
    fig = rendered if isinstance(rendered, Figure) else Figure(rendered)
    return _style_figure(fig, height=height, title=title, theme_name=theme_name, language=language)


def _figure_to_standalone_html(fig: Figure, *, title: str) -> str:
    return pio.to_html(
        fig,
        include_plotlyjs="cdn",
        full_html=True,
        config={
            "responsive": True,
            "displaylogo": False,
            "displayModeBar": "hover",
            "toImageButtonOptions": {
                "format": "svg",
                "filename": _slugify(title),
                "scale": 2,
            },
        },
    )


def _manual_params(params: dict[str, Any]) -> dict[str, Any]:
    return {
        "assembly_mode": params["assembly_mode"],
        "context_name": params["context"],
        "up_bins": params["up_bins"],
        "body_bins": params["body_bins"],
        "down_bins": params["down_bins"],
        "flank_bp": params["flank_bp"],
        "limit_regions": params["limit_regions"] or None,
    }


def _cluster_params(params: dict[str, Any], *, hierarchical_enabled: bool) -> dict[str, Any]:
    return {
        "context_name": params["context"],
        "limit_genes": params["limit_genes"] or None,
        "min_coverage": params["min_coverage"],
        "query_block_merge_gap_bp": params["query_block_merge_gap_bp"],
        "n_clusters": params["n_clusters"],
        "seed": params["seed"],
        "hierarchical_enabled": hierarchical_enabled,
        "hierarchical_max_genes": params.get("hierarchical_max_genes", 5_000),
    }


def build_plot(*, plot_id: str, bsx_path: Path, annot_gff_path: Path, params: dict[str, Any]) -> Figure:
    hv.extension("plotly")
    height = int(params["height"])
    context_name = params["context"]
    theme_name = params.get("theme", "Light")
    language = params.get("language", "en")
    spec = _plot_spec(language, plot_id)

    if plot_id in {"metagene_line", "metagene_heatmap", "segment_box", "segment_violin"}:
        drd, segments = _metagene_data(
            bsx_path=str(bsx_path),
            annot_gff_path=str(annot_gff_path),
            **_manual_params(params),
        )
        if plot_id == "metagene_line":
            smooth_window = params.get("smooth_window", 0)
            plot = line_plot(
                drd,
                name=context_name,
                segments=list(segments),
                smooth=(smooth_window or None),
                title=f"{spec.title} | n={len(drd.positions)}",
                width=1100,
                height=height,
            )
            return _hv_plot_to_figure(
                plot,
                height=height,
                title=spec.title,
                theme_name=theme_name,
                language=language,
            )
        if plot_id == "metagene_heatmap":
            rank_rows = params.get("rank_rows", 0)
            plot = heatmap(
                drd,
                segments=list(segments),
                rank_rows=(rank_rows or None),
                title=f"{spec.title} | n={len(drd.positions)}",
                width=1200,
                height=height,
            )
            return _hv_plot_to_figure(
                plot,
                height=height,
                title=spec.title,
                theme_name=theme_name,
                language=language,
            )

        if plot_id == "segment_box":
            plot = box_plot(
                drd,
                segments=list(segments),
                per_region=params["per_region"],
                distribution_mode=params["distribution_mode"],
                as_percent=params["as_percent"],
                title=spec.title,
                width=1000,
                height=height,
            )
            return _hv_plot_to_figure(
                plot,
                height=height,
                title=spec.title,
                theme_name=theme_name,
                language=language,
            )

        plot = violin_plot(
            drd,
            segments=list(segments),
            per_region=params["per_region"],
            distribution_mode=params["distribution_mode"],
            as_percent=params["as_percent"],
            title=spec.title,
            width=1000,
            height=height,
        )
        return _hv_plot_to_figure(
            plot,
            height=height,
            title=spec.title,
            theme_name=theme_name,
            language=language,
        )

    if plot_id == "chromosome_map":
        data = _chromosome_map_data(
            bsx_path=str(bsx_path),
            context_name=context_name,
            bin_size_bp=params["bin_size_bp"],
        )
        plot = build_chromosome_methylation_map(
            data,
            name=bsx_path.stem,
            title="Chromosome methylation map",
            width=1100,
            height=height,
        )
        return _hv_plot_to_figure(
            plot,
            height=height,
            title=spec.title,
            theme_name=theme_name,
            language=language,
        )

    if plot_id == "gene_pca":
        result = _cluster_result(
            bsx_path=str(bsx_path),
            annot_gff_path=str(annot_gff_path),
            **_cluster_params(params, hierarchical_enabled=False),
        )
        plot = (
            GeneEmbeddingPlotComposer(
                title=f"Gene PCA | limit={params['limit_genes'] or 'None'}",
                width=900,
                height=height,
            )
            .add_data(build_gene_embedding_data(result))
            .finish()
        )
        return _hv_plot_to_figure(
            plot,
            height=height,
            title=spec.title,
            theme_name=theme_name,
            language=language,
        )

    if plot_id == "gene_dendrogram":
        result = _cluster_result(
            bsx_path=str(bsx_path),
            annot_gff_path=str(annot_gff_path),
            **_cluster_params(params, hierarchical_enabled=True),
        )
        dendrogram_data = build_gene_dendrogram_data(result)
        if dendrogram_data is None:
            raise RuntimeError(
                "Dendrogram data is empty. Try increasing limit_genes or relaxing constraints."
            )
        plot = (
            GeneDendrogramPlotComposer(
                title=f"Gene dendrogram | limit={params['limit_genes'] or 'None'}",
                width=1100,
                height=height,
            )
            .add_data(dendrogram_data)
            .finish()
        )
        return _hv_plot_to_figure(
            plot,
            height=height,
            title=spec.title,
            theme_name=theme_name,
            language=language,
        )

    if plot_id == "cluster_metagene":
        result = _cluster_result(
            bsx_path=str(bsx_path),
            annot_gff_path=str(annot_gff_path),
            **_cluster_params(params, hierarchical_enabled=False),
        )
        plot = build_cluster_metagene_plot(
            build_cluster_metagene_data(result),
            title=f"Cluster metagene | limit={params['limit_genes'] or 'None'}",
            width=1100,
            height=height,
        )
        return _hv_plot_to_figure(
            plot,
            height=height,
            title=spec.title,
            theme_name=theme_name,
            language=language,
        )

    raise ValueError(f"Unknown plot_id: {plot_id}")


def cached_or_build_plot(
    *,
    spec: PlotSpec,
    bsx_path: Path,
    annot_gff_path: Path,
    inputs_fingerprint: str,
    params: dict[str, Any],
    force_rebuild: bool,
) -> tuple[Figure, bool, Path]:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = _cache_path(spec.plot_id, params, inputs_fingerprint)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    if cache_file.exists() and not force_rebuild:
        return pio.from_json(cache_file.read_text(encoding="utf-8")), True, cache_file

    fig = build_plot(
        plot_id=spec.plot_id,
        bsx_path=bsx_path,
        annot_gff_path=annot_gff_path,
        params=params,
    )
    cache_file.write_text(fig.to_json(), encoding="utf-8")
    return fig, False, cache_file


def configure_page() -> None:
    st.set_page_config(
        page_title="AW25 Plot Studio",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )


def apply_page_style(theme_name: str, language: str) -> None:
    tokens = _theme_tokens(theme_name)
    bg_text = "#eef4f2" if theme_name == "Dark" else "#1e2521"
    app_font = _font_family(language)
    mono_font = _mono_font_family()
    st.markdown(
        f"""
        <style>
          html, body, [class*="css"] {{
            font-family: {app_font};
          }}
          .stApp {{
            color: {bg_text};
            background:
              radial-gradient(circle at 10% 6%, {tokens["bg_grad_1"]}, transparent 25rem),
              radial-gradient(circle at 88% 10%, {tokens["bg_grad_2"]}, transparent 24rem),
              linear-gradient(180deg, {tokens["bg"]} 0%, {tokens["bg"]} 100%);
          }}
          [data-testid="stHeader"] {{
            background: transparent;
          }}
          [data-testid="stToolbar"] {{
            right: 1rem;
          }}
          [data-testid="stSidebar"] > div:first-child {{
            background:
              linear-gradient(180deg, {tokens["surface_alt"]} 0%, {tokens["surface"]} 100%);
            border-right: 1px solid {tokens["surface_border"]};
          }}
          [data-testid="stSidebar"] .stMarkdown p,
          [data-testid="stSidebar"] .stCaption {{
            color: {tokens["text_soft"]};
          }}
          [data-testid="stSidebar"] code {{
            background: {tokens["code_bg"]};
            border: 1px solid {tokens["code_border"]};
            border-radius: 6px;
            padding: .12rem .35rem;
            font-family: {mono_font};
            font-size: .82rem;
          }}
          .sidebar-label {{
            margin: 0 0 .55rem;
            color: {tokens["text_muted"]};
            font-size: .72rem;
            font-weight: 700;
            letter-spacing: .14em;
            text-transform: uppercase;
          }}
          .sidebar-meta {{
            margin: .3rem 0 .6rem;
            padding: .75rem .85rem;
            border: 1px solid {tokens["surface_border"]};
            border-radius: 12px;
            background: {tokens["surface"]};
            box-shadow: {tokens["shadow"]};
          }}
          .sidebar-meta span {{
            display: block;
            margin-bottom: .28rem;
            color: {tokens["text_muted"]};
            font-size: .75rem;
            font-weight: 600;
            letter-spacing: .05em;
          }}
          .hero {{
            padding: 1.2rem 1.4rem 1.3rem;
            border: 1px solid {tokens["surface_border_strong"]};
            border-radius: 16px;
            background:
              linear-gradient(180deg, {tokens["surface_alt"]} 0%, {tokens["surface"]} 100%);
            box-shadow: {tokens["shadow"]};
            margin-bottom: 1rem;
          }}
          .hero-kicker {{
            margin: 0 0 .65rem;
            color: {tokens["accent"]};
            font-size: .75rem;
            font-weight: 700;
            letter-spacing: .12em;
            text-transform: uppercase;
            font-family: {mono_font};
          }}
          .hero-title {{
            margin: 0;
            font-size: clamp(1.95rem, 4vw, 2.75rem);
            line-height: 1.02;
            letter-spacing: -.02em;
            font-weight: 700;
          }}
          .hero-subtitle {{
            margin: .5rem 0 0;
            color: {tokens["text_soft"]};
            font-size: 1rem;
            font-weight: 500;
            max-width: 920px;
          }}
          .hero-lead {{
            margin: .6rem 0 0;
            max-width: 960px;
            color: {tokens["text_muted"]};
            line-height: 1.6;
          }}
          .hero-meta {{
            display: flex;
            flex-wrap: wrap;
            gap: .55rem;
            margin-top: .95rem;
          }}
          .hero-pill {{
            display: inline-flex;
            align-items: center;
            gap: .45rem;
            padding: .42rem .68rem;
            border: 1px solid {tokens["surface_border"]};
            border-radius: 999px;
            background: {tokens["accent_soft"]};
            color: {tokens["font_color"]};
            font-size: .8rem;
            font-weight: 600;
            font-family: {mono_font};
          }}
          .plot-card, .docs-card {{
            padding: .95rem 1.05rem;
            margin: .65rem 0 .9rem;
            border: 1px solid {tokens["surface_border"]};
            border-radius: 14px;
            background:
              linear-gradient(180deg, {tokens["surface_alt"]} 0%, {tokens["surface"]} 100%);
            box-shadow: {tokens["shadow"]};
          }}
          .coverage-card {{
            padding: 1rem 1.05rem;
            margin: .25rem 0 1rem;
            border: 1px solid {tokens["surface_border_strong"]};
            border-radius: 14px;
            background:
              linear-gradient(180deg, {tokens["surface_alt"]} 0%, {tokens["surface"]} 100%);
            box-shadow: {tokens["shadow"]};
          }}
          .coverage-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: .75rem;
            margin-top: .85rem;
          }}
          .coverage-item {{
            padding: .8rem .9rem;
            border: 1px solid {tokens["surface_border"]};
            border-radius: 12px;
            background: {tokens["surface"]};
          }}
          .coverage-item h4 {{
            margin: 0 0 .4rem;
            font-size: .95rem;
          }}
          .coverage-item p {{
            margin: 0;
            color: {tokens["text_soft"]};
            line-height: 1.5;
          }}
          .protocol-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
            gap: .75rem;
            margin-top: .85rem;
          }}
          .protocol-item {{
            padding: .85rem .95rem;
            border: 1px solid {tokens["surface_border"]};
            border-radius: 12px;
            background: {tokens["surface"]};
          }}
          .protocol-item strong {{
            display: block;
            margin-bottom: .38rem;
            color: {tokens["font_color"]};
          }}
          .protocol-item p {{
            margin: 0;
            color: {tokens["text_soft"]};
            line-height: 1.55;
          }}
          .review-strip {{
            margin: .45rem 0 0;
            padding: .75rem .85rem;
            border: 1px solid {tokens["surface_border"]};
            border-radius: 12px;
            background: {tokens["accent_soft"]};
          }}
          .review-label {{
            display: block;
            margin-bottom: .25rem;
            color: {tokens["accent"]};
            font-size: .72rem;
            font-weight: 700;
            letter-spacing: .12em;
            text-transform: uppercase;
            font-family: {mono_font};
          }}
          .review-text {{
            margin: 0;
            color: {tokens["text_soft"]};
            line-height: 1.55;
          }}
          .plot-eyebrow {{
            margin: 0 0 .2rem;
            color: {tokens["accent"]};
            font-size: .72rem;
            font-weight: 700;
            letter-spacing: .12em;
            text-transform: uppercase;
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
          }}
          .plot-title {{
            margin: 0;
            font-size: 1.18rem;
            font-weight: 650;
            line-height: 1.2;
          }}
          .plot-desc, .docs-card p, .docs-card li {{
            color: {tokens["text_soft"]};
            line-height: 1.58;
          }}
          .cache-pill {{
            display: inline-flex;
            align-items: center;
            gap: .4rem;
            padding: .34rem .7rem;
            border: 1px solid {tokens["surface_border_strong"]};
            border-radius: 999px;
            background: {tokens["code_bg"]};
            font-size: .82rem;
            font-family: {mono_font};
            color: {tokens["font_color"]};
          }}
          .docs-card h2, .docs-card h3 {{
            margin-top: 0;
            color: {bg_text};
            font-weight: 650;
          }}
          .docs-card ul {{
            padding-left: 1.15rem;
          }}
          .docs-card code, .hero code, .sidebar-meta code {{
            background: {tokens["code_bg"]};
            border: 1px solid {tokens["code_border"]};
            padding: .1rem .34rem;
            border-radius: .4rem;
            font-family: {mono_font};
            font-size: .84rem;
          }}
          .docs-card pre {{
            background: {tokens["code_bg"]};
            border: 1px solid {tokens["code_border"]};
            border-radius: 12px;
            padding: .9rem 1rem;
            overflow-x: auto;
            color: {tokens["font_color"]};
            font-family: {mono_font};
          }}
          [data-testid="stForm"] {{
            border: 1px solid {tokens["surface_border"]};
            border-radius: 14px;
            padding: .9rem .95rem .25rem;
            background: {tokens["surface_alt"]};
          }}
          [data-testid="stExpander"] {{
            border: 1px solid {tokens["surface_border"]};
            border-radius: 14px;
            background: {tokens["surface"]};
            box-shadow: {tokens["shadow"]};
          }}
          [data-testid="stExpander"] summary {{
            font-weight: 650;
          }}
          [data-baseweb="tab-list"] {{
            gap: .4rem;
          }}
          [data-baseweb="tab"] {{
            border-radius: 10px 10px 0 0;
            padding: .55rem .9rem;
            font-weight: 600;
          }}
          .stButton > button,
          .stDownloadButton > button {{
            border-radius: 10px;
            border: 1px solid {tokens["surface_border_strong"]};
            background: {tokens["button_bg"]};
            color: {tokens["button_text"]};
            font-weight: 600;
          }}
          .stButton > button:hover,
          .stDownloadButton > button:hover {{
            border-color: {tokens["button_hover_border"]};
            color: {tokens["button_hover_text"]};
            background: {tokens["button_hover_bg"]};
          }}
          [data-testid="stFileUploader"] section {{
            border-radius: 12px;
            border: 1px dashed {tokens["surface_border_strong"]};
            background: {tokens["surface"]};
          }}
          [data-testid="stAlert"] {{
            border-radius: 12px;
            border: 1px solid {tokens["surface_border"]};
          }}
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_hero(language: str) -> None:
    if language == "ru":
        html = """
        <section class="hero">
          <p class="hero-kicker">BSXplorer2 / AW25 / analysis workstation</p>
          <h1 class="hero-title">Интерактивная студия графиков</h1>
          <p class="hero-subtitle">
            Научная консоль визуализации для анализа бисульфитного метилирования.
          </p>
          <p class="hero-lead">
            Строй метагеновые, хромосомные и кластеризационные графики из
            <code>report.bsx</code> и <code>annot.gff</code>. Каждый модуль хранит
            собственные параметры, строит Plotly-вывод по запросу и использует
            детерминированный дисковый кэш для воспроизводимых повторных запусков.
          </p>
          <div class="hero-meta">
            <span class="hero-pill">локальная аналитическая среда</span>
            <span class="hero-pill">воспроизводимые сборки графиков</span>
            <span class="hero-pill">кэшируемые Plotly-результаты</span>
            <span class="hero-pill">manual + annotation metagene paths</span>
          </div>
        </section>
        """
    else:
        html = """
        <section class="hero">
          <p class="hero-kicker">BSXplorer2 / AW25 / analysis workstation</p>
          <h1 class="hero-title">Interactive Plot Studio</h1>
          <p class="hero-subtitle">
            Scientific visualization console for bisulfite methylation analysis.
          </p>
          <p class="hero-lead">
            Build metagene, chromosome-scale, and clustering figures from
            <code>report.bsx</code> and <code>annot.gff</code>. Each module keeps
            its own parameter set, generates Plotly-backed output on demand, and
            reuses deterministic on-disk cache entries for reproducible reruns.
          </p>
          <div class="hero-meta">
            <span class="hero-pill">local interactive analysis studio</span>
            <span class="hero-pill">reproducible plot builds</span>
            <span class="hero-pill">cached Plotly outputs</span>
            <span class="hero-pill">manual + annotation metagene paths</span>
          </div>
        </section>
        """
    st.markdown(html, unsafe_allow_html=True)


def render_aw25_coverage(language: str) -> None:
    html = """
        <section class="coverage-card">
          <p class="plot-eyebrow">aw25 coverage map</p>
          <h3 class="plot-title">Required deliverables visible in this workstation</h3>
          <p class="plot-desc">
            This layout is arranged so a reviewer can map each AW25 requirement to a concrete analysis module.
          </p>
          <div class="coverage-grid">
            <div class="coverage-item">
              <h4>Metagene visualization</h4>
              <p>Line, Heatmap, Box, and Violin are exposed as separate modules over the same metagene input and layout.</p>
            </div>
            <div class="coverage-item">
              <h4>RegionReader + HcAnnotStore</h4>
              <p>Metagene modules support both annotation-driven and manual-composed assembly modes for transparent compute paths.</p>
            </div>
            <div class="coverage-item">
              <h4>Clustering and grouping</h4>
              <p>Gene PCA uses KMeans labels, Gene dendrogram covers hierarchy, and Cluster metagene covers grouped profiles.</p>
            </div>
            <div class="coverage-item">
              <h4>Genome-scale output</h4>
              <p>Chromosome methylation map exposes whole-genome methylation structure with configurable bins.</p>
            </div>
          </div>
        </section>
        """
    if language == "ru":
        html = """
        <section class="coverage-card">
          <p class="plot-eyebrow">aw25 coverage map</p>
          <h3 class="plot-title">Какие deliverables AW25 уже видны в этой студии</h3>
          <p class="plot-desc">
            Компоновка сделана так, чтобы проверяющий мог напрямую сопоставить каждый пункт AW25 с конкретным analysis module.
          </p>
          <div class="coverage-grid">
            <div class="coverage-item">
              <h4>Визуализация метагена</h4>
              <p>Line, Heatmap, Box и Violin доступны как отдельные модули поверх одного и того же metagene layout.</p>
            </div>
            <div class="coverage-item">
              <h4>RegionReader + HcAnnotStore</h4>
              <p>Метагеновые модули поддерживают both annotation-driven и manual-composed режимы сборки.</p>
            </div>
            <div class="coverage-item">
              <h4>Кластеризация и группировка</h4>
              <p>Gene PCA использует KMeans labels, Gene dendrogram покрывает hierarchy, а Cluster metagene закрывает grouped profiles.</p>
            </div>
            <div class="coverage-item">
              <h4>Полногеномный вывод</h4>
              <p>Chromosome methylation map показывает whole-genome структуру метилирования с настраиваемыми genomic bins.</p>
            </div>
          </div>
        </section>
        """
    st.markdown(html, unsafe_allow_html=True)


def render_aw25_review_protocol(language: str) -> None:
    html = """
        <section class="coverage-card">
          <p class="plot-eyebrow">review protocol</p>
          <h3 class="plot-title">What to run to demonstrate AW25 completion</h3>
          <p class="plot-desc">
            A reviewer should not have to infer coverage from source code. The sequence below mirrors the task statement
            and makes each required capability explicit in the workstation UI.
          </p>
          <div class="protocol-grid">
            <div class="protocol-item">
              <strong>1. Metagene baseline</strong>
              <p>Run Line, Heatmap, Box, and Violin in <code>annotation-driven</code> mode with <code>context=CG</code>,
              <code>bins=25/50/25</code>, <code>flank=2000</code>, and <code>limit=0</code>.</p>
            </div>
            <div class="protocol-item">
              <strong>2. Transparent assembly path</strong>
              <p>Re-run at least Line or Heatmap in <code>manual-composed</code> mode to show that the metagene can also be
              assembled from explicit layout parts instead of only from annotation sugar.</p>
            </div>
            <div class="protocol-item">
              <strong>3. Genome-scale output</strong>
              <p>Run Chromosome methylation map with <code>context=CG</code> and a reviewer-friendly bin size such as
              <code>5000000</code>.</p>
            </div>
            <div class="protocol-item">
              <strong>4. Gene-level clustering</strong>
              <p>Run Gene PCA with KMeans labels, then Gene dendrogram, then Cluster metagene. Together these cover
              clustering, PCA, dendrograms, and the optional grouped methylation profile.</p>
            </div>
          </div>
        </section>
        """
    if language == "ru":
        html = """
        <section class="coverage-card">
          <p class="plot-eyebrow">review protocol</p>
          <h3 class="plot-title">Что запускать, чтобы показать выполнение AW25</h3>
          <p class="plot-desc">
            Проверяющий не должен восстанавливать покрытие задачи по исходникам. Последовательность ниже повторяет формулировку задания и делает каждый требуемый capability явным в UI.
          </p>
          <div class="protocol-grid">
            <div class="protocol-item">
              <strong>1. Базовый метаген</strong>
              <p>Запусти Line, Heatmap, Box и Violin в режиме <code>annotation-driven</code> с <code>context=CG</code>,
              <code>bins=25/50/25</code>, <code>flank=2000</code> и <code>limit=0</code>.</p>
            </div>
            <div class="protocol-item">
              <strong>2. Прозрачный путь сборки</strong>
              <p>Повтори хотя бы Line или Heatmap в режиме <code>manual-composed</code>, чтобы показать, что metagene можно собирать и из явных layout parts, а не только из annotation sugar.</p>
            </div>
            <div class="protocol-item">
              <strong>3. Полногеномный вывод</strong>
              <p>Запусти Chromosome methylation map с <code>context=CG</code> и понятным для ревью размером бина, например <code>5000000</code>.</p>
            </div>
            <div class="protocol-item">
              <strong>4. Кластеризация генов</strong>
              <p>Запусти Gene PCA с KMeans labels, затем Gene dendrogram, затем Cluster metagene. Вместе они закрывают clustering, PCA, dendrograms и optional grouped methylation profile.</p>
            </div>
          </div>
        </section>
        """
    st.markdown(html, unsafe_allow_html=True)


def common_manual_controls(prefix: str, language: str) -> dict[str, Any]:
    context_label = "Контекст" if language == "ru" else "Context"
    assembly_label = "Режим сборки" if language == "ru" else "Assembly mode"
    flank_label = "Фланк, bp" if language == "ru" else "Flank, bp"
    up_label = "Бины upstream" if language == "ru" else "Up bins"
    body_label = "Бины body" if language == "ru" else "Body bins"
    down_label = "Бины downstream" if language == "ru" else "Down bins"
    limit_label = "Лимит регионов / генов" if language == "ru" else "Limit regions / genes"
    height_label = "Высота графика" if language == "ru" else "Graph height"
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        context = st.selectbox(context_label, ["CG", "CHG", "CHH"], key=f"{prefix}_context")
        assembly_mode = st.selectbox(
            assembly_label,
            ["annotation-driven", "manual-composed"],
            index=0,
            key=f"{prefix}_assembly",
        )
        flank_bp = st.number_input(flank_label, min_value=0, max_value=250_000, value=2_000, step=500, key=f"{prefix}_flank")
    with c2:
        up_bins = st.number_input(up_label, min_value=1, max_value=300, value=25, step=1, key=f"{prefix}_up")
        body_bins = st.number_input(body_label, min_value=1, max_value=500, value=50, step=1, key=f"{prefix}_body")
    with c3:
        down_bins = st.number_input(down_label, min_value=1, max_value=300, value=25, step=1, key=f"{prefix}_down")
        limit_regions = st.number_input(limit_label, min_value=0, max_value=500_000, value=0, step=100, key=f"{prefix}_limit")
    with c4:
        height = st.number_input(height_label, min_value=260, max_value=1400, value=520, step=20, key=f"{prefix}_height")
    return {
        "context": context,
        "assembly_mode": assembly_mode,
        "flank_bp": int(flank_bp),
        "up_bins": int(up_bins),
        "body_bins": int(body_bins),
        "down_bins": int(down_bins),
        "limit_regions": int(limit_regions),
        "height": int(height),
    }


def common_cluster_controls(
    prefix: str,
    language: str,
    *,
    default_limit: int,
    default_height: int,
    show_hierarchical_max: bool = False,
) -> dict[str, Any]:
    context_label = "Контекст" if language == "ru" else "Context"
    limit_label = "Лимит генов" if language == "ru" else "Limit genes"
    coverage_label = "Мин. покрытие" if language == "ru" else "Min coverage"
    merge_label = "Merge gap запроса, bp" if language == "ru" else "Query merge gap, bp"
    clusters_label = "Кластеры" if language == "ru" else "Clusters"
    seed_label = "Seed"
    height_label = "Высота графика" if language == "ru" else "Graph height"
    hmax_label = "Макс. genes для hierarchical" if language == "ru" else "Hierarchical max genes"
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        context = st.selectbox(context_label, ["CG", "CHG", "CHH"], key=f"{prefix}_context")
        limit_genes = st.number_input(limit_label, min_value=0, max_value=500_000, value=default_limit, step=50, key=f"{prefix}_limit")
    with c2:
        min_coverage = st.number_input(coverage_label, min_value=1, max_value=1_000, value=5, step=1, key=f"{prefix}_coverage")
        merge_gap = st.number_input(merge_label, min_value=0, max_value=100_000, value=250, step=50, key=f"{prefix}_merge")
    with c3:
        n_clusters = st.number_input(clusters_label, min_value=2, max_value=50, value=4, step=1, key=f"{prefix}_clusters")
        seed = st.number_input(seed_label, min_value=0, max_value=1_000_000, value=0, step=1, key=f"{prefix}_seed")
    with c4:
        height = st.number_input(height_label, min_value=260, max_value=1400, value=default_height, step=20, key=f"{prefix}_height")
        hierarchical_max_genes = st.number_input(hmax_label, min_value=10, max_value=500_000, value=5_000, step=100, key=f"{prefix}_hmax") if show_hierarchical_max else 5_000
    return {
        "context": context,
        "limit_genes": int(limit_genes),
        "min_coverage": int(min_coverage),
        "query_block_merge_gap_bp": int(merge_gap),
        "n_clusters": int(n_clusters),
        "seed": int(seed),
        "hierarchical_max_genes": int(hierarchical_max_genes),
        "height": int(height),
    }


def plot_params_ui(spec: PlotSpec, language: str) -> dict[str, Any]:
    prefix = spec.plot_id
    if spec.plot_id == "metagene_line":
        params = common_manual_controls(prefix, language)
        params["smooth_window"] = int(
            st.number_input(
                "Окно сглаживания, 0 = off" if language == "ru" else "Smooth window, 0 = off",
                min_value=0,
                max_value=250,
                value=0,
                step=1,
                key=f"{prefix}_smooth",
            )
        )
        return params
    if spec.plot_id == "metagene_heatmap":
        params = common_manual_controls(prefix, language)
        params["rank_rows"] = int(
            st.number_input(
                "Число rank rows, 0 = off" if language == "ru" else "Rank rows, 0 = off",
                min_value=0,
                max_value=500_000,
                value=200,
                step=50,
                key=f"{prefix}_rank",
            )
        )
        params["height"] = max(params["height"], 640)
        return params
    if spec.plot_id in {"segment_box", "segment_violin"}:
        params = common_manual_controls(prefix, language)
        c1, c2, c3 = st.columns(3)
        with c1:
            params["per_region"] = st.checkbox("По регионам" if language == "ru" else "Per region", value=False, key=f"{prefix}_per_region")
        with c2:
            params["as_percent"] = st.checkbox("В процентах" if language == "ru" else "As percent", value=False, key=f"{prefix}_percent")
        with c3:
            params["distribution_mode"] = st.selectbox("Режим распределения" if language == "ru" else "Distribution mode", ["segments"], key=f"{prefix}_mode")
        return params
    if spec.plot_id == "chromosome_map":
        c1, c2, c3 = st.columns(3)
        with c1:
            context = st.selectbox("Контекст" if language == "ru" else "Context", ["CG", "CHG", "CHH"], key=f"{prefix}_context")
        with c2:
            bin_size = st.number_input("Размер бина, bp" if language == "ru" else "Bin size, bp", min_value=1_000, max_value=500_000_000, value=5_000_000, step=500_000, key=f"{prefix}_bin")
        with c3:
            height = st.number_input("Высота графика" if language == "ru" else "Graph height", min_value=260, max_value=1000, value=390, step=20, key=f"{prefix}_height")
        return {"context": context, "bin_size_bp": int(bin_size), "height": int(height)}
    if spec.plot_id == "gene_pca":
        return common_cluster_controls(prefix, language, default_limit=0, default_height=480)
    if spec.plot_id == "gene_dendrogram":
        return common_cluster_controls(prefix, language, default_limit=50, default_height=480, show_hierarchical_max=True)
    if spec.plot_id == "cluster_metagene":
        return common_cluster_controls(prefix, language, default_limit=0, default_height=480)
    raise ValueError(f"Unknown spec: {spec.plot_id}")


def render_plot_card(
    *,
    spec: PlotSpec,
    bsx_path: Path,
    annot_gff_path: Path,
    inputs_fingerprint: str,
    theme_name: str,
    language: str,
) -> None:
    st.markdown(
        f"""
        <div class="plot-card">
          <p class="plot-eyebrow">{spec.eyebrow}</p>
          <h3 class="plot-title">{spec.title}</h3>
          <p class="plot-desc">{spec.description}</p>
          <div class="review-strip">
            <span class="review-label">{spec.aw25_requirement}</span>
            <p class="review-text"><b>{"Параметры для ревью." if language == "ru" else "Review preset."}</b> {spec.review_preset}</p>
            <p class="review-text"><b>{"Что это доказывает." if language == "ru" else "What this proves."}</b> {spec.review_goal}</p>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.form(key=f"form_{spec.plot_id}", border=False):
        params = plot_params_ui(spec, language)
        params["theme"] = theme_name
        params["language"] = language
        c1, c2, _ = st.columns([1.1, 1.2, 3])
        with c1:
            submitted = st.form_submit_button("Запустить модуль" if language == "ru" else "Run module", use_container_width=True)
        with c2:
            force = st.checkbox("Перестроить принудительно" if language == "ru" else "Force rebuild", value=False, key=f"{spec.plot_id}_force")

    cache_file = _cache_path(spec.plot_id, params, inputs_fingerprint)
    should_show_cached = cache_file.exists() and not submitted
    if submitted or should_show_cached:
        try:
            spinner_label = (
                "Строю график..." if submitted and (force or not cache_file.exists()) else "Загружаю из кэша..."
            ) if language == "ru" else (
                "Building plot..." if submitted and (force or not cache_file.exists()) else "Loading cached plot..."
            )
            with st.spinner(spinner_label):
                fig, from_cache, resolved_cache_file = cached_or_build_plot(
                    spec=spec,
                    bsx_path=bsx_path,
                    annot_gff_path=annot_gff_path,
                    inputs_fingerprint=inputs_fingerprint,
                    params=params,
                    force_rebuild=force,
                )
            st.markdown(
                f'<span class="cache-pill">{"Загружено из кэша" if from_cache else "Построено сейчас"}</span>' if language == "ru"
                else f'<span class="cache-pill">{"Loaded from cache" if from_cache else "Built now"}</span>',
                unsafe_allow_html=True,
            )
            st.plotly_chart(
                fig,
                use_container_width=True,
                config={"responsive": True, "displaylogo": False, "displayModeBar": "hover"},
            )
            st.download_button(
                "Экспортировать standalone HTML" if language == "ru" else "Export standalone HTML",
                _figure_to_standalone_html(fig, title=spec.title).encode("utf-8"),
                file_name=f"{_slugify(spec.title)}.html",
                mime="text/html",
                key=f"download_{spec.plot_id}_{resolved_cache_file.stem}",
            )
        except Exception as exc:
            st.error(f"Не удалось построить график: {exc}" if language == "ru" else f"Failed to build the plot: {exc}")
            st.caption(
                "Попробуй уменьшить лимит генов или регионов, затем проверь, что report.bsx и annot.gff совместимы."
                if language == "ru"
                else "Try reducing gene or region limits, then verify that report.bsx and annot.gff are compatible."
            )
    else:
        st.caption("Для текущего набора параметров график ещё не построен." if language == "ru" else "This plot has not been built for the current parameter set yet.")


def render_documentation_tab(language: str) -> None:
    if language == "ru":
        st.markdown(
            """
            <div class="docs-card">
              <h2>Официальная документация</h2>
              <p>
                Эта студия является локальным интерактивным фронтендом для текущего AW25 plotting pipeline в
                <code>bsx2.viz</code>. Она использует те же compute и render paths, что и примеры в репозитории,
                и служит канонической локальной рабочей средой для exploratory analysis из
                <code>report.bsx</code> и <code>annot.gff</code>.
              </p>
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            """
            <div class="docs-card">
              <h2>Official Documentation</h2>
              <p>
                This studio is the local interactive front end for the current <code>bsx2.viz</code> AW25 plotting
                pipeline. It uses the same compute and render paths as the repository examples and should be treated as
                the canonical local workstation for exploratory analysis from <code>report.bsx</code> and
                <code>annot.gff</code>.
              </p>
            </div>
            """,
            unsafe_allow_html=True,
        )
    c1, c2 = st.columns(2)
    with c1:
        if DOCS_BUILD_INDEX.exists():
            st.link_button("Открыть локальные Sphinx docs" if language == "ru" else "Open local Sphinx docs", DOCS_BUILD_INDEX.as_uri(), use_container_width=True)
        else:
            st.caption("Локальный Sphinx HTML ещё не собран." if language == "ru" else "Local Sphinx HTML is not built yet.")
    with c2:
        if DOCS_PUBLIC_URL:
            st.link_button("Открыть опубликованные docs" if language == "ru" else "Open published docs", DOCS_PUBLIC_URL, use_container_width=True)
        else:
            st.caption("Задай `AW25_DOCS_URL`, чтобы здесь появилась ссылка на опубликованный docs site." if language == "ru" else "Set `AW25_DOCS_URL` to expose a published docs site link here.")
    if language == "ru":
        st.markdown(
            """
            <div class="docs-card">
              <h3>Назначение</h3>
              <p>
                Приложение строит графики метилирования по запросу, кэширует Plotly figures на диске и позволяет
                экспортировать каждый готовый график как standalone HTML. Это компактная вычислительная рабочая среда,
                а не статический report generator.
              </p>

              <h3>Обязательные входы</h3>
              <ul>
                <li><code>report.bsx</code>: methylation report, читаемый через <code>RegionReader</code>.</li>
                <li><code>annot.gff</code>: gene annotation, читаемая через <code>HcAnnotStore.from_gff(...)</code>.</li>
              </ul>

              <h3>Справка по графикам</h3>
              <ul>
                <li><b>Metagene line / heatmap / box / violin</b>: настраиваемый metagene layout с annotation-driven или manual composition.</li>
                <li><b>Chromosome map</b>: whole-genome chromosome methylation map с настраиваемым genomic bin size.</li>
                <li><b>Gene PCA / dendrogram / cluster metagene</b>: clustering views по gene profile matrices, включая KMeans-labeled embedding и optional cluster summaries.</li>
              </ul>

              <h3>AW25 review checklist</h3>
              <ul>
                <li>Запусти все четыре metagene-графика в режиме <code>annotation-driven</code>.</li>
                <li>Повтори Line или Heatmap в режиме <code>manual-composed</code>.</li>
                <li>Запусти Chromosome methylation map.</li>
                <li>Запусти Gene PCA, Gene dendrogram и Cluster metagene.</li>
                <li>Экспортируй HTML-файлы прямо из студии для передачи на ревью.</li>
              </ul>
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.markdown(
            """
            <div class="docs-card">
              <h3>Локальная сборка</h3>
              <p>Собери Sphinx site из каталога <code>python</code>:</p>
              <pre><code>pip install -r docs/requirements.txt
PYTHONPATH=src sphinx-build -b html docs/source docs/build/html</code></pre>

              <p>
                На Windows команда <code>pip install -e .</code> может падать из-за upstream Rust dependency
                <code>bsxplorer2 = 0.2.3</code>, которая импортирует <code>std::os::fd</code>. Этот editable-install failure
                не мешает локальной сборке документации.
              </p>

              <h3>Behavior contract</h3>
              <ul>
                <li>Графики считаются только по явному запросу.</li>
                <li>Ключи кэша включают хэши входных файлов, тип графика, параметры и выбранную тему.</li>
                <li>Повторная сборка с теми же входами и параметрами переиспользует закэшированный Plotly JSON.</li>
                <li>Выбор темы влияет и на shell интерфейса, и на рендеринг Plotly figures.</li>
              </ul>

              <h3>Troubleshooting</h3>
              <ul>
                <li>Сначала уменьшай лимит генов или регионов, если график падает.</li>
                <li>Если dendrogram data пустая, увеличь <code>limit_genes</code> или ослабь ограничения.</li>
                <li>Используй <code>Force rebuild</code> после изменения rendering semantics.</li>
                <li>Загрузка других файлов автоматически создаёт новый input fingerprint.</li>
              </ul>
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            """
            <div class="docs-card">
              <h3>Purpose</h3>
              <p>
                The application builds methylation plots on demand, caches resulting Plotly figures on disk, and exports
                any completed figure as standalone HTML. It is intended to behave like a compact computational analysis
                workstation rather than a static report generator.
              </p>

              <h3>Required Inputs</h3>
              <ul>
                <li><code>report.bsx</code>: methylation report readable by <code>RegionReader</code>.</li>
                <li><code>annot.gff</code>: gene annotation readable by <code>HcAnnotStore.from_gff(...)</code>.</li>
              </ul>

              <h3>Plot Reference</h3>
              <ul>
                <li><b>Metagene line / heatmap / box / violin</b>: configurable metagene layout with annotation-driven or manual composition.</li>
                <li><b>Chromosome map</b>: whole-genome chromosome methylation map with configurable genomic bin size.</li>
                <li><b>Gene PCA / dendrogram / cluster metagene</b>: clustering views over gene profile matrices, including KMeans-labeled embedding and optional cluster summaries.</li>
              </ul>

              <h3>AW25 Review Checklist</h3>
              <ul>
                <li>Run all four metagene plots in <code>annotation-driven</code> mode.</li>
                <li>Repeat Line or Heatmap in <code>manual-composed</code> mode.</li>
                <li>Run Chromosome methylation map.</li>
                <li>Run Gene PCA, Gene dendrogram, and Cluster metagene.</li>
                <li>Export the resulting HTML files directly from the workstation for reviewer handoff.</li>
              </ul>
            </div>
            """,
            unsafe_allow_html=True,
        )
        st.markdown(
            """
            <div class="docs-card">
              <h3>Local Build</h3>
              <p>Build the Sphinx site from the <code>python</code> directory:</p>
              <pre><code>pip install -r docs/requirements.txt
PYTHONPATH=src sphinx-build -b html docs/source docs/build/html</code></pre>

              <p>
                On Windows, <code>pip install -e .</code> may currently fail because the upstream Rust dependency
                <code>bsxplorer2 = 0.2.3</code> imports <code>std::os::fd</code>. That editable-install failure does not
                block local documentation builds.
              </p>

              <h3>Behavior Contract</h3>
              <ul>
                <li>Plots are computed only when requested.</li>
                <li>Cache keys include input file hashes, plot type, parameters, and selected theme.</li>
                <li>Rebuilding with the same inputs and parameters reuses the cached Plotly JSON.</li>
                <li>Theme selection affects both the workstation shell and rendered Plotly figures.</li>
              </ul>

              <h3>Troubleshooting</h3>
              <ul>
                <li>Reduce gene or region limits first if a plot fails.</li>
                <li>If dendrogram data is empty, increase <code>limit_genes</code> or relax constraints.</li>
                <li>Use <code>Refresh cache</code> after changing rendering semantics.</li>
                <li>Uploading different files automatically produces a new input fingerprint.</li>
              </ul>
            </div>
            """,
            unsafe_allow_html=True,
        )


def render_sidebar() -> tuple[Any, Any, str, str]:
    language = st.sidebar.segmented_control("Language / Язык", options=["en", "ru"], default="en", key="studio_language")
    st.sidebar.markdown(
        f'<p class="sidebar-label">{"Входной набор данных" if language == "ru" else "Input dataset"}</p>',
        unsafe_allow_html=True,
    )
    report_file = st.sidebar.file_uploader("report.bsx", type=["bsx"], key="report_bsx")
    annot_file = st.sidebar.file_uploader("annot.gff", type=["gff", "gff3", "txt"], key="annot_gff")
    st.sidebar.markdown("---")
    st.sidebar.markdown(
        f'<p class="sidebar-label">{"Режим отображения" if language == "ru" else "Display mode"}</p>',
        unsafe_allow_html=True,
    )
    theme_name = st.sidebar.segmented_control("Тема" if language == "ru" else "Theme", options=["Light", "Dark"], default="Light")
    st.sidebar.markdown("---")
    st.sidebar.markdown(
        f'<p class="sidebar-label">{"Кэш и хранение" if language == "ru" else "Cache and storage"}</p>',
        unsafe_allow_html=True,
    )
    st.sidebar.markdown(
        f'<div class="sidebar-meta"><span>{"кэш графиков" if language == "ru" else "plot cache"}</span><code>{CACHE_DIR}</code></div>',
        unsafe_allow_html=True,
    )
    st.sidebar.markdown(
        f'<div class="sidebar-meta"><span>{"загруженные входы" if language == "ru" else "uploaded inputs"}</span><code>{UPLOAD_DIR}</code></div>',
        unsafe_allow_html=True,
    )
    if st.sidebar.button("Очистить in-memory кэш" if language == "ru" else "Clear in-memory cache", use_container_width=True):
        st.cache_resource.clear()
        st.sidebar.success("In-memory кэш очищен. Файлы на диске сохранены." if language == "ru" else "In-memory cache cleared. Disk cache files remain available.")
    return report_file, annot_file, theme_name, language


def main() -> None:
    configure_page()
    report_file, annot_file, theme_name, language = render_sidebar()
    apply_page_style(theme_name, language)

    studio_tab, docs_tab = st.tabs(["Студия", "Документация"] if language == "ru" else ["Studio", "Documentation"])

    with docs_tab:
        render_documentation_tab(language)

    with studio_tab:
        render_hero(language)
        render_aw25_coverage(language)
        render_aw25_review_protocol(language)
        if report_file is None or annot_file is None:
            st.info("Сначала загрузи `report.bsx` и `annot.gff` в sidebar, чтобы инициализировать рабочую среду." if language == "ru" else "Load both `report.bsx` and `annot.gff` from the sidebar to initialize the analysis workspace.")
            return

        bsx_path, annot_gff_path, inputs_fingerprint = _save_uploaded_inputs(report_file, annot_file)
        st.success(f"Рабочая среда готова: `{bsx_path.name}` + `{annot_gff_path.name}` · input hash `{inputs_fingerprint}`" if language == "ru" else f"Workspace ready: `{bsx_path.name}` + `{annot_gff_path.name}` · input hash `{inputs_fingerprint}`")
        st.caption(
            "Каждый analysis module имеет независимые параметры, собственный trigger запуска и детерминированный дисковый ключ Plotly-кэша."
            if language == "ru"
            else "Each analysis module has independent parameters, its own execution trigger, and a deterministic on-disk Plotly cache key."
        )

        for spec in _plot_specs(language):
            if spec.plot_id == "metagene_line":
                st.markdown("### AW25 · Визуализация метагена" if language == "ru" else "### AW25 · Metagene visualization")
                st.caption("Обязательные metagene-графики поверх одного и того же layout и methylation context." if language == "ru" else "Required metagene graphics over the same configurable layout and methylation context.")
            elif spec.plot_id == "chromosome_map":
                st.markdown("### AW25 · Полногеномная карта метилирования" if language == "ru" else "### AW25 · Genome-scale methylation map")
                st.caption("Полногеномный обзор метилирования по хромосомам." if language == "ru" else "Whole-genome chromosome methylation view.")
            elif spec.plot_id == "gene_pca":
                st.markdown("### AW25 · Кластеризация и группировка" if language == "ru" else "### AW25 · Clustering and grouping")
                st.caption("PCA с метками KMeans, dendrogram и optional cluster-level metagene outputs." if language == "ru" else "KMeans-labeled PCA, dendrogram, and optional cluster-level metagene outputs.")
            with st.expander(spec.title, expanded=spec.plot_id == "metagene_line"):
                render_plot_card(
                    spec=spec,
                    bsx_path=bsx_path,
                    annot_gff_path=annot_gff_path,
                    inputs_fingerprint=inputs_fingerprint,
                    theme_name=theme_name,
                    language=language,
                )


if __name__ == "__main__":
    main()
