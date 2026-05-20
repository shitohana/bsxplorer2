"""Reusable DMR visualization curve specifications.

This module treats DMR outputs as interval-centric visualization inputs and
builds lightweight curve specifications plus optional summary tables/figures.
It intentionally does not implement a new plotting framework.  The spec/result
objects are dashboard-friendly wrappers; rendering delegates to simple existing
Python plotting backends when available, and meta-region line rendering can
reuse the existing ``DiscreteRegionData`` line plot layer.

Input assumptions are deliberately permissive: DMR tables, region-count tables,
and support matrices may use common column aliases that are normalized by
``bsx2.viz.dmr_curve_data``.  Missing optional columns produce warnings rather
than hard failures.

Limitations:
- DMR summary plots can be built from a DMR table alone.
- PCA and heatmap curves require region-counts plus design metadata.
- Line/meta-DMR profiles require positional methylation signal or a precomputed
  ``DiscreteRegionData`` object; a DMR table alone is not enough.
- Curve specs describe visualization intent and do not replace statistical
  testing, Regional Evidence, beta-binomial validation, or caller-specific
  assumptions.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .dmr_curve_data import (
    COUNT_COLUMNS,
    ANNOTATION_COLUMNS,
    classify_annotation_table,
    classify_design_table,
    classify_dmr_table,
    classify_region_counts_table,
    coerce_table,
    find_best_dmr_table,
    normalize_dmr_curve_columns,
    read_annotation_table,
    read_beta_binom_table,
    read_caller_support_table,
    read_design_table,
    read_dmr_table,
    read_region_counts,
)


@dataclass
class DmrCurveSpec:
    """Portable description of a DMR visualization curve."""

    curve_id: str
    curve_type: str
    title: str
    description: str
    input_tables: dict[str, Any]
    data_mapping: dict[str, Any]
    plot_params: dict[str, Any]
    filters: dict[str, Any]
    renderer: str
    created_by: str = "bsx2.dmr_curves"

    def validate(self) -> list[str]:
        warnings: list[str] = []
        for key in ("curve_id", "curve_type", "title", "renderer"):
            if not getattr(self, key):
                warnings.append(f"missing_required_spec_field:{key}")
        return warnings

    def short_label(self) -> str:
        return f"{self.curve_type}:{self.curve_id}"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, path: str | Path) -> None:
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8")

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "DmrCurveSpec":
        return cls(**value)

    @classmethod
    def from_json(cls, path: str | Path) -> "DmrCurveSpec":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


@dataclass
class DmrCurveResult:
    """Curve spec plus optional computed table, figure object, and warnings."""

    spec: DmrCurveSpec
    figure: object | None = None
    table: pd.DataFrame | None = None
    summary: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    quality_status: str = "exploratory"
    quality_reasons: list[str] = field(default_factory=list)

    def save_table(self, path: str | Path) -> bool:
        if self.table is None:
            return False
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        self.table.to_csv(out_path, sep="\t", index=False)
        return True

    def save_summary(self, path: str | Path) -> None:
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "spec": self.spec.to_dict(),
            "summary": _json_safe(self.summary),
            "warnings": list(self.warnings),
            "quality_status": self.quality_status,
            "quality_reasons": list(self.quality_reasons),
        }
        out_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    def save_figure(self, path: str | Path) -> bool:
        if self.figure is None:
            return False
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if hasattr(self.figure, "savefig"):
            self.figure.savefig(out_path, bbox_inches="tight", dpi=150)
            return True
        if hasattr(self.figure, "write_html"):
            self.figure.write_html(str(out_path))
            return True
        try:
            import holoviews as hv  # type: ignore

            hv.save(self.figure, out_path)
            return True
        except Exception:
            return False

    def to_dashboard_payload(self) -> dict[str, Any]:
        preview: list[dict[str, Any]] = []
        if self.table is not None:
            preview = self.table.head(20).replace({np.nan: None}).to_dict(orient="records")
        return {
            "curve_id": self.spec.curve_id,
            "curve_type": self.spec.curve_type,
            "title": self.spec.title,
            "description": self.spec.description,
            "summary": _json_safe(self.summary),
            "warnings": list(self.warnings),
            "quality_status": self.quality_status,
            "quality_reasons": list(self.quality_reasons),
            "table_preview": _json_safe(preview),
            "figure_type": _figure_type(self.figure),
            "spec": self.spec.to_dict(),
        }


@dataclass
class DmrCurveBundle:
    """Collection of portable DMR curve specs for dashboards or reports."""

    bundle_id: str
    curves: list[DmrCurveSpec]
    inputs: dict[str, Any]
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "bundle_id": self.bundle_id,
            "curves": [curve.to_dict() for curve in self.curves],
            "inputs": self.inputs,
            "notes": self.notes,
        }

    def to_json(self, path: str | Path) -> None:
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8")

    @classmethod
    def from_json(cls, path: str | Path) -> "DmrCurveBundle":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(
            bundle_id=payload["bundle_id"],
            curves=[DmrCurveSpec.from_dict(item) for item in payload.get("curves", [])],
            inputs=payload.get("inputs", {}),
            notes=payload.get("notes", []),
        )

    def list_curves(self) -> list[str]:
        return [curve.short_label() for curve in self.curves]

    def render_all(self, results: Iterable[DmrCurveResult] | None = None) -> list[DmrCurveResult]:
        return list(results or [])

    def save_manifest(self, path: str | Path, results: Iterable[DmrCurveResult] | None = None) -> None:
        rows = []
        result_by_id = {result.spec.curve_id: result for result in (results or [])}
        for spec in self.curves:
            result = result_by_id.get(spec.curve_id)
            rows.append(
                {
                    "curve_id": spec.curve_id,
                    "curve_type": spec.curve_type,
                    "title": spec.title,
                    "renderer": spec.renderer,
                    "has_table": bool(result and result.table is not None),
                    "has_figure": bool(result and result.figure is not None),
                    "quality_status": result.quality_status if result else "unknown",
                    "quality_reasons": ";".join(result.quality_reasons) if result else "",
                    "n_warnings": len(result.warnings) if result else 0,
                    "warnings": ";".join(result.warnings) if result else "",
                }
            )
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def make_dmr_line_curve(
    dmr_table: Any | None = None,
    discrete_region_data: Any | None = None,
    mode: str = "scaled",
    curve_id: str = "dmr_line_profile",
    render: bool = True,
) -> DmrCurveResult:
    """Create a DMR-centered/scaled line profile spec.

    A DMR table alone cannot define a line profile.  The function returns a
    warning-only result unless positional signal or ``DiscreteRegionData`` is
    supplied.
    """

    spec = _spec(
        curve_id,
        "dmr_line",
        "DMR methylation profile",
        "DMR-centered or scaled methylation profile using positional signal.",
        {"dmr_table": _input_name(dmr_table), "discrete_region_data": bool(discrete_region_data)},
        {"x": "relative_position", "y": "methylation_density"},
        {"mode": mode},
        {},
        "bsx2.viz.render.line.line_plot",
    )
    warnings: list[str] = []
    figure = None
    if discrete_region_data is None:
        warnings.append(
            "line_profile_requires_positional_signal_or_precomputed_DiscreteRegionData"
        )
        quality_status = "skipped"
        quality_reasons = ["requires positional methylation signal or precomputed DiscreteRegionData"]
    elif render:
        quality_status = "thesis_ready"
        quality_reasons = ["precomputed DiscreteRegionData supplied"]
        try:
            from .render.line import line_plot

            figure = line_plot(discrete_region_data, title=spec.title)
        except Exception as exc:
            warnings.append(f"line_render_failed:{exc}")
            quality_status = "warning"
            quality_reasons = ["DiscreteRegionData supplied but rendering failed"]
    else:
        quality_status = "thesis_ready"
        quality_reasons = ["precomputed DiscreteRegionData supplied"]
    return DmrCurveResult(
        spec=spec,
        figure=figure,
        table=None,
        summary={"mode": mode},
        warnings=warnings,
        quality_status=quality_status,
        quality_reasons=quality_reasons,
    )


def make_dmr_chromosome_curve(
    dmr_table: Any,
    group_by: str | None = None,
    curve_id: str = "dmr_chromosome_distribution",
    render: bool = True,
    input_scope: str = "input",
) -> DmrCurveResult:
    dmr, warnings = coerce_table(dmr_table, read_dmr_table)
    semantic_type, semantic_warnings = classify_dmr_table(dmr)
    warnings.extend(semantic_warnings)
    n_unique_chromosomes = int(dmr["chrom"].nunique(dropna=True)) if dmr is not None and "chrom" in dmr else 0
    title_suffix = ""
    quality_status = "thesis_ready"
    quality_reasons = ["region-level DMR table with chromosome column"]
    if input_scope.startswith("top") or semantic_type == "top_region_level_dmr":
        title_suffix = f", top {len(dmr)} input"
        quality_reasons.append("input is a top/subset table, not genome-wide distribution")
    if n_unique_chromosomes == 1:
        warnings.append("chromosome_distribution_single_chromosome_input_check_subset_or_column_mapping")
        title_suffix = " (input has 1 chromosome)"
        quality_status = "warning"
        quality_reasons = ["single chromosome input may indicate a subset/top table or chrom column mapping issue"]
    spec = _spec(
        curve_id,
        "dmr_chromosome",
        "DMR distribution by chromosome" + title_suffix,
        "Counts DMR intervals per chromosome and optional group.",
        {"dmr_table": _input_name(dmr_table)},
        {"x": "chrom", "y": "n_dmrs_in_input", "group": group_by},
        {"input_scope": input_scope},
        {},
        "bsx2.viz.compute.chrline_or_table",
    )
    if dmr.empty or "chrom" not in dmr:
        warnings.append("dmr_chromosome_curve_requires_chrom_column")
        table = pd.DataFrame(columns=["chrom", "n_dmrs_in_input"])
        quality_status = "skipped"
        quality_reasons = ["chromosome distribution requires a chrom column"]
    elif group_by and group_by in dmr:
        table = dmr.groupby(["chrom", group_by], dropna=False).size().reset_index(name="n_dmrs_in_input")
    else:
        table = dmr.groupby("chrom", dropna=False).size().reset_index(name="n_dmrs_in_input")
    figure = _bar_figure(table, x="chrom", y="n_dmrs_in_input", title=spec.title, y_label="n_DMRs_in_input") if render else None
    n_rows = int(table["n_dmrs_in_input"].sum()) if "n_dmrs_in_input" in table else 0
    top_fraction = float(table["n_dmrs_in_input"].max() / n_rows) if n_rows else 0.0
    summary = {
        "semantic_type": semantic_type,
        "n_chromosomes": int(table["chrom"].nunique()) if "chrom" in table else 0,
        "n_dmrs_in_input": n_rows,
        "top_chromosome_fraction": top_fraction,
    }
    if n_unique_chromosomes == 1:
        summary["interpretation_note"] = "This may indicate a subset/top table or chrom column mapping issue."
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_box_curve(
    dmr_table: Any | None = None,
    region_counts: Any | None = None,
    design: Any | None = None,
    mode: str = "delta_by_context",
    curve_id: str = "dmr_box",
    render: bool = True,
) -> DmrCurveResult:
    title = "DMR methylation by condition" if mode == "methylation_by_condition" else "DMR delta methylation by context"
    y_label = "methylation = mC / total" if mode == "methylation_by_condition" else "delta methylation"
    spec = _spec(
        curve_id,
        "dmr_box",
        title,
        "Box plot for DMR delta/q values or region methylation by condition.",
        {"dmr_table": _input_name(dmr_table), "region_counts": _input_name(region_counts), "design": _input_name(design)},
        {"mode": mode},
        {},
        {},
        "bsx2.viz.render.box_or_matplotlib_boxplot",
    )
    table, warnings = _distribution_table(dmr_table, region_counts, design, mode)
    figure = _box_figure(table, group="group", value="value", title=spec.title, y_label=y_label) if render else None
    summary = {"mode": mode, "n_values": int(len(table))}
    quality_status = "thesis_ready" if not table.empty and not _has_required_warning(warnings) else ("skipped" if table.empty else "warning")
    quality_reasons = [_mode_quality_reason(mode)] if quality_status == "thesis_ready" else warnings
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_violin_curve(
    dmr_table: Any | None = None,
    region_counts: Any | None = None,
    design: Any | None = None,
    mode: str = "delta_by_context",
    curve_id: str = "dmr_violin",
    render: bool = True,
) -> DmrCurveResult:
    title = "DMR methylation by condition" if mode == "methylation_by_condition" else "DMR delta methylation by context"
    y_label = "methylation = mC / total" if mode == "methylation_by_condition" else "delta methylation"
    spec = _spec(
        curve_id,
        "dmr_violin",
        title,
        "Violin plot for DMR delta/q values or region methylation by condition.",
        {"dmr_table": _input_name(dmr_table), "region_counts": _input_name(region_counts), "design": _input_name(design)},
        {"mode": mode},
        {},
        {},
        "bsx2.viz.render.violin_or_matplotlib_violinplot",
    )
    table, warnings = _distribution_table(dmr_table, region_counts, design, mode)
    figure = _violin_figure(table, group="group", value="value", title=spec.title, y_label=y_label) if render else None
    summary = {"mode": mode, "n_values": int(len(table))}
    quality_status = "thesis_ready" if not table.empty and not _has_required_warning(warnings) else ("skipped" if table.empty else "warning")
    quality_reasons = [_mode_quality_reason(mode)] if quality_status == "thesis_ready" else warnings
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_pca_curve(
    region_counts: Any,
    design: Any,
    min_total: int = 5,
    top_n: int = 1000,
    impute: str = "median",
    curve_id: str = "dmr_pca",
    render: bool = True,
) -> DmrCurveResult:
    spec = _spec(
        curve_id,
        "dmr_pca",
        "DMR methylation PCA",
        "Sample PCA from region-by-sample DMR methylation matrix.",
        {"region_counts": _input_name(region_counts), "design": _input_name(design)},
        {"rows": "regions", "columns": "samples", "value": "Y/m"},
        {"min_total": min_total, "top_n": top_n, "impute": impute},
        {},
        "numpy_svd_compatible_with_dimred_layer",
    )
    matrix, matrix_warnings = _region_methylation_matrix(region_counts, design, min_total, top_n, impute)
    warnings = list(matrix_warnings)
    if matrix.empty or matrix.shape[1] < 2:
        warnings.append("pca_requires_at_least_two_samples")
        table = pd.DataFrame(columns=["sample_id", "PC1", "PC2"])
        summary = {"explained_variance": []}
        figure = None
    else:
        filled = matrix.T.to_numpy(dtype=float)
        centered = filled - filled.mean(axis=0, keepdims=True)
        try:
            u, s, _ = np.linalg.svd(centered, full_matrices=False)
            coords = u[:, :2] * s[:2]
            if coords.shape[1] == 1:
                coords = np.column_stack([coords[:, 0], np.zeros(coords.shape[0])])
            variance = (s**2) / max(centered.shape[0] - 1, 1)
            explained = (variance / variance.sum()).tolist() if variance.sum() > 0 else [0.0, 0.0]
            table = pd.DataFrame({"sample_id": matrix.columns, "PC1": coords[:, 0], "PC2": coords[:, 1]})
            design_df, design_warnings = coerce_table(design, read_design_table)
            warnings.extend(design_warnings)
            if "sample_id" in design_df:
                table = table.merge(design_df.drop_duplicates("sample_id"), on="sample_id", how="left")
            summary = {"explained_variance": explained[:2], "n_regions": int(matrix.shape[0]), "n_samples": int(matrix.shape[1])}
            x_label = f"PC1 ({explained[0] * 100:.1f}% variance)" if explained else "PC1"
            y_label = f"PC2 ({explained[1] * 100:.1f}% variance)" if len(explained) > 1 else "PC2"
            figure = _scatter_figure(table, "PC1", "PC2", "condition", spec.title, x_label=x_label, y_label=y_label) if render else None
        except Exception as exc:
            warnings.append(f"pca_failed:{exc}")
            table = pd.DataFrame(columns=["sample_id", "PC1", "PC2"])
            summary = {"explained_variance": []}
            figure = None
    required_warnings = _has_required_warning(warnings)
    quality_status = "thesis_ready" if not table.empty and not required_warnings else ("skipped" if table.empty else "warning")
    quality_reasons = ["valid region_counts/design sample mapping; PCA scores computed"] if quality_status == "thesis_ready" else warnings
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_heatmap_curve(
    region_counts: Any,
    design: Any | None = None,
    dmr_table: Any | None = None,
    top_n: int = 1000,
    min_total: int = 5,
    sort_by: str = "variance",
    curve_id: str = "dmr_heatmap",
    render: bool = True,
) -> DmrCurveResult:
    spec = _spec(
        curve_id,
        "dmr_heatmap",
        f"DMR methylation heatmap, top {top_n} variable regions",
        "Region-by-sample methylation matrix for DMRs.",
        {"region_counts": _input_name(region_counts), "design": _input_name(design), "dmr_table": _input_name(dmr_table)},
        {"rows": "regions", "columns": "samples", "value": "Y/m"},
        {"top_n": top_n, "min_total": min_total, "sort_by": sort_by},
        {},
        "bsx2.viz.render.heatmap_or_matrix_table",
    )
    matrix, warnings = _region_methylation_matrix(region_counts, design, min_total, top_n, "median")
    if sort_by in {"abs_delta", "q_value", "evidence_class"} and dmr_table is not None:
        dmr, dmr_warnings = coerce_table(dmr_table, read_dmr_table)
        warnings.extend(dmr_warnings)
        matrix = _sort_matrix_by_dmr_table(matrix, dmr, sort_by, top_n)
    table = matrix.reset_index().rename(columns={"index": "region_id"})
    figure = _heatmap_figure(matrix, spec.title) if render else None
    summary = {"n_regions": int(matrix.shape[0]), "n_samples": int(matrix.shape[1]), "selection": sort_by}
    required_warnings = _has_required_warning([w for w in warnings if w != "matrix_imputation_median_used"])
    quality_status = "thesis_ready" if not table.empty and not required_warnings else ("skipped" if table.empty else "warning")
    quality_reasons = ["valid region_counts/design sample mapping; DMR methylation matrix created"] if quality_status == "thesis_ready" else warnings
    if quality_status == "thesis_ready" and "matrix_imputation_median_used" in warnings:
        quality_reasons.append("median imputation used for missing methylation values")
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_volcano_curve(
    dmr_table: Any,
    curve_id: str = "dmr_volcano",
    render: bool = True,
    q_cap: float = 50.0,
) -> DmrCurveResult:
    dmr, warnings = coerce_table(dmr_table, read_dmr_table)
    spec = _spec(
        curve_id,
        "dmr_volcano",
        "DMR volcano plot",
        "Delta methylation versus capped -log10(q-value).",
        {"dmr_table": _input_name(dmr_table)},
        {"x": "delta", "y": f"-log10(q_value), capped at {q_cap:g}"},
        {"q_cap": q_cap},
        {},
        "matplotlib_scatter_fallback",
    )
    if "delta" not in dmr or "q_value" not in dmr:
        warnings.append("volcano_requires_delta_and_q_value")
        table = pd.DataFrame(columns=["region_id", "delta", "q_value", "neg_log10_q_capped", "q_capped"])
        figure = None
        quality_status = "skipped"
        quality_reasons = ["volcano requires delta and q_value columns"]
    else:
        table = dmr[[col for col in ["region_id", "chrom", "start", "end", "context", "delta", "q_value", "evidence_class"] if col in dmr]].copy()
        transformed = safe_neg_log10_q(table["q_value"], cap=q_cap)
        table["neg_log10_q_capped"] = transformed["neg_log10_q_capped"]
        table["q_capped"] = transformed["q_capped"]
        capped_count = int(table["q_capped"].sum())
        if capped_count:
            warnings.append(f"volcano_q_values_capped_count={capped_count}")
        figure = _scatter_figure(
            table,
            "delta",
            "neg_log10_q_capped",
            "evidence_class",
            spec.title,
            y_label=f"-log10(q), capped at {q_cap:g}",
            hline=-math.log10(0.05),
            alpha=0.5,
        ) if render else None
        quality_status = "warning" if capped_count > max(10, len(table) * 0.1) else "thesis_ready"
        quality_reasons = ["q-values transformed with finite cap for plotting"]
        if quality_status == "warning":
            quality_reasons.append("many q-values were zero or below plotting cap")
    summary = {"n_points": int(len(table))}
    if "q_capped" in table:
        summary["q_capped_count"] = int(table["q_capped"].sum())
        summary["q_cap"] = q_cap
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_evidence_bar_curve(
    dmr_table: Any,
    curve_id: str = "dmr_evidence_class_counts",
    render: bool = True,
) -> DmrCurveResult:
    dmr, warnings = coerce_table(dmr_table, read_dmr_table)
    spec = _spec(
        curve_id,
        "dmr_evidence_bar",
        "DMR evidence class counts",
        "Counts DMRs by evidence class.",
        {"dmr_table": _input_name(dmr_table)},
        {"x": "evidence_class", "y": "n_dmrs"},
        {},
        {},
        "bar_table",
    )
    if "evidence_class" not in dmr:
        warnings.append("evidence_bar_requires_evidence_class")
        table = pd.DataFrame(columns=["evidence_class", "n_dmrs"])
        quality_status = "skipped"
        quality_reasons = ["evidence_class column is required"]
    else:
        table = dmr.groupby("evidence_class", dropna=False).size().reset_index(name="n_dmrs")
        quality_status = "thesis_ready"
        quality_reasons = ["region-level DMR table includes evidence_class"]
    figure = _bar_figure(table, "evidence_class", "n_dmrs", spec.title) if render else None
    summary = {"n_classes": int(len(table)), "n_dmrs": int(table["n_dmrs"].sum()) if "n_dmrs" in table else 0}
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_caller_support_curve(
    caller_support: Any,
    curve_id: str = "dmr_caller_support",
    render: bool = True,
) -> DmrCurveResult:
    support, warnings = coerce_table(caller_support, read_caller_support_table)
    spec = _spec(
        curve_id,
        "dmr_caller_support",
        "External caller support",
        "Distribution of DMR support across external callers.",
        {"caller_support": _input_name(caller_support)},
        {"x": "n_callers_supporting", "y": "n_regions"},
        {},
        {},
        "bar_table",
    )
    if "n_callers_supporting" in support:
        table = support.groupby("n_callers_supporting", dropna=False).size().reset_index(name="n_regions")
        quality_status = "thesis_ready"
        quality_reasons = ["caller support matrix includes n_callers_supporting"]
    else:
        support_cols = [col for col in support.columns if col.endswith("_support")]
        if support_cols:
            tmp = support.copy()
            tmp["n_callers_supporting"] = tmp[support_cols].apply(lambda row: sum(_truthy(value) for value in row), axis=1)
            table = tmp.groupby("n_callers_supporting", dropna=False).size().reset_index(name="n_regions")
            quality_status = "thesis_ready"
            quality_reasons = ["caller support matrix includes caller-specific support columns"]
        else:
            warnings.append("caller_support_curve_requires_n_callers_supporting_or_support_columns")
            table = pd.DataFrame(columns=["n_callers_supporting", "n_regions"])
            quality_status = "skipped"
            quality_reasons = ["caller support columns are missing"]
    figure = _bar_figure(table, "n_callers_supporting", "n_regions", spec.title) if render else None
    summary = {"n_regions": int(table["n_regions"].sum()) if "n_regions" in table else 0}
    return DmrCurveResult(spec, figure, table, summary, warnings, quality_status, quality_reasons)


def make_dmr_annotation_summary_curve(
    annotation: pd.DataFrame,
    warnings: list[str] | None = None,
    input_name: str = "DataFrame",
    render: bool = True,
) -> DmrCurveResult:
    """Create a semantically checked DMR annotation/enrichment summary curve."""

    return _annotation_summary_curve(annotation, list(warnings or []), input_name, render)


def build_default_dmr_curve_bundle(
    dmr_table: Any,
    region_counts: Any | None = None,
    design: Any | None = None,
    beta_binom: Any | None = None,
    caller_support: Any | None = None,
    annotation: Any | None = None,
    out_dir: str | Path | None = None,
    top_n: int = 1000,
    min_total: int = 5,
    render: bool = True,
    thesis_ready_only: bool = False,
) -> tuple[DmrCurveBundle, list[DmrCurveResult]]:
    """Build and optionally save a default reusable DMR curve bundle."""

    dmr, dmr_warnings = coerce_table(dmr_table, read_dmr_table)
    results: list[DmrCurveResult] = [
        make_dmr_chromosome_curve(dmr, render=render, input_scope="full_input"),
        make_dmr_evidence_bar_curve(dmr, render=render),
        make_dmr_box_curve(dmr_table=dmr, mode="delta_by_context", curve_id="dmr_box_delta_by_context", render=render),
        make_dmr_violin_curve(dmr_table=dmr, mode="delta_by_context", curve_id="dmr_violin_delta_by_context", render=render),
        make_dmr_volcano_curve(dmr, render=render),
        make_dmr_line_curve(dmr_table=dmr, render=render),
    ]
    if dmr_warnings:
        results[0].warnings = dmr_warnings + results[0].warnings

    if region_counts is not None and design is not None:
        results.append(make_dmr_heatmap_curve(region_counts, design, dmr, top_n=top_n, min_total=min_total, render=render))
        results.append(make_dmr_pca_curve(region_counts, design, min_total=min_total, top_n=top_n, render=render))
        results.append(
            make_dmr_box_curve(region_counts=region_counts, design=design, mode="methylation_by_condition", curve_id="dmr_box_methylation_by_condition", render=render)
        )
        results.append(
            make_dmr_violin_curve(region_counts=region_counts, design=design, mode="methylation_by_condition", curve_id="dmr_violin_methylation_by_condition", render=render)
        )
    if caller_support is not None:
        results.append(make_dmr_caller_support_curve(caller_support, render=render))
    if beta_binom is not None:
        beta, beta_warnings = coerce_table(beta_binom, read_beta_binom_table)
        results.append(_beta_binom_summary_curve(beta, beta_warnings, _input_name(beta_binom), render=render))
    if annotation is not None:
        ann, ann_warnings = coerce_table(annotation, read_annotation_table)
        results.append(_annotation_summary_curve(ann, ann_warnings, _input_name(annotation), render=render))

    bundle = DmrCurveBundle(
        bundle_id=f"dmr_curve_bundle_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}",
        curves=[result.spec for result in results],
        inputs={
            "dmr_table": _input_name(dmr_table),
            "region_counts": _input_name(region_counts),
            "design": _input_name(design),
            "beta_binom": _input_name(beta_binom),
            "caller_support": _input_name(caller_support),
            "annotation": _input_name(annotation),
            "top_n": top_n,
            "min_total": min_total,
            "thesis_ready_only": thesis_ready_only,
        },
        notes=[
            "Curve specs are reusable dashboard objects and do not replace statistical testing.",
            "Line/meta-DMR curve is warning-only unless positional signal or DiscreteRegionData is supplied.",
        ],
    )
    if out_dir is not None:
        save_dmr_curve_bundle_outputs(bundle, results, out_dir, thesis_ready_only=thesis_ready_only)
    return bundle, results


def save_dmr_curve_bundle_outputs(
    bundle: DmrCurveBundle,
    results: list[DmrCurveResult],
    out_dir: str | Path,
    formats: Iterable[str] = ("json", "tsv", "png"),
    thesis_ready_only: bool = False,
) -> None:
    root = Path(out_dir)
    specs_dir = root / "specs"
    tables_dir = root / "tables"
    figures_dir = root / "figures"
    logs_dir = root / "logs"
    for path in (specs_dir, tables_dir, figures_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=True)
    bundle.to_json(root / "dmr_curve_bundle.json")
    bundle.save_manifest(root / "dmr_curve_manifest.tsv", results)
    warning_rows = []
    for result in results:
        result.spec.to_json(specs_dir / f"{result.spec.curve_id}.json")
        save_payload = (not thesis_ready_only) or result.quality_status == "thesis_ready"
        if save_payload and result.table is not None and "tsv" in formats:
            result.save_table(tables_dir / f"{result.spec.curve_id}.tsv")
        if save_payload and result.figure is not None and "png" in formats:
            result.save_figure(figures_dir / f"{result.spec.curve_id}.png")
        result.save_summary(logs_dir / f"{result.spec.curve_id}.summary.json")
        for warning in _dedupe_strings(result.warnings):
            warning_rows.append({"curve_id": result.spec.curve_id, "warning": warning})
    pd.DataFrame(warning_rows, columns=["curve_id", "warning"]).to_csv(root / "warnings.tsv", sep="\t", index=False)
    _write_bundle_summary(root / "dmr_curve_bundle_summary.md", bundle, results)


def _distribution_table(
    dmr_table: Any | None,
    region_counts: Any | None,
    design: Any | None,
    mode: str,
) -> tuple[pd.DataFrame, list[str]]:
    warnings: list[str] = []
    if mode == "methylation_by_condition":
        if region_counts is None or design is None:
            return pd.DataFrame(columns=["group", "value"]), ["methylation_by_condition_requires_region_counts_and_design"]
        counts, count_warnings = coerce_table(region_counts, read_region_counts)
        design_df, design_warnings = coerce_table(design, read_design_table)
        warnings.extend(count_warnings + design_warnings)
        _, counts_class_warnings = classify_region_counts_table(counts)
        _, design_class_warnings = classify_design_table(design_df, counts)
        warnings.extend(counts_class_warnings + design_class_warnings)
        if not {"region_id", "sample_id", "Y", "m"}.issubset(counts.columns) or "sample_id" not in design_df:
            warnings.append("region_counts_missing_required_columns_for_methylation_distribution")
            return pd.DataFrame(columns=["group", "value"]), warnings
        data = counts.copy()
        data["value"] = np.where(data["m"] > 0, data["Y"] / data["m"], np.nan)
        if "condition" not in data.columns and "condition" in design_df.columns:
            data = data.merge(design_df[["sample_id", "condition"]].drop_duplicates("sample_id"), on="sample_id", how="left")
        elif "condition" in data.columns and "condition" in design_df.columns:
            data = data.merge(
                design_df[["sample_id", "condition"]].drop_duplicates("sample_id"),
                on="sample_id",
                how="left",
                suffixes=("", "_design"),
            )
        condition_col = "condition" if "condition" in data.columns else "condition_design"
        if condition_col not in data.columns:
            warnings.append("methylation_by_condition_missing_condition_column;using_unknown")
            data["group"] = "unknown"
        else:
            data["group"] = data[condition_col].fillna("unknown")
        return data[["region_id", "sample_id", "group", "value"]].dropna(subset=["value"]), warnings

    dmr, dmr_warnings = coerce_table(dmr_table, read_dmr_table)
    warnings.extend(dmr_warnings)
    if mode == "q_by_evidence_class":
        value_col, group_col = "q_value", "evidence_class"
    else:
        value_col, group_col = "delta", "context"
    if value_col not in dmr:
        warnings.append(f"{mode}_requires_{value_col}")
        return pd.DataFrame(columns=["group", "value"]), warnings
    if group_col not in dmr:
        warnings.append(f"{mode}_requires_{group_col}")
        return pd.DataFrame(columns=["group", "value"]), warnings
    table = dmr[[col for col in ["region_id", group_col, value_col] if col in dmr]].copy()
    table = table.rename(columns={group_col: "group", value_col: "value"})
    return table.dropna(subset=["value"]), warnings


def _region_methylation_matrix(
    region_counts: Any,
    design: Any | None,
    min_total: int,
    top_n: int,
    impute: str,
) -> tuple[pd.DataFrame, list[str]]:
    counts, warnings = coerce_table(region_counts, read_region_counts)
    _, count_class_warnings = classify_region_counts_table(counts)
    warnings.extend(count_class_warnings)
    design_df = None
    if design is not None:
        design_df, design_warnings = coerce_table(design, read_design_table)
        _, design_class_warnings = classify_design_table(design_df, counts)
        warnings.extend(design_warnings + design_class_warnings)
    if not {"region_id", "sample_id", "Y", "m"}.issubset(counts.columns):
        warnings.append("region_counts_matrix_requires_region_id_sample_id_Y_m")
        return pd.DataFrame(), warnings
    data = counts.copy()
    data.loc[data["m"] < min_total, "Y"] = np.nan
    data["methylation"] = np.where(data["m"] >= min_total, data["Y"] / data["m"], np.nan)
    matrix = data.pivot_table(index="region_id", columns="sample_id", values="methylation", aggfunc="mean")
    if design_df is not None and "sample_id" in design_df.columns:
        sample_order = [sample for sample in design_df["sample_id"].astype(str).tolist() if sample in set(matrix.columns.astype(str))]
        if sample_order:
            matrix.columns = matrix.columns.astype(str)
            matrix = matrix.loc[:, sample_order]
    if matrix.empty:
        warnings.append("region_methylation_matrix_empty_after_filtering")
        return matrix, warnings
    variances = matrix.var(axis=1, skipna=True).sort_values(ascending=False)
    keep = variances.head(top_n).index
    matrix = matrix.loc[keep]
    if impute == "drop":
        matrix = matrix.dropna(axis=0)
    else:
        if matrix.isna().any().any():
            warnings.append("matrix_imputation_median_used")
        matrix = matrix.apply(lambda row: row.fillna(row.median()), axis=1)
        matrix = matrix.fillna(matrix.stack().median() if not matrix.stack().empty else 0.0)
    return matrix, warnings


def _sort_matrix_by_dmr_table(matrix: pd.DataFrame, dmr: pd.DataFrame, sort_by: str, top_n: int) -> pd.DataFrame:
    if matrix.empty or "region_id" not in dmr:
        return matrix.head(top_n)
    dmr_indexed = dmr.drop_duplicates("region_id").set_index("region_id")
    common = [idx for idx in matrix.index if idx in dmr_indexed.index]
    if not common:
        return matrix.head(top_n)
    if sort_by == "abs_delta" and "delta" in dmr_indexed:
        order = dmr_indexed.loc[common, "delta"].abs().sort_values(ascending=False).index
    elif sort_by == "q_value" and "q_value" in dmr_indexed:
        order = dmr_indexed.loc[common, "q_value"].sort_values(ascending=True).index
    elif sort_by == "evidence_class" and "evidence_class" in dmr_indexed:
        ranks = {"strong": 0, "moderate": 1, "weak": 2, "candidate_only": 3, "qc_limited": 4}
        order = dmr_indexed.loc[common, "evidence_class"].map(ranks).fillna(99).sort_values().index
    else:
        return matrix.head(top_n)
    return matrix.loc[list(order)[:top_n]]


def _beta_binom_summary_curve(beta: pd.DataFrame, warnings: list[str], input_name: str, render: bool) -> DmrCurveResult:
    spec = _spec(
        "dmr_beta_binomial_summary",
        "dmr_beta_binomial_bar",
        "Beta-binomial validation summary",
        "Counts beta-binomial validation status classes.",
        {"beta_binom": input_name},
        {"x": "qc_flag_or_significance", "y": "n_regions"},
        {},
        {},
        "bar_table",
    )
    if beta.empty:
        table = pd.DataFrame(columns=["status", "n_regions"])
        quality_status = "skipped"
        quality_reasons = ["beta-binomial table empty"]
    elif "qc_flag" in beta:
        table = beta.groupby("qc_flag", dropna=False).size().reset_index(name="n_regions").rename(columns={"qc_flag": "status"})
        quality_status = "thesis_ready"
        quality_reasons = ["beta-binomial QC flags available"]
    elif "significant_beta_binom" in beta:
        table = beta.groupby("significant_beta_binom", dropna=False).size().reset_index(name="n_regions").rename(columns={"significant_beta_binom": "status"})
        quality_status = "thesis_ready"
        quality_reasons = ["beta-binomial significance column available"]
    else:
        warnings.append("beta_binom_summary_missing_qc_or_significance_columns")
        table = pd.DataFrame(columns=["status", "n_regions"])
        quality_status = "skipped"
        quality_reasons = ["beta-binomial summary requires qc_flag or significant_beta_binom"]
    figure = _bar_figure(table, "status", "n_regions", spec.title) if render else None
    return DmrCurveResult(spec, figure, table, {"n_rows": int(len(beta))}, warnings, quality_status, quality_reasons)


def _annotation_summary_curve(annotation: pd.DataFrame, warnings: list[str], input_name: str, render: bool) -> DmrCurveResult:
    semantic_type, semantic_warnings = classify_annotation_table(annotation)
    warnings.extend(semantic_warnings)
    warnings = _dedupe_strings(warnings)
    lower = {str(col).strip().lower().replace("-", "_").replace(".", "_").replace(" ", "_"): col for col in annotation.columns}
    ann_col = next((lower[key] for key in [c.lower() for c in ANNOTATION_COLUMNS] if key in lower), None)
    count_col = next((lower[key] for key in [c.lower() for c in COUNT_COLUMNS] if key in lower), None)
    title = "DMR annotation composition"
    y_col = "n_DMRs"
    table = pd.DataFrame(columns=["annotation", y_col])
    quality_status = "skipped"
    quality_reasons = ["annotation summary requires region-level annotation or a numeric count column"]
    if semantic_type == "region_level_annotation" and ann_col:
        table = annotation.groupby(ann_col, dropna=False).size().reset_index(name=y_col).rename(columns={ann_col: "annotation"})
        quality_status = "thesis_ready"
        quality_reasons = ["region-level annotation input with DMR identifiers/coordinates"]
    elif semantic_type == "enrichment_summary" and ann_col and count_col:
        title = "Annotation enrichment summary"
        y_col = count_col
        table = annotation[[ann_col, count_col]].copy().rename(columns={ann_col: "annotation"})
        table[count_col] = pd.to_numeric(table[count_col], errors="coerce")
        warnings = _dedupe_strings(warnings + ["annotation_input_is_summary_not_region_level"])
        quality_status = "exploratory"
        quality_reasons = ["annotation input is an enrichment/summary table, not one row per DMR"]
    else:
        warnings = _dedupe_strings(warnings + ["annotation_summary_not_plotted_requires_region_level_or_count_column"])
    spec = _spec(
        "dmr_annotation_summary",
        "dmr_annotation_bar",
        title,
        "Annotation composition for DMRs or explicit annotation enrichment counts.",
        {"annotation": input_name},
        {"x": "annotation", "y": y_col},
        {"semantic_type": semantic_type},
        {},
        "bar_table",
    )
    figure = _bar_figure(table, "annotation", y_col, spec.title, y_label=y_col) if render and not table.empty else None
    return DmrCurveResult(spec, figure, table, {"n_rows": int(len(annotation)), "semantic_type": semantic_type}, warnings, quality_status, quality_reasons)


def _spec(
    curve_id: str,
    curve_type: str,
    title: str,
    description: str,
    input_tables: dict[str, Any],
    data_mapping: dict[str, Any],
    plot_params: dict[str, Any],
    filters: dict[str, Any],
    renderer: str,
) -> DmrCurveSpec:
    return DmrCurveSpec(curve_id, curve_type, title, description, input_tables, data_mapping, plot_params, filters, renderer)


def _has_required_warning(warnings: Iterable[str]) -> bool:
    markers = (
        "requires_",
        "_requires_",
        "missing_required",
        "missing_columns",
        "missing_condition",
        "missing_group",
        "missing_coordinate",
        "samples_in_counts_missing_from_design",
        "samples_in_design_missing_from_counts",
    )
    return any(any(marker in warning for marker in markers) for warning in warnings)


def _mode_quality_reason(mode: str) -> str:
    if mode == "methylation_by_condition":
        return "valid region_counts/design input supports methylation by condition"
    if mode == "q_by_evidence_class":
        return "valid DMR table supports q-value distribution by evidence class"
    return "valid DMR table supports delta methylation by context"


def _dedupe_strings(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _bar_figure(table: pd.DataFrame, x: str, y: str, title: str, y_label: str | None = None) -> object | None:
    if table.empty or x not in table or y not in table:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(table[x].astype(str), table[y].astype(float), color="#38BDF8")
        ax.set_title(title)
        ax.set_xlabel(x)
        ax.set_ylabel(y_label or y)
        ax.tick_params(axis="x", rotation=45)
        fig.tight_layout()
        return fig
    except Exception:
        return None


def _box_figure(table: pd.DataFrame, group: str, value: str, title: str, y_label: str | None = None) -> object | None:
    if table.empty or group not in table or value not in table:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        groups = [name for name, _ in table.groupby(group, dropna=False)]
        data = [table.loc[table[group] == name, value].dropna().astype(float).to_numpy() for name in groups]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.boxplot(data, labels=[str(item) for item in groups], patch_artist=True)
        ax.set_title(title)
        ax.set_xlabel(group)
        ax.set_ylabel(y_label or value)
        ax.tick_params(axis="x", rotation=45)
        fig.tight_layout()
        return fig
    except Exception:
        return None


def _violin_figure(table: pd.DataFrame, group: str, value: str, title: str, y_label: str | None = None) -> object | None:
    if table.empty or group not in table or value not in table:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        groups = [name for name, _ in table.groupby(group, dropna=False)]
        data = [table.loc[table[group] == name, value].dropna().astype(float).to_numpy() for name in groups]
        data = [values for values in data if len(values) > 0]
        if not data:
            return None
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.violinplot(data, showmeans=True)
        ax.set_xticks(range(1, len(groups) + 1), [str(item) for item in groups])
        ax.set_title(title)
        ax.set_xlabel(group)
        ax.set_ylabel(y_label or value)
        ax.tick_params(axis="x", rotation=45)
        fig.tight_layout()
        return fig
    except Exception:
        return None


def _scatter_figure(
    table: pd.DataFrame,
    x: str,
    y: str,
    color: str | None,
    title: str,
    x_label: str | None = None,
    y_label: str | None = None,
    hline: float | None = None,
    alpha: float = 0.8,
) -> object | None:
    if table.empty or x not in table or y not in table:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 5))
        if color and color in table:
            for label, group_df in table.groupby(color, dropna=False):
                ax.scatter(group_df[x], group_df[y], label=str(label), alpha=alpha, s=20)
            ax.legend(fontsize=8)
        else:
            ax.scatter(table[x], table[y], alpha=alpha, s=20, color="#38BDF8")
        if hline is not None:
            ax.axhline(hline, color="#EF4444", linestyle="--", linewidth=1, alpha=0.8)
        ax.set_title(title)
        ax.set_xlabel(x_label or x)
        ax.set_ylabel(y_label or y)
        fig.tight_layout()
        return fig
    except Exception:
        return None


def _heatmap_figure(matrix: pd.DataFrame, title: str) -> object | None:
    if matrix.empty:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        display = matrix.head(200)
        fig, ax = plt.subplots(figsize=(8, 6))
        image = ax.imshow(display.to_numpy(dtype=float), aspect="auto", interpolation="nearest", cmap="viridis")
        ax.set_title(title)
        ax.set_xlabel("sample")
        ax.set_ylabel("DMR")
        ax.set_xticks(range(display.shape[1]), [str(col) for col in display.columns], rotation=45, ha="right")
        ax.set_yticks([])
        fig.colorbar(image, ax=ax, label="methylation")
        fig.tight_layout()
        return fig
    except Exception:
        return None


def _write_bundle_summary(path: Path, bundle: DmrCurveBundle, results: list[DmrCurveResult]) -> None:
    by_quality: dict[str, list[DmrCurveResult]] = {}
    for result in results:
        by_quality.setdefault(result.quality_status, []).append(result)
    lines = [
        "# DMR Curve Bundle Summary",
        "",
        f"Bundle id: `{bundle.bundle_id}`",
        "",
        "This bundle contains reusable DMR visualization specs and optional summary tables/figures. It does not run DMR calling or replace statistical testing.",
        "",
        "## Thesis-Ready Curves",
        "",
    ]
    for result in by_quality.get("thesis_ready", []):
        lines.append(f"- `{result.spec.curve_id}`: {result.spec.title}.")
    if not by_quality.get("thesis_ready"):
        lines.append("- None.")
    lines.extend(["", "## Exploratory Curves", ""])
    for result in by_quality.get("exploratory", []):
        lines.append(f"- `{result.spec.curve_id}`: {result.spec.title}; reasons: {'; '.join(result.quality_reasons)}")
    if not by_quality.get("exploratory"):
        lines.append("- None.")
    lines.extend(["", "## Warning / Skipped Curves", ""])
    for key in ("warning", "skipped"):
        for result in by_quality.get(key, []):
            lines.append(f"- `{result.spec.curve_id}` ({key}): {'; '.join(result.quality_reasons or result.warnings)}")
    if not by_quality.get("warning") and not by_quality.get("skipped"):
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Curve Manifest",
            "",
            "| curve_id | type | quality | table | figure | warnings |",
            "|---|---|---|---:|---:|---|",
        ]
    )
    for result in results:
        lines.append(
            f"| `{result.spec.curve_id}` | `{result.spec.curve_type}` | `{result.quality_status}` | {result.table is not None} | {result.figure is not None} | {'; '.join(result.warnings)} |"
        )
    lines.extend(
        [
            "",
            "## Input Table Semantic QC",
            "",
            "- DMR table alone supports summary plots such as chromosome distribution, evidence class counts, delta distributions, and volcano plots.",
            "- Region-counts plus design support methylation distributions, PCA, and heatmap plots.",
            "- Positional signal or precomputed `DiscreteRegionData` is required for line/meta-DMR profiles.",
            "- Annotation plots require region-level annotation; enrichment-only tables are shown only when they include explicit numeric counts.",
            "",
            "## Notes For Interpretation",
            "",
            "- Quality status `thesis_ready` means the input schema matches the curve's biological interpretation.",
            "- Quality status `exploratory` means the curve can be useful but should not be described as direct DMR-level evidence without caveats.",
            "- Quality status `skipped` means the required semantic input was absent.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def safe_neg_log10_q(values: pd.Series, cap: float = 50.0) -> pd.DataFrame:
    """Return capped ``-log10(q)`` values plus a cap indicator."""

    numeric = pd.to_numeric(values, errors="coerce")
    transformed = []
    capped = []
    for value in numeric:
        if pd.isna(value):
            transformed.append(np.nan)
            capped.append(False)
        elif value <= 0:
            transformed.append(float(cap))
            capped.append(True)
        else:
            score = -math.log10(float(value))
            transformed.append(float(min(score, cap)))
            capped.append(score > cap)
    return pd.DataFrame({"neg_log10_q_capped": transformed, "q_capped": capped})


def _neg_log10(values: pd.Series) -> pd.Series:
    return safe_neg_log10_q(values, cap=50.0)["neg_log10_q_capped"]


def _figure_type(figure: object | None) -> str:
    if figure is None:
        return "none"
    module = type(figure).__module__.lower()
    if "matplotlib" in module:
        return "matplotlib"
    if "plotly" in module:
        return "plotly"
    if "holoviews" in module:
        return "holoviews"
    return type(figure).__name__


def _input_name(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, pd.DataFrame):
        return "DataFrame"
    return str(value)


def _truthy(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "support", "supported"}
    return bool(value)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if np.isnan(value) else float(value)
    if isinstance(value, float):
        return None if math.isnan(value) else value
    if not isinstance(value, (dict, list, tuple)):
        try:
            missing = pd.isna(value)
            if isinstance(missing, (bool, np.bool_)) and missing:
                return None
        except Exception:
            pass
    return value


__all__ = [
    "DmrCurveBundle",
    "DmrCurveResult",
    "DmrCurveSpec",
    "build_default_dmr_curve_bundle",
    "make_dmr_box_curve",
    "make_dmr_caller_support_curve",
    "make_dmr_annotation_summary_curve",
    "make_dmr_chromosome_curve",
    "make_dmr_evidence_bar_curve",
    "make_dmr_heatmap_curve",
    "make_dmr_line_curve",
    "make_dmr_pca_curve",
    "make_dmr_violin_curve",
    "make_dmr_volcano_curve",
    "safe_neg_log10_q",
    "save_dmr_curve_bundle_outputs",
]
