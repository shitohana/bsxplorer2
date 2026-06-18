"""Genome-wide DMR distribution (Manhattan-style) shared layout + renderer.

Lays DMRs out along concatenated chromosomes and plots significance
(``-log10 q``, or ``|delta|`` when q is unavailable) at each DMR midpoint,
colored by methylation context. Used by the matplotlib bundle and the
interactive report.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd


def _natural_key(chrom: str):
    parts = re.split(r"(\d+)", str(chrom))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def compute_manhattan_layout(rows: pd.DataFrame) -> dict:
    """Return concatenated-genome layout for DMR rows.

    Result keys: ``data`` (with ``_gpos``, ``_y``), ``ticks`` (list of
    ``(position, chrom)``), ``y_label``, ``shading`` (list of ``(x0, x1)`` per
    alternate chromosome), or an empty dict if the input is unusable.
    """
    if rows.empty or not {"chrom", "start", "end"}.issubset(rows.columns):
        return {}
    df = rows.copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    if df.empty:
        return {}
    df["_mid"] = (df["start"] + df["end"]) / 2.0

    q = pd.to_numeric(df.get("q_value", pd.Series(index=df.index, dtype=float)), errors="coerce")
    if q.notna().any():
        df["_y"] = -np.log10(q.clip(lower=1e-300))
        y_label = "-log10 q"
    else:
        df["_y"] = pd.to_numeric(df.get("delta", pd.Series(index=df.index, dtype=float)), errors="coerce").abs()
        y_label = "|delta|"
    df = df.dropna(subset=["_y"])
    if df.empty:
        return {}

    chroms = sorted(df["chrom"].unique(), key=_natural_key)
    gap = 0.0
    offset = 0.0
    offsets: dict[str, float] = {}
    ticks: list[tuple[float, str]] = []
    shading: list[tuple[float, float]] = []
    for i, chrom in enumerate(chroms):
        chrom_len = float(df.loc[df["chrom"] == chrom, "end"].max())
        offsets[chrom] = offset
        ticks.append((offset + chrom_len / 2.0, chrom))
        if i % 2 == 1:
            shading.append((offset, offset + chrom_len))
        offset += chrom_len + gap
    df["_gpos"] = df.apply(lambda r: offsets[r["chrom"]] + r["_mid"], axis=1)
    return {"data": df, "ticks": ticks, "y_label": y_label, "shading": shading, "total": offset}


def build_manhattan_figure(rows: pd.DataFrame, *, title: str = "Genome-wide DMR distribution"):
    """Build (do not save) a matplotlib Manhattan figure colored by context."""
    layout = compute_manhattan_layout(rows)
    if not layout:
        return None
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from dmr_validation_framework.reports.palette import CONTEXT_COLORS, color_for, legend_outside

    df = layout["data"]
    ctx_col = "_context" if "_context" in df.columns else "context" if "context" in df.columns else None
    fig, ax = plt.subplots(figsize=(12, 4.6))
    for x0, x1 in layout["shading"]:
        ax.axvspan(x0, x1, color="#F3F4F6", zorder=0)
    if ctx_col:
        for ctx, sub in df.groupby(ctx_col):
            ax.scatter(sub["_gpos"], sub["_y"], s=10, alpha=0.7, color=color_for(ctx, CONTEXT_COLORS), label=str(ctx))
        legend_outside(ax, title="context")
    else:
        ax.scatter(df["_gpos"], df["_y"], s=10, alpha=0.7)
    ax.set_xticks([pos for pos, _ in layout["ticks"]])
    ax.set_xticklabels([label for _, label in layout["ticks"]], rotation=60, ha="right", fontsize=8)
    ax.set_ylabel(layout["y_label"])
    ax.set_xlabel("chromosome")
    ax.set_title(title)
    fig.tight_layout()
    return fig
