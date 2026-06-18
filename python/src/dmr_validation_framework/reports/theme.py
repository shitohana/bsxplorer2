"""Global look-and-feel for all DMR validation report figures.

A single place that styles both rendering stacks so the static matplotlib
bundle and the interactive HoloViews/Plotly report share one consistent visual
identity: fonts, title weight, gridlines, no top/right spines, frameless
legends, high-resolution export.
"""

from __future__ import annotations

MPL_RC = {
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.axisbelow": True,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": "#E5E7EB",
    "grid.linewidth": 0.8,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
}


def apply_matplotlib_theme() -> None:
    """Apply the shared matplotlib rcParams (idempotent)."""
    try:
        import matplotlib

        matplotlib.rcParams.update(MPL_RC)
    except Exception:
        pass


def apply_holoviews_theme() -> None:
    """Apply shared HoloViews defaults. Must run after ``hv.extension(...)``."""
    try:
        import holoviews as hv

        hv.opts.defaults(
            hv.opts.Points(fontscale=1.1),
            hv.opts.Bars(fontscale=1.1),
            hv.opts.BoxWhisker(fontscale=1.1),
            hv.opts.Curve(fontscale=1.1),
            hv.opts.Scatter(fontscale=1.1),
            hv.opts.HeatMap(fontscale=1.0),
        )
    except Exception:
        # Backend-specific option validation can fail; theming is best-effort.
        pass
