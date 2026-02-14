from __future__ import annotations

# Compatibility module: keep legacy import path `bsx2.plots.polars_html`.

from .box import box_html, box_html_from_annot
from .heatmap import heatmap_html, heatmap_html_from_annot
from .line import line_html, line_html_from_annot
from .violin import violin_html, violin_html_from_annot

__all__ = [
    "line_html",
    "heatmap_html",
    "box_html",
    "violin_html",
    "line_html_from_annot",
    "heatmap_html_from_annot",
    "box_html_from_annot",
    "violin_html_from_annot",
]
