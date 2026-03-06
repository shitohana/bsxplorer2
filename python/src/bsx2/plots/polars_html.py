from __future__ import annotations

# Compatibility module: keep legacy import path `bsx2.plots.polars_html`.

from .box import box_html
from .heatmap import heatmap_html
from .line import line_html
from .violin import violin_html

__all__ = [
    "line_html",
    "heatmap_html",
    "box_html",
    "violin_html",
]
