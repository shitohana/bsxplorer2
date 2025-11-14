import os
import numpy as np

# Рекомендуемая защита для headless-режима (WSL без дисплея)
os.environ.setdefault("MPLBACKEND", "Agg")

# 1) Параметры путей (используем Linux-копию; fallback на /mnt/c если не найдены)
bsx_path = os.path.expanduser('~/bsxplorer2-wsl/bsxplorer2/tests/data/report.bsx')
gff_path = os.path.expanduser('~/bsxplorer2-wsl/bsxplorer2/tests/data/annot.gff')
if not os.path.exists(bsx_path):
    bsx_path = '/mnt/c/Users/sysue/bsx/bsxplorer2/bsxplorer2/tests/data/report.bsx'
if not os.path.exists(gff_path):
    gff_path = '/mnt/c/Users/sysue/bsx/bsxplorer2/bsxplorer2/tests/data/annot.gff'

print("Using bsx:", bsx_path)
print("Using gff:", gff_path)

# 2) Импорт API
from bsx2.io import RegionReader
from bsx2._bsx2 import HcAnnotStore, AggMethod
from bsx2.plots import Segment
from bsx2.plots.polars_html import (
    line_html_from_annot, heatmap_html_from_annot, box_html_from_annot, violin_html_from_annot
)
from bsx2.plots import line_plot, heatmap, box_plot, violin_plot  # HoloViews

# 3) Инициализация источников
rr = RegionReader(bsx_path)
annot = HcAnnotStore.from_gff(gff_path)

# 4) Настройка фильтров (опционально; важно для чистоты данных)
try:
    rr.filter_coverage_gt(5)
except (AttributeError, TypeError, ValueError):
    pass

# 5) Задаём произвольные сегменты (не привязаны к promoter/body/terminator)
segs = [Segment("up", 25), Segment("body", 50), Segment("down", 25)]

# 6) Генерация интерактивных HTML (Plotly) — standalone
print("Generating Plotly HTML ...")
html_line = line_html_from_annot(
    rr, annot, segments=segs, agg="mean",
    agg_method=AggMethod.Mean, feature_type=None, limit=100  # ограничим 100 регионов для скорости
)
html_heat = heatmap_html_from_annot(rr, annot, segments=segs, agg_method=AggMethod.Mean, feature_type=None, limit=100)
html_box  = box_html_from_annot(rr, annot, segments=segs, agg_method=AggMethod.Mean, feature_type=None, limit=100)
html_violin = violin_html_from_annot(rr, annot, segments=segs, agg_method=AggMethod.Mean, feature_type=None, limit=100)

# 7) Сохранение HTML
open("line.html", "w", encoding="utf-8").write(html_line)
open("heatmap.html", "w", encoding="utf-8").write(html_heat)
open("box.html", "w", encoding="utf-8").write(html_box)
open("violin.html", "w", encoding="utf-8").write(html_violin)
print("Saved: line.html, heatmap.html, box.html, violin.html")

# 8) HoloViews smoke (типы объектов), без сохранения
print("HoloViews smoke ...")
from bsx2.plots.metagene import compute_from_annot

# --- ВАЖНО: Инициализация backend ДО создания hv-объектов ---
import holoviews as hv
# быстрый статичный вариант (ничего ставить не нужно)
hv.extension('matplotlib')

# Если хочешь интерактивный HTML, раскомментируй две строки ниже и поставь bokeh:
#   pip install bokeh
# hv.extension('bokeh')

# Вычислим список контуров
contigs, _ = [], []
try:
    from bsx2.plots.metagene import collect_contigs_from_hcannot
    contigs, _ = collect_contigs_from_hcannot(annot, feature_type=None, limit=20)
except (AttributeError, TypeError, ValueError, RuntimeError) as e:
    print("collect_contigs_from_hcannot failed:", e)

if contigs:
    hv_curve = line_plot(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
    hv_hm = heatmap(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
    hv_box = box_plot(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
    hv_violin = violin_plot(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)

    if not isinstance(hv_curve, hv.Curve):
        raise TypeError("line_plot must return a holoviews Curve")
    if not isinstance(hv_hm, hv.HeatMap):
        raise TypeError("heatmap must return a holoviews HeatMap")
    if not isinstance(hv_box, hv.BoxWhisker):
        raise TypeError("box_plot must return a holoviews BoxWhisker")
    if not isinstance(hv_violin, hv.Violin):
        raise TypeError("violin_plot must return a holoviews Violin")
    print("HoloViews types OK")
else:
    print("No contigs from annot; skipping HoloViews smoke")
