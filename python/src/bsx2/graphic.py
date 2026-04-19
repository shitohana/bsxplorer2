from __future__ import annotations

from pathlib import Path

import holoviews as hv
import plotly.io as pio

from bsx2 import Context, HcAnnotStore, RegionReader
from bsx2.viz import (
    AnnotProfileLayout,
    AnnotProfilePart,
    compute_from_annot,
    heatmap,
    line_plot,
)

repo = Path.cwd()

bsx_path = repo / "bsxplorer2" / "tests" / "data" / "report.bsx"
gff_path = repo / "bsxplorer2" / "tests" / "data" / "annot.gff"

if not bsx_path.exists():
    raise FileNotFoundError(bsx_path)
if not gff_path.exists():
    raise FileNotFoundError(gff_path)

reader = RegionReader(str(bsx_path))
reader.clear_filters()
reader.filter_context(Context.CG)

annot = HcAnnotStore.from_gff(str(gff_path))

layout = AnnotProfileLayout(
    (
        AnnotProfilePart("up", 25, source="flank5", flank_bp=2000),
        AnnotProfilePart("body", 50, source="gene"),
        AnnotProfilePart("down", 25, source="flank3", flank_bp=2000),
    )
)
segments = list(layout.segments)

limit = None

drd = compute_from_annot(
    reader,
    annot,
    layout=layout,
    segments=segments,
    limit=limit,
)

print("profiles built:", len(drd.positions))

hv.extension("plotly")

line_hv = line_plot(
    drd,
    name="CG",
    segments=segments,
    smooth=None,
    title=f"All genes up/body/down (CG), n={len(drd.positions)}",
    width=1000,
    height=450,
)
heatmap_hv = heatmap(
    drd,
    segments=segments,
    rank_rows=200,
    title=f"All genes up/body/down heatmap (CG), n={len(drd.positions)}",
    width=1000,
    height=800,
)

line_fig = hv.render(line_hv, backend="plotly")
heatmap_fig = hv.render(heatmap_hv, backend="plotly")

output_dir = repo / "python" / "examples"
output_dir.mkdir(parents=True, exist_ok=True)
line_output = output_dir / "_tmp_metagene_profile_line_graphic_match.html"
heatmap_output = output_dir / "_tmp_metagene_profile_heatmap_graphic_match.html"

pio.write_html(line_fig, line_output, include_plotlyjs="cdn", full_html=True)
pio.write_html(heatmap_fig, heatmap_output, include_plotlyjs="cdn", full_html=True)

print("written:")
print(f"  {line_output}")
print(f"  {heatmap_output}")
