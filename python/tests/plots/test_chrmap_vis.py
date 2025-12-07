import numpy as np
import pytest

import holoviews as hv
hv.extension("bokeh", logo=False)

from bsx2.plots.chrmap import ChrLineData, ChrBoxData
from bsx2.plots.chrmap_vis import chr_line_hv, chr_box_hv


def test_chr_line_hv_smoke():
    line = ChrLineData(
        x=np.arange(5),
        y=np.linspace(0, 100, 5),
        x_ticks=[1, 3],
        x_labels=["chr1", "chr2"],
        borders=np.array([0, 3, 5]),
        lower=np.linspace(0, 50, 5),
        upper=np.linspace(50, 100, 5),
    )
    obj = chr_line_hv(line, label="sample")
    assert isinstance(obj, hv.Overlay)
    # ensure children exist
    assert len(obj) >= 1


def test_chr_box_hv_smoke():
    box = ChrBoxData(labels=["chr1", "chr2"], values=[np.array([10.0, 20.0]), np.array([30.0])])
    obj_box = chr_box_hv(box, kind="box")
    assert isinstance(obj_box, hv.BoxWhisker)
    obj_violin = chr_box_hv(box, kind="violin")
    assert isinstance(obj_violin, hv.Violin)
