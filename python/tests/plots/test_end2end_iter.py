import os
import pytest


@pytest.mark.end2end
def test_end2end_iter_contigs_html_and_hv_smoke() -> None:
    """End-to-end on repo test data using RegionReader.iter_contigs.

    Skips if repo test data are not present. Verifies:
    - HTML generation from annot (Plotly)
    - HoloViews object types for line/heat/box/violin
    """
    # Resolve test data paths in common locations
    candidates_bsx = [
        os.path.expanduser("~/bsxplorer2-wsl/bsxplorer2/tests/data/report.bsx"),
        "/mnt/c/Users/sysue/bsx/bsxplorer2/bsxplorer2/tests/data/report.bsx",
        os.path.join(os.getcwd(), "..", "bsxplorer2", "tests", "data", "report.bsx"),
    ]
    candidates_gff = [
        os.path.expanduser("~/bsxplorer2-wsl/bsxplorer2/tests/data/annot.gff"),
        "/mnt/c/Users/sysue/bsx/bsxplorer2/bsxplorer2/tests/data/annot.gff",
        os.path.join(os.getcwd(), "..", "bsxplorer2", "tests", "data", "annot.gff"),
    ]

    bsx_path = next((p for p in candidates_bsx if os.path.exists(p)), None)
    gff_path = next((p for p in candidates_gff if os.path.exists(p)), None)
    if not bsx_path or not gff_path:
        pytest.skip("Test data (.bsx/.gff) not found; skipping end-to-end iter_contigs test")

    # Imports (within test to allow skip above)
    import holoviews as hv
    hv.extension("matplotlib")  # choose available backend for smoke

    from bsx2.io import RegionReader
    from bsx2._bsx2 import HcAnnotStore, AggMethod
    from bsx2.plots import (
        Segment,
        line_plot,
        heatmap,
        box_plot,
        violin_plot,
    )
    from bsx2.plots.metagene import collect_contigs_from_hcannot
    from bsx2.plots.polars_html import (
        line_html_from_annot,
        heatmap_html_from_annot,
        box_html_from_annot,
        violin_html_from_annot,
    )

    rr = RegionReader(bsx_path)
    annot = HcAnnotStore.from_gff(gff_path)

    # Collect a limited set of contigs for fast test
    contigs, _labels = collect_contigs_from_hcannot(annot, feature_type=None, limit=20)
    assert contigs, "No contigs collected from annot for end-to-end test"

    segs = [Segment("up", 10), Segment("body", 20), Segment("down", 10)]

    # HTML generation from annot wrappers (Plotly)
    html_line = line_html_from_annot(
        rr,
        annot,
        segments=segs,
        agg="mean",
        agg_method=AggMethod.Mean,
        feature_type=None,
        limit=20,
        full_html=False,
    )
    html_heat = heatmap_html_from_annot(rr, annot, segments=segs, agg_method=AggMethod.Mean, feature_type=None, limit=20)
    html_box = box_html_from_annot(rr, annot, segments=segs, agg_method=AggMethod.Mean, feature_type=None, limit=20)
    html_violin = violin_html_from_annot(rr, annot, segments=segs, agg_method=AggMethod.Mean, feature_type=None, limit=20)

    for doc in (html_line, html_heat, html_box, html_violin):
        assert isinstance(doc, str) and ("<div" in doc or "plotly" in doc.lower())

    # HoloViews smoke (types)
    hv_curve = line_plot(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
    hv_hm = heatmap(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
    hv_box = box_plot(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)
    hv_violin = violin_plot(rr, contigs=contigs, segments=segs, agg_method=AggMethod.Mean)

    assert isinstance(hv_curve, hv.Curve)
    assert isinstance(hv_hm, hv.HeatMap)
    assert isinstance(hv_box, hv.BoxWhisker)
    assert isinstance(hv_violin, hv.Violin)

