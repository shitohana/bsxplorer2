from __future__ import annotations

import numpy as np
import pytest

bsx2 = pytest.importorskip("bsx2")
if not hasattr(bsx2, "Context"):
    pytest.skip("bsx2 extension is unavailable", allow_module_level=True)

viz = pytest.importorskip("bsx2.viz")

AnnotProfileLayout = viz.AnnotProfileLayout
AnnotProfilePart = viz.AnnotProfilePart
MetageneProfileSegment = viz.MetageneProfileSegment
box_plot = viz.box_plot
build_annotation_metagene = viz.build_annotation_metagene
build_manual_metagene = viz.build_manual_metagene
collect_layout_parts_from_hcannot = viz.collect_layout_parts_from_hcannot
compute_from_annot = viz.compute_from_annot
heatmap = viz.heatmap
line_plot = viz.line_plot
violin_plot = viz.violin_plot


class FakeSeries:
    def __init__(self, values):
        self._values = np.asarray(values, dtype=float)

    def to_numpy(self):
        return self._values


class FakeBatch:
    def __init__(self, positions, density):
        self._positions = FakeSeries(positions)
        self._density = FakeSeries(density)

    def position(self):
        return self._positions

    def density(self):
        return self._density


class QueryReader:
    def __init__(self, batches):
        self._batches = dict(batches)

    def _key(self, contig):
        return (str(contig.seqname), int(contig.start), int(contig.end))

    def query(self, contig):
        return self._batches.get(self._key(contig))

    def reset(self):
        return None


class FakeAttributes:
    def __init__(self, parent=None):
        self.parent = [] if parent is None else [parent]


class FakeEntry:
    def __init__(self, feature_type, contig, entry_id=None, parent=None):
        self.feature_type = feature_type
        self.contig = contig
        self.id = entry_id
        self.attributes = FakeAttributes(parent=parent)


class FakeAnnot:
    def __init__(self, entries):
        self._entries = list(entries)

    def iter(self):
        return iter(self._entries)

    def __iter__(self):
        return iter(self._entries)


def _is_holoviews_object(obj) -> bool:
    return obj.__class__.__module__.startswith("holoviews")


def _make_contig(seqname: str, start: int, end: int, strand: str = "+"):
    strand_value = bsx2.Strand.Forward if strand == "+" else bsx2.Strand.Reverse
    return bsx2.Contig(seqname, start, end, strand_value)


def _sample_layout() -> AnnotProfileLayout:
    return AnnotProfileLayout(
        (
            AnnotProfilePart("enhancer", 1, source="feature", feature_type="enhancer"),
            AnnotProfilePart("promoter", 1, source="feature", feature_type="promoter"),
            AnnotProfilePart("gene", 1, source="gene", required=True),
            AnnotProfilePart("terminator", 1, source="feature", feature_type="terminator"),
        )
    )


def _sample_annot():
    gene = _make_contig("chr1", 100, 200)
    enhancer = _make_contig("chr1", 20, 40)
    promoter = _make_contig("chr1", 50, 100)
    terminator = _make_contig("chr1", 200, 240)
    gene2 = _make_contig("chr1", 300, 380)
    promoter2 = _make_contig("chr1", 260, 300)

    annot = FakeAnnot(
        [
            FakeEntry("gene", gene, entry_id="gene_1"),
            FakeEntry("enhancer", enhancer, parent="gene_1"),
            FakeEntry("promoter", promoter, parent="gene_1"),
            FakeEntry("terminator", terminator, parent="gene_1"),
            FakeEntry("gene", gene2, entry_id="gene_2"),
            FakeEntry("promoter", promoter2, parent="gene_2"),
        ]
    )

    batches = {
        ("chr1", 20, 40): FakeBatch([30], [0.10]),
        ("chr1", 50, 100): FakeBatch([75], [0.20]),
        ("chr1", 100, 200): FakeBatch([150], [0.30]),
        ("chr1", 200, 240): FakeBatch([220], [0.40]),
        ("chr1", 260, 300): FakeBatch([280], [0.25]),
        ("chr1", 300, 380): FakeBatch([340], [0.35]),
    }
    return annot, QueryReader(batches)


def test_collect_layout_parts_from_hcannot_skips_missing_required_gene_parts() -> None:
    annot, _ = _sample_annot()

    part_map = collect_layout_parts_from_hcannot(
        annot,
        layout=_sample_layout(),
    )

    assert [label for label in part_map["gene"][1]] == ["gene_1", "gene_2"]
    assert [label for label in part_map["terminator"][1]] == ["gene_1"]


def test_layout_composer_is_public_and_legacy_name_is_removed() -> None:
    assert hasattr(viz, "compose_layout_drd")
    assert not hasattr(viz, "combine_parts_drd")
    assert hasattr(viz, "build_annotation_metagene")
    assert hasattr(viz, "build_manual_metagene")


def test_collect_layout_parts_from_hcannot_builds_flanks_for_negative_strand() -> None:
    layout = AnnotProfileLayout(
        (
            AnnotProfilePart("promoter", 1, source="flank5", flank_bp=40),
            AnnotProfilePart("gene", 1, source="gene", required=True),
            AnnotProfilePart("terminator", 1, source="flank3", flank_bp=30),
        )
    )
    gene = _make_contig("chr1", 100, 200, strand="-")
    annot = FakeAnnot([FakeEntry("gene", gene, entry_id="gene_neg")])

    part_map = collect_layout_parts_from_hcannot(annot, layout=layout)

    promoter = part_map["promoter"][0][0]
    body = part_map["gene"][0][0]
    terminator = part_map["terminator"][0][0]
    assert (promoter.start, promoter.end, promoter.strand) == (
        200,
        240,
        bsx2.Strand.Reverse,
    )
    assert (body.start, body.end, body.strand) == (100, 200, bsx2.Strand.Reverse)
    assert (terminator.start, terminator.end, terminator.strand) == (
        70,
        100,
        bsx2.Strand.Reverse,
    )


def test_collect_layout_parts_from_hcannot_skips_orphan_feature_parents() -> None:
    layout = AnnotProfileLayout(
        (AnnotProfilePart("promoter", 1, source="feature", feature_type="promoter"),)
    )
    annot = FakeAnnot(
        [FakeEntry("promoter", _make_contig("chr1", 50, 80), parent="missing_gene")]
    )

    part_map = collect_layout_parts_from_hcannot(annot, layout=layout)

    assert part_map["promoter"] == ([], [])


def test_collect_layout_parts_from_hcannot_merges_duplicate_feature_spans() -> None:
    layout = AnnotProfileLayout(
        (
            AnnotProfilePart("promoter", 1, source="feature", feature_type="promoter"),
            AnnotProfilePart("gene", 1, source="gene", required=True),
        )
    )
    annot = FakeAnnot(
        [
            FakeEntry("gene", _make_contig("chr1", 100, 200), entry_id="gene_1"),
            FakeEntry("promoter", _make_contig("chr1", 40, 60), parent="gene_1"),
            FakeEntry("promoter", _make_contig("chr1", 60, 90), parent="gene_1"),
        ]
    )

    part_map = collect_layout_parts_from_hcannot(annot, layout=layout)

    promoter = part_map["promoter"][0][0]
    assert (promoter.start, promoter.end) == (40, 90)
    assert part_map["promoter"][1] == ["gene_1"]


def test_collect_layout_parts_from_hcannot_drops_inconsistent_optional_spans() -> None:
    layout = AnnotProfileLayout(
        (
            AnnotProfilePart("promoter", 1, source="feature", feature_type="promoter"),
            AnnotProfilePart("gene", 1, source="gene", required=True),
        )
    )
    annot = FakeAnnot(
        [
            FakeEntry("gene", _make_contig("chr1", 100, 200), entry_id="gene_1"),
            FakeEntry("promoter", _make_contig("chr1", 40, 60), parent="gene_1"),
            FakeEntry("promoter", _make_contig("chr2", 60, 90), parent="gene_1"),
        ]
    )

    part_map = collect_layout_parts_from_hcannot(annot, layout=layout)

    assert part_map["promoter"] == ([], [])
    assert part_map["gene"][1] == ["gene_1"]


def test_compute_from_annot_layout_supports_arbitrary_parts() -> None:
    annot, reader = _sample_annot()

    data = compute_from_annot(
        reader,
        annot,
        layout=_sample_layout(),
    )

    assert data.labels == ["gene_1", "gene_2"]
    np.testing.assert_allclose(data.positions[0], np.array([0.125, 0.375, 0.625, 0.875]))
    np.testing.assert_allclose(data.densities[0], np.array([0.10, 0.20, 0.30, 0.40]))
    np.testing.assert_allclose(data.positions[1], np.array([0.375, 0.625]))
    np.testing.assert_allclose(data.densities[1], np.array([0.25, 0.35]))


def test_annotation_and_manual_metagene_builders_match() -> None:
    annot, reader = _sample_annot()
    layout = _sample_layout()

    by_annotation = build_annotation_metagene(
        reader,
        annot,
        layout=layout,
    )
    part_map = collect_layout_parts_from_hcannot(annot, layout=layout)
    by_manual = build_manual_metagene(
        reader,
        part_map=part_map,
        layout=layout,
    )

    assert by_annotation.labels == by_manual.labels == ["gene_1", "gene_2"]
    for left, right in zip(by_annotation.positions, by_manual.positions):
        np.testing.assert_allclose(left, right)
    for left, right in zip(by_annotation.densities, by_manual.densities):
        np.testing.assert_allclose(left, right)


def test_compute_from_annot_layout_accepts_explicit_segments() -> None:
    annot, reader = _sample_annot()

    data = compute_from_annot(
        reader,
        annot,
        layout=_sample_layout(),
        segments=[
            MetageneProfileSegment("enh", 1),
            MetageneProfileSegment("prom", 1),
            MetageneProfileSegment("body", 1),
            MetageneProfileSegment("term", 1),
        ],
    )

    assert data.labels == ["gene_1", "gene_2"]
    np.testing.assert_allclose(data.positions[0], np.array([0.125, 0.375, 0.625, 0.875]))
    np.testing.assert_allclose(data.densities[0], np.array([0.10, 0.20, 0.30, 0.40]))


@pytest.mark.parametrize(
    ("layout", "message"),
    [
        (AnnotProfileLayout(tuple()), "must not be empty"),
        (
            AnnotProfileLayout(
                (
                    AnnotProfilePart("dup", 1, source="gene"),
                    AnnotProfilePart("dup", 1, source="feature", feature_type="promoter"),
                )
            ),
            "names must be unique",
        ),
        (
            AnnotProfileLayout((AnnotProfilePart("bad", 1, source="unknown"),)),
            "source must be one of",
        ),
        (
            AnnotProfileLayout((AnnotProfilePart("promoter", 1, source="flank5"),)),
            "flank_bp must be set",
        ),
    ],
)
def test_compute_from_annot_rejects_invalid_layouts(layout, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        compute_from_annot(QueryReader({}), FakeAnnot([]), layout=layout)


def test_compute_from_annot_layout_smoke_renders_all_plot_types() -> None:
    annot, reader = _sample_annot()
    data = compute_from_annot(reader, annot, layout=_sample_layout())

    line = line_plot(data, name="sample", segments=list(_sample_layout().segments))
    hm = heatmap(data, segments=list(_sample_layout().segments))
    box = box_plot(data, segments=list(_sample_layout().segments), as_percent=False)
    violin = violin_plot(data, segments=list(_sample_layout().segments), as_percent=False)

    assert _is_holoviews_object(line)
    assert _is_holoviews_object(hm)
    assert _is_holoviews_object(box)
    assert _is_holoviews_object(violin)
