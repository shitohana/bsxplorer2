import math

import polars as pl
import pytest
from bsx2 import _bsx2
from bsx2.metagene import Metagene
from bsx2.types import (
    AggMethod,
    BsxBatch,
    BsxColumns,
    Context,
    ContextData,
    Contig,
    GenomicPosition,
    GffEntry,
    GffEntryAttributes,
    HcAnnotStore,
    ReportTypeSchema,
    Strand,
)


def test_column_constructor_and_dtype():
    batch = BsxBatch(
        "chr1", None, [10, 20], [True, False], [True, None], [1, 0], [2, 0]
    )
    assert batch.height() == 2
    assert batch.position().to_list() == [10, 20]
    density = batch.density().to_list()
    assert density[0] == 0.5
    assert math.isnan(density[1])  # Preserve the core's zero-coverage convention.
    assert BsxBatch.empty(pl.Categorical).is_empty()
    assert BsxBatch.schema() == BsxColumns.schema()
    with pytest.raises(ValueError, match="categorical"):
        BsxBatch.empty(pl.Int32)
    with pytest.raises(ValueError, match="same length"):
        BsxBatch("chr1", None, [10], [], [], [], [])
    with pytest.raises(ValueError, match="sorted"):
        BsxBatch("chr1", None, [20, 10], [True] * 2, [True] * 2, [1] * 2, [2] * 2)
    with pytest.raises(OverflowError):
        BsxBatch("chr1", None, [-1], [True], [True], [1], [2])
    with pytest.raises(OverflowError):
        BsxBatch("chr1", None, [10], [True], [True], [65536], [65536])


@pytest.mark.parametrize(
    "change",
    [
        {"position": [20, 10, 30]},
        {"position": [10, 10, 30]},
        {"position": [10, None, 30]},
        {"chr": ["chr1", "chr2", "chr1"]},
    ],
)
def test_dataframe_validation_rejects_invalid_batch(sample_batch, change):
    df = sample_batch.data().with_columns(
        [pl.Series(name, values) for name, values in change.items()]
    )
    with pytest.raises(ValueError):
        BsxBatch.from_dataframe(df)


def test_dataframe_missing_columns(sample_batch):
    with pytest.raises(ValueError):
        BsxBatch.from_dataframe(sample_batch.data().drop("count_total"))


def test_batch_partition_bounds_and_known_results(sample_batch):
    positions, densities = sample_batch.partition([1, 3], AggMethod.Mean)
    assert positions == pytest.approx([10 / 21, 1.0])
    assert densities == pytest.approx([0.2, 0.5])
    for breakpoints in [[0], [4], [3, 1], [1, 1], [3, 3]]:
        with pytest.raises((ValueError, RuntimeError)):
            sample_batch.partition(breakpoints, AggMethod.Mean)
    with pytest.raises((ValueError, RuntimeError)):
        sample_batch.discretise(0, AggMethod.Mean)
    positions, densities = sample_batch.normalized()
    assert positions == pytest.approx([0, 10 / 21, 20 / 21])
    assert densities == pytest.approx([0.2, 0.4, 0.6])
    assert BsxBatch.empty().normalized() == ([], [])
    maximal = BsxBatch(
        "chr1", None, [0, 2**32 - 1], [True] * 2, [True] * 2, [1] * 2, [2] * 2
    )
    assert maximal.normalized()[0] == pytest.approx([0, (2**32 - 1) / 2**32])
    with pytest.raises(OverflowError):
        maximal.partition([1], AggMethod.Mean)


def test_batch_segmentation_and_probabilities(sample_batch):
    positions, densities = sample_batch.shrink(min_size=1)
    assert positions[-1] == 1
    assert len(positions) == len(densities)
    assert BsxBatch.empty().shrink(1) == ([], [])
    for kwargs in [
        {"min_size": 0},
        {"min_size": 1, "beta": float("nan")},
        {"min_size": 1, "beta": -1},
    ]:
        with pytest.raises(ValueError):
            sample_batch.shrink(**kwargs)
    for mean, pvalue in [(-1, 0.05), (0.5, float("nan")), (0.5, 2)]:
        with pytest.raises(ValueError):
            sample_batch.as_binom(mean, pvalue)


def test_unchecked_extension_cannot_corrupt_python_batch(sample_batch):
    first = sample_batch.slice(0, 2)
    assert first.extend_unchecked(sample_batch.slice(2, 1)) is None
    assert first.height() == 3
    with pytest.raises(ValueError):
        first.extend_unchecked(sample_batch)
    assert first.height() == 3
    assert first == sample_batch
    assert first != BsxBatch.empty()
    with pytest.raises(NotImplementedError):
        _ = first < sample_batch
    assert "BsxBatch" in repr(first)


def test_context_data_and_batch_alignment():
    context = ContextData(b"CGACG")
    positions, strands, contexts = context.take()
    assert len(positions) == len(strands) == len(contexts) == len(context)
    assert all(isinstance(s, Strand) for s in strands)
    assert all(isinstance(c, Context) for c in contexts)
    batch = BsxBatch("chr1", None, [0], [True], [True], [1], [2])
    aligned = batch.add_context_data(context)
    assert aligned.height() >= batch.height()
    assert context.take()[0] == positions  # take is a copied view, not consumption.
    assert GffEntryAttributes().other == {}


def test_mutable_coordinate_properties_and_comparisons():
    gpos = GenomicPosition("chr1", 10)
    gpos.position = 20
    gpos.seqname = "chr2"
    assert gpos == GenomicPosition("chr2", 20)
    assert gpos > GenomicPosition("chr2", 10)
    assert gpos >= GenomicPosition("chr2", 20)
    assert gpos <= GenomicPosition("chr2", 20)
    with pytest.raises(ValueError):
        _ = gpos < GenomicPosition("chr1", 20)
    with pytest.raises(OverflowError):
        _ = GenomicPosition("chr1", 2**32 - 1) + GenomicPosition("chr1", 1)
    contig = Contig("chr1", 10, 30, Strand.Forward)
    contig.start = 15
    contig.end = 25
    contig.seqname = "chr2"
    contig.strand = Strand.Reverse
    assert contig == Contig("chr2", 15, 25, Strand.Reverse)
    for field, value in [("start", 26), ("end", 14)]:
        with pytest.raises(ValueError):
            setattr(contig, field, value)
    assert contig.length() == 10
    assert not Contig("chr1", 10, 10, Strand.Null).is_empty()
    assert Contig("", 0, 0, Strand.Null).is_empty()
    assert "Contig" in repr(contig)


def test_annotations_and_id_contract(tmp_path):
    gff = tmp_path / "genes.gff"
    gff.write_text(
        "chr1\ttest\tgene\t10\t30\t.\t+\t.\tID=gene1\n"
        "chr1\ttest\texon\t15\t20\t.\t+\t.\tID=exon1;Parent=gene1\n"
    )
    store = HcAnnotStore.from_gff(gff)
    entries = {entry.id: (key, entry) for key, entry in store}
    gene_id = entries["gene1"][0]
    exon_id = entries["exon1"][0]
    store.init_tree()
    store.init_imap()
    assert store.get_parent(exon_id) == gene_id
    assert store.get_children(gene_id) == [exon_id]
    gene = store.get_entry(gene_id)
    assert gene is not None
    assert gene.id == "gene1"
    assert store.get_entries_regex("^gene1$")[0].id == "gene1"
    assert gene_id in store.genomic_query(Contig("chr1", 11, 29, Strand.Null))
    assert store.get_feature_types()["gene"] == [gene_id]
    with pytest.raises(ValueError):
        store.get_entries_regex("[")
    store.add_flanks([gene_id], 5, "downstream_")
    assert store.len() == 3
    store.init_tree()
    children = store.get_children(gene_id)
    assert children is not None
    assert len(children) == 2
    iterator = store.iter()
    assert iter(iterator) is iterator
    assert len(list(iterator)) == 3
    with pytest.raises(StopIteration):
        next(iterator)
    inserted = store.insert(GffEntry(Contig("chr2", 0, 10, Strand.Null), id="new"))
    entry = store.get_entry(inserted)
    assert entry is not None
    assert entry.id == "new"
    with pytest.raises(ValueError, match="coordinate"):
        store.add_flanks([inserted], -20, "upstream_")
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t10\t30\tgene1\t0\t+\n")
    assert HcAnnotStore.from_bed(bed).len() == 1
    with pytest.raises(FileNotFoundError):
        HcAnnotStore.from_gff(tmp_path / "missing")


def test_all_report_schema_members():
    for schema in [
        ReportTypeSchema.Bismark,
        ReportTypeSchema.CgMap,
        ReportTypeSchema.BedGraph,
        ReportTypeSchema.Coverage,
    ]:
        assert list(schema.schema()) == schema.col_names()
        assert schema.chr_col() in schema.col_names()
        assert schema.position_col() in schema.col_names()
        assert isinstance(schema.need_align(), bool)
        assert schema.strand_col() is None or schema.strand_col() in schema.col_names()
        assert (
            schema.context_col() is None or schema.context_col() in schema.col_names()
        )


def test_native_merge_validates_and_preserves_pairs():
    positions, densities = _bsx2.merge_metagene_values(
        [[0.5], [0.1, 0.9]], [[0.7], [0.2, 0.8]]
    )
    assert positions == [0.1, 0.5, 0.9]
    assert densities == [0.2, 0.7, 0.8]
    assert _bsx2.merge_metagene_values([], []) == ([], [])
    for positions, densities in [
        ([[0.1]], []),
        ([[0.1]], [[0.1, 0.2]]),
        ([[float("nan")]], [[0.5]]),
    ]:
        with pytest.raises(ValueError):
            _bsx2.merge_metagene_values(positions, densities)
    assert math.isnan(_bsx2.merge_metagene_values([[0.1]], [[float("nan")]])[1][0])


def test_metagene_protocol_and_empty_data():
    gene = Metagene()
    assert gene.gather() == ([], [])
    gene["one"] = ([0.1, 0.5], [0.2, 0.8])
    assert len(gene) == 1
    assert gene["one"] == ([0.1, 0.5], [0.2, 0.8])
    assert list(gene.densities()) == [0.2, 0.8]
    other = Metagene()
    other["two"] = ([0.9], [0.5])
    gene |= other
    assert set(gene.keys()) == {"one", "two"}
    assert gene.gather()[0] == [0.1, 0.5, 0.9]
    assert gene.remove("two") == ([0.9], [0.5])
    assert gene.get("missing") is None
    with pytest.raises(KeyError):
        _ = gene["missing"]
    with pytest.raises(ValueError):
        gene.insert("bad", [0.1], [])
    with pytest.raises(ValueError):
        gene["bad"] = ([0.1],)
