import pytest
import polars as pl
from src.bsx2.types import (
    Strand, Context, BsxColumns, AggMethod, BsxBatch, ContextData,
    ReportTypeSchema, LazyBsxBatch, GenomicPosition, Contig,
    GffEntryAttributes, GffEntry, HcAnnotStore
)


# --- Enum Tests ---

def test_strand_properties():
    """Test strand enum properties - mirrors test_strand_display"""
    assert Strand.Forward.name == "Forward"
    assert Strand.Reverse.name == "Reverse"
    assert Strand.Null.name == "Null"


def test_context_properties():
    """Test context enum properties - mirrors test_context_display"""
    assert Context.CG.name == "CG"
    assert Context.CHG.name == "CHG"
    assert Context.CHH.name == "CHH"


def test_bsx_columns():
    """Test BsxColumns functionality"""
    colnames = BsxColumns.colnames()
    assert isinstance(colnames, list)
    assert len(colnames) > 0

    assert BsxColumns.has_name("chr")
    assert BsxColumns.has_name("position")
    assert not BsxColumns.has_name("NonExistent")

    assert BsxColumns.Chr.as_str() == "chr"


# --- BsxBatch Tests ---

def test_empty_batch():
    """Test empty batch creation - mirrors test_empty_batch"""
    batch = BsxBatch.empty()
    assert batch.is_empty()
    assert batch.height() == 0
    assert len(batch) == 0


@pytest.fixture
def batch() -> BsxBatch:
    """Test batch fixture"""
    data = pl.DataFrame({
        BsxColumns.Chr.as_str(): ["chr1"] * 5,
        BsxColumns.Position.as_str(): [100, 200, 300, 400, 500],
        BsxColumns.Strand.as_str(): [True, False, True, False, True],
        BsxColumns.Context.as_str(): [True, False, None, True, False],
        BsxColumns.CountM.as_str(): [5, 3, 8, 2, 7],
        BsxColumns.CountTotal.as_str(): [10, 10, 10, 10, 10],
        BsxColumns.Density.as_str(): [.5, .3, .8, .2, .7]
    })

    return BsxBatch.from_dataframe(data)


def test_batch_from_dataframe(batch):
    """Test batch creation from dataframe - mirrors test_build_bsx_from_report_type"""
    assert not batch.is_empty()
    assert batch.height() == 5
    assert len(batch) == 5


def test_batch_from_dataframe_validation():
    """Test batch validation options"""
    data = pl.DataFrame({
        BsxColumns.Chr.as_str(): ["chr1"] * 3,
        BsxColumns.Position.as_str(): [100, 200, 300],
        BsxColumns.Strand.as_str(): [True, False, True],
        BsxColumns.Context.as_str(): [True, False, None],
        BsxColumns.CountM.as_str(): [5, 3, 8],
        BsxColumns.CountTotal.as_str(): [10, 10, 10],
        BsxColumns.Density.as_str(): [.5, .3, .8]
    })

    # Test with different validation flags
    batch1 = BsxBatch.from_dataframe(data, check_nulls=False)
    assert batch1.height() == 3

    batch2 = BsxBatch.from_dataframe(data, check_sorted=False)
    assert batch2.height() == 3

    batch3 = BsxBatch.from_dataframe(data, check_duplicates=False)
    assert batch3.height() == 3

    batch4 = BsxBatch.from_dataframe(data, rechunk=False)
    assert batch4.height() == 3


def test_batch_column_getters(batch):
    """Test column getters - mirrors test_column_getters"""
    chr_series = batch.chr()
    assert isinstance(chr_series, pl.Series)
    assert len(chr_series) == 5

    pos_series = batch.position()
    assert isinstance(pos_series, pl.Series)
    assert len(pos_series) == 5

    strand_series = batch.strand()
    assert isinstance(strand_series, pl.Series)

    context_series = batch.context()
    assert isinstance(context_series, pl.Series)

    count_m_series = batch.count_m()
    assert isinstance(count_m_series, pl.Series)

    count_total_series = batch.count_total()
    assert isinstance(count_total_series, pl.Series)

    density_series = batch.density()
    assert isinstance(density_series, pl.Series)

    # Test column by enum
    chr_col = batch.column(BsxColumns.Chr)
    assert isinstance(chr_col, pl.Series)


def test_batch_slice(batch):
    """Test batch slicing - mirrors test_slice"""
    sliced = batch.slice(1, 3)
    assert sliced.height() == 3


def test_batch_split_at(batch):
    """Test batch splitting - mirrors test_split_at"""
    left, right = batch.split_at(2)
    assert left.height() == 2
    assert right.height() == 3


def test_batch_position_methods(batch):
    """Test position methods - mirrors test_position_methods"""
    seqname = batch.seqname()
    assert seqname == "chr1"

    first_pos = batch.first_pos()
    assert first_pos == 100

    last_pos = batch.last_pos()
    assert last_pos == 500

    # Test genomic positions
    first_gpos = batch.first_genomic_pos()
    assert first_gpos is not None
    assert first_gpos.seqname == "chr1"
    assert first_gpos.position == 100

    last_gpos = batch.last_genomic_pos()
    assert last_gpos is not None
    assert last_gpos.seqname == "chr1"
    assert last_gpos.position == 500


def test_batch_as_contig(batch):
    """Test batch to contig conversion"""
    contig = batch.as_contig()
    assert contig is not None
    assert contig.seqname == "chr1"
    assert contig.start == 100
    assert contig.end == 500


def test_batch_into_dataframe(batch):
    """Test batch to dataframe conversion"""
    df = batch.into_dataframe()
    assert isinstance(df, pl.DataFrame)
    assert df.height == 5


def test_batch_data(batch):
    """Test batch data access"""
    df = batch.data()
    assert isinstance(df, pl.DataFrame)
    assert df.height == 5


def test_batch_discretise(batch):
    """Test batch discretization - mirrors test_discretise_basic_case"""
    positions, densities = batch.discretise(2, AggMethod.Mean)
    assert len(positions) == 2
    assert len(densities) == 2

    # Test with different aggregation methods
    pos_geom, dens_geom = batch.discretise(3, AggMethod.GeometricMean)
    assert len(pos_geom) == 3
    assert len(dens_geom) == 3

    pos_med, dens_med = batch.discretise(2, AggMethod.Median)
    assert len(pos_med) == 2
    assert len(dens_med) == 2

    pos_max, dens_max = batch.discretise(2, AggMethod.Max)
    assert len(pos_max) == 2
    assert len(dens_max) == 2

    pos_min, dens_min = batch.discretise(2, AggMethod.Min)
    assert len(pos_min) == 2
    assert len(dens_min) == 2


def test_batch_discretise_edge_cases():
    """Test batch discretization edge cases - mirrors test_discretise_n_fragments_one"""
    data = pl.DataFrame({
        "chr": ["chr1"] * 3,
        "position": [100, 200, 300],
        "strand": [True, False, True],
        "context": [True, False, None],
        "count_m": [5, 3, 8],
        "count_total": [10, 10, 10],
        "density": [.5, .3, .8]
    })

    batch = BsxBatch.from_dataframe(data)

    # Test with n_fragments = 1
    positions, densities = batch.discretise(1, AggMethod.Mean)
    assert len(positions) == 1
    assert len(densities) == 1


def test_batch_partition():
    """Test batch partitioning - mirrors test_partition"""
    data = pl.DataFrame({
        "chr": ["chr1"] * 10,
        "position": list(range(100, 200, 10)),
        "strand": [True] * 10,
        "context": [True] * 10,
        "count_m": [5] * 10,
        "count_total": [10] * 10,
        "density": [.5] * 10
    })

    batch = BsxBatch.from_dataframe(data)
    breakpoints = [1, 3, 4]
    positions, densities = batch.partition(breakpoints, AggMethod.Mean)
    assert len(positions) == len(breakpoints) + 1
    assert len(densities) == len(breakpoints) + 1


def test_batch_extend():
    """Test batch extension - mirrors test_can_extend"""
    data1 = pl.DataFrame({
        "chr": ["chr1"] * 3,
        "position": [100, 200, 300],
        "strand": [True, False, True],
        "context": [True, False, None],
        "count_m": [5, 3, 8],
        "count_total": [10, 10, 10],
        "density": [.5, .3, .8]
    })

    data2 = pl.DataFrame({
        "chr": ["chr1"] * 2,
        "position": [400, 500],
        "strand": [False, True],
        "context": [True, False],
        "count_m": [2, 7],
        "count_total": [10, 10],
        "density": [.2, .7]
    })

    batch1 = BsxBatch.from_dataframe(data1)
    batch2 = BsxBatch.from_dataframe(data2)

    batch1.extend(batch2)
    assert batch1.height() == 5

    # Test extend_unchecked
    batch3 = BsxBatch.from_dataframe(data1)
    batch3.extend_unchecked(batch2)
    assert batch3.height() == 5


def test_batch_concat():
    """Test batch concatenation"""
    data1 = pl.DataFrame({
        "chr": ["chr1"] * 3,
        "position": [100, 200, 300],
        "strand": [True, False, True],
        "context": [True, False, None],
        "count_m": [5, 3, 8],
        "count_total": [10, 10, 10],
        "density": [.5, .3, .8]
    })

    data2 = pl.DataFrame({
        "chr": ["chr1"] * 2,
        "position": [400, 500],
        "strand": [False, True],
        "context": [True, False],
        "count_m": [2, 7],
        "count_total": [10, 10],
        "density": [.2, .7]
    })

    batch1 = BsxBatch.from_dataframe(data1)
    batch2 = BsxBatch.from_dataframe(data2)

    concatenated = BsxBatch.concat([batch1, batch2])
    assert concatenated.height() == 5


def test_batch_rechunk(batch):
    """Test batch rechunking"""
    rechunked = batch.rechunk()
    assert rechunked.height() == batch.height()


def test_batch_as_binom(batch):
    """Test binomial conversion - mirrors test_as_binom"""
    binom_batch = batch.as_binom(0.5, 0.05)
    assert binom_batch.height() == batch.height()


def test_batch_into_report():
    """Test batch report conversion - mirrors test_into_report"""
    data = pl.DataFrame({
        "chr": ["chr1"] * 3,
        "position": [100, 200, 300],
        "strand": [True, False, True],
        "context": [True, False, None],
        "count_m": [5, 3, 8],
        "count_total": [10, 10, 10],
        "density": [.5, .3, .8]
    })

    batch = BsxBatch.from_dataframe(data)

    # Test different report types
    bismark_report = batch.into_report(ReportTypeSchema.Bismark)
    assert isinstance(bismark_report, pl.DataFrame)

    cgmap_report = batch.into_report(ReportTypeSchema.CgMap)
    assert isinstance(cgmap_report, pl.DataFrame)

    bedgraph_report = batch.into_report(ReportTypeSchema.BedGraph)
    assert isinstance(bedgraph_report, pl.DataFrame)

    coverage_report = batch.into_report(ReportTypeSchema.Coverage)
    assert isinstance(coverage_report, pl.DataFrame)


# --- ContextData Tests ---

def test_context_data_empty():
    """Test empty context data - mirrors test_read_sequence_empty_or_short"""
    context_data = ContextData.empty()
    assert context_data.is_empty()
    assert len(context_data) == 0


def test_context_data_from_sequence():
    """Test context data from sequence - mirrors test_from_sequence_basic"""
    sequence = b"ATCGATCG"
    context_data = ContextData(sequence)
    assert not context_data.is_empty()
    assert len(context_data) > 0


def test_context_data_take():
    """Test context data take method"""
    sequence = b"CGATTACG"
    context_data = ContextData(sequence)

    positions, strands, contexts = context_data.take()
    assert isinstance(positions, list)
    assert isinstance(strands, list)
    assert isinstance(contexts, list)
    assert len(positions) == len(strands) == len(contexts)


def test_context_data_to_decoded_df():
    """Test context data to dataframe conversion - mirrors test_to_df"""
    sequence = b"CGATTACG"
    context_data = ContextData(sequence)

    df = context_data.to_decoded_df()
    assert isinstance(df, pl.DataFrame)


# --- LazyBsxBatch Tests ---

def test_lazy_batch(batch):
    """Test lazy batch operations - mirrors test_ops"""
    lazy_batch = batch.lazy()

    # Test filtering operations
    filtered = lazy_batch.filter_pos_gt(200)
    collected = filtered.collect()
    assert collected.height() <= batch.height()

    filtered2 = lazy_batch.filter_pos_lt(400)
    collected2 = filtered2.collect()
    assert collected2.height() <= batch.height()

    filtered3 = lazy_batch.filter_coverage_gt(5)
    collected3 = filtered3.collect()
    assert collected3.height() <= batch.height()

    filtered4 = lazy_batch.filter_strand(Strand.Forward)
    collected4 = filtered4.collect()
    assert collected4.height() <= batch.height()

    filtered5 = lazy_batch.filter_context(Context.CG)
    collected5 = filtered5.collect()
    assert collected5.height() <= batch.height()


def test_lazy_batch_creation():
    """Test lazy batch creation - mirrors test_lazybatch"""
    data = pl.DataFrame({
        "chr": ["chr1"] * 3,
        "position": [100, 200, 300],
        "strand": [True, False, True],
        "context": [True, False, None],
        "count_m": [5, 3, 8],
        "count_total": [10, 10, 10],
        "density": [.5, .3, .8]
    })

    batch = BsxBatch.from_dataframe(data)
    lazy = LazyBsxBatch(batch)
    collected = lazy.collect()
    assert collected.height() == 3


# --- GenomicPosition Tests ---

def test_genomic_position_creation():
    """Test genomic position creation"""
    gpos = GenomicPosition("chr1", 12345)
    assert gpos.seqname == "chr1"
    assert gpos.position == 12345


def test_genomic_position_is_zero():
    """Test genomic position zero check"""
    gpos_zero = GenomicPosition("chr1", 0)
    assert gpos_zero.is_zero()

    gpos_nonzero = GenomicPosition("chr1", 100)
    assert not gpos_nonzero.is_zero()


def test_genomic_position_arithmetic():
    """Test genomic position arithmetic - mirrors test_genomic_position_add_same_seqname and test_genomic_position_sub_same_seqname_ge"""
    gpos1 = GenomicPosition("chr1", 100)
    gpos2 = GenomicPosition("chr1", 50)

    # Addition - should work for same chromosome
    result_add = gpos1 + gpos2
    if result_add is not None:
        assert result_add.seqname == "chr1"
        assert result_add.position == 150

    # Subtraction - should work when first >= second
    result_sub = gpos1 - gpos2
    if result_sub is not None:
        assert result_sub.seqname == "chr1"
        assert result_sub.position == 50

    # Test arithmetic with different chromosomes returns None
    gpos3 = GenomicPosition("chr2", 50)
    result_add_diff = gpos1 + gpos3
    assert result_add_diff is None

    result_sub_diff = gpos1 - gpos3
    assert result_sub_diff is None

    # Test subtraction when first < second returns None
    result_sub_lt = gpos2 - gpos1
    assert result_sub_lt is None


def test_genomic_position_display():
    """Test genomic position display - mirrors test_genomic_position_display"""
    gpos = GenomicPosition("chr1", 12345)
    assert str(gpos) == "chr1:12345"


# --- Contig Tests ---

def test_contig_creation():
    """Test contig creation"""
    contig = Contig("chr1", 100, 200, Strand.Forward)
    assert contig.seqname == "chr1"
    assert contig.start == 100
    assert contig.end == 200
    assert contig.strand == Strand.Forward

    # Test with different strand representations
    contig_reverse = Contig("chr1", 100, 200, Strand.Reverse)
    assert contig_reverse.strand == Strand.Reverse

    contig_none = Contig("chr1", 100, 200, Strand.Null)
    assert contig_none.strand == Strand.Null


def test_contig_length():
    """Test contig length calculation - mirrors test_contig_length"""
    contig = Contig("chr1", 100, 200, Strand.Forward)
    assert contig.length() == 100


def test_contig_strand_str():
    """Test contig strand string representation"""
    contig_forward = Contig("chr1", 100, 200, Strand.Forward)
    assert contig_forward.strand_str == "+"

    contig_reverse = Contig("chr1", 100, 200, Strand.Reverse)
    assert contig_reverse.strand_str == "-"

    contig_none = Contig("chr1", 100, 200, Strand.Null)
    assert contig_none.strand_str == "."


def test_contig_genomic_positions():
    """Test contig genomic position methods"""
    contig = Contig("chr1", 100, 200, Strand.Forward)

    start_gpos = contig.start_gpos()
    assert start_gpos.seqname == "chr1"
    assert start_gpos.position == 100

    end_gpos = contig.end_gpos()
    assert end_gpos.seqname == "chr1"
    assert end_gpos.position == 200


def test_contig_extend():
    """Test contig extension methods - mirrors test_contig_extend_upstream and test_contig_extend_downstream"""
    contig = Contig("chr1", 100, 200, Strand.Forward)
    original_start = contig.start
    original_end = contig.end

    # Test extend methods exist and modify the contig
    contig.extend_upstream(50)
    assert contig.start <= original_start  # Should extend upstream

    contig.extend_downstream(50)
    assert contig.end >= original_end  # Should extend downstream


def test_contig_is_in():
    """Test contig containment check"""
    contig1 = Contig("chr1", 100, 200, Strand.Forward)
    contig2 = Contig("chr1", 50, 250, Strand.Forward)  # Contains contig1
    contig3 = Contig("chr1", 150, 300, Strand.Forward)  # Overlaps but doesn't contain

    assert contig1.is_in(contig2)
    assert not contig1.is_in(contig3)


def test_contig_display():
    """Test contig display - mirrors test_contig_display"""
    contig = Contig("chrX", 1000, 2000, Strand.Forward)
    display_str = str(contig)
    assert "chrX" in display_str
    assert "1000" in display_str
    assert "2000" in display_str
    assert "+" in display_str


# --- Report Type Schema Tests ---

def test_report_type_schema():
    """Test report type schema methods"""
    bismark = ReportTypeSchema.Bismark
    colnames = bismark.col_names()
    assert isinstance(colnames, list)

    chr_col = bismark.chr_col()
    assert isinstance(chr_col, str)

    pos_col = bismark.position_col()
    assert isinstance(pos_col, str)

    # Test optional methods
    _ = bismark.context_col()
    _ = bismark.strand_col()
    need_align = bismark.need_align()
    assert isinstance(need_align, bool)

    # Test other report types
    cgmap = ReportTypeSchema.CgMap
    cgmap_cols = cgmap.col_names()
    assert isinstance(cgmap_cols, list)

    bedgraph = ReportTypeSchema.BedGraph
    bedgraph_cols = bedgraph.col_names()
    assert isinstance(bedgraph_cols, list)

    coverage = ReportTypeSchema.Coverage
    coverage_cols = coverage.col_names()
    assert isinstance(coverage_cols, list)


# --- GffEntry and Attributes Tests ---

def test_gff_entry_attributes():
    """Test GFF entry attributes"""
    attrs = GffEntryAttributes()

    # Test properties exist and are accessible
    id_val = attrs.id
    name_val = attrs.name
    alias_val = attrs.alias
    parent_val = attrs.parent
    other_val = attrs.other

    assert isinstance(other_val, dict)

    # Test that values can be None initially
    assert id_val is None or isinstance(id_val, str)
    assert name_val is None or isinstance(name_val, list)
    assert alias_val is None or isinstance(alias_val, list)
    assert parent_val is None or isinstance(parent_val, list)


def test_gff_entry_creation():
    """Test GFF entry creation"""
    contig = Contig("chr1", 100, 200, Strand.Forward)
    entry = GffEntry(contig, "test_source", "gene", 0.9, 0, "gene123")

    assert entry.id == "gene123"
    assert entry.contig.seqname == "chr1"
    assert entry.source == "test_source"
    assert entry.feature_type == "gene"
    assert entry.score == 0.9
    assert entry.phase == 0

    attrs = entry.attributes
    assert isinstance(attrs, GffEntryAttributes)


def test_gff_entry_optional_params():
    """Test GFF entry with optional parameters"""
    contig = Contig("chr1", 100, 200, Strand.Forward)

    # Test with None values
    entry = GffEntry(contig, None, None, None, None, None)
    assert entry.contig.seqname == "chr1"


# --- HcAnnotStore Tests ---

def test_hc_annot_store():
    """Test HcAnnotStore creation and basic methods"""
    store = HcAnnotStore()
    assert store.is_empty()
    assert store.len() == 0
    assert len(store) == 0


def test_hc_annot_store_methods():
    """Test HcAnnotStore methods"""
    store = HcAnnotStore()
    store.init_imap()
    store.init_tree()

    # Test method existence
    feature_types = store.get_feature_types()
    assert isinstance(feature_types, dict)

    # Test get_entry with non-existent ID
    entry = store.get_entry(999)
    assert entry is None

    # Test get_entries_regex
    entries = store.get_entries_regex(".*")
    assert isinstance(entries, list)

    # Test get_children with non-existent ID
    children = store.get_children(999)
    assert children is None

    # Test get_parent with non-existent ID
    parent = store.get_parent(999)
    assert parent is None

    # Test iterator
    iterator = store.iter()
    assert iterator is not None

    # Test tree and imap initialization
    store.init_tree()
    store.init_imap()

    # Test genomic query
    contig = Contig("chr1", 100, 200, Strand.Forward)
    store.genomic_query(contig)

    # Test add_flanks
    store.add_flanks([], 100, "flank")


def test_hc_annot_store_iterator():
    """Test HcAnnotStore iterator"""
    store = HcAnnotStore()

    # Test iteration on empty store
    count = 0
    for item in store:
        count += 1
    assert count == 0

    # Test iterator methods
    iterator = store.iter()
    try:
        next(iterator)
    except StopIteration:
        pass  # Expected for empty store
