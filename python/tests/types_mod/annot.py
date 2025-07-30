from bsx2.types import GffEntry, Contig, Strand, BatchIndex

def create_test_entry(id: str, seqname: str, start: int, end: int) -> GffEntry:
    """Helper function to create a test GffEntry.

    This function is not meant to be used directly but is used by the
    doctests.
    """
    contig = Contig(seqname, start, end, "Forward")
    return GffEntry(contig, source="test", feature_type="gene", id=id)


def create_test_entry_with_parent(
    id: str,
    seqname: str,
    start: int,
    end: int,
    parent: list[str] | None = None
) -> GffEntry:
    """Helper function to create a test GffEntry with parent information.

    This function is not meant to be used directly but is used by the
    doctests.
    """
    contig = Contig(seqname, start, end, "Forward")
    # Note: Parent functionality would need to be exposed in Python bindings
    return GffEntry(contig, source="test", feature_type="gene", id=id)


def create_mock_index() -> BatchIndex:
    """Create a mock BatchIndex for testing"""
    index = BatchIndex()

    # Insert contigs for chr1 and chr2 with batch indices
    index.insert("chr1", 0, 10000, 0)
    index.insert("chr2", 0, 20000, 1)

    return index


def test_batch_index_operations():
    """Test basic BatchIndex operations"""
    index = create_mock_index()

    # Test chromosome order
    chr_order = index.get_chr_order()
    assert "chr1" in chr_order
    assert "chr2" in chr_order

    # Test chromosome index lookup
    chr1_idx = index.get_chr_index("chr1")
    chr2_idx = index.get_chr_index("chr2")
    assert chr1_idx is not None
    assert chr2_idx is not None

    # Test finding overlapping batches
    overlapping = index.find("chr1", 5000, 8000)
    assert overlapping is not None
    assert 0 in overlapping


def test_contig_operations():
    """Test basic Contig operations"""
    contig = Contig("chr1", 100, 200, "Forward")

    assert contig.seqname == "chr1"
    assert contig.start == 100
    assert contig.end == 200
    assert contig.strand == Strand.Forward
    assert contig.length() == 100

    # Test extensions
    contig.extend_upstream(50)
    assert contig.start == 50

    contig.extend_downstream(50)
    assert contig.end == 250

    # Test containment
    smaller_contig = Contig("chr1", 60, 70, "Forward")
    assert smaller_contig.is_in(contig)


def test_gff_entry_creation():
    """Test GffEntry creation and basic operations"""
    contig = Contig("chr1", 100, 200, "Forward")
    entry = GffEntry(contig, source="test", feature_type="gene", id="gene123")

    assert entry.id == "gene123"
    assert entry.source == "test"
    assert entry.feature_type == "gene"
    assert entry.contig.seqname == "chr1"
    assert entry.contig.start == 100
    assert entry.contig.end == 200
