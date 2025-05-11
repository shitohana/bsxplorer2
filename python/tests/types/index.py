import pytest
from bsx2.types import BatchIndex, Contig

def test_new():
    """Test the constructor."""
    index = BatchIndex()
    assert isinstance(index, BatchIndex)

def test_insert_and_find():
    """Test inserting contigs and finding overlaps."""
    index = BatchIndex()

    contig1 = Contig("chr1", 1, 100)
    contig2 = Contig("chr1", 50, 150)
    contig3 = Contig("chr2", 1, 100)

    index.insert(contig1.seqname, contig1.start, contig1.end, 1)
    index.insert(contig2.seqname, contig2.start, contig2.end, 2)
    index.insert(contig3.seqname, contig3.start, contig3.end, 3)

    # Test finding overlapping contigs
    query_contig1 = Contig("chr1", 20, 60)
    result1 = index.find(query_contig1.seqname, query_contig1.start, query_contig1.end)
    assert result1 is not None
    assert result1 == [1, 2]

    # Test finding non-overlapping contigs
    query_contig2 = Contig("chr1", 200, 300)
    result2 = index.find(query_contig2.seqname, query_contig2.start, query_contig2.end)
    assert result2 is None

    # Test finding contigs on different chromosomes
    query_contig3 = Contig("chr2", 50, 60)
    result3 = index.find(query_contig3.seqname, query_contig3.start, query_contig3.end)
    assert result3 is not None
    assert result3 == [3] # Single element list order is fixed

    # Test finding on a non-existent chromosome
    query_contig4 = Contig("chrX", 10, 20)
    result4 = index.find(query_contig4.seqname, query_contig4.start, query_contig4.end)
    assert result4 is None

def test_get_chr_order():
    """Test retrieving chromosome order."""
    index = BatchIndex()
    index.insert("chrA", 1, 100, 1)
    index.insert("chrC", 1, 100, 2)
    index.insert("chrB", 1, 100, 3)

    order = index.get_chr_order()
    # IndexSet preserves insertion order
    assert order == ["chrA", "chrC", "chrB"]

def test_get_chr_index():
    """Test retrieving the index of a chromosome in the order."""
    index = BatchIndex()
    index.insert("chrA", 1, 100, 1)
    index.insert("chrC", 1, 100, 2)
    index.insert("chrB", 1, 100, 3)

    assert index.get_chr_index("chrA") == 0
    assert index.get_chr_index("chrC") == 1
    assert index.get_chr_index("chrB") == 2
    assert index.get_chr_index("chrX") is None

def test_sort():
    """Test sorting a list of contigs."""
    index = BatchIndex()
    # Establish a specific chromosome order
    index.insert("chr2", 1, 10, 1)
    index.insert("chr1", 1, 10, 2)
    index.insert("chr3", 1, 10, 3)

    contig1 = Contig("chr1", 50, 150)
    contig2 = Contig("chr1", 1, 100)
    contig3 = Contig("chr2", 1, 100)
    contig4 = Contig("chr2", 50, 150)
    contig5 = Contig("chr3", 10, 20)

    contigs_to_sort = [contig1, contig2, contig3, contig4, contig5]

    # The index uses the chromosome order established by insertions
    # Expected order: chr2 (by start), then chr1 (by start), then chr3
    expected_sorted = [contig3, contig4, contig2, contig1, contig5]

    sorted_contigs = index.sort(contigs_to_sort)

    assert len(sorted_contigs) == len(expected_sorted)
    for i in range(len(sorted_contigs)):
        assert sorted_contigs[i] == expected_sorted[i]

@pytest.fixture
def temp_index_file(tmp_path):
    """Fixture to provide a temporary file path for saving/loading."""
    return tmp_path / "batch_index.bin"

def test_save_and_load(temp_index_file):
    """Test saving the index to a file and loading it back."""
    original_index = BatchIndex()
    original_index.insert("chrA", 1, 100, 1)
    original_index.insert("chrB", 50, 150, 2)
    original_index.insert("chrA", 200, 300, 3)

    # Save the index
    original_index.save(str(temp_index_file))

    # Load the index
    loaded_index = BatchIndex.load(str(temp_index_file))

    # Test if the loaded index has the same content
    # Check chromosome order
    assert loaded_index.get_chr_order() == original_index.get_chr_order()

    # Check content by performing queries
    query1 = Contig("chrA", 60, 120)
    result1_orig = original_index.find(query1.seqname, query1.start, query1.end)
    result1_loaded = loaded_index.find(query1.seqname, query1.start, query1.end)
    assert result1_orig == result1_loaded

    query2 = Contig("chrB", 100, 200)
    result2_orig = original_index.find(query2.seqname, query2.start, query2.end)
    result2_loaded = loaded_index.find(query2.seqname, query2.start, query2.end)
    assert result2_orig == result2_loaded

    query3 = Contig("chrC", 10, 20) # Non-existent
    result3_orig = original_index.find(query3.seqname, query3.start, query3.end)
    result3_loaded = loaded_index.find(query3.seqname, query3.start, query3.end)
    assert result3_orig is None and result3_loaded is None
