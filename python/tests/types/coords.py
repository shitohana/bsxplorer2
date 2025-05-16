from bsx2.types import Contig, GenomicPosition
import pytest

@pytest.fixture
def dummy_genomic_position() -> GenomicPosition:
    return GenomicPosition("chr1", 123)

@pytest.fixture
def dummy_contig() -> Contig:
    return Contig("chr1", 123, 456)

# --- GenomicPosition Tests ---

def test_genomic_position_creation(dummy_genomic_position: GenomicPosition):
    gp = dummy_genomic_position
    assert gp.seqname == "chr1"
    assert gp.position == 123

def test_genomic_position_str_repr(dummy_genomic_position: GenomicPosition):
    gp = dummy_genomic_position
    assert str(gp) == "chr1:123"
    assert repr(gp) == "PyGenomicPosition(seqname='chr1', position=123)"

def test_genomic_position_equality():
    gp1a = GenomicPosition("chr1", 100)
    gp1b = GenomicPosition("chr1", 100)
    gp2_pos = GenomicPosition("chr1", 200)
    gp2_chr = GenomicPosition("chr2", 100)

    assert gp1a == gp1b
    assert not (gp1a != gp1b)

    assert gp1a != gp2_pos
    assert not (gp1a == gp2_pos)

    assert gp1a != gp2_chr # Equality check works across different seqnames
    assert not (gp1a == gp2_chr)

def test_genomic_position_ordering_same_chr():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr1", 200)
    gp3 = GenomicPosition("chr1", 100)

    assert gp1 < gp2
    assert gp1 <= gp2
    assert gp1 <= gp3
    assert gp2 > gp1
    assert gp2 >= gp1
    assert gp3 >= gp1

    assert not (gp1 > gp2)
    assert not (gp2 < gp1)

def test_genomic_position_ordering_diff_chr():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr2", 50) # Smaller position, different chromosome

    with pytest.raises(ValueError, match="Cannot order genomic positions on different chromosomes/seqnames"):
        _ = gp1 < gp2
    with pytest.raises(ValueError, match="Cannot order genomic positions on different chromosomes/seqnames"):
        _ = gp1 <= gp2
    with pytest.raises(ValueError, match="Cannot order genomic positions on different chromosomes/seqnames"):
        _ = gp1 > gp2
    with pytest.raises(ValueError, match="Cannot order genomic positions on different chromosomes/seqnames"):
        _ = gp1 >= gp2

def test_genomic_position_addition_same_chr():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr1", 50)
    result = gp1 + gp2
    assert result is not None
    assert result.seqname == "chr1"
    assert result.position == 150

def test_genomic_position_addition_diff_chr():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr2", 50)
    result = gp1 + gp2
    assert result is None

def test_genomic_position_subtraction_same_chr_valid():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr1", 50)
    result = gp1 - gp2
    assert result is not None
    assert result.seqname == "chr1"
    assert result.position == 50

def test_genomic_position_subtraction_same_chr_invalid_result():
    gp1 = GenomicPosition("chr1", 50)
    gp2 = GenomicPosition("chr1", 100)
    result = gp1 - gp2 # 50 - 100 would be negative
    assert result is None

def test_genomic_position_subtraction_diff_chr():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr2", 50)
    result = gp1 - gp2
    assert result is None

def test_genomic_position_is_zero():
    gp1 = GenomicPosition("chr1", 0)
    assert gp1.is_zero()
    gp2 = GenomicPosition("chr1", 1)
    assert not gp2.is_zero()
    gp3 = GenomicPosition("", 0) # seqname not checked by is_zero in python
    assert gp3.is_zero()


# --- Contig Tests ---

def test_contig_creation_valid():
    c = Contig("chr1", 100, 200, "+")
    assert c.seqname == "chr1"
    assert c.start == 100
    assert c.end == 200
    assert c.strand_str == "+"

    c_rev = Contig("chrX", 0, 10, "REVERSE")
    assert c_rev.strand_str == "-"

    c_none = Contig("chrM", 5, 15, ".")
    assert c_none.strand_str == "."

    c_none_word = Contig("chrM", 5, 15, "None")
    assert c_none_word.strand_str == "."


def test_contig_creation_invalid_start_end():
    with pytest.raises(ValueError, match="Start position must be less than or equal to end position"):
        Contig("chr1", 200, 100)

def test_contig_creation_invalid_strand():
    with pytest.raises(ValueError, match="Invalid strand value: 'x'"):
        Contig("chr1", 100, 200, "x")

def test_contig_properties_default_strand():
    c = Contig("chr1", 100, 200) # Default strand is "."
    assert c.seqname == "chr1"
    assert c.start == 100
    assert c.end == 200
    assert c.strand_str == "."

def test_contig_properties_custom_strand(dummy_contig: Contig):
    c = dummy_contig
    assert c.seqname == "chr1"
    assert c.start == 123
    assert c.end == 456
    assert c.strand_str == "." # dummy_contig has default strand

    c_fwd = Contig("chr2", 10, 20, "Forward")
    assert c_fwd.strand_str == "+"

def test_contig_str_repr(dummy_contig: Contig):
    c = dummy_contig
    assert str(c) == "chr1:123-456 (.)"
    assert repr(c) == "PyContig(seqname='chr1', start=123, end=456, strand='.')"

    c_fwd = Contig("chr2", 10, 20, "+")
    assert str(c_fwd) == "chr2:10-20 (+)"
    assert repr(c_fwd) == "PyContig(seqname='chr2', start=10, end=20, strand='+')"

def test_contig_length(dummy_contig: Contig):
    c = dummy_contig # 123-456
    assert c.length() == 456 - 123

    c_zero_len = Contig("chr1", 100, 100)
    assert c_zero_len.length() == 0

def test_contig_start_end_gpos(dummy_contig: Contig):
    c = dummy_contig
    start_gp = c.start_gpos()
    assert start_gp.seqname == "chr1"
    assert start_gp.position == 123

    end_gp = c.end_gpos()
    assert end_gp.seqname == "chr1"
    assert end_gp.position == 456

def test_contig_extend_upstream(dummy_contig: Contig):
    c = Contig(dummy_contig.seqname, dummy_contig.start, dummy_contig.end, dummy_contig.strand_str)
    c.extend_upstream(23)
    assert c.start == 100
    assert c.end == 456

def test_contig_extend_upstream_saturating():
    c = Contig("chr1", 10, 20)
    c.extend_upstream(50)
    assert c.start == 0 # Saturates at 0
    assert c.end == 20

def test_contig_extend_downstream(dummy_contig: Contig):
    c = Contig(dummy_contig.seqname, dummy_contig.start, dummy_contig.end, dummy_contig.strand_str)
    c.extend_downstream(44)
    assert c.start == 123
    assert c.end == 500

def test_contig_is_in():
    outer = Contig("chr1", 100, 500)
    inner_same_strand = Contig("chr1", 200, 400, ".")
    inner_diff_strand = Contig("chr1", 200, 400, "+") # is_in ignores strand

    partial_overlap_left = Contig("chr1", 50, 150)
    partial_overlap_right = Contig("chr1", 450, 550)
    outside_left = Contig("chr1", 0, 50)
    outside_right = Contig("chr1", 550, 600)
    diff_chr = Contig("chr2", 200, 400)
    identical = Contig("chr1", 100, 500)

    assert inner_same_strand.is_in(outer)
    assert inner_diff_strand.is_in(outer) # Strand is not checked by is_in in Rust version
    assert identical.is_in(outer)

    assert not outer.is_in(inner_same_strand)
    assert not partial_overlap_left.is_in(outer)
    assert not partial_overlap_right.is_in(outer)
    assert not outside_left.is_in(outer)
    assert not outside_right.is_in(outer)
    assert not diff_chr.is_in(outer)

def test_contig_is_empty():
    c1 = Contig("chr1", 0, 0)
    assert c1.is_empty()
    c2 = Contig("chr1", 0, 1)
    assert not c2.is_empty()
    c3 = Contig("chr1", 1, 1) # start != 0
    assert not c3.is_empty()
    c4 = Contig("", 0, 0) # Empty seqname would not match Rust, but Python only checks start/end
    assert c4.is_empty()


def test_contig_equality():
    c1a = Contig("chr1", 100, 200, "+")
    c1b = Contig("chr1", 100, 200, "+")
    c_diff_strand = Contig("chr1", 100, 200, "-")
    c_diff_start = Contig("chr1", 101, 200, "+")
    c_diff_end = Contig("chr1", 100, 201, "+")
    c_diff_chr = Contig("chr2", 100, 200, "+")

    assert c1a == c1b
    assert not (c1a != c1b)

    assert c1a != c_diff_strand
    assert c1a != c_diff_start
    assert c1a != c_diff_end
    assert c1a != c_diff_chr

def test_contig_ordering_same_chr_no_intersect():
    c1 = Contig("chr1", 100, 200)
    c2 = Contig("chr1", 250, 300) # c1 < c2
    c3 = Contig("chr1", 50, 90)   # c3 < c1

    assert c1 < c2
    assert c1 <= c2
    assert c2 > c1
    assert c2 >= c1

    assert c3 < c1
    assert c3 <= c1
    assert c1 > c3
    assert c1 >= c3

    # Ordering ignores strand
    c1_plus = Contig("chr1", 100, 200, "+")
    c2_minus = Contig("chr1", 250, 300, "-")
    assert c1_plus < c2_minus

def test_contig_ordering_same_chr_intersect():
    c1 = Contig("chr1", 100, 200)
    c_intersect = Contig("chr1", 150, 250)

    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 < c_intersect
    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 <= c_intersect # Even for identical contigs, strict ordering might be false
    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 > c_intersect
    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 >= c_intersect

    # Test identical contigs for <= and >=
    c_identical = Contig("chr1", 100, 200)
    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 <= c_identical
    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 >= c_identical


def test_contig_ordering_diff_chr():
    c1 = Contig("chr1", 100, 200)
    c_diff_chr = Contig("chr2", 50, 150)

    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 < c_diff_chr
    with pytest.raises(NotImplementedError, match="Contigs on different seqnames or with intersecting regions are not strictly comparable."):
        _ = c1 > c_diff_chr
