from bsx2.types import Contig, GenomicPosition, Strand
import pytest

# --- GenomicPosition Tests ---

def test_genomic_position_creation():
    gp = GenomicPosition("chr1", 100)
    assert gp.seqname == "chr1"
    assert gp.position == 100

def test_genomic_position_eq_and_ne():
    gp1a = GenomicPosition("chr1", 100)
    gp1b = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr1", 200)
    gp3 = GenomicPosition("chr2", 100)

    assert gp1a == gp1b
    assert not (gp1a != gp1b)

    assert gp1a != gp2
    assert not (gp1a == gp2)

    assert gp1a != gp3
    assert not (gp1a == gp3)

def test_genomic_position_add_same_seqname():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr1", 50)
    result = gp1 + gp2
    assert result is not None
    assert result.seqname == "chr1"
    assert result.position == 150

def test_genomic_position_add_different_seqnames():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr2", 50)
    result = gp1 + gp2
    assert result is None

def test_genomic_position_sub_same_seqname_ge():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr1", 50)
    result = gp1 - gp2
    assert result is not None
    assert result.seqname == "chr1"
    assert result.position == 50

def test_genomic_position_sub_same_seqname_lt():
    gp1 = GenomicPosition("chr1", 50)
    gp2 = GenomicPosition("chr1", 100)
    result = gp1 - gp2
    assert result is None

def test_genomic_position_sub_different_seqnames():
    gp1 = GenomicPosition("chr1", 100)
    gp2 = GenomicPosition("chr2", 50)
    result = gp1 - gp2
    assert result is None

def test_genomic_position_display():
    gp = GenomicPosition("chr1", 12345)
    assert str(gp) == "chr1:12345"

def test_genomic_position_is_zero():
    gp1 = GenomicPosition("chr1", 0)
    assert gp1.is_zero()
    gp2 = GenomicPosition("chr1", 1)
    assert not gp2.is_zero()

# --- Contig Tests ---

def test_contig_new_invalid_range_panics():
    with pytest.raises(ValueError, match="Start position must be less than or equal to end position"):
        Contig("chr1", 100, 50)

def test_contig_strand():
    contig = Contig("chr1", 1, 10, "+")
    assert contig.strand == Strand.Forward

def test_contig_length():
    contig = Contig("chr1", 10, 100)
    assert contig.length() == 90

def test_contig_extend_upstream():
    contig = Contig("chr1", 100, 200)
    contig.extend_upstream(50)
    assert contig.start == 50
    assert contig.end == 200

    # Test saturating behavior
    contig_saturate = Contig("chr1", 10, 20)
    contig_saturate.extend_upstream(50)
    assert contig_saturate.start == 0  # Start is 0 due to saturating_sub
    assert contig_saturate.end == 20

def test_contig_extend_downstream():
    contig = Contig("chr1", 100, 200)
    contig.extend_downstream(50)
    assert contig.start == 100
    assert contig.end == 250

def test_contig_start_end_gpos():
    contig = Contig("chr1", 100, 200)
    start_gp = contig.start_gpos()
    assert start_gp.seqname == "chr1"
    assert start_gp.position == 100

    end_gp = contig.end_gpos()
    assert end_gp.seqname == "chr1"
    assert end_gp.position == 200

def test_contig_display():
    contig = Contig("chrX", 1000, 2000, "+")
    assert str(contig) == "chrX:1000-2000 (+)"

def test_contig_is_empty():
    contig1 = Contig("chr1", 0, 0)
    assert contig1.is_empty()

    contig2 = Contig("chr1", 0, 1)
    assert not contig2.is_empty()

def test_contig_is_in():
    outer = Contig("chr1", 100, 500)
    inner = Contig("chr1", 200, 400)
    partial = Contig("chr1", 50, 150)
    outside = Contig("chr1", 600, 700)
    diff_chr = Contig("chr2", 200, 400)

    assert inner.is_in(outer)
    assert not partial.is_in(outer)
    assert not outside.is_in(outer)
    assert not diff_chr.is_in(outer)
    assert not outer.is_in(inner)

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
