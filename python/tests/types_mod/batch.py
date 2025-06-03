from bsx2.types import BsxBatch, AggMethod, ReportTypeSchema, ContextData
import polars as pl

def test_batch():
    """Test basic batch creation and properties"""
    # Create a simple batch
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True, True, False]  # True for Forward, False for Reverse
    context = [True, False, None]  # True for CG, False for CHG, None for CHH
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    assert not batch.is_empty()
    assert batch.height() == 3
    assert batch.seqname() == "chr1"
    assert batch.first_pos() == 100
    assert batch.last_pos() == 300

def test_empty_batch():
    """Test empty batch creation"""
    empty_batch = BsxBatch.empty()

    assert empty_batch.is_empty()
    assert empty_batch.height() == 0

def test_partition():
    """Test partition functionality"""
    chr_name = "chr1"
    positions = [10, 20, 30, 40, 50]
    strand = [True] * 5
    context = [True] * 5
    count_m = [1, 2, 3, 4, 5]
    count_total = [10, 20, 30, 40, 50]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    # Test partition with breakpoints
    breakpoints = [1, 3]
    densities, coverages = batch.partition(breakpoints, AggMethod.Mean)

    assert len(densities) == 3  # Should create 3 partitions
    assert len(coverages) == 3

def test_discretise_basic_case():
    """Test basic discretisation"""
    chr_name = "chr1"
    positions = [10, 20, 30, 40, 50]
    strand = [True] * 5
    context = [True] * 5
    count_m = [1, 2, 3, 4, 5]
    count_total = [10, 20, 30, 40, 50]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    # Test discretisation into 3 fragments
    densities, coverages = batch.discretise(3, AggMethod.Mean)

    assert len(densities) == 3
    assert len(coverages) == 3

def test_discretise_n_fragments_one():
    """Test discretisation with n_fragments = 1"""
    chr_name = "chr1"
    positions = [10, 20, 30]
    strand = [True] * 3
    context = [True] * 3
    count_m = [1, 2, 3]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    densities, coverages = batch.discretise(1, AggMethod.Mean)

    assert len(densities) == 1
    assert len(coverages) == 1

def test_discretise_different_agg_methods():
    """Test different aggregation methods"""
    chr_name = "chr1"
    positions = [10, 20, 30, 40]
    strand = [True] * 4
    context = [True] * 4
    count_m = [1, 2, 3, 4]
    count_total = [10, 20, 30, 40]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    # Test different aggregation methods
    for method in [AggMethod.Mean, AggMethod.Max, AggMethod.Min, AggMethod.Median]:
        densities, coverages = batch.discretise(2, method)
        assert len(densities) == 2
        assert len(coverages) == 2

def test_into_report():
    """Test conversion to different report formats"""
    chr_name = "chr1"
    positions = [100, 200]
    strand = [True, False]
    context = [True, False]
    count_m = [5, 10]
    count_total = [10, 20]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    # Test conversion to Bismark format
    bismark_df = batch.into_report(ReportTypeSchema.Bismark)
    assert isinstance(bismark_df, pl.DataFrame)
    assert bismark_df.height == 2

def test_add_context_data():
    """Test adding context data"""
    chr_name = "chr1"
    positions = [100, 200]
    strand = [True, False]
    context = [None, None]  # No context initially
    count_m = [5, 10]
    count_total = [10, 20]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    # Create context data
    context_data = ContextData([65, 84, 67, 71])  # ATCG sequence

    # Add context data
    updated_batch = batch.add_context_data(context_data)
    assert updated_batch.height() == batch.height()

def test_column_getters():
    """Test column getter methods"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True, True, False]
    context = [True, False, None]
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    # Test individual column getters
    chr_col = batch.chr()
    pos_col = batch.position()
    strand_col = batch.strand()
    context_col = batch.context()
    count_m_col = batch.count_m()
    count_total_col = batch.count_total()
    density_col = batch.density()

    assert chr_col.len() == 3
    assert pos_col.len() == 3
    assert strand_col.len() == 3
    assert context_col.len() == 3
    assert count_m_col.len() == 3
    assert count_total_col.len() == 3
    assert density_col.len() == 3

def test_split_at():
    """Test split_at method"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True, True, False]
    context = [True, False, None]
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    left, right = batch.split_at(1)

    assert left.height() == 1
    assert right.height() == 2
    assert left.first_pos() == 100
    assert right.first_pos() == 200

def test_slice():
    """Test slice method"""
    chr_name = "chr1"
    positions = [100, 200, 300, 400]
    strand = [True] * 4
    context = [True] * 4
    count_m = [5, 10, 15, 20]
    count_total = [10, 20, 30, 40]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    sliced = batch.slice(1, 2)

    assert sliced.height() == 2
    assert sliced.first_pos() == 200
    assert sliced.last_pos() == 300

def test_position_methods():
    """Test position-related methods"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True] * 3
    context = [True] * 3
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    assert batch.first_pos() == 100
    assert batch.last_pos() == 300

    first_gpos = batch.first_genomic_pos()
    last_gpos = batch.last_genomic_pos()

    assert first_gpos is not None
    assert last_gpos is not None
    assert first_gpos.seqname == "chr1"
    assert first_gpos.position == 100
    assert last_gpos.position == 300

def test_get_methylation_stats():
    """Test methylation statistics"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True] * 3
    context = [True] * 3
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    stats = batch.get_methylation_stats()

    assert stats.total_coverage() > 0
    assert 0.0 <= stats.mean_methylation() <= 1.0

def test_get_coverage_dist():
    """Test coverage distribution"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True] * 3
    context = [True] * 3
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    coverage_dist = batch.get_coverage_dist()

    assert isinstance(coverage_dist, dict)
    assert len(coverage_dist) > 0

def test_get_context_stats():
    """Test context statistics"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True] * 3
    context = [True, False, None]  # CG, CHG, CHH
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    context_stats = batch.get_context_stats()

    assert isinstance(context_stats, dict)

def test_get_strand_stats():
    """Test strand statistics"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True, True, False]  # Forward, Forward, Reverse
    context = [True] * 3
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    strand_stats = batch.get_strand_stats()

    assert isinstance(strand_stats, dict)

def test_as_binom():
    """Test binomial filtering"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True] * 3
    context = [True] * 3
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    filtered = batch.as_binom(0.5, 0.05)

    assert filtered.height() <= batch.height()

def test_can_extend():
    """Test extend functionality"""
    chr_name = "chr1"
    positions1 = [100, 200]
    positions2 = [300, 400]
    strand = [True] * 2
    context = [True] * 2
    count_m = [5, 10]
    count_total = [10, 20]

    batch1 = BsxBatch(chr_name, None, positions1, strand, context, count_m, count_total)
    batch2 = BsxBatch(chr_name, None, positions2, strand, context, count_m, count_total)

    # Test extend
    batch1.extend(batch2)

    assert batch1.height() == 4
    assert batch1.last_pos() == 400

def test_concat():
    """Test concatenation of batches"""
    chr_name = "chr1"
    positions1 = [100, 200]
    positions2 = [300, 400]
    strand = [True] * 2
    context = [True] * 2
    count_m = [5, 10]
    count_total = [10, 20]

    batch1 = BsxBatch(chr_name, None, positions1, strand, context, count_m, count_total)
    batch2 = BsxBatch(chr_name, None, positions2, strand, context, count_m, count_total)

    concatenated = BsxBatch.concat([batch1, batch2])

    assert concatenated.height() == 4
    assert concatenated.first_pos() == 100
    assert concatenated.last_pos() == 400

def test_lazy():
    """Test lazy batch operations"""
    chr_name = "chr1"
    positions = [100, 200, 300]
    strand = [True] * 3
    context = [True] * 3
    count_m = [5, 10, 15]
    count_total = [10, 20, 30]

    batch = BsxBatch(chr_name, None, positions, strand, context, count_m, count_total)

    lazy_batch = batch.lazy()

    # Test lazy operations
    filtered_lazy = lazy_batch.filter_pos_gt(150)
    filtered_batch = filtered_lazy.collect()

    assert filtered_batch.height() == 2
    assert filtered_batch.first_pos() == 200
