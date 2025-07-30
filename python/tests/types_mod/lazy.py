import pytest
import polars as pl
from bsx2.types import (
    LazyBsxBatch,
    BsxBatch,
    Strand,
    Context,
    ReportTypeSchema,
    AggMethod,
    MethylationStats,
)


def create_test_batch() -> BsxBatch:
    """Create a test batch similar to the Rust implementation"""
    return BsxBatch(
        chr="chr1",
        chr_dtype=None,
        positions=[10, 20, 30],
        strand=[True, True, False],  # Forward, Forward, Reverse
        context=[True, False, None],  # CG, CHG, CHH
        count_m=[5, 10, 15],
        count_total=[10, 20, 30],
    )


def create_empty_batch() -> BsxBatch:
    """Create an empty test batch"""
    return BsxBatch.empty()


@pytest.fixture
def test_batch() -> BsxBatch:
    return create_test_batch()


@pytest.fixture
def empty_batch() -> BsxBatch:
    return create_empty_batch()


@pytest.fixture
def test_lazy_batch(test_batch: BsxBatch) -> LazyBsxBatch:
    return test_batch.lazy()


def test_lazy_batch_creation(test_lazy_batch: LazyBsxBatch):
    """Test basic creation of LazyBsxBatch"""
    assert isinstance(test_lazy_batch, LazyBsxBatch)


def test_lazy_collect(test_batch: BsxBatch, test_lazy_batch: LazyBsxBatch):
    """Test collect method for LazyBsxBatch"""
    collected_batch = test_lazy_batch.collect()

    assert isinstance(collected_batch, BsxBatch)
    assert test_batch.data().equals(collected_batch.data())


@pytest.mark.parametrize(
    "filter_func, pos, expected_positions",
    [
        ("filter_pos_lt", 25, [10, 20]),
        ("filter_pos_gt", 15, [20, 30]),
    ],
)
def test_lazy_filter_pos(
    test_lazy_batch: LazyBsxBatch,
    filter_func: str,
    pos: int,
    expected_positions: list[int],
):
    """Test position filtering for LazyBsxBatch"""
    filtered_lazy_batch = getattr(test_lazy_batch, filter_func)(pos)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.position().to_list() == expected_positions


def test_lazy_filter_coverage_gt(test_lazy_batch: LazyBsxBatch):
    """Test filter_coverage_gt for LazyBsxBatch"""
    # Filter coverage > 15 (should keep 20 and 30, filter out 10)
    filtered_lazy_batch = test_lazy_batch.filter_coverage_gt(15)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.count_total().to_list() == [20, 30]


@pytest.mark.parametrize(
    "strand, expected_positions",
    [
        (Strand.Forward, [10, 20]),
        (Strand.Reverse, [30]),
        (Strand.Null, []),
    ],
)
def test_lazy_filter_strand(
    test_lazy_batch: LazyBsxBatch, strand: Strand, expected_positions: list[int]
):
    """Test filter_strand for LazyBsxBatch"""
    filtered_lazy_batch = test_lazy_batch.filter_strand(strand)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.position().to_list() == expected_positions


@pytest.mark.parametrize(
    "context, expected_positions",
    [
        (Context.CG, [10]),
        (Context.CHG, [20]),
    ],
)
def test_lazy_filter_context(
    test_lazy_batch: LazyBsxBatch, context: Context, expected_positions: list[int]
):
    """Test filter_context for LazyBsxBatch"""
    filtered_lazy_batch = test_lazy_batch.filter_context(context)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.position().to_list() == expected_positions


def test_lazy_batch_operations(test_lazy_batch: LazyBsxBatch):
    """Test chaining multiple operations on LazyBsxBatch"""
    # Chain filters: position > 15 AND strand Forward
    result = (test_lazy_batch
              .filter_pos_gt(15)
              .filter_strand(Strand.Forward)
              .collect())

    assert result.position().to_list() == [20]
    assert result.height() == 1


def test_lazy_batch_empty(empty_batch: BsxBatch):
    """Test LazyBsxBatch with empty batch"""
    lazy_empty = empty_batch.lazy()
    collected = lazy_empty.collect()

    assert collected.is_empty()
    assert collected.height() == 0


def test_lazy_batch_multiple_filters(test_lazy_batch: LazyBsxBatch):
    """Test multiple filters applied to LazyBsxBatch"""
    result = (test_lazy_batch
              .filter_pos_gt(5)
              .filter_pos_lt(25)
              .filter_coverage_gt(15)
              .collect())

    # Should keep position 20 (pos > 5, pos < 25, coverage > 15)
    assert result.position().to_list() == [20]
    assert result.count_total().to_list() == [20]


def test_lazy_batch_no_matches(test_lazy_batch: LazyBsxBatch):
    """Test LazyBsxBatch filters that match no rows"""
    # Filter for positions > 100 (no matches)
    result = test_lazy_batch.filter_pos_gt(100).collect()

    assert result.is_empty()
    assert result.height() == 0


def test_lazy_batch_preserve_schema(test_lazy_batch: LazyBsxBatch):
    """Test that LazyBsxBatch preserves schema after operations"""
    filtered = test_lazy_batch.filter_pos_gt(15).collect()
    original = test_lazy_batch.collect()

    # Schema should be the same
    assert filtered.data().schema == original.data().schema

    # Column types should be preserved
    assert filtered.position().dtype == original.position().dtype
    assert filtered.count_m().dtype == original.count_m().dtype


@pytest.mark.parametrize(
    "report_type",
    [
        ReportTypeSchema.Bismark,
        ReportTypeSchema.CgMap,
        ReportTypeSchema.BedGraph,
        ReportTypeSchema.Coverage,
    ],
)
def test_lazy_into_report(test_lazy_batch: LazyBsxBatch, report_type: ReportTypeSchema):
    """Test conversion to different report types from LazyBsxBatch"""
    # Collect first, then convert to report
    batch = test_lazy_batch.collect()
    report_df = batch.into_report(report_type)

    assert isinstance(report_df, pl.DataFrame)
    assert report_df.height == 3  # Should have same number of rows

    # Check that required columns exist based on report type
    expected_columns = report_type.col_names()
    for col in expected_columns:
        assert col in report_df.columns


def test_lazy_batch_integration_with_batch_methods(test_lazy_batch: LazyBsxBatch):
    """Test that LazyBsxBatch can be used with regular BsxBatch methods after collect"""
    # Filter and collect
    filtered_batch = test_lazy_batch.filter_pos_gt(15).collect()

    # Test batch methods work on the result
    assert filtered_batch.seqname() == "chr1"
    assert filtered_batch.first_pos() == 20
    assert filtered_batch.last_pos() == 30

    # Test statistical methods
    stats = filtered_batch.get_methylation_stats()
    assert isinstance(stats, MethylationStats)

    # Test partitioning
    densities, coverages = filtered_batch.partition([1], AggMethod.Mean)
    assert len(densities) == 2  # Two partitions
    assert len(coverages) == 2
