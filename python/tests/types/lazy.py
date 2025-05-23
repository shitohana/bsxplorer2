import pytest
import polars as pl
from bsx2.types import (
    LazyBsxBatch,
    BsxBatch,
    Strand,
    Context,
    ReportTypeSchema,
)
from tests.types import (  # type: ignore
    create_dummy_batch,
    create_dummy_df,
)


@pytest.fixture
def create_dummy_lazy(create_dummy_batch) -> LazyBsxBatch:
    return create_dummy_batch.lazy()


def test_lazy_batch_creation(create_dummy_lazy: LazyBsxBatch):
    """Test basic creation of LazyBsxBatch"""
    lazy_batch = create_dummy_lazy
    assert isinstance(lazy_batch, LazyBsxBatch)


def test_lazy_collect(
    create_dummy_batch: BsxBatch, create_dummy_lazy: LazyBsxBatch
):
    """Test collect method for LazyBsxBatch"""
    original_batch = create_dummy_batch
    lazy_batch = create_dummy_lazy

    collected_batch = lazy_batch.collect()

    assert isinstance(collected_batch, BsxBatch)
    assert original_batch.data().equals(collected_batch.data())


@pytest.mark.parametrize(
    "filter_func, pos, expected_positions",
    [
        ("filter_pos_lt", 25, [10, 20]),
        ("filter_pos_gt", 15, [20, 30]),
    ],
)
def test_lazy_filter_pos(
    create_dummy_lazy: LazyBsxBatch,
    filter_func: str,
    pos: int,
    expected_positions: list[int],
):
    """Test position filtering for LazyBsxBatch"""
    lazy_batch = create_dummy_lazy
    filtered_lazy_batch = getattr(lazy_batch, filter_func)(pos)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.position().to_list() == expected_positions


def test_lazy_filter_coverage_lt(
    create_dummy_lazy: LazyBsxBatch,
):
    """Test filter_coverage_lt for LazyBsxBatch"""
    lazy_batch = create_dummy_lazy
    # Filter coverage < 25 (should keep 10 and 20, filter out 30)
    filtered_lazy_batch = lazy_batch.filter_coverage_lt(25)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.count_total().to_list() == [10, 20]


@pytest.mark.parametrize(
    "strand, expected_positions",
    [
        (Strand.Forward, [10, 20]),
        (Strand.Reverse, [30]),
        (Strand.Null, []),
    ],
)
def test_lazy_filter_strand(
    create_dummy_lazy: LazyBsxBatch, strand: Strand, expected_positions: list[int]
):
    """Test filter_strand for LazyBsxBatch"""
    lazy_batch = create_dummy_lazy
    filtered_lazy_batch = lazy_batch.filter_strand(strand)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.position().to_list() == expected_positions


@pytest.mark.parametrize(
    "context, expected_positions",
    [
        (Context.CG, [10]),
        (Context.CHG, [20]),
        (Context.CHH, [30]),
    ],
)
def test_lazy_filter_context(
    create_dummy_lazy: LazyBsxBatch, context: Context, expected_positions: list[int]
):
    """Test filter_context for LazyBsxBatch"""
    lazy_batch = create_dummy_lazy
    filtered_lazy_batch = lazy_batch.filter_context(context)
    collected_batch = filtered_lazy_batch.collect()

    assert collected_batch.position().to_list() == expected_positions


def test_lazy_mark_low_coverage(create_dummy_lazy: LazyBsxBatch):
    """Test mark_low_coverage for LazyBsxBatch"""
    lazy_batch = create_dummy_lazy
    # Mark rows with coverage < 25 as low coverage
    marked_lazy_batch = lazy_batch.mark_low_coverage(25)
    collected_batch = marked_lazy_batch.collect()

    assert collected_batch.count_m().to_list() == [0, 0, 15]
    assert collected_batch.count_total().to_list() == [0, 0, 30]


def test_lazy_stats(create_dummy_lazy: LazyBsxBatch):
    """Test stats for LazyBsxBatch"""
    lazy_batch = create_dummy_lazy
    stats = lazy_batch.stats()

    # Based on create_dummy_df: density=[0.5, 0.5, 0.5], count_m=[5, 10, 15], count_total=[10, 20, 30]
    # Total coverage: 10 + 20 + 30 = 60
    # Mean density: (0.5 + 0.5 + 0.5) / 3 = 0.5
    # Variance density: should be close to 0 for constant density

    assert stats.total_coverage() == 60
    assert stats.mean_methylation() == pytest.approx(0.5)
    assert stats.methylation_var() == pytest.approx(0.0)


@pytest.mark.parametrize(
    "report_type",
    [
        ReportTypeSchema.Bismark,
        ReportTypeSchema.CgMap,
        ReportTypeSchema.BedGraph,
        ReportTypeSchema.Coverage,
    ],
)
def test_lazy_into_report(
    create_dummy_lazy: LazyBsxBatch, report_type: ReportTypeSchema
):
    # TODO
    pass


def test_lazy_align_with_contexts(create_dummy_lazy: LazyBsxBatch):
    # TODO
    pass
