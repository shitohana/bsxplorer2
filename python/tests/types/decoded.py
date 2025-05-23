import pytest
from bsx2.types import BsxBatch
from tests.types import create_dummy_batch, create_dummy_df #type: ignore
import polars as pl

def test_batch_creation(create_dummy_batch):
    """Test basic creation of BsxBatch"""
    batch = create_dummy_batch

    assert not batch.is_empty()
    assert batch.chr_val() == "chr1"
    assert batch.height() == 3
    assert len(batch) == 3
    assert batch.start_pos() == 10
    assert batch.end_pos() == 30

@pytest.mark.parametrize("field, correct", [
    ("chr", ["chr1", "chr1", "chr1"]),
    ("position", [10, 20, 30]),
    ("strand", ["+", "+", "-"]),
    ("context", ["CG", "CHG", "CHH"]),
    ("count_m", [5, 10, 15]),
    ("count_total", [10, 20, 30]),
    ("density", [0.5, 0.5, 0.5])
])
def test_batch_accessors(create_dummy_batch, field, correct):
    """Test BsxBatch accessors"""
    batch = create_dummy_batch

    method = getattr(batch, field)
    assert method().to_list() == correct

def test_batch_take_and_data(create_dummy_batch):
    """Test take and data methods"""
    batch = create_dummy_batch

    taken_df = batch.take()
    data_df = batch.data()

    assert taken_df.equals(data_df, null_equal=True)

def test_batch_split_at(create_dummy_batch):
    """Test split_at method"""
    batch = create_dummy_batch

    batch1, batch2 = batch.split_at(1)

    assert batch1.height() == 1
    assert batch2.height() == 2

    assert batch1.chr_val() == "chr1"
    assert batch2.chr_val() == "chr1"

    assert batch1.position().to_list() == [10]
    assert batch2.position().to_list() == [20, 30]

def test_batch_filter_mask(create_dummy_batch):
    """Test filter_mask method"""
    batch = create_dummy_batch

    # Create a mask for positions > 15
    mask = pl.Series([pos > 15 for pos in batch.position().to_list()])
    filtered = batch.filter_mask(mask)

    assert filtered.height() == 2
    assert filtered.position().to_list() == [20, 30]

@pytest.mark.xfail
def test_batch_vstack(create_dummy_batch):
    """Test vstack method"""
    batch1 = create_dummy_batch

    # Create a second batch with only the first row
    batch2 = BsxBatch.from_dataframe_unchecked(batch1.data().slice(0, 1))

    combined = batch1.vstack(batch2)

    assert combined.height() == 4
    assert combined.position().to_list() == [10, 20, 30, 10]

def test_batch_empty():
    """Test empty batch creation"""
    empty_batch = BsxBatch.empty()

    assert empty_batch.is_empty()
    assert empty_batch.height() == 0
    assert len(empty_batch) == 0

def test_batch_equality(create_dummy_df):
    """Test batch equality comparison"""
    df = create_dummy_df
    batch1 = BsxBatch(df)
    batch2 = BsxBatch(df)
    batch3 = BsxBatch(df.slice(0, 2))  # Only first two rows

    assert batch1 == batch2
    assert batch1 != batch3