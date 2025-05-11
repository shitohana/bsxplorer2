import pytest
from bsx2.types import BsxBatch, EncodedBsxBatch, encode
from tests.types import create_dummy_decoded, create_dummy_df, create_dummy_encoded #type: ignore
import polars as pl

def test_encoded_batch_creation(create_dummy_encoded):
    """Test basic creation of EncodedBsxBatch"""
    batch = create_dummy_encoded

    assert not batch.is_empty()
    assert batch.chr_val() == "chr1"
    assert batch.height() == 3
    assert len(batch) == 3
    assert batch.start_pos() == 10
    assert batch.end_pos() == 30

@pytest.mark.parametrize("field, correct", [
    ("chr", ["chr1", "chr1", "chr1"]),
    ("position", [10, 20, 30]),
    ("strand", [True, True, False]),
    ("context", [True, False, None]),
    ("count_m", [5, 10, 15]),
    ("count_total", [10, 20, 30]),
    ("density", [0.5, 0.5, 0.5])
])
def test_encoded_batch_accessors(create_dummy_encoded, field, correct):
    """Test EncodedBsxBatch accessors"""
    batch = create_dummy_encoded

    method = getattr(batch, field)
    assert method().to_list() == correct

def test_encoded_batch_take_and_data(create_dummy_encoded):
    """Test take and data methods for EncodedBsxBatch"""
    batch = create_dummy_encoded


    taken_df = batch.take()
    data_df = batch.data()

    assert taken_df.equals(data_df, null_equal=True)

def test_encoded_batch_split_at(create_dummy_encoded):
    """Test split_at method for EncodedBsxBatch"""
    batch = create_dummy_encoded

    batch1, batch2 = batch.split_at(1)

    assert batch1.height() == 1
    assert batch2.height() == 2

    assert batch1.chr_val() == "chr1"
    assert batch2.chr_val() == "chr1"

    assert batch1.position().to_list() == [10]
    assert batch2.position().to_list() == [20, 30]

def test_encoded_batch_filter_mask(create_dummy_encoded):
    """Test filter_mask method for EncodedBsxBatch"""
    batch = create_dummy_encoded

    # Create a mask for positions > 15
    mask = pl.Series([pos > 15 for pos in batch.position().to_list()])
    filtered = batch.filter_mask(mask)

    assert filtered.height() == 2
    assert filtered.position().to_list() == [20, 30]

@pytest.mark.xfail
def test_encoded_batch_vstack(create_dummy_encoded):
    """Test vstack method for EncodedBsxBatch"""
    batch1 =create_dummy_encoded

    # Create a second batch with only the first row
    batch2 = create_dummy_encoded

    combined = batch1.vstack(batch2)

    assert combined.height() == 4
    assert combined.position().to_list() == [10, 20, 30, 10]

def test_encoded_batch_empty():
    """Test empty encoded batch creation"""
    empty_batch = EncodedBsxBatch.empty()

    assert empty_batch.is_empty()
    assert empty_batch.height() == 0
    assert len(empty_batch) == 0

def test_encoded_batch_equality(create_dummy_encoded, create_dummy_df):
    """Test encoded batch equality comparison"""
    batch1 = create_dummy_encoded
    batch2 = create_dummy_encoded
    batch3 = encode(BsxBatch(create_dummy_df.slice(0, 2)), chr_values=["chr1"])  # Only first two rows

    assert batch1.data().equals(batch2.data())
    assert not batch1.data().equals(batch3.data())
