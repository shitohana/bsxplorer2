import polars as pl
import pytest
from bsx2.types import BsxBatch


@pytest.fixture
def create_dummy_df() -> pl.DataFrame:
    return pl.DataFrame(dict(
        chr=["chr1", "chr1", "chr1"],
        position=[10, 20, 30],
        strand=["+", "+", "-"],
        context=["CG", "CHG", "CHH"],
        count_m=[5, 10, 15],
        count_total=[10, 20, 30],
        density=[0.5, 0.5, 0.5]
    ))

@pytest.fixture
def create_dummy_batch(create_dummy_df) -> BsxBatch:
    return BsxBatch.from_dataframe(create_dummy_df)
