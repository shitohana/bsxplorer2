from pathlib import Path

import polars as pl
import pytest
from bsx2.io import BsxFileWriter
from bsx2.types import BsxBatch


@pytest.fixture
def sample_batch() -> BsxBatch:
    return BsxBatch.from_dataframe(
        pl.DataFrame(
            {
                "chr": ["chr1"] * 3,
                "position": [10, 20, 30],
                "strand": [True, False, True],
                "context": [True, False, None],
                "count_m": [2, 4, 6],
                "count_total": [10, 10, 10],
                "density": [0.2, 0.4, 0.6],
            }
        ),
        chr_values=["chr1", "chr2"],
    )


@pytest.fixture
def bsx_path(tmp_path: Path, sample_batch: BsxBatch) -> Path:
    path = tmp_path / "tiny.bsx"
    with BsxFileWriter(path, ["chr1", "chr2"]) as writer:
        writer.write_batch(sample_batch.slice(0, 2))
        writer.write_batch(sample_batch.slice(2, 1))
    return path


@pytest.fixture
def fasta_path(tmp_path: Path) -> Path:
    path = tmp_path / "reference.fa"
    # Reference is only used for alignment-dependent report formats.
    path.write_text(">chr1\n" + "ACG" * 20 + "\n>chr2\n" + "ACG" * 20 + "\n")
    path.with_suffix(".fa.fai").write_text(
        "chr1\t60\t6\t60\t61\nchr2\t60\t73\t60\t61\n"
    )
    return path
