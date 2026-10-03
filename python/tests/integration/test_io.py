import io

import polars as pl
import pytest
from bsx2.io import (
    BsxFileReader,
    BsxFileWriter,
    Compression,
    FilterOperation,
    IpcCompression,
    RegionReader,
    ReportReader,
    ReportWriter,
)
from bsx2.types import BatchIndex, BsxBatch, Context, Contig, ReportTypeSchema, Strand
from polars.testing import assert_frame_equal


def require_batch(batch: BsxBatch | None) -> BsxBatch:
    assert batch is not None
    return batch


@pytest.mark.parametrize("compression", [None, IpcCompression.LZ4, IpcCompression.ZSTD])
def test_bsx_roundtrip(tmp_path, sample_batch, compression):
    path = tmp_path / "roundtrip.bsx"
    with BsxFileWriter(path, ["chr1", "chr2"], compression) as writer:
        writer.write_batch(sample_batch)
    reader = BsxFileReader(path)
    assert reader.blocks_total == 1
    assert reader.n_threads >= 1
    actual = reader.get_batch(0)
    assert actual is not None
    assert_frame_equal(actual.data(), sample_batch.data())
    assert reader.get_batch(1) is None
    assert reader.get_batches([0, 1])[1] is None
    reader.cache_batches([0])
    assert [b.position().to_list() for b in reader] == [[10, 20, 30]]
    assert [b.position().to_list() for b in reader] == [[10, 20, 30]]
    with pytest.raises(StopIteration):
        next(reader)


def test_reader_restart_after_partial_iteration(bsx_path):
    reader = BsxFileReader(bsx_path)
    assert next(iter(reader)).position().to_list() == [10, 20]
    assert [b.position().to_list() for b in reader] == [[10, 20], [30]]
    assert BsxFileReader(str(bsx_path)).blocks_total == 2
    with bsx_path.open("rb") as handle:
        assert require_batch(
            BsxFileReader(handle).get_batch(1)
        ).position().to_list() == [30]


@pytest.mark.parametrize("reader_class", [BsxFileReader, RegionReader])
def test_reader_failures(tmp_path, reader_class):
    with pytest.raises(FileNotFoundError):
        reader_class(tmp_path / "missing")
    with pytest.raises(OSError):
        reader_class(io.BytesIO(b"invalid"))
    malformed = tmp_path / "invalid.bsx"
    malformed.write_bytes(b"not an IPC file")
    with pytest.raises(RuntimeError):
        reader_class(malformed)


def test_bsx_filelike_writer_and_lifecycle(tmp_path, sample_batch):
    stream = io.BytesIO()
    writer = BsxFileWriter(stream, ["chr1", "chr2"])
    writer.write_batch(sample_batch)
    writer.close()
    writer.close()
    with pytest.raises(OSError, match="closed"):
        writer.write_batch(sample_batch)
    path = tmp_path / "bytes.bsx"
    path.write_bytes(stream.getvalue())
    assert require_batch(BsxFileReader(path).get_batch(0)).height() == 3
    with (
        pytest.raises(RuntimeError, match="sentinel"),
        BsxFileWriter(tmp_path / "context.bsx", ["chr1", "chr2"]) as w,
    ):
        w.write_batch(sample_batch)
        raise RuntimeError("sentinel")
    assert BsxFileReader(tmp_path / "context.bsx").blocks_total == 1


@pytest.mark.parametrize("factory", ["from_sink_and_fai", "from_sink_and_fasta"])
def test_bsx_writer_factories(tmp_path, sample_batch, fasta_path, factory):
    source = (
        fasta_path if factory.endswith("fasta") else fasta_path.with_suffix(".fa.fai")
    )
    path = tmp_path / "factory.bsx"
    with getattr(BsxFileWriter, factory)(path, source) as writer:
        writer.write_batch(sample_batch)
    assert require_batch(BsxFileReader(path).get_batch(0)).height() == 3


def test_batch_index_persistence_and_bounds(tmp_path):
    index = BatchIndex()
    index.insert("chr1", 10, 30, 0)
    index.insert("chr2", 10, 30, 1)
    assert index.find("chr1", 15, 25) == [0]
    assert index.find("unknown", 15, 25) is None
    assert index.find("chr1", 50, 60) == []
    assert index.get_chr_order() == ["chr1", "chr2"]
    assert index.get_chr_index("chr2") == 1
    assert index.get_chr_index("missing") is None
    contigs = [Contig("chr2", 10, 30, Strand.Null), Contig("chr1", 10, 30, Strand.Null)]
    assert [c.seqname for c in index.sort(contigs)] == ["chr1", "chr2"]
    path = tmp_path / "index.bci"
    index.save(path)
    assert BatchIndex.load(path).find("chr1", 15, 25) == [0]
    index.save(str(path))
    assert BatchIndex.load(str(path)).get_chr_order() == ["chr1", "chr2"]
    with pytest.raises(ValueError):
        index.insert("chr1", 30, 10, 0)
    with pytest.raises(ValueError):
        index.find("chr1", 30, 10)


@pytest.mark.parametrize(
    "reader_factory",
    [RegionReader, lambda p: RegionReader.from_reader(BsxFileReader(p))],
)
def test_region_queries_and_iterators(bsx_path, reader_factory):
    reader = reader_factory(bsx_path)
    region = Contig("chr1", 10, 31, Strand.Null)
    assert reader.chr_order() == ["chr1"]
    assert reader.index().get_chr_index("chr1") == 0
    assert require_batch(reader.query(region)).position().to_list() == [10, 20, 30]
    with pytest.raises(RuntimeError, match="seqname"):
        reader.query(Contig("chr2", 10, 31, Strand.Null))
    absent = Contig("chr1", 40, 50, Strand.Null)
    assert reader.query(absent) is None
    reader.reset()
    iterator = reader.iter_contigs([absent] * 2000 + [region])
    assert iter(iterator) is iterator
    assert next(iterator).position().to_list() == [10, 20, 30]
    with pytest.raises(StopIteration):
        next(iterator)
    with pytest.raises(StopIteration):
        next(iterator)


@pytest.mark.parametrize(
    "method,value,expected",
    [
        ("filter_pos_lt", 25, [10, 20]),
        ("filter_pos_gt", 15, [20, 30]),
        ("filter_coverage_gt", 11, []),
        ("filter_strand", Strand.Forward, [10, 30]),
        ("filter_context", Context.CHG, [20]),
    ],
)
def test_region_filters(bsx_path, method, value, expected):
    reader = RegionReader(bsx_path)
    region = Contig("chr1", 10, 31, Strand.Null)
    getattr(reader, method)(value)
    assert require_batch(reader.query(region)).position().to_list() == expected
    reader.clear_filters()
    assert require_batch(reader.query(region)).height() == 3


@pytest.mark.parametrize(
    "variant,value",
    [
        ("PosLt", 25),
        ("PosGt", 15),
        ("CoverageGt", 5),
        ("Strand", Strand.Forward),
        ("Context", Context.CG),
    ],
)
def test_filter_operation_variants(bsx_path, variant, value):
    operation = getattr(FilterOperation, variant)(value)
    assert isinstance(operation, FilterOperation)
    assert operation.value == value
    reader = RegionReader(bsx_path)
    reader.add_filter(operation)
    assert reader.query(Contig("chr1", 10, 31, Strand.Null)) is not None


@pytest.mark.parametrize("schema", [ReportTypeSchema.Bismark, ReportTypeSchema.CgMap])
@pytest.mark.parametrize(
    "compression",
    [
        Compression.No,
        Compression.Gz,
        Compression.Zstd,
        Compression.Lz4,
        Compression.Xz2,
        Compression.Bzip2,
        Compression.Zip,
    ],
)
def test_report_compression_roundtrip(
    tmp_path, sample_batch, fasta_path, schema, compression
):
    path = tmp_path / "report.txt"
    writer = ReportWriter(path, schema, compression=compression)
    writer.write_batch(sample_batch)
    writer.close()
    writer.close()
    reader = ReportReader(
        path,
        schema,
        compression=compression,
        batch_size=2,
        fai_path=fasta_path.with_suffix(".fa.fai"),
    )
    assert iter(reader) is reader
    batches = list(reader)
    assert [p for b in batches for p in b.position().to_list()] == [10, 20, 30]
    assert [m for b in batches for m in b.count_m().to_list()] == [2, 4, 6]
    with pytest.raises(StopIteration):
        next(reader)
    with pytest.raises(RuntimeError, match="closed"):
        writer.write_batch(sample_batch)
    with pytest.raises(RuntimeError, match="closed"):
        writer.write_df(sample_batch.into_report(schema))


@pytest.mark.parametrize(
    "schema",
    [
        ReportTypeSchema.Bismark,
        ReportTypeSchema.CgMap,
        ReportTypeSchema.BedGraph,
        ReportTypeSchema.Coverage,
    ],
)
def test_report_dataframe_writer_and_reader(tmp_path, sample_batch, fasta_path, schema):
    path = tmp_path / "frame.tsv"
    with ReportWriter(path, schema) as writer:
        writer.write_df(sample_batch.into_report(schema))
    reader = ReportReader(path, schema, fasta_path=fasta_path)
    batches = list(reader)
    assert sum(b.height() for b in batches) > 0
    assert all(isinstance(b, BsxBatch) for b in batches)


def test_report_filelike_sink(sample_batch):
    stream = io.BytesIO()
    with ReportWriter(stream, ReportTypeSchema.Bismark) as writer:
        writer.write_batch(sample_batch)
    assert len(stream.getvalue().splitlines()) == 3


def test_report_failures(tmp_path):
    with pytest.raises(FileNotFoundError):
        ReportReader(tmp_path / "missing", ReportTypeSchema.Bismark)
    with pytest.raises(ValueError):
        ReportReader(tmp_path / "missing", ReportTypeSchema.Bismark, chunk_size=0)
    with pytest.raises(ValueError):
        ReportReader(tmp_path / "missing", ReportTypeSchema.Bismark, batch_size=0)
    with pytest.raises(ValueError):
        ReportWriter(tmp_path / "out", ReportTypeSchema.Bismark, n_threads=0)


def test_report_close_propagates_sink_flush_errors(sample_batch):
    class FailingFlush(io.BytesIO):
        def flush(self):
            raise OSError("flush failed")

    writer = ReportWriter(FailingFlush(), ReportTypeSchema.Bismark)
    writer.write_batch(sample_batch)
    with pytest.raises(OSError, match="flush failed"):
        writer.close()
    writer.close()


def test_bsx_chromosome_transitions(tmp_path, sample_batch):
    second = BsxBatch.from_dataframe(
        sample_batch.data().with_columns(pl.lit("chr2").alias("chr")),
        chr_values=["chr1", "chr2"],
    )
    path = tmp_path / "chromosomes.bsx"
    with BsxFileWriter(path, ["chr1", "chr2"]) as writer:
        writer.write_batch(sample_batch)
        writer.write_batch(second)
    assert [batch.seqname() for batch in BsxFileReader(path)] == ["chr1", "chr2"]
    reader = RegionReader(path)
    assert reader.chr_order() == ["chr1", "chr2"]
    queried = require_batch(reader.query(Contig("chr2", 10, 31, Strand.Null)))
    assert queried.seqname() == "chr2"
    assert queried.position().to_list() == [10, 20, 30]
