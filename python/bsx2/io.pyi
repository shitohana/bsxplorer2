from enum import Enum
from typing import Optional, Union, List, BinaryIO
from pathlib import Path

from .types import ReportTypeSchema, BsxBatch, Contig, BatchIndex, Strand, Context
import polars as pl

__all__ = [
    "BsxFileReader",
    "BsxFileWriter",
    "IpcCompression",
    "Compression",
    "ReportReader",
    "ReportWriter",
    "RegionReader",
    "FilterOperation",
    "RegionReaderIterator",
]

class FilterOperation(Enum):
    """Filter operations for RegionReader."""
    PosGt = "PosGt"
    PosLt = "PosLt"
    CoverageGt = "CoverageGt"
    Strand = "Strand"
    Context = "Context"


class RegionReaderIterator:
    """Iterator for reading BSX data for multiple contigs using a RegionReader."""

    def __iter__(self) -> "RegionReaderIterator":
        """Return self as iterator."""
        pass

    def __next__(self) -> Optional[BsxBatch]:
        """Get the next contig's data."""
        pass

class RegionReader:
    """
    A reader for BSX files that operates on specific regions of the genome.

    RegionReader is a reader for BSX files that operates on a specific region of
    the genome. It uses a BatchIndex to efficiently find and cache only the batches
    that contain data for queried genomic regions, minimizing memory usage and I/O.
    """

    def __init__(self, file: Union[str, BinaryIO]) -> None:
        """
        Creates a new RegionReader instance.

        Initializes the reader and builds an index for the BSX file.
        The index is built upon creation, which might take time for large files.

        Parameters
        ----------
        file : str or file-like object
            Path to the BSX file or a file-like object.

        Raises
        ------
        FileNotFoundError
            If the input BSX file cannot be found.
        RuntimeError
            If there's an error reading the file or building the index.
        """
        pass

    @classmethod
    def from_reader(cls, reader: "BsxFileReader") -> "RegionReader":
        """
        Create a RegionReader from an existing BsxFileReader.

        Parameters
        ----------
        reader : BsxFileReader
            An existing BsxFileReader instance.

        Returns
        -------
        RegionReader
            A new RegionReader instance.
        """
        pass

    def query(self, contig: Contig) -> Optional[pl.DataFrame]:
        """
        Queries the BSX file for data within a specific genomic region.

        This method uses the internal index and cache to efficiently retrieve
        data for the specified contig. It will read and cache batches as needed
        to fulfill the query.

        Parameters
        ----------
        contig : Contig
            The genomic region to query.

        Returns
        -------
        DataFrame or None
            A Polars DataFrame containing the data within the specified region,
            or None if no data is found for the chromosome or within the range.

        Raises
        ------
        ValueError
            If the contig seqname is not found in the index.
        RuntimeError
            If an internal error occurs during batch retrieval, caching, or assembly.
        """
        pass

    def reset(self) -> None:
        """
        Resets the internal cache.

        Clears all cached batches, forcing subsequent queries to re-read data
        from the file. Use this if memory usage becomes high or if you need
        to ensure data freshness across different query patterns.
        """
        pass

    def chr_order(self) -> List[str]:
        """
        Returns the chromosome order defined in the BSX file's index.

        Returns
        -------
        list[str]
            A list of chromosome names in the order they appear in the BSX file's index.
        """
        pass

    def add_filter(self, filter_op: FilterOperation) -> None:
        """
        Add a filter to be applied to query results.

        Parameters
        ----------
        filter_op : FilterOperation
            The filter operation to add.
        """
        pass

    def clear_filters(self) -> None:
        """Clear all filters."""
        pass

    def filter_pos_lt(self, value: int) -> None:
        """Add a filter for positions less than value."""
        pass

    def filter_pos_gt(self, value: int) -> None:
        """Add a filter for positions greater than value."""
        pass

    def filter_coverage_gt(self, value: int) -> None:
        """Add a filter for coverage greater than value."""
        pass

    def filter_strand(self, value: Strand) -> None:
        """Add a filter for specific strand."""
        pass

    def filter_context(self, value: Context) -> None:
        """Add a filter for specific context."""
        pass

    def index(self) -> BatchIndex:
        """
        Returns the BatchIndex of the BSX file.

        Returns
        -------
        BatchIndex
            The index of batches in the BSX file
        """
        pass

    def iter_contigs(self, contigs: List[Contig]) -> RegionReaderIterator:
        """
        Create an iterator over a list of contigs.

        Parameters
        ----------
        contigs : list[Contig]
            A list of Contig regions to iterate over.

        Returns
        -------
        RegionReaderIterator
            An iterator over the contigs.
        """
        pass

class Compression:
    """Compression algorithms supported for reading/writing report files."""

    pass

class ReportReader:
    """
    Reads methylation report files batch by batch.

    Reads report data and yields batches. This class provides functionality for
    reading various methylation report file formats, such as Bismark, CgMap,
    BedGraph, and Coverage reports. It handles file parsing, optional alignment
    with FASTA context, and yields data in standardized BsxBatch objects.
    """
    def __init__(
        self,
        path: Union[str, Path],
        report_type: ReportTypeSchema,
        chunk_size: int = 10000,
        fasta_path: Optional[Union[str, Path]] = None,
        fai_path: Optional[Union[str, Path]] = None,
        batch_size: int = 100000,
        n_threads: Optional[int] = None,
        low_memory: bool = False,
        compression: Optional[Compression] = None,
    ) -> None:
        """Creates a new ReportReader instance.

        Parameters
        ----------
        path : str or Path
            Path to the input report file.
        report_type : ReportTypeSchema
            The format of the report file (e.g., ReportTypeSchema.Bismark).
        chunk_size : int, optional
            The target number of records per yielded batch. Defaults to 10000.
        fasta_path : str or Path, optional
            Path to the reference FASTA file. Required for BedGraph/Coverage formats
            or when alignment is needed.
        fai_path : str or Path, optional
            Path to the FASTA index file (.fai). Can be used instead of `fasta_path`
            for determining chromosome order/names if FASTA content is not needed
            for alignment.
        batch_size : int, optional
            The number of lines read from the file at once internally. Defaults to 100000.
        n_threads : int, optional
            Number of threads to use for parsing. Defaults to using available cores.
        low_memory : bool, optional
            Whether to use a low-memory parsing mode. Defaults to False.
        compression : Compression, optional
            Compression type. Defaults to None (autodetect).

        Returns
        -------
        ReportReader
            A new instance of the ReportReader.

        Raises
        ------
        FileNotFoundError
            If the input report file or FASTA/FAI file cannot be found.
        ValueError
            If required parameters (like FASTA for certain types) are missing.
        RuntimeError
            If there's an error during reader initialization or parsing.
        """
        pass

    def __iter__(self) -> "ReportReader":
        pass

    def __next__(self) -> Optional[BsxBatch]:
        """Retrieves the next batch of methylation data.

        Returns
        -------
        BsxBatch or None
            The next batch of data, or None when iteration is complete.

        Raises
        ------
        StopIteration
            When there are no more batches to read.
        RuntimeError
            If an error occurs during reading or processing a batch.
        """
        pass

class ReportWriter:
    """
    Writes methylation report data to a file.

    Writes report data to a sink in CSV format based on a specified schema.
    This class handles writing BsxBatch objects or Polars DataFrames
    to a specified file path in a chosen report format (Bismark, CGmap, etc.).
    It manages file handling, formatting, and optional compression.
    """
    def __init__(
        self,
        file: Union[str, BinaryIO],
        schema: ReportTypeSchema,
        n_threads: int = 1,
        compression: Optional[Compression] = None,
        compression_level: Optional[int] = None,
    ) -> None:
        """Creates a new ReportWriter instance.

        Parameters
        ----------
        file : str or file-like object
            Path to the output file or a file-like object.
        schema : ReportTypeSchema
            The desired output format schema (e.g., ReportTypeSchema.Bismark).
        n_threads : int, optional
            Number of threads to use for writing. Defaults to 1.
        compression : Compression, optional
            Compression type. Defaults to None.
        compression_level : int, optional
            Compression level if applicable. Defaults depend on the compression type.

        Returns
        -------
        ReportWriter
            A new instance of the ReportWriter.

        Raises
        ------
        ValueError
            If the compression type is invalid.
        RuntimeError
            If the file cannot be created or the writer fails to initialize.
        """
        pass

    def write_batch(self, batch: BsxBatch) -> None:
        """Writes a BsxBatch to the output file.

        Parameters
        ----------
        batch : BsxBatch
            The batch of data to write.

        Raises
        ------
        RuntimeError
            If the writer is already closed or if writing fails.
        """
        pass

    def write_df(self, df: pl.DataFrame) -> None:
        """Writes a Polars DataFrame to the output file.

        Note: The DataFrame schema should ideally match the writer's schema,
        though the underlying writer might perform some conversions.

        Parameters
        ----------
        df : polars.DataFrame
            The DataFrame to write.

        Raises
        ------
        RuntimeError
            If the writer is already closed or if writing fails.
        """
        pass

    def close(self) -> None:
        """Closes the writer and finalizes the output file.

        This method should be called explicitly when done writing,
        especially if not using a `with` statement context manager
        (which is not directly implemented here but recommended in Python usage).
        It ensures all buffered data is flushed to the file.
        """
        pass

class BsxFileReader:
    """
    A reader for .bsx files based on Apache Arrow IPC format, optimized for
    reading batches in parallel.

    This reader is optimized for parallel access to data batches within a .bsx file,
    leveraging memory mapping and multithreading for efficient data access.
    """
    def __init__(self, file: Union[str, BinaryIO]) -> None:
        """
        Creates a new BsxFileReader from a file handle.

        This will memory map the file and read its metadata and dictionaries.
        It also initializes a thread pool and thread-local handles for parallel
        reading.

        Parameters
        ----------
        file : str or file-like object
            Path to the BSX file or a file-like object.
        """
        pass

    def get_batch(self, batch_idx: int) -> Optional[BsxBatch]:
        """
        Reads a single batch from the file at the given index.

        Parameters
        ----------
        batch_idx : int
            The index of the batch to retrieve.

        Returns
        -------
        BsxBatch or None
            The requested batch, or None if the index is out of bounds.
        """
        pass

    def get_batches(self, batch_indices: List[int]) -> List[Optional[BsxBatch]]:
        """
        Reads multiple batches from the file at the given indices in parallel.

        Parameters
        ----------
        batch_indices : list[int]
            List of batch indices to retrieve.

        Returns
        -------
        list[BsxBatch or None]
            List of batches, with None for indices that are out of bounds.
        """
        pass

    def cache_batches(self, batch_indices: List[int]) -> None:
        """
        Reads multiple batches from the file at the given indices in parallel
        and adds them to the internal cache.

        Parameters
        ----------
        batch_indices : list[int]
            List of batch indices to cache.
        """
        pass

    def n_threads(self) -> int:
        """
        Returns the number of threads used by the internal thread pool.

        Returns
        -------
        int
            Number of threads.
        """
        pass

    def blocks_total(self) -> int:
        """
        Returns the total number of data batches (blocks) in the file.

        Returns
        -------
        int
            The total number of batches.
        """
        pass

    def __iter__(self) -> "BsxFileReader":
        pass

    def __next__(self) -> Optional[BsxBatch]:
        pass

class IpcCompression:
    """Compression algorithms for IPC files."""

    LZ4: "IpcCompression"
    ZSTD: "IpcCompression"

class BsxFileWriter:
    """
    Writer for BSX data in Arrow IPC format with optional compression.

    A writer for creating .bsx files from data batches, supporting different
    data sources for schema information (e.g., FASTA index) and optional compression.
    """
    def __init__(
        self,
        sink: Union[str, BinaryIO],
        chr_names: List[str],
        compression: Optional[IpcCompression] = None,
    ) -> None:
        """
        Creates a new BSX IPC writer.

        Parameters
        ----------
        sink : str or file-like object
            Path or file-like object to write to.
        chr_names : list[str]
            List of chromosome names.
        compression : IpcCompression, optional
            Compression algorithm to use (LZ4, ZSTD). Defaults to None (uncompressed).
        """
        pass

    @staticmethod
    def from_sink_and_fai(
        sink: Union[str, BinaryIO],
        fai_path: Union[str, Path],
        compression: Optional[IpcCompression] = None,
    ) -> "BsxFileWriter":
        """
        Creates a writer from a sink and FASTA index file.

        Parameters
        ----------
        sink : str or file-like object
            Path or file-like object to write to.
        fai_path : str or Path
            Path to the FASTA index file.
        compression : IpcCompression, optional
            Compression algorithm (LZ4, ZSTD). Defaults to None.

        Returns
        -------
        BsxFileWriter
            A new writer instance.
        """
        pass

    @staticmethod
    def from_sink_and_fasta(
        sink: Union[str, BinaryIO],
        fasta_path: Union[str, Path],
        compression: Optional[IpcCompression] = None,
    ) -> "BsxFileWriter":
        """
        Creates a writer from a sink and FASTA file.

        Parameters
        ----------
        sink : str or file-like object
            Path or file-like object to write to.
        fasta_path : str or Path
            Path to the FASTA file.
        compression : IpcCompression, optional
            Compression algorithm (LZ4, ZSTD). Defaults to None.

        Returns
        -------
        BsxFileWriter
            A new writer instance.
        """
        pass

    def write_batch(self, batch: BsxBatch) -> None:
        """
        Writes an encoded BSX batch.

        Parameters
        ----------
        batch : BsxBatch
            The BSX batch to encode and write.
        """
        pass

    def close(self) -> None:
        """
        Finalizes the IPC file and closes the writer.

        Ensures resources are properly cleaned up when dropped.
        """
        pass

    def __enter__(self) -> "BsxFileWriter":
        pass

    def __exit__(
        self,
        exc_type: Optional[type],
        exc_value: Optional[BaseException],
        traceback: Optional[object],
    ) -> None:
        pass
