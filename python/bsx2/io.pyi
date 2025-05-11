from typing import Optional, Union, List, BinaryIO, Callable, Any
from pathlib import Path

from .types import ReportTypeSchema, BsxBatch, Contig, EncodedBsxBatch, BatchIndex
import polars as pl

__all__ = [
    "RegionReader",
    "Compression",
    "ReportReader",
    "ReportWriter",
    "BsxFileReader",
    "IpcCompression",
    "BsxIpcWriter"
]

class RegionReader:
    """
    Reads BSX files efficiently by caching relevant batches for a genomic region.

    This reader is designed for quickly retrieving data within specific
    genomic coordinates, suitable for interactive exploration or querying.
    It maintains an internal cache of decoded batches that overlap the
    most recently queried regions.
    """
    def __init__(self, file: Union[str, BinaryIO]) -> None:
        """
        Creates a new RegionReader instance.

        Initializes the reader and builds an index for the BSX file.
        The index is built upon creation, which might take time for large files.

        Parameters
        ----------
        path : str or Path
            Path to the BSX file.

        Raises
        ------
        FileNotFoundError
            If the input BSX file cannot be found.
        RuntimeError
            If there's an error reading the file or building the index.
        """
        pass

    def query(self, contig: Contig, postprocess_fn: Callable[[EncodedBsxBatch], EncodedBsxBatch]) -> Optional[pl.DataFrame]:
        """
        Queries the BSX file for data within a specific genomic region.

        This method uses the internal index and cache to efficiently retrieve
        data for the specified chromosome and coordinate range. It will read
        and cache batches as needed to fulfill the query.

        Parameters
        ----------
        seqname : str
            The name of the chromosome.
        start : int
            The start position of the region (inclusive, 0-based).
        end : int
            The end position of the region (exclusive, 0-based).

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

    def set_preprocess_fn(self, preprocess_fn: Callable[[Any], Any]) -> None:
        """
        Sets the preprocessing function applied to each batch before it is cached.

        Parameters
        ----------
        preprocess_fn : callable
            A function that takes an EncodedBsxBatch and returns a processed EncodedBsxBatch

        Raises
        ------
        ValueError
            If preprocess_fn is not callable
        """
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

class Compression:
    """Compression algorithms."""
    GZIP: str
    ZSTD: str
    BGZIP: str
    XZ: str
    NONE: str

class ReportReader:
    """Reads methylation report files batch by batch.

    This class provides an iterator interface to read various methylation report
    formats (like Bismark, CGmap, BedGraph, Coverage) efficiently.
    It handles file parsing, optional alignment with FASTA context, and
    yields data in standardized BsxBatch objects.
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
        queue_len: int = 1000,
        compression: Optional[str] = None,
        compression_level: Optional[int] = None,
    ) -> None:
        """Creates a new ReportReader instance.

        Parameters
        ----------
        path : str
            Path to the input report file.
        report_type : ReportTypeSchema
            The format of the report file (e.g., ReportTypeSchema.Bismark).
        chunk_size : int, optional
            The target number of records per yielded batch. Defaults to 10000.
        fasta_path : str, optional
            Path to the reference FASTA file. Required for BedGraph/Coverage formats
            or when alignment is needed.
        fai_path : str, optional
            Path to the FASTA index file (.fai). Can be used instead of `fasta_path`
            for determining chromosome order/names if FASTA content is not needed
            for alignment.
        batch_size : int, optional
            The number of lines read from the file at once internally. Defaults to 100000.
        n_threads : int, optional
            Number of threads to use for parsing. Defaults to using available cores.
        low_memory : bool, optional
            Whether to use a low-memory parsing mode. Defaults to False.
        queue_len : int, optional
            Internal buffer size for parsed batches. Defaults to 1000.
        compression : str, optional
            Compression type ('gzip', 'zstd', 'bgzip', 'xz'). Defaults to None (autodetect).
        compression_level : int, optional
            Compression level if applicable. Defaults depend on the compression type.

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
        BsxBatch
            The next batch of data.

        Raises
        ------
        StopIteration
            When there are no more batches to read.
        RuntimeError
            If an error occurs during reading or processing a batch.
        """
        pass

class ReportWriter:
    """Writes methylation report data to a file.

    This class handles writing BsxBatch objects or Polars DataFrames
    to a specified file path in a chosen report format (Bismark, CGmap, etc.).
    It manages file handling, formatting, and optional compression.
    """
    def __init__(
        self,
        file: Union[str, BinaryIO],
        schema: ReportTypeSchema,
        n_threads: int = 1,
        compression: Optional[str] = None,
        compression_level: Optional[int] = None,
    ) -> None:
        """Creates a new ReportWriter instance.

        Parameters
        ----------
        path : str
            Path to the output file.
        schema : ReportTypeSchema
            The desired output format schema (e.g., ReportTypeSchema.Bismark).
        n_threads : int, optional
            Number of threads to use for writing. Defaults to 1.
        compression : str, optional
            Compression type ('gzip', 'zstd', 'bgzip', 'xz'). Defaults to None.
        compression_level : int, optional
            Compression level (1-21 for zstd, 1-9 for gzip/bgzip, 0-9 for xz).
            Defaults vary by type (e.g., 3 for zstd, 6 for gzip).

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
    """Reader for BSX files.

    Parameters
    ----------
    file : str or file-like object
        Path to the BSX file or a file-like object.
    """
    def __init__(self, file: Union[str, BinaryIO]) -> None:
        pass

    @staticmethod
    def from_file_and_index(file: Union[str, BinaryIO], index: Union[str, BinaryIO]) -> "BsxFileReader":
        """Create a reader from a file and an index.

        Parameters
        ----------
        file : str or file-like object
            Path to the BSX file or a file-like object.
        index : str or file-like object
            Path to the index file or a file-like object.

        Returns
        -------
        BsxFileReader
            A new reader instance.
        """
        pass

    def get_batch(self, batch_idx: int) -> Optional[pl.DataFrame]:
        """Get a specific batch by index.

        Parameters
        ----------
        batch_idx : int
            The index of the batch to retrieve.

        Returns
        -------
        DataFrame or None
            The requested batch as a Polars DataFrame, or None if the index is out of bounds.
        """
        pass

    def blocks_total(self) -> int:
        """Get the total number of batches in the file.

        Returns
        -------
        int
            The total number of batches.
        """
        pass

    def index(self) -> BatchIndex:
        """Get the index of the BSX file.

        Returns
        -------
        BatchIndex
            The index of the BSX file

        Raises
        ------
        RuntimeError
            If there's an error reading the index.
        """
        pass

    def set_index(self, index: BatchIndex) -> None:
        """Set the index of the BSX file.

        Parameters
        ----------
        index : BatchIndex
            The index to set.
        """
        pass

    def query(self, query: Contig) -> Optional[EncodedBsxBatch]:
        """Query the BSX file for data within a specific genomic region.

        Parameters
        ----------
        query : Contig
            The genomic region to query.

        Returns
        -------
        Optional[EncodedBsxBatch]
            The batch data for the specified region, or None if not found.

        Raises
        ------
        RuntimeError
            If there's an error during the query.
        """
        pass

    def __iter__(self) -> "BsxFileReader":
        pass

    def __next__(self) -> Optional[pl.DataFrame]:
        pass

class IpcCompression:
    """Compression algorithms for IPC."""
    LZ4: str
    ZSTD: str

class BsxIpcWriter:
    """Writer for BSX data in Arrow IPC format.

    Parameters
    ----------
    sink : str or file-like object
        Path or file-like object to write to.
    chr_names : list[str]
        List of chromosome names.
    compression : str, optional
        Compression algorithm to use ('uncompressed', 'lz4', 'zstd'). Defaults to None (uncompressed).
    custom_metadata : bytes, optional
        Custom metadata to embed in the IPC file.
    """
    def __init__(
        self,
        sink: Union[str, BinaryIO],
        chr_names: List[str],
        compression: Optional[IpcCompression] = None,
    ) -> None:
        pass

    @staticmethod
    def from_sink_and_fai(
        sink: Union[str, BinaryIO],
        fai_path: Union[str, Path],
        compression: Optional[IpcCompression] = None,
    ) -> "BsxIpcWriter":
        """Create a writer using a FASTA index file (.fai) to get chromosome names.

        Parameters
        ----------
        sink : str or file-like object
            Path or file-like object to write to.
        fai_path : str
            Path to the FASTA index file.
        compression : str, optional
            Compression algorithm ('uncompressed', 'lz4', 'zstd'). Defaults to None.
        custom_metadata : bytes, optional
            Custom metadata.

        Returns
        -------
        BsxIpcWriter
            A new writer instance.
        """
        pass

    @staticmethod
    def from_sink_and_fasta(
        sink: Union[str, BinaryIO],
        fasta_path: Union[str, Path],
        compression: Optional[IpcCompression] = None,
    ) -> "BsxIpcWriter":
        """Create a writer using a FASTA file to get chromosome names (will index if needed).

        Parameters
        ----------
        sink : str or file-like object
            Path or file-like object to write to.
        fasta_path : str
            Path to the FASTA file.
        compression : str, optional
            Compression algorithm ('uncompressed', 'lz4', 'zstd'). Defaults to None.
        custom_metadata : bytes, optional
            Custom metadata.

        Returns
        -------
        BsxIpcWriter
            A new writer instance.
        """
        pass

    def write_encoded_batch(self, batch: pl.DataFrame) -> None:
        """Write an already encoded BSX batch (DataFrame).

        Parameters
        ----------
        batch : DataFrame
            The encoded BSX batch (Polars DataFrame) to write.
        """
        pass

    def write_batch(self, batch: pl.DataFrame) -> None:
        """Encode and write a BSX batch (DataFrame).

        Parameters
        ----------
        batch : DataFrame
            The BSX batch (Polars DataFrame) to encode and write.
        """
        pass

    def close(self) -> None:
        """Finalize the IPC file and close the writer."""
        pass

    def __enter__(self) -> "BsxIpcWriter":
        pass

    def __exit__(self, exc_type: Optional[type], exc_value: Optional[BaseException], traceback: Optional[object]) -> None:
        pass
