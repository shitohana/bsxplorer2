from typing import Optional, Union, List
from pathlib import Path

from bsx2.types import ReportTypeSchema, BsxBatch
import polars as pl

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
        ...

    def __iter__(self) -> "ReportReader":
        ...

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
        ...

class ReportWriter:
    """Writes methylation report data to a file.

    This class handles writing BsxBatch objects or Polars DataFrames
    to a specified file path in a chosen report format (Bismark, CGmap, etc.).
    It manages file handling, formatting, and optional compression.
    """
    def __init__(
        self,
        path: Union[str, Path],
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
        ...

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
        ...

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
        ...

    def close(self) -> None:
        """Closes the writer and finalizes the output file.

        This method should be called explicitly when done writing,
        especially if not using a `with` statement context manager
        (which is not directly implemented here but recommended in Python usage).
        It ensures all buffered data is flushed to the file.
        """
        ...

class BsxFileReader:
    """Reader for BSX files.

    Parameters
    ----------
    file : str or file-like object
        Path to the BSX file or a file-like object.
    """
    def __init__(self, file: Union[str, Path]) -> None:
        ...

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
        ...

    def blocks_total(self) -> int:
        """Get the total number of batches in the file.

        Returns
        -------
        int
            The total number of batches.
        """
        ...

    def __iter__(self) -> "BsxFileReader":
        ...

    def __next__(self) -> Optional[pl.DataFrame]:
        ...

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
        sink: Union[str, Path],
        chr_names: List[str],
        compression: Optional[IpcCompression] = None,
    ) -> None:
        ...

    @staticmethod
    def from_sink_and_fai(
        sink: Union[str, Path],
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
        ...

    @staticmethod
    def from_sink_and_fasta(
        sink: Union[str, Path],
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
        ...

    def write_encoded_batch(self, batch: pl.DataFrame) -> None:
        """Write an already encoded BSX batch (DataFrame).

        Parameters
        ----------
        batch : DataFrame
            The encoded BSX batch (Polars DataFrame) to write.
        """
        ...

    def write_batch(self, batch: pl.DataFrame) -> None:
        """Encode and write a BSX batch (DataFrame).

        Parameters
        ----------
        batch : DataFrame
            The BSX batch (Polars DataFrame) to encode and write.
        """
        ...

    def close(self) -> None:
        """Finalize the IPC file and close the writer."""
        ...

    def __enter__(self) -> "BsxIpcWriter":
        ...

    def __exit__(self, exc_type: Optional[type], exc_value: Optional[BaseException], traceback: Optional[object]) -> None:
        ...


