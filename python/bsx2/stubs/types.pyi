import polars as pl
from typing import Optional, List

class EncodedBsxBatch:
    """
    Represents a batch of BS-Seq data in an encoded format.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        check_nulls: bool = True,
        check_sorted: bool = True,
        check_duplicates: bool = True,
        rechunk: bool = True,
        check_single_chr: bool = True,
        context_data: ContextData = None,
        report_schema: ReportTypeSchema = None,
    ) -> None:
        """
        Create a new EncodedBsxBatch from a DataFrame.

        Performs various checks to ensure data integrity unless disabled.

        Parameters
        ----------
        data : pl.DataFrame
            Input Polars DataFrame containing BS-Seq data.
        check_nulls : bool, optional
            If True (default), check for null values in essential columns.
        check_sorted : bool, optional
            If True (default), check if the data is sorted by position.
        check_duplicates : bool, optional
            If True (default), check for duplicate entries (chr, pos, strand).
        rechunk : bool, optional
            If True (default), rechunk the DataFrame for optimal performance.
        check_single_chr : bool, optional
            If True (default), ensure the batch contains data for only one chromosome.
        context_data : ContextData, optional
            Context data associated with the batch.
        report_schema : ReportTypeSchema, optional
            Schema defining the report type.

        Raises
        ------
        ValueError
            If any of the enabled checks fail.
        """
        ...

    @classmethod
    def from_dataframe_unchecked(cls, data: pl.DataFrame) -> "EncodedBsxBatch":
        """
        Create a new EncodedBsxBatch from a DataFrame without performing any checks.

        Warning: Use with caution, assumes the input DataFrame is valid and correctly formatted.

        Parameters
        ----------
        data : pl.DataFrame
            Input Polars DataFrame.

        Returns
        -------
        EncodedBsxBatch
            A new instance of EncodedBsxBatch.
        """
        ...

    @staticmethod
    def schema() -> pl.Schema:
        """
        Get the Polars schema for an EncodedBsxBatch.

        Returns
        -------
        pl.Schema
            The Polars schema definition.
        """
        ...

    @staticmethod
    def empty() -> "EncodedBsxBatch":
        """
        Create an empty EncodedBsxBatch with the correct schema.

        Returns
        -------
        EncodedBsxBatch
            An empty batch instance.
        """
        ...

    @property
    def chr(self) -> pl.Series:
        """
        Access the chromosome column as a Polars Series.

        Returns
        -------
        pl.Series
            The chromosome data (Categorical).
        """
        ...

    @property
    def chr_categorical(self) -> pl.Series:
        """
        Access the chromosome column as a Categorical Series.

        Returns
        -------
        pl.Series
            The chromosome data as a Polars Categorical Series.

        Raises
        ------
        ValueError
            If the 'chr' column is not found or is not categorical.
        """
        ...

    @property
    def position(self) -> pl.Series:
        """
        Access the position column as a Polars Series.

        Returns
        -------
        pl.Series
            The position data (UInt32).
        """
        ...

    @property
    def strand(self) -> pl.Series:
        """
        Access the strand column as a Polars Series.

        Returns
        -------
        pl.Series
            The strand data (Categorical).
        """
        ...

    @property
    def context(self) -> pl.Series:
        """
        Access the context column as a Polars Series.

        Returns
        -------
        pl.Series
            The context data (Categorical).
        """
        ...

    @property
    def count_m(self) -> pl.Series:
        """
        Access the methylated counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The methylated count data (UInt32).
        """
        ...

    @property
    def count_total(self) -> pl.Series:
        """
        Access the total counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The total count data (UInt32).
        """
        ...

    @property
    def density(self) -> pl.Series:
        """
        Access the density column as a Polars Series.

        Returns
        -------
        pl.Series
            The density data (Float32).
        """
        ...

    def is_empty(self) -> bool:
        """
        Check if the batch is empty.

        Returns
        -------
        bool
            True if the batch contains no rows, False otherwise.
        """
        ...

    def split_at(self, index: int) -> tuple["EncodedBsxBatch", "EncodedBsxBatch"]:
        """
        Split the batch into two at the given index.

        Parameters
        ----------
        index : int
            The row index at which to split the batch.

        Returns
        -------
        tuple[EncodedBsxBatch, EncodedBsxBatch]
            A tuple containing two new batches, the first with rows up to `index`,
            the second with rows from `index` onwards.
        """
        ...

    @property
    def data(self) -> pl.DataFrame:
        """
        Get a reference to the underlying Polars DataFrame.

        Returns
        -------
        pl.DataFrame
            A clone of the internal DataFrame.
        """
        ...

    @property
    def take(self) -> pl.DataFrame:
        """
        Consume the batch and return the underlying Polars DataFrame.

        Returns
        -------
        pl.DataFrame
            The internal DataFrame.
        """
        ...

    def chr_val(self) -> str:
        """
        Get the unique chromosome value present in the batch.

        Assumes the batch contains data for only a single chromosome.

        Returns
        -------
        str
            The chromosome name.

        Raises
        ------
        ValueError
            If the batch is empty or contains multiple chromosomes.
        """
        ...

    def start_pos(self) -> Optional[int]:
        """
        Get the starting genomic position in the batch.

        Returns
        -------
        int or None
            The minimum position value, or None if the batch is empty.
        """
        ...

    def end_pos(self) -> Optional[int]:
        """
        Get the ending genomic position in the batch.

        Returns
        -------
        int or None
            The maximum position value, or None if the batch is empty.
        """
        ...

    def vstack(self, other: "EncodedBsxBatch") -> "EncodedBsxBatch":
        """
        Vertically stack this batch with another batch.

        Creates a new batch containing rows from both batches. Performs checks
        on the combined data.

        Parameters
        ----------
        other : EncodedBsxBatch
            The batch to stack underneath this one.

        Returns
        -------
        EncodedBsxBatch
            A new batch containing the combined data.

        Raises
        ------
        ValueError
            If the stacking operation results in invalid data (e.g., unsorted, duplicates).
        """
        ...

    def extend(self, other: "EncodedBsxBatch") -> None:
        """
        Extend this batch by appending rows from another batch in-place.

        Modifies the current batch. Performs checks after extending.

        Parameters
        ----------
        other : EncodedBsxBatch
            The batch whose rows will be appended.

        Raises
        ------
        ValueError
            If the extend operation results in invalid data (e.g., unsorted, duplicates).
        """
        ...

    def filter_mask(self, mask: pl.Series) -> "EncodedBsxBatch":
        """
        Filter the batch using a boolean mask Series.

        Creates a new batch containing only the rows where the mask is True.

        Parameters
        ----------
        mask : pl.Series
            A Polars Series of boolean values with the same length as the batch.

        Returns
        -------
        EncodedBsxBatch
            A new batch containing the filtered rows.

        Raises
        ------
        TypeError
            If the mask is not a boolean Series.
        ValueError
            If the mask length does not match the batch height or other Polars errors occur.
        """
        ...

    def height(self) -> int:
        """
        Get the number of rows in the batch.

        Returns
        -------
        int
            The height (number of rows) of the batch.
        """
        ...

    def __len__(self) -> int:
        """
        Get the number of rows in the batch (implements Python's `len()`).

        Returns
        -------
        int
            The height (number of rows) of the batch.
        """
        ...

    def __repr__(self) -> str:
        """
        Get a string representation of the batch.

        Returns
        -------
        str
            A string showing shape and chromosome.
        """
        ...

    def __richcmp__(self, other: "EncodedBsxBatch", op: int) -> bool:
        """
        Compare two batches for equality (implements `==` and `!=`).

        Parameters
        ----------
        other : EncodedBsxBatch
            The batch to compare against.
        op : int
            The comparison operation (Eq or Ne).

        Returns
        -------
        bool
            True if the batches are equal/unequal based on `op`, False otherwise.

        Raises
        ------
        NotImplementedError
            If the comparison operation is not `==` or `!=`.
        """
        ...

class BsxBatch:
    """
    Represents a batch of BS-Seq data.
    """

    def __init__(
        self,
        data: pl.DataFrame,
        check_nulls: bool = True,
        check_sorted: bool = True,
        check_duplicates: bool = True,
        rechunk: bool = True,
        check_single_chr: bool = True,
        context_data: ContextData = None,
        report_schema: ReportTypeSchema = None,
    ) -> None:
        """
        Create a new BsxBatch from a DataFrame.

        Performs various checks to ensure data integrity unless disabled.

        Parameters
        ----------
        data : pl.DataFrame
            Input Polars DataFrame containing BS-Seq data.
        check_nulls : bool, optional
            If True (default), check for null values in essential columns.
        check_sorted : bool, optional
            If True (default), check if the data is sorted by position.
        check_duplicates : bool, optional
            If True (default), check for duplicate entries (chr, pos, strand).
        rechunk : bool, optional
            If True (default), rechunk the DataFrame for optimal performance.
        check_single_chr : bool, optional
            If True (default), ensure the batch contains data for only one chromosome.
        context_data : ContextData, optional
            Context data associated with the batch.
        report_schema : ReportTypeSchema, optional
            Schema defining the report type.

        Raises
        ------
        ValueError
            If any of the enabled checks fail.
        """
        ...

    @classmethod
    def from_dataframe_unchecked(cls, data: pl.DataFrame) -> "BsxBatch":
        """
        Create a new BsxBatch from a DataFrame without performing any checks.

        Warning: Use with caution, assumes the input DataFrame is valid and correctly formatted.

        Parameters
        ----------
        data : pl.DataFrame
            Input Polars DataFrame.

        Returns
        -------
        BsxBatch
            A new instance of BsxBatch.
        """
        ...

    @staticmethod
    def schema() -> pl.Schema:
        """
        Get the Polars schema for a BsxBatch.

        Returns
        -------
        pl.Schema
            The Polars schema definition.
        """
        ...

    @staticmethod
    def empty() -> "BsxBatch":
        """
        Create an empty BsxBatch with the correct schema.

        Returns
        -------
        BsxBatch
            An empty batch instance.
        """
        ...

    @property
    def chr(self) -> pl.Series:
        """
        Access the chromosome column as a Polars Series.

        Returns
        -------
        pl.Series
            The chromosome data (Utf8).
        """
        ...

    @property
    def position(self) -> pl.Series:
        """
        Access the position column as a Polars Series.

        Returns
        -------
        pl.Series
            The position data (UInt32).
        """
        ...

    @property
    def strand(self) -> pl.Series:
        """
        Access the strand column as a Polars Series.

        Returns
        -------
        pl.Series
            The strand data (Utf8).
        """
        ...

    @property
    def context(self) -> pl.Series:
        """
        Access the context column as a Polars Series.

        Returns
        -------
        pl.Series
            The context data (Utf8).
        """
        ...

    @property
    def count_m(self) -> pl.Series:
        """
        Access the methylated counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The methylated count data (UInt32).
        """
        ...

    @property
    def count_total(self) -> pl.Series:
        """
        Access the total counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The total count data (UInt32).
        """
        ...

    @property
    def density(self) -> pl.Series:
        """
        Access the density column as a Polars Series.

        Returns
        -------
        pl.Series
            The density data (Float32).
        """
        ...

    def is_empty(self) -> bool:
        """
        Check if the batch is empty.

        Returns
        -------
        bool
            True if the batch contains no rows, False otherwise.
        """
        ...

    def split_at(self, index: int) -> tuple["BsxBatch", "BsxBatch"]:
        """
        Split the batch into two at the given index.

        Parameters
        ----------
        index : int
            The row index at which to split the batch.

        Returns
        -------
        tuple[BsxBatch, BsxBatch]
            A tuple containing two new batches, the first with rows up to `index`,
            the second with rows from `index` onwards.
        """
        ...

    @property
    def data(self) -> pl.DataFrame:
        """
        Get a reference to the underlying Polars DataFrame.

        Returns
        -------
        pl.DataFrame
            A clone of the internal DataFrame.
        """
        ...

    def take(self) -> pl.DataFrame:
        """
        Consume the batch and return the underlying Polars DataFrame.

        Returns
        -------
        pl.DataFrame
            The internal DataFrame.
        """
        ...

    def chr_val(self) -> str:
        """
        Get the unique chromosome value present in the batch.

        Assumes the batch contains data for only a single chromosome.

        Returns
        -------
        str
            The chromosome name.

        Raises
        ------
        ValueError
            If the batch is empty or contains multiple chromosomes.
        """
        ...

    def start_pos(self) -> Optional[int]:
        """
        Get the starting genomic position in the batch.

        Returns
        -------
        int or None
            The minimum position value, or None if the batch is empty.
        """
        ...

    def end_pos(self) -> Optional[int]:
        """
        Get the ending genomic position in the batch.

        Returns
        -------
        int or None
            The maximum position value, or None if the batch is empty.
        """
        ...

    def vstack(self, other: "BsxBatch") -> "BsxBatch":
        """
        Vertically stack this batch with another batch.

        Creates a new batch containing rows from both batches. Performs checks
        on the combined data.

        Parameters
        ----------
        other : BsxBatch
            The batch to stack underneath this one.

        Returns
        -------
        BsxBatch
            A new batch containing the combined data.

        Raises
        ------
        ValueError
            If the stacking operation results in invalid data (e.g., unsorted, duplicates).
        """
        ...

    def extend(self, other: "BsxBatch") -> None:
        """
        Extend this batch by appending rows from another batch in-place.

        Modifies the current batch. Performs checks after extending.

        Parameters
        ----------
        other : BsxBatch
            The batch whose rows will be appended.

        Raises
        ------
        ValueError
            If the extend operation results in invalid data (e.g., unsorted, duplicates).
        """
        ...

    def filter_mask(self, mask: pl.Series) -> "BsxBatch":
        """
        Filter the batch using a boolean mask Series.

        Creates a new batch containing only the rows where the mask is True.

        Parameters
        ----------
        mask : pl.Series
            A Polars Series of boolean values with the same length as the batch.

        Returns
        -------
        BsxBatch
            A new batch containing the filtered rows.

        Raises
        ------
        TypeError
            If the mask is not a boolean Series.
        ValueError
            If the mask length does not match the batch height or other Polars errors occur.
        """
        ...

    def height(self) -> int:
        """
        Get the number of rows in the batch.

        Returns
        -------
        int
            The height (number of rows) of the batch.
        """
        ...

    def __len__(self) -> int:
        """
        Get the number of rows in the batch (implements Python's `len()`).

        Returns
        -------
        int
            The height (number of rows) of the batch.
        """
        ...

    def __repr__(self) -> str:
        """
        Get a string representation of the batch.

        Returns
        -------
        str
            A string showing shape and chromosome.
        """
        ...

    def __richcmp__(self, other: "BsxBatch", op: int) -> bool:
        """
        Compare two batches for equality (implements `==` and `!=`).

        Parameters
        ----------
        other : BsxBatch
            The batch to compare against.
        op : int
            The comparison operation (Eq or Ne).

        Returns
        -------
        bool
            True if the batches are equal/unequal based on `op`, False otherwise.

        Raises
        ------
        NotImplementedError
            If the comparison operation is not `==` or `!=`.
        """
        ...

class ReportTypeSchema:
    """
    Represents different input/output file formats for methylation data.
    """

    Bismark: "ReportTypeSchema"
    CgMap: "ReportTypeSchema"
    BedGraph: "ReportTypeSchema"
    Coverage: "ReportTypeSchema"

    def col_names(self) -> List[str]:
        """
        Get the list of column names for this report format.

        Returns
        -------
        list[str]
            A list of column names.
        """
        ...

    def schema(self) -> pl.Schema:
        """
        Get the Polars schema for this report format.

        Returns
        -------
        pl.Schema
            The Polars schema definition.
        """
        ...

    def chr_col(self) -> str:
        """
        Get the name of the chromosome column for this format.

        Returns
        -------
        str
            The chromosome column name.
        """
        ...

    def position_col(self) -> str:
        """
        Get the name of the position column for this format.

        Notes
        -----
        For some formats like BedGraph, this might be 'start'.

        Returns
        -------
        str
            The position column name.
        """
        ...

    def context_col(self) -> Optional[str]:
        """
        Get the name of the context column, if available for this format.

        Returns
        -------
        Optional[str]
            The context column name or None if not applicable.
        """
        ...

    def strand_col(self) -> Optional[str]:
        """
        Get the name of the strand column, if available for this format.

        Returns
        -------
        Optional[str]
            The strand column name or None if not applicable.
        """
        ...

    def need_align(self) -> bool:
        """
        Check if this report format typically requires alignment with context data.

        Notes
        -----
        Formats like BedGraph or Coverage often lack explicit strand/context
        information and benefit from alignment with genomic context.

        Returns
        -------
        bool
            True if alignment is generally needed, False otherwise.
        """
        ...

class Strand:
    """
    Represents DNA strand information.
    """

    Forward: "Strand"
    Reverse: "Strand"
    None_: "Strand"  # Renamed to avoid conflict with Python's None

class Context:
    """
    Represents methylation context (sequence context).
    """

    CG: "Context"
    CHG: "Context"
    CHH: "Context"

class ContextData:
    """
    Stores methylation context information derived from a DNA sequence.

    This class pre-calculates the positions, strands, and contexts (CG, CHG, CHH)
    for cytosines within a given sequence, allowing for efficient alignment
    with methylation data.
    """

    def __init__(self, sequence: bytes) -> None:
        """
        Creates a new ContextData instance from a DNA sequence.

        Parameters
        ----------
        sequence : bytes
            The DNA sequence as bytes (e.g., b"ACGTACGT").

        Returns
        -------
        ContextData
            A new instance containing context information.
        """
        ...

    def to_decoded_df(self) -> pl.DataFrame:
        """
        Converts the context information into a Polars DataFrame compatible with BsxBatch (decoded format).

        The resulting DataFrame contains 'chr', 'position', 'strand', and 'context' columns
        suitable for joining or aligning with decoded methylation data.

        Returns
        -------
        pl.DataFrame
            A Polars DataFrame with context information.
        """
        ...

    def to_encoded_df(self) -> pl.DataFrame:
        """
        Converts the context information into a Polars DataFrame compatible with EncodedBsxBatch.

        The resulting DataFrame contains 'chr', 'position', 'strand', and 'context' columns
        with types suitable for joining or aligning with encoded methylation data (e.g., boolean strand/context).

        Returns
        -------
        pl.DataFrame
            A Polars DataFrame with encoded context information.
        """
        ...
