from enum import Enum
from typing import Callable, Dict, Optional, List, Tuple

import polars as pl

__all__ = [
    "ReportTypeSchema",
    "Strand",
    "Context",
    "ContextData",
    "GenomicPosition",
    "Contig",
    "BatchIndex",
    "GffEntry",
    "AnnotStore",
    "MethylationStats",
    "EncodedBsxBatch",
    "BsxBatch",
    "LazyBsxBatch",
    "LazyEncodedBsxBatch"
]

class ReportTypeSchema(Enum):
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
    Null: "Strand"  # Renamed to avoid conflict with Python's None

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

    @staticmethod
    def empty() -> "ContextData":
        """
        Creates an empty `ContextData`.
        """
        ...

    def is_empty(self) -> bool:
        """
        Checks if the `ContextData` is empty.
        """
        ...

    def take(self) -> Tuple[List[int], List[Strand], List[Context]]:
        """
        Takes the data, returning vectors of positions, strands, and contexts.
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

    def __len__(self) -> int:
        """
        Returns the number of entries in the `ContextData`.
        """
        ...


class GenomicPosition:
    """
    Represents a genomic position on a chromosome.

    Attributes
    ----------
    seqname : str
        The name of the chromosome or sequence.
    position : int
        The 0-based or 1-based genomic position. The interpretation depends
        on the convention used elsewhere (often 0-based in bioinformatics).
    """

    def __init__(self, seqname: str, position: int) -> None:
        """
        Create a new GenomicPosition.

        Parameters
        ----------
        seqname : str
            The name of the chromosome or sequence.
        position : int
            The genomic position.
        """
        ...

    @property
    def seqname(self) -> str:
        """
        Get the sequence name.
        """
        ...

    @property
    def position(self) -> int:
        """
        Get the position.
        """
        ...

    def __str__(self) -> str:
        """
        Get the string representation (seqname:position).
        """
        ...

    def __repr__(self) -> str:
        """
        Get the developer representation.
        """
        ...

    def __lt__(self, other: "GenomicPosition") -> bool:
        ...

    def __le__(self, other: "GenomicPosition") -> bool:
        ...

    def __richcmp__(self, other: "GenomicPosition", op: int) -> bool:
        ...

    def __gt__(self, other: "GenomicPosition") -> bool:
        ...

    def __ge__(self, other: "GenomicPosition") -> bool:
        ...

    def is_zero(self) -> bool:
        """
        Checks if the `GenomicPosition` is zero (empty seqname and zero position).
        """
        ...

    def __add__(self, other: "GenomicPosition") -> Optional["GenomicPosition"]:
        """
        Add positions. Only possible if seqnames are identical.

        Parameters
        ----------
        other : GenomicPosition
            The position to add.

        Returns
        -------
        GenomicPosition or None
            A new GenomicPosition representing the sum, or None if seqnames differ.
        """
        ...

    def __sub__(self, other: "GenomicPosition") -> Optional["GenomicPosition"]:
        """
        Subtract positions. Only possible if seqnames are identical and self >= other.

        Parameters
        ----------
        other : GenomicPosition
            The position to subtract.

        Returns
        -------
        GenomicPosition or None
            A new GenomicPosition representing the difference, or None if seqnames
            differ or the other position is greater.
        """
        ...


class Contig:
    """
    Represents a contiguous genomic region with an optional strand.

    Attributes
    ----------
    seqname : str
        The name of the chromosome or sequence.
    start : int
        The start position of the region (inclusive).
    end : int
        The end position of the region (exclusive).
    """

    def __init__(self, seqname: str, start: int, end: int, strand: str = ".") -> None:
        """
        Create a new Contig.

        Parameters
        ----------
        seqname : str
            The name of the chromosome or sequence.
        start : int
            The start position (inclusive).
        end : int
            The end position (exclusive).
        strand : str, optional
            The strand ('+', '-', or '.', case-insensitive). Defaults to '.'.

        Raises
        ------
        ValueError
            If start > end or the strand value is invalid.
        """
        ...

    @property
    def seqname(self) -> str:
        """
        Get the sequence name.
        """
        ...

    @property
    def start(self) -> int:
        """
        Get the start position (inclusive).
        """
        ...

    @property
    def end(self) -> int:
        """
        Get the end position (exclusive).
        """
        ...

    @property
    def strand(self) -> Strand: # In Rust PyO3, it's PyStrand which is an enum.
        """
        Get the strand as a Strand enum value.
        """
        ...

    @property
    def strand_str(self) -> str:
        """
        Get the strand as a string ('+', '-', or '.').
        """
        ...

    def length(self) -> int:
        """
        Calculate the length of the contig (end - start).
        """
        ...

    def start_gpos(self) -> GenomicPosition:
        """
        Get the start position as a GenomicPosition.
        """
        ...

    def end_gpos(self) -> GenomicPosition:
        """
        Get the end position as a GenomicPosition.
        """
        ...

    def extend_upstream(self, length: int) -> None:
        """
        Extend the contig upstream by `length` bases in-place.

        The start position will not go below 0.

        Parameters
        ----------
        length : int
            The number of bases to extend upstream.
        """
        ...

    def extend_downstream(self, length: int) -> None:
        """
        Extend the contig downstream by `length` bases in-place.

        The end position will be increased by `length`.

        Parameters
        ----------
        length : int
            The number of bases to extend downstream.
        """
        ...

    def is_in(self, other: "Contig") -> bool:
        """
        Check if this contig is completely contained within another contig.

        Requires both contigs to be on the same chromosome. Strand is not considered.

        Parameters
        ----------
        other : Contig
            The contig to check for containment within.

        Returns
        -------
        bool
            True if this contig is within the other, False otherwise.
        """
        ...

    def is_empty(self) -> bool:
        """
        Checks if the contig is empty (zero start, zero end, and empty seqname).
        """
        ...

    def __str__(self) -> str:
        """
        Get the string representation (seqname:start-end(strand)).
        """
        ...

    def __repr__(self) -> str:
        """
        Get the developer representation.
        """
        ...

    def __richcmp__(self, other: "Contig", op: int) -> bool:
        ...

    def __lt__(self, other: "Contig") -> bool:
        ...

    def __le__(self, other: "Contig") -> bool:
        ...

    def __gt__(self, other: "Contig") -> bool:
        ...

    def __ge__(self, other: "Contig") -> bool:
        ...

class BatchIndex:
    """
    Index for batches in a BSX file.

    Stores an interval tree for each chromosome to quickly find batches
    that overlap a given genomic region. Also maintains the order of
    chromosomes as they were inserted.

    Attributes
    ----------
    map : dict
        A dictionary mapping chromosome names (str) to interval trees.
        Each interval tree maps genomic ranges to batch indices (int).
    chr_order : list
        A list of chromosome names (str) in the order they were first encountered.
    """
    def __init__(self):
        """
        Initialize a new, empty BatchIndex.
        """
        ...

    def insert(self, seqname: str, start: int, end: int, batch_idx: int):
        """
        Insert a contig and its corresponding batch index.

        The contig is defined by its sequence name, start, and end positions.
        The chromosome order is updated if the sequence name is new.

        Parameters
        ----------
        seqname : str
            The name of the sequence/chromosome.
        start : int
            The start position of the contig.
        end : int
            The end position of the contig.
        batch_idx : int
            The index of the batch corresponding to this contig.
        """
        ...

    def find(self, seqname: str, start: int, end: int) -> Optional[List[int]]:
        """
        Find the batch indices that overlap with a given contig.

        Parameters
        ----------
        seqname : str
            The name of the sequence/chromosome.
        start : int
            The start position of the query contig.
        end : int
            The end position of the query contig.

        Returns
        -------
        Optional[List[int]]
            A list of batch indices that overlap the query contig.
            Returns None if the chromosome is not found in the index.
            Returns an empty list if the chromosome is found but no batches overlap.
        """
        ...

    def get_chr_order(self) -> List[str]:
        """
        Return the chromosome order.

        Returns
        -------
        List[str]
            A list of chromosome names in the order they were inserted.
        """
        ...

    def get_chr_index(self, chr_name: str) -> Optional[int]:
        """
        Get the index of a chromosome in the established chromosome order.

        Parameters
        ----------
        chr_name : str
            The name of the chromosome.

        Returns
        -------
        Optional[int]
            The 0-based index of the chromosome in the order, or None if not found.
        """
        ...

    def save(self, filename: str) -> None:
        """
        Serialize the BatchIndex to a file using bincode.

        Parameters
        ----------
        filename : str
            The path to the file where the index will be saved.
        """
        ...

    @staticmethod
    def load(filename: str) -> "BatchIndex":
        """
        Deserialize a BatchIndex from a file using bincode.

        Parameters
        ----------
        filename : str
            The path to the file from which to load the index.

        Returns
        -------
        BatchIndex
            The deserialized BatchIndex object.
        """
        ...

    def sort(self, contigs: List[Contig]) -> List[Contig]:
        """
        Sorts a list of contigs based on the BatchIndex's chromosome order
        and then by start position within each chromosome.
        """
        ...

class GffEntry:
    """
    Represents a GFF (General Feature Format) entry.
    """

    def __init__(self,
                 contig: Contig,
                 source: Optional[str] = None,
                 feature_type: Optional[str] = None,
                 score: Optional[float] = None,
                 phase: Optional[int] = None, # u8 in Rust
                 id: Optional[str] = None
                 ) -> None:
        """
        Create a new GffEntry.

        Parameters
        ----------
        contig : Contig
            The genomic coordinates of the feature.
        source : str, optional
            The source of the feature (e.g., "HAVANA", "Ensembl").
        feature_type : str, optional
            The type of the feature (e.g., "gene", "exon", "CDS").
        score : float, optional
            The score of the feature, if applicable.
        phase : int, optional
            For CDS features, the phase (0, 1, or 2).
        id : str, optional
            A unique identifier for this GFF entry. If not provided, one may be generated.
        """
        ...

    @property
    def id(self) -> str:
        """
        Get the ID of the GFF entry.
        """
        ...

    @property
    def contig(self) -> Contig:
        """
        Get the Contig object representing the coordinates of this GFF entry.
        """
        ...

    @property
    def source(self) -> str:
        """
        Get the source of the GFF entry.
        """
        ...

    @property
    def feature_type(self) -> str:
        """
        Get the feature type of the GFF entry.
        """
        ...

    @property
    def score(self) -> Optional[float]:
        """
        Get the score of the GFF entry, if available.
        """
        ...

    @property
    def phase(self) -> Optional[int]:
        """
        Get the phase of the GFF entry (typically for CDS features), if available.
        """
        ...

    def __repr__(self) -> str:
        """
        Return a developer-friendly string representation of the GffEntry.
        """
        ...


class AnnotStore:
    """
    A store for GFF/BED annotation entries.

    Manages GFF entries, allowing for insertion, removal, and querying
    of parent-child relationships. It can also generate flanking regions
    (upstream/downstream) for specified entries.
    """
    def __init__(self):
        """
        Initialize a new, empty AnnotStore.
        """
        ...

    @staticmethod
    def from_gff(path: str) -> "AnnotStore":
        """
        Create an AnnotStore from a GFF file.
        (Note: This is a high-level constructor; specific parsing logic not shown here.)

        Parameters
        ----------
        path : str
            Path to the GFF file.

        Returns
        -------
        AnnotStore
            An AnnotStore populated with entries from the GFF file.
        """
        ...

    @staticmethod
    def from_bed(path: str) -> "AnnotStore":
        """
        Create an AnnotStore from a BED file.
        (Note: This is a high-level constructor; specific parsing logic not shown here.)

        Parameters
        ----------
        path : str
            Path to the BED file.

        Returns
        -------
        AnnotStore
            An AnnotStore populated with entries from the BED file.
        """
        ...

    def len(self) -> int:
        """
        Return the number of entries in the store.

        Returns
        -------
        int
            The total number of annotation entries.
        """
        ...

    def is_empty(self) -> bool:
        """
        Check if the annotation store is empty.

        Returns
        -------
        bool
            True if the store contains no entries, False otherwise.
        """
        ...

    def insert(self, entry: GffEntry) -> None:
        """
        Inserts a GFF entry into the store.

        If an entry with the same ID already exists, the behavior might depend
        on the underlying Rust implementation (e.g., it might be ignored or overwrite).
        The Python wrapper typically returns None on success or raises an error.

        Parameters
        ----------
        entry : GffEntry
            The GFF entry to insert.
        """
        ...

    def remove(self, id: str) -> Optional[GffEntry]:
        """
        Removes a GFF entry from the store by ID.

        Parameters
        ----------
        id : str
            The ID of the GFF entry to remove.

        Returns
        -------
        Optional[GffEntry]
            The removed entry if found, or None if no entry with the given ID exists.
        """
        ...

    def get_children(self, id: str) -> Optional[List[str]]:
        """
        Gets the children of an entry by its ID.

        Parameters
        ----------
        id : str
            The ID of the parent entry.

        Returns
        -------
        Optional[List[str]]
            A list of IDs of child entries, or None if the parent ID is not found
            or has no children.
        """
        ...

    def get_parents(self, id: str) -> Optional[List[str]]:
        """
        Gets the parents of an entry by its ID.

        Parameters
        ----------
        id : str
            The ID of the child entry.

        Returns
        -------
        Optional[List[str]]
            A list of IDs of parent entries, or None if the child ID is not found
            or has no parents.
        """
        ...

    def id_map(self) -> Dict[str, GffEntry]:
        """
        Returns a copy of the internal ID map, mapping entry IDs to GffEntry objects.
        """
        ...

    def parent_map(self) -> Dict[str, List[str]]:
        """
        Returns a copy of the internal parent map.
        Maps child entry IDs to a list of their parent entry IDs.
        """
        ...

    def children_map(self) -> Dict[str, List[str]]:
        """
        Returns a copy of the internal children map.
        Maps parent entry IDs to a list of their child entry IDs.
        """
        ...

    def add_upstream(self, length: int, selector: Optional[Callable[[GffEntry], bool]] = None) -> None:
        """
        Add upstream regions to selected entries.

        For each entry matching the selector (or all entries if no selector),
        an upstream region of the specified length is created and added to the store.
        The new upstream entry will have its parent attribute set to the ID of the
        original entry.

        Parameters
        ----------
        length : int
            The length of the upstream region to add.
        selector : callable, optional
            A function that takes a GffEntry and returns True if an upstream
            region should be added for it. If None, applies to all entries.
            (Note: current Rust wrapper applies to all).
        """
        ...

    def add_downstream(self, length: int, selector: Optional[Callable[[GffEntry], bool]] = None) -> None:
        """
        Add downstream regions to selected entries.

        For each entry matching the selector (or all entries if no selector),
        a downstream region of the specified length is created and added to the store.
        The new downstream entry will have its parent attribute set to the ID of the
        original entry.

        Parameters
        ----------
        length : int
            The length of the downstream region to add.
        selector : callable, optional
            A function that takes a GffEntry and returns True if a downstream
            region should be added for it. If None, applies to all entries.
            (Note: current Rust wrapper applies to all).
        """
        ...

    def add_flanks(self, length: int, selector: Optional[Callable[[GffEntry], bool]] = None) -> None:
        """
        Add both upstream and downstream regions (flanks) to selected entries.

        This is a convenience method that calls `add_upstream` and `add_downstream`.

        Parameters
        ----------
        length : int
            The length of the upstream and downstream regions to add.
        selector : callable, optional
            A function that takes a GffEntry and returns True if flanks
            should be added for it. If None, applies to all entries.
            (Note: current Rust wrapper applies to all).
        """
        ...

    def iter_sorted(self, index: BatchIndex) -> List[Contig]: # Python stub returns List[Contig], Rust returns AnnotIterator
        """
        Iterate over annotation entries, sorted according to a BatchIndex.

        The entries are grouped by chromosome index (from the BatchIndex) and
        then sorted by their start position within each chromosome group.

        Parameters
        ----------
        index : BatchIndex
            A BatchIndex object providing the chromosome order and indices.

        Returns
        -------
        List[Contig]
            A list of Contig objects representing the annotation entries in sorted order.
            (Note: The Python stub specifies List[Contig], derived from PyContig from Rust which might represent a GffEntry's coordinates.)
        """
        ...

    def get_entry_ids(self) -> List[str]:
        """
        Get a list of all entry IDs present in the store.

        Returns
        -------
        List[str]
            A list of unique GFF entry IDs.
        """
        ...

    def __len__(self) -> int:
        """
        Return the number of entries in the store.

        Returns
        -------
        int
            The total number of annotation entries.
        """
        ...

    def __repr__(self) -> str:
        """
        Return a string representation of the AnnotStore.

        Returns
        -------
        str
            A string indicating the type and number of entries.
        """
        ...


class MethylationStats:
    """
    Represents comprehensive methylation statistics for DNA sequencing data.

    Stores information about methylation levels, variance, coverage
    distribution, and context-specific (CG, CHG, CHH) and strand-specific
    methylation data.
    """
    def __init__(self):
        """
        Initialize a new, empty MethylationStats instance.

        All statistical values are initialized to zero or empty collections.
        """
        ...

    @staticmethod
    def from_data(
        mean_methylation: float,
        variance_methylation: float,
        coverage_distribution: Dict[int, int], # Rust is u16, u32
        context_methylation: Dict[str, Tuple[float, int]], # Python uses str for Context/Strand keys
        strand_methylation: Dict[str, Tuple[float, int]]
    ) -> "MethylationStats":
        """
        Create a MethylationStats instance from pre-calculated statistics.

        Parameters
        ----------
        mean_methylation : float
            The overall mean methylation level.
        variance_methylation : float
            The variance in methylation levels.
        coverage_distribution : Dict[int, int]
            Mapping of coverage depths (int) to their frequencies (int).
        context_methylation : Dict[str, Tuple[float, int]]
            Mapping of sequence contexts (str, e.g., "CG", "CHG", "CHH")
            to a tuple of (sum of methylation levels, count of positions).
        strand_methylation : Dict[str, Tuple[float, int]]
            Mapping of DNA strands (str, e.g., "+", "-")
            to a tuple of (sum of methylation levels, count of positions).

        Returns
        -------
        MethylationStats
            A new MethylationStats instance populated with the provided data.
            NaN values in mean_methylation and variance_methylation are replaced with 0.0.
        """
        ...

    def merge(self, other: "MethylationStats") -> None:
        """
        Merges another MethylationStats instance into this one.

        This performs a weighted merge of statistics, with weights based on
        coverage. The method updates mean methylation, variance, coverage
        distribution, and context/strand-specific statistics.
        """
        ...

    def finalize_methylation(self) -> None:
        """
        Finalizes methylation statistics by converting accumulated sums to mean
        values.

        This should be called after all merging operations are complete and
        before retrieving final statistics.
        """
        ...

    def total_coverage(self) -> int: # Rust is u32
        """
        Computes the total sequencing coverage across all positions.
        """
        ...

    def mean_methylation(self) -> float:
        """
        Computes the genome-wide mean methylation level.
        """
        ...

    @staticmethod
    def merge_multiple(stats_list: List["MethylationStats"]) -> "MethylationStats":
        """
        Merges multiple MethylationStats instances to create a genome-wide
        summary.
        """
        ...

    def coverage_distribution(self) -> Dict[int, int]: # Rust is u16, u32
        """
        Returns the coverage distribution map.
        Key: coverage depth (int), Value: number of positions with that coverage (int).
        """
        ...

    def methylation_var(self) -> float:
        """
        Returns the variance of methylation levels.
        """
        ...

    def context_methylation(self) -> Dict[str, Tuple[float, int]]:
        """
        Returns methylation statistics per sequence context.
        Key: context (str), Value: (mean methylation level, count of positions).
        Note: Rust version might store sum before finalization. Python stub implies mean.
        """
        ...

    def strand_methylation(self) -> Dict[str, Tuple[float, int]]:
        """
        Returns methylation statistics per DNA strand.
        Key: strand (str), Value: (mean methylation level, count of positions).
        Note: Rust version might store sum before finalization. Python stub implies mean.
        """
        ...


class LazyBsxBatch:
    """
    A lazy representation of a BSX batch for efficient query operations on decoded data.

    Allows for chaining of filter and transformation operations that are only
    executed when the `collect` method is called or when statistics are computed.
    """
    def __init__(self, batch: "BsxBatch"):
        """
        Initialize a LazyBsxBatch from a BsxBatch.

        Parameters
        ----------
        batch : BsxBatch
            The BsxBatch to operate on lazily.
        """
        ...

    def into_report(self, report_type: ReportTypeSchema) -> pl.DataFrame:
        """
        Convert the batch to a specified report type format.

        Transforms the internal DataFrame to match the column names and structure
        of the target report type (e.g., Bismark, CgMap).

        Parameters
        ----------
        report_type : ReportTypeSchema
            The target report format.

        Returns
        -------
        pl.DataFrame
            A Polars DataFrame formatted according to the specified report type.
        """
        ...

    def collect(self) -> "BsxBatch":
        """
        Execute all queued lazy operations and return a concrete BsxBatch.

        Returns
        -------
        BsxBatch
            A new BsxBatch containing the results of the lazy computations.
        """
        ...

    def filter_pos_lt(self, pos: int) -> "LazyBsxBatch": # Rust is u32
        """
        Filter positions less than the specified value.

        Parameters
        ----------
        pos : int
            The position value. Rows with 'position' < `pos` will be kept.

        Returns
        -------
        LazyBsxBatch
            A new LazyBsxBatch with the filter applied.
        """
        ...

    def filter_pos_gt(self, pos: int) -> "LazyBsxBatch": # Rust is u32
        """
        Filter positions greater than the specified value.

        Parameters
        ----------
        pos : int
            The position value. Rows with 'position' > `pos` will be kept.

        Returns
        -------
        LazyBsxBatch
            A new LazyBsxBatch with the filter applied.
        """
        ...

    def filter_coverage_lt(self, coverage: int) -> "LazyBsxBatch": # Rust is u32
        """
        Filter entries with total coverage less than the specified value.

        Parameters
        ----------
        coverage : int
            The coverage threshold. Rows with 'count_total' < `coverage` will be kept.

    Returns
        -------
        LazyBsxBatch
            A new LazyBsxBatch with the filter applied.
        """
        ...

    def filter_strand(self, strand: Strand) -> "LazyBsxBatch":
        """
        Filter entries by strand value.

        Parameters
        ----------
        strand : Strand
            The strand value to filter by (e.g., Strand.Forward, Strand.Reverse).

        Returns
        -------
        LazyBsxBatch
            A new LazyBsxBatch with the filter applied.
        """
        ...

    def filter_context(self, context: Context) -> "LazyBsxBatch":
        """
        Filter entries by context value.

        Parameters
        ----------
        context : Context
            The context value to filter by (e.g., Context.CG, Context.CHG).

        Returns
        -------
        LazyBsxBatch
            A new LazyBsxBatch with the filter applied.
        """
        ...

    def mark_low_coverage(self, threshold: int) -> "LazyBsxBatch": # Rust is u32
        """
        Mark entries with coverage below a threshold.

        Rows where 'count_total' is less than `threshold` will have their
        'count_total' and 'count_m' set to 0, and 'density' set to NaN.

        Parameters
        ----------
        threshold : int
            The coverage threshold.

        Returns
        -------
        LazyBsxBatch
            A new LazyBsxBatch with low coverage entries marked.
        """
        ...

    def align_with_contexts(self, context_data: ContextData, chr_val: str) -> "LazyBsxBatch":
        """
        Align the batch with provided context data for a specified chromosome.

        Performs a left join of the context data DataFrame with the batch's data,
        using 'position' as the join key. Drops existing 'context' and 'strand'
        columns from the batch before joining. Fills nulls introduced by the join:
        'chr' with `chr_val`, 'count_m' and 'count_total' with 0, 'density' with NaN.

        Parameters
        ----------
        context_data : ContextData
            The context data to align with.
        chr_val : str
            The chromosome value to fill in for new rows from `context_data`.

        Returns
        -------
        LazyBsxBatch
            A new LazyBsxBatch with aligned context information.
        """
        ...

    def stats(self) -> MethylationStats:
        """
        Calculate methylation statistics for the batch.

        Computes overall mean methylation, variance, coverage distribution,
        and context-specific and strand-specific methylation means.

        Returns
        -------
        MethylationStats
            An object containing the calculated statistics.
        """
        ...

class LazyEncodedBsxBatch:
    """
    A lazy representation of an EncodedBsxBatch for efficient query operations.

    Allows for chaining of filter and transformation operations that are only
    executed when the `collect` method is called or when statistics are computed.
    Operates on data where 'strand' and 'context' may be boolean encoded.
    """
    def __init__(self, batch: "EncodedBsxBatch"):
        """
        Initialize a LazyEncodedBsxBatch from an EncodedBsxBatch.

        Parameters
        ----------
        batch : EncodedBsxBatch
            The EncodedBsxBatch to operate on lazily.
        """
        ...

    def into_report(self, report_type: ReportTypeSchema) -> pl.DataFrame:
        """
        Convert the batch to a specified report type format.

        Transforms the internal DataFrame to match the column names and structure
        of the target report type. Handles decoding of boolean 'strand' and 'context'
        columns if necessary for the report type.

        Parameters
        ----------
        report_type : ReportTypeSchema
            The target report format.

        Returns
        -------
        pl.DataFrame
            A Polars DataFrame formatted according to the specified report type.
        """
        ...

    def collect(self) -> "EncodedBsxBatch":
        """
        Execute all queued lazy operations and return a concrete EncodedBsxBatch.

        Returns
        -------
        EncodedBsxBatch
            A new EncodedBsxBatch containing the results of the lazy computations.
        """
        ...

    def filter_pos_lt(self, pos: int) -> "LazyEncodedBsxBatch": # Rust is u32
        """
        Filter positions less than the specified value.

        Parameters
        ----------
        pos : int
            The position value. Rows with 'position' < `pos` will be kept.

        Returns
        -------
        LazyEncodedBsxBatch
            A new LazyEncodedBsxBatch with the filter applied.
        """
        ...

    def filter_pos_gt(self, pos: int) -> "LazyEncodedBsxBatch": # Rust is u32
        """
        Filter positions greater than the specified value.

        Parameters
        ----------
        pos : int
            The position value. Rows with 'position' > `pos` will be kept.

        Returns
        -------
        LazyEncodedBsxBatch
            A new LazyEncodedBsxBatch with the filter applied.
        """
        ...

    def filter_coverage_lt(self, coverage: int) -> "LazyEncodedBsxBatch": # Rust is u32
        """
        Filter entries with total coverage less than the specified value.

        Parameters
        ----------
        coverage : int
            The coverage threshold. Rows with 'count_total' < `coverage` will be kept.

        Returns
        -------
        LazyEncodedBsxBatch
            A new LazyEncodedBsxBatch with the filter applied.
        """
        ...

    def filter_strand(self, strand: Strand) -> "LazyEncodedBsxBatch":
        """
        Filter entries by strand value.

        Compares with the 'strand' column, which may be boolean encoded
        (True for '+', False for '-', None for '.').

        Parameters
        ----------
        strand : Strand
            The strand value to filter by.

        Returns
        -------
        LazyEncodedBsxBatch
            A new LazyEncodedBsxBatch with the filter applied.
        """
        ...

    def filter_context(self, context: Context) -> "LazyEncodedBsxBatch":
        """
        Filter entries by context value.

        Compares with the 'context' column, which may be boolean encoded
        (True for 'CG', False for 'CHG', None for 'CHH').

        Parameters
        ----------
        context : Context
            The context value to filter by.

        Returns
        -------
        LazyEncodedBsxBatch
            A new LazyEncodedBsxBatch with the filter applied.
        """
        ...

    def mark_low_coverage(self, threshold: int) -> "LazyEncodedBsxBatch": # Rust is u32
        """
        Mark entries with coverage below a threshold.

        Rows where 'count_total' is less than `threshold` will have their
        'count_total' and 'count_m' set to 0, and 'density' set to NaN.

        Parameters
        ----------
        threshold : int
            The coverage threshold.

        Returns
        -------
        LazyEncodedBsxBatch
            A new LazyEncodedBsxBatch with low coverage entries marked.
        """
        ...

    def align_with_contexts(self, context_data: ContextData, chr_val: str) -> "LazyEncodedBsxBatch":
        """
        Align the batch with provided encoded context data for a specified chromosome.

        Performs a left join of the context data DataFrame (in encoded form)
        with the batch's data. See `LazyBsxBatch.align_with_contexts` for more details.

        Parameters
        ----------
        context_data : ContextData
            The context data to align with.
        chr_val : str
            The chromosome value to fill in for new rows.

        Returns
        -------
        LazyEncodedBsxBatch
            A new LazyEncodedBsxBatch with aligned context information.
        """
        ...

    def stats(self) -> MethylationStats:
        """
        Calculate methylation statistics for the encoded batch.

        Computes statistics similar to `LazyBsxBatch.stats`, accounting for
        potentially encoded 'strand' and 'context' values.

        Returns
        -------
        MethylationStats
            An object containing the calculated statistics.
        """
        ...


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
        context_data: Optional[ContextData] = None,
        report_schema: Optional[ReportTypeSchema] = None,
        chr_values: Optional[List[str]] = None
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
        chr_values : List[str], optional
            A list of unique chromosome string values to establish a categorical type.
            Required if converting from formats that don't inherently define categorical order.

        Raises
        ------
        ValueError
            If any of the enabled checks fail.
        """
        ...

    def lazy(self) -> LazyEncodedBsxBatch:
        """
        Create a lazy EncodedBsxBatch from the current batch.

        Returns
        -------
        LazyEncodedBsxBatch
            A lazy EncodedBsxBatch.
        """

    @staticmethod
    def from_dataframe_unchecked(data: pl.DataFrame) -> "EncodedBsxBatch":
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

    def chr(self) -> pl.Series:
        """
        Access the chromosome column as a Polars Series.

        Returns
        -------
        pl.Series
            The chromosome data (Categorical).
        """
        ...

    def position(self) -> pl.Series:
        """
        Access the position column as a Polars Series.

        Returns
        -------
        pl.Series
            The position data (UInt32).
        """
        ...

    def strand(self) -> pl.Series:
        """
        Access the strand column as a Polars Series.

        Returns
        -------
        pl.Series
            The strand data (Boolean, where True is '+', False is '-', Null is '.').
        """
        ...

    def context(self) -> pl.Series:
        """
        Access the context column as a Polars Series.

        Returns
        -------
        pl.Series
            The context data (Boolean, where True is 'CG', False is 'CHG', Null is 'CHH').
        """
        ...

    def count_m(self) -> pl.Series:
        """
        Access the methylated counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The methylated count data (Int16).
        """
        ...

    def count_total(self) -> pl.Series:
        """
        Access the total counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The total count data (Int16).
        """
        ...

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

    def start_pos(self) -> int: # Rust is u32
        """
        Get the starting genomic position in the batch.

        Returns
        -------
        int
            The minimum position value. Returns 0 if the batch is empty.
        """
        ...

    def end_pos(self) -> int: # Rust is u32
        """
        Get the ending genomic position in the batch.

        Returns
        -------
        int
            The maximum position value. Returns 0 if the batch is empty.
        """
        ...

    def start_gpos(self) -> GenomicPosition:
        """
        Returns the genomic position at the start of the batch.
        """
        ...

    def end_gpos(self) -> GenomicPosition:
        """
        Returns the genomic position at the end of the batch.
        """
        ...

    def as_contig(self) -> Optional[Contig]:
        """
        Represents the batch as a Contig object.

        The contig spans from the first position to the last position in the batch,
        inclusive. The end of the contig will be `self.end_pos() + 1` to
        adhere to the common exclusive-end convention for contigs.
        The strand will be set to '.' (None).

        Returns
        -------
        Contig or None
            A Contig object representing the genomic span of the batch,
            or None if the batch is empty.
        """
        ...

    def vstack(self, other: "EncodedBsxBatch") -> "EncodedBsxBatch": # Rust takes &[<Py $wrap_type>]
        """
        Vertically stack this batch with other batches.

        Creates a new batch containing rows from this batch and all other batches.
        Performs checks on the combined data.

        Parameters
        ----------
        other : List[EncodedBsxBatch]
            A list of batches to stack underneath this one.

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

    def extend(self, other: "EncodedBsxBatch") -> None: # Rust takes &[<Py $wrap_type>]
        """
        Extend this batch by appending rows from other batches in-place.

        Modifies the current batch. Performs checks after extending.

        Parameters
        ----------
        other : List[EncodedBsxBatch]
            A list of batches whose rows will be appended.

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
        context_data: Optional[ContextData] = None,
        report_schema: Optional[ReportTypeSchema] = None,
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
        chr_values : List[str], optional
            A list of unique chromosome string values. Used if `report_schema` implies
            conversion to categorical chromosome types, like in `EncodedBsxBatch`.

        Raises
        ------
        ValueError
            If any of the enabled checks fail.
        """
        ...

    def lazy(self) -> LazyBsxBatch:
        """
        Create a lazy EncodedBsxBatch from the current batch.

        Returns
        -------
        LazyBsxBatch
            A lazy BsxBatch.
        """

    @staticmethod
    def from_dataframe_unchecked(data: pl.DataFrame) -> "BsxBatch":
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

    def chr(self) -> pl.Series:
        """
        Access the chromosome column as a Polars Series.

        Returns
        -------
        pl.Series
            The chromosome data (Utf8).
        """
        ...

    def position(self) -> pl.Series:
        """
        Access the position column as a Polars Series.

        Returns
        -------
        pl.Series
            The position data (UInt64).
        """
        ...

    def strand(self) -> pl.Series:
        """
        Access the strand column as a Polars Series.

        Returns
        -------
        pl.Series
            The strand data (Utf8).
        """
        ...

    def context(self) -> pl.Series:
        """
        Access the context column as a Polars Series.

        Returns
        -------
        pl.Series
            The context data (Utf8).
        """
        ...

    def count_m(self) -> pl.Series:
        """
        Access the methylated counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The methylated count data (UInt32).
        """
        ...

    def count_total(self) -> pl.Series:
        """
        Access the total counts column as a Polars Series.

        Returns
        -------
        pl.Series
            The total count data (UInt32).
        """
        ...

    def density(self) -> pl.Series:
        """
        Access the density column as a Polars Series.

        Returns
        -------
        pl.Series
            The density data (Float64).
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

    def start_pos(self) -> int: # Rust is u32
        """
        Get the starting genomic position in the batch.

        Returns
        -------
        int
            The minimum position value. Returns 0 if the batch is empty.
        """
        ...

    def end_pos(self) -> int: # Rust is u32
        """
        Get the ending genomic position in the batch.

        Returns
        -------
        int
            The maximum position value. Returns 0 if the batch is empty.
        """
        ...

    def start_gpos(self) -> GenomicPosition:
        """
        Returns the genomic position at the start of the batch.
        """
        ...

    def end_gpos(self) -> GenomicPosition:
        """
        Returns the genomic position at the end of the batch.
        """
        ...

    def as_contig(self) -> Optional[Contig]: # Rust is PyResult<Option<PyContig>>
        """
        Represents the batch as a Contig object.

        The contig spans from the first position to the last position in the batch,
        inclusive. The end of the contig will be `self.end_pos() + 1` to
        adhere to the common exclusive-end convention for contigs.
        The strand will be set to '.' (None).

        Returns
        -------
        Contig or None
            A Contig object representing the genomic span of the batch,
            or None if the batch is empty.
        """
        ...

    def vstack(self, other: "BsxBatch") -> "BsxBatch":
        """
        Vertically stack this batch with other batches.

        Creates a new batch containing rows from this batch and all other batches.
        Performs checks on the combined data.

        Parameters
        ----------
        other : List[BsxBatch]
            A list of batches to stack underneath this one.

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
        Extend this batch by appending rows from other batches in-place.

        Modifies the current batch. Performs checks after extending.

        Parameters
        ----------
        other : List[BsxBatch]
            A list of batches whose rows will be appended.

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

def encode(decoded: BsxBatch, chr_values: List[str]) -> EncodedBsxBatch:
    """
    Converts a decoded batch to an encoded batch format.
    """
    ...

def decode(encoded: EncodedBsxBatch) -> BsxBatch:
    """
    Converts an encoded batch to a decoded batch format.
    """
    ...
