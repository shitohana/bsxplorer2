from enum import Enum
from typing import Callable, Dict, Optional, List, Tuple, Sequence
import polars as pl

__all__ = [
    "Strand",
    "Context",
    "BsxBatch",
    "AggMethod",
    "BsxColumns",
    "ContextData",
    "ReportTypeSchema",
    "LazyBsxBatch",
    "Contig",
    "GenomicPosition",
    "MethylationStats",
    "AnnotStore",
    "GffEntry",
    "GffEntryAttributes",
    "BatchIndex",
]

class Strand(Enum):
    """Represents genomic strand."""

    Forward = "Forward"
    Reverse = "Reverse"
    Null = "Null"

    @property
    def name(self) -> str:
        """Returns the name of the strand."""
        ...

class Context(Enum):
    """Represents methylation context."""

    CG = "CG"
    CHG = "CHG"
    CHH = "CHH"

    @property
    def name(self) -> str:
        """Returns the name of the context."""
        ...

class BsxColumns(Enum):
    """Represents the columns expected in BSX data."""

    Chr = "Chr"
    Position = "Position"
    Strand = "Strand"
    Context = "Context"
    CountM = "CountM"
    CountTotal = "CountTotal"
    Density = "Density"

    @staticmethod
    def schema() -> pl.Schema:
        """Returns the Polars Schema for the BSX columns."""
        ...

    def as_str(self) -> str:
        """Returns the string representation of the column name."""
        ...

    @staticmethod
    def colnames() -> List[str]:
        """Returns an array containing all BSX column names as strings."""
        ...

    @staticmethod
    def has_name(name: str) -> bool:
        """Checks if the given string matches any of the BSX column names."""
        ...

class AggMethod(Enum):
    """Enumerates different aggregation methods for methylation density values."""

    Mean = "Mean"
    GeometricMean = "GeometricMean"
    Median = "Median"
    Max = "Max"
    Min = "Min"

class BsxBatch:
    """Represents a batch of methylation data as a Polars DataFrame.

    This struct wraps a Polars DataFrame with specific columns expected for
    methylation data (chromosome, position, strand, context, counts, density).
    It provides methods for accessing columns, manipulating data, and performing
    specific methylation-related operations.
    """

    def __init__(
        self,
        chr: str,
        chr_dtype: Optional[pl.DataType],
        positions: Sequence[int],
        strand: Sequence[bool],
        context: Sequence[Optional[bool]],
        count_m: Sequence[int],
        count_total: Sequence[int],
    ) -> None:
        """Creates a new BsxBatch from columnar data. All vectors must have the same length."""
        ...

    @staticmethod
    def from_dataframe(
        data: pl.DataFrame,
        check_nulls: bool = True,
        check_sorted: bool = True,
        check_duplicates: bool = True,
        rechunk: bool = True,
        check_single_chr: bool = True,
        chr_values: Optional[Sequence[str]] = None,
    ) -> "BsxBatch":
        """Creates a BsxBatch from a Polars DataFrame after validation and casting."""
        ...

    @staticmethod
    def schema() -> pl.Schema:
        """Returns the Polars Schema for the BSX columns."""
        ...

    @staticmethod
    def empty(chr_dtype: Optional[pl.DataType] = None) -> "BsxBatch":
        """Creates an empty BsxBatch with the correct schema."""
        ...

    @staticmethod
    def concat(batches: List["BsxBatch"]) -> "BsxBatch":
        """Concatenates multiple BsxBatch instances into a single batch."""
        ...

    def chr(self) -> pl.Series:
        """Returns the chromosome column."""
        ...

    def position(self) -> pl.Series:
        """Returns the position column."""
        ...

    def strand(self) -> pl.Series:
        """Returns the strand column."""
        ...

    def context(self) -> pl.Series:
        """Returns the context column."""
        ...

    def count_m(self) -> pl.Series:
        """Returns the methylated count column."""
        ...

    def count_total(self) -> pl.Series:
        """Returns the total count column."""
        ...

    def density(self) -> pl.Series:
        """Returns the density column."""
        ...

    def column(self, name: BsxColumns) -> pl.Series:
        """Gets a specific column by its enum identifier."""
        ...

    def is_empty(self) -> bool:
        """Checks if the BsxBatch is empty (contains no rows)."""
        ...

    def split_at(self, index: int) -> Tuple["BsxBatch", "BsxBatch"]:
        """Splits the BsxBatch into two BsxBatches at the given index."""
        ...

    def rechunk(self) -> "BsxBatch":
        """Rechunks the underlying DataFrame in place."""
        ...

    def data(self) -> pl.DataFrame:
        """Returns a reference to the underlying Polars DataFrame."""
        ...

    def into_dataframe(self) -> pl.DataFrame:
        """Consumes the batch and returns the inner DataFrame."""
        ...

    def slice(self, start: int, length: int) -> "BsxBatch":
        """Creates a slice of the BsxBatch."""
        ...

    def add_context_data(self, context_data: "ContextData") -> "BsxBatch":
        """Adds context data (strand and context) to the BsxBatch."""
        ...

    def extend(self, other: "BsxBatch") -> None:
        """Appends another BsxBatch to this one, performing validation."""
        ...

    def extend_unchecked(self, other: "BsxBatch") -> None:
        """Appends another BsxBatch
        to this one without performing validation."""
        ...

    def shrink(
        self, min_size: int, beta: Optional[float] = None
    ) -> Tuple[List[float], List[float]]:
        """Segments the methylation data using a specified algorithm and aggregates densities."""
        ...

    def discretise(
        self, n_fragments: int, method: AggMethod
    ) -> Tuple[List[float], List[float]]:
        """
        Discretises the batch into a fixed number of fragments based on genomic position and aggregates densities."""
        ...

    def partition(
        self, breakpoints: List[int], method: AggMethod
    ) -> Tuple[List[float], List[float]]:
        """Partitions the batch based on the provided breakpoints and aggregates densities within each partition."""
        ...

    def seqname(self) -> Optional[str]:
        """Gets the sequence name (chromosome) for the batch."""
        ...

    def first_pos(self) -> Optional[int]:
        """Gets the position of the first site in the batch."""
        ...

    def last_pos(self) -> Optional[int]:
        """Gets the position of the last site in the batch."""
        ...

    def first_genomic_pos(self) -> Optional["GenomicPosition"]:
        """Gets the genomic position of the first site in the batch."""
        ...

    def last_genomic_pos(self) -> Optional["GenomicPosition"]:
        """Gets the genomic position of the last site in the batch."""
        ...

    def as_contig(self) -> Optional["Contig"]:
        """Represents the batch as a Contig if possible."""
        ...

    def get_methylation_stats(self) -> "MethylationStats":
        """Calculates and returns various methylation statistics for the batch."""
        ...

    def get_coverage_dist(self) -> Dict[int, int]:
        """Gets the distribution of coverage counts across the batch."""
        ...

    def get_context_stats(self) -> Dict[Context, Tuple[float, int]]:
        """Gets methylation statistics grouped by context (CG/CHG/CHH)."""
        ...

    def get_strand_stats(self) -> Dict[Strand, Tuple[float, int]]:
        """Gets methylation statistics grouped by strand (+/-)."""
        ...

    def as_binom(self, mean: float, pvalue: float) -> "BsxBatch":
        """Filters sites based on a binomial test p-value and replaces counts and density."""
        ...

    def into_report(self, report_type: "ReportTypeSchema") -> pl.DataFrame:
        """Converts the BsxBatch into a DataFrame formatted according to the specified ReportType."""
        ...

    def lazy(self) -> "LazyBsxBatch":
        """Converts the BsxBatch into a LazyBsxBatch."""
        ...

    def height(self) -> int:
        """Returns the number of rows (sites) in the batch."""
        ...

class ContextData:
    """Stores genomic context (like CG, CHG, CHH sites) derived from sequence information."""

    def __init__(self, sequence: List[int]) -> None:
        """Creates a ContextData from a sequence of bytes."""
        ...

    @staticmethod
    def empty() -> "ContextData":
        """Creates an empty ContextData."""
        ...

    def is_empty(self) -> bool:
        """Checks if the ContextData is empty."""
        ...

    def take(self) -> Tuple[List[int], List[Strand], List[Context]]:
        """Takes the data, returning vectors of positions, strands, and contexts."""
        ...

    def to_decoded_df(self) -> pl.DataFrame:
        """Converts the ContextData to a Polars DataFrame."""
        ...

class ReportTypeSchema(Enum):
    """Enumerates different report type schemas."""

    Bismark = "Bismark"
    CgMap = "CgMap"
    BedGraph = "BedGraph"
    Coverage = "Coverage"

    def col_names(self) -> List[str]:
        """Returns the column names for this report type."""
        ...

    def schema(self) -> pl.Schema:
        """Returns the schema for this report type."""
        ...

    def chr_col(self) -> str:
        """Returns the chromosome column name."""
        ...

    def position_col(self) -> str:
        """Returns the position column name."""
        ...

    def context_col(self) -> Optional[str]:
        """Returns the context column name if present."""
        ...

    def strand_col(self) -> Optional[str]:
        """Returns the strand column name if present."""
        ...

    def need_align(self) -> bool:
        """Returns whether alignment is needed for this report type."""
        ...

class LazyBsxBatch:
    """A wrapper around a Polars LazyFrame for BsxBatch, allowing for lazy evaluation of operations."""

    def __init__(self, batch: BsxBatch) -> None:
        """Creates a LazyBsxBatch from a BsxBatch."""
        ...

    def collect(self) -> BsxBatch:
        """Collects the lazy batch into a BsxBatch."""
        ...

    def filter_pos_lt(self, pos: int) -> "LazyBsxBatch":
        """Filters positions less than the specified value."""
        ...

    def filter_pos_gt(self, pos: int) -> "LazyBsxBatch":
        """Filters positions greater than the specified value."""
        ...

    def filter_coverage_gt(self, coverage: int) -> "LazyBsxBatch":
        """Filters entries with coverage greater than the specified value."""
        ...

    def filter_strand(self, strand: Strand) -> "LazyBsxBatch":
        """Filters entries by strand value."""
        ...

    def filter_context(self, context: Context) -> "LazyBsxBatch":
        """Filters entries by context value."""
        ...

class GenomicPosition:
    """Represents a genomic position with a sequence name and a position."""

    def __init__(self, seqname: str, position: int) -> None:
        """Creates a new GenomicPosition."""
        ...

    @property
    def seqname(self) -> str:
        """Returns the sequence name."""
        ...

    @property
    def position(self) -> int:
        """Returns the position."""
        ...

    def is_zero(self) -> bool:
        """Checks if the position is zero."""
        ...

    def __add__(self, other: "GenomicPosition") -> Optional["GenomicPosition"]:
        """Adds two GenomicPositions."""
        ...

    def __sub__(self, other: "GenomicPosition") -> Optional["GenomicPosition"]:
        """Subtracts two GenomicPositions."""
        ...

class Contig:
    """Represents a contig with a sequence name, start position, end position, and strand."""

    def __init__(self, seqname: str, start: int, end: int, strand: str = ".") -> None:
        """Creates a new Contig."""
        ...

    @property
    def seqname(self) -> str:
        """Returns the sequence name."""
        ...

    @property
    def start(self) -> int:
        """Returns the start position."""
        ...

    @property
    def end(self) -> int:
        """Returns the end position."""
        ...

    @property
    def strand(self) -> Strand:
        """Returns the strand."""
        ...

    def strand_str(self) -> str:
        """Returns the strand as a string."""
        ...

    def length(self) -> int:
        """Returns the length of the contig."""
        ...

    def start_gpos(self) -> GenomicPosition:
        """Returns the start position as a GenomicPosition."""
        ...

    def end_gpos(self) -> GenomicPosition:
        """Returns the end position as a GenomicPosition."""
        ...

    def extend_upstream(self, length: int) -> None:
        """Extends the contig upstream by a given length."""
        ...

    def extend_downstream(self, length: int) -> None:
        """Extends the contig downstream by a given length."""
        ...

    def is_in(self, other: "Contig") -> bool:
        """Checks if this contig is fully contained within another contig."""
        ...

    def is_empty(self) -> bool:
        """Checks if the contig is empty."""
        ...

class MethylationStats:
    """Comprehensive methylation statistics."""

    def __init__(self) -> None:
        """Creates an empty instance."""
        ...

    @staticmethod
    def from_data(
        mean_methylation: float,
        variance_methylation: float,
        coverage_distribution: Dict[int, int],
        context_methylation: Dict[str, Tuple[float, int]],
        strand_methylation: Dict[str, Tuple[float, int]],
    ) -> "MethylationStats":
        """Creates an instance from data."""
        ...

    def merge(self, other: "MethylationStats") -> None:
        """Merges another instance into this one using weighted averages based on coverage."""
        ...

    def finalize_methylation(self) -> None:
        """Finalizes statistics by converting sums to means."""
        ...

    def total_coverage(self) -> int:
        """Computes the total sequencing coverage."""
        ...

    def mean_methylation(self) -> float:
        """Computes the genome-wide mean methylation."""
        ...

    @staticmethod
    def merge_multiple(stats_list: List["MethylationStats"]) -> "MethylationStats":
        """Merges multiple instances."""
        ...

    def coverage_distribution(self) -> Dict[int, int]:
        """Returns the coverage distribution."""
        ...

    def methylation_var(self) -> float:
        """Returns the methylation variance."""
        ...

    def context_methylation(self) -> Dict[str, Tuple[float, int]]:
        """Returns methylation statistics per context."""
        ...

    def strand_methylation(self) -> Dict[str, Tuple[float, int]]:
        """Returns methylation statistics per strand."""
        ...

class GffEntryAttributes:
    """A structured representation of the key-value pairs found in the GFF attributes column."""

    def __init__(self) -> None:
        """Creates a new GffEntryAttributes."""
        ...

    @property
    def get_id(self) -> Optional[str]:
        """Returns the ID attribute."""
        ...

    @property
    def get_name(self) -> Optional[List[str]]:
        """Returns the Name attribute."""
        ...

    @property
    def get_alias(self) -> Optional[List[str]]:
        """Returns the Alias attribute."""
        ...

    @property
    def get_parent(self) -> Optional[List[str]]:
        """Returns the Parent attribute."""
        ...

class GffEntry:
    """Represents a single feature/annotation entry from a file like GFF or BED."""

    def __init__(
        self,
        contig: Contig,
        source: Optional[str] = None,
        feature_type: Optional[str] = None,
        score: Optional[float] = None,
        phase: Optional[int] = None,
        id: Optional[str] = None,
    ) -> None:
        """Creates a new GffEntry from a contig."""
        ...

    @property
    def id(self) -> str:
        """Returns the ID of the entry."""
        ...

    @property
    def contig(self) -> Contig:
        """Returns the contig of the entry."""
        ...

    @property
    def source(self) -> str:
        """Returns the source of the entry."""
        ...

    @property
    def feature_type(self) -> str:
        """Returns the feature type of the entry."""
        ...

    @property
    def score(self) -> Optional[float]:
        """Returns the score of the entry."""
        ...

    @property
    def phase(self) -> Optional[int]:
        """Returns the phase of the entry."""
        ...

    @property
    def attributes(self) -> GffEntryAttributes:
        """Returns the attributes of the entry."""
        ...

class AnnotStoreIterator:
    """Iterator for AnnotStore entries."""

    def __iter__(self) -> "AnnotStoreIterator":
        """Returns the iterator."""
        ...

    def __next__(self) -> Optional[GffEntry]:
        """Returns the next entry."""
        ...

class AnnotStore:
    """A data structure to store and manage genomic annotations (GffEntry).

    It uses multiple internal structures for efficient lookup:
    - A HashMap for ID-based lookup of GffEntry.
    - An id_tree::Tree to represent parent-child relationships between entries.
    - A HashMap mapping sequence names to IntervalTree for genomic range queries.
    """

    def __init__(self) -> None:
        """Creates a new empty AnnotStore."""
        ...

    @staticmethod
    def from_gff(path: str) -> "AnnotStore":
        """Creates a new AnnotStore by reading and parsing a GFF file."""
        ...

    @staticmethod
    def from_bed(path: str) -> "AnnotStore":
        """Creates a new AnnotStore by reading and parsing a BED file."""
        ...

    @staticmethod
    def with_capacity(capacity: int) -> "AnnotStore":
        """Creates a new empty AnnotStore with a pre-allocated capacity."""
        ...

    def len(self) -> int:
        """Returns the number of entries in the store."""
        ...

    def is_empty(self) -> bool:
        """Returns True if the store contains no entries."""
        ...

    def insert(self, entry: GffEntry) -> None:
        """Inserts a GffEntry into the store."""
        ...

    def get_entry(self, id: str) -> Optional[GffEntry]:
        """Retrieves a reference to a GffEntry by its ID."""
        ...

    def get_entries_regex(self, pattern: str) -> List[GffEntry]:
        """Retrieves a vector of references to GffEntry objects whose IDs match a regex pattern."""
        ...

    def get_children(self, id: str) -> Optional[List[str]]:
        """Retrieves the IDs of the direct children of a given entry in the tree structure."""
        ...

    def get_parent(self, id: str) -> Optional[str]:
        """Retrieves the ID of the direct parent of a given entry in the tree structure."""
        ...

    def genomic_query(self, contig: Contig) -> Optional[List[str]]:
        """Queries the store for entries overlapping a given genomic contig/range."""
        ...

    def get_feature_types(self) -> List[str]:
        """Returns a vector of unique feature types present in the store."""
        ...

    def add_flanks(
        self, parents: list[str], flank: int, prefix: str
    ) -> None:
        """Adds flanking regions to annotation map."""
        ...

    def sort_self(self) -> None:
        """Sorts the children of each node in the internal tree based on the start position."""
        ...

    def iter(self) -> AnnotStoreIterator:
        """Returns an iterator over all entries in the store."""
        ...

class BatchIndex:
    """Python wrapper for BatchIndex."""

    def __init__(self) -> None:
        """Creates a new BatchIndex."""
        ...

    def insert(self, seqname: str, start: int, end: int, batch_idx: int) -> None:
        """Insert a contig and its corresponding batch index."""
        ...

    def find(self, seqname: str, start: int, end: int) -> Optional[List[int]]:
        """Find batch indices that overlap with a given contig."""
        ...

    def get_chr_order(self) -> List[str]:
        """Get chromosome order."""
        ...

    def get_chr_index(self, chr: str) -> Optional[int]:
        """Get index of a chromosome in the order."""
        ...

    def save(self, filename: str) -> None:
        """Save index to a file."""
        ...

    @staticmethod
    def load(filename: str) -> "BatchIndex":
        """Load index from a file."""
        ...

    def sort(self, contigs: List[Contig]) -> List[Contig]:
        """Sort contigs using the index."""
        ...
