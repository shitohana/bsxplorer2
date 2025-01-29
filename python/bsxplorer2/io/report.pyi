from enum import Enum
from polars import Schema, DataFrame
from ..misc import BsxBatch, EncodedBsxBatch

class ReportTypeSchema(Enum):
    Bismark = "Bismark"
    CgMap = "CgMap"
    BedGraph = "BedGraph"
    Coverage = "Coverage"

    def get_schema(self) -> Schema: ...

class ReportReader:
    def __init__(
            self,
            report_path: str,
            report_type: ReportTypeSchema,
            /,
            rechunk: bool = false,
            n_threads: int | None = None,
            low_memory: bool = false,
            n_rows: int | None = None,
            chunk_size: int = 10_000,
            skip: int = 0,
            fasta_path: str | None = None,
            fai_path: str | None = None,
            batch_per_read: int = 16,
            batch_size: int = 2 << 20
    ) -> 'ReportReader': ...

    def __iter__(self) -> 'ReportReader': ...
    def __next__(self) -> 'BsxBatch': ...


class ReportWriter:
    def __init__(self, sink: str, schema: ReportTypeSchema, n_threads: int) -> ReportWriter: ...
    def write_batch(self, batch: BsxBatch): ...
    def write_df(self, df: DataFrame): ...