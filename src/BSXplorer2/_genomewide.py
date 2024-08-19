import uuid
from dataclasses import dataclass
from pathlib import Path

import pyarrow.dataset as pads

from ._batches import ReportTypes
from ._readers import UniversalReader, UniversalBatch
from ._report import ReportStats


@dataclass
class Chromosome:
        uuid: uuid.UUID
        stats: dict
        source: str | Path

        @property
        def dataset(self):
            return pads.dataset(source=self.source)


class GenomeWide:
    __slots__ = {
        "chromosomes": list[Chromosome]
    }

    @classmethod
    def from_report(
            cls,
            file: str | Path,
            report_type: ReportTypes,
            use_threads: bool = True,
            **reader_kwargs
    ):
        report_stats = ReportStats.from_report(file, report_type, use_threads, **reader_kwargs)
