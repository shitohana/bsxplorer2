import dataclasses
from collections import OrderedDict, UserDict
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ._batches import ReportTypes, UniversalBatch
from ._readers import UniversalReader


class ChromStats(UserDict):

    __init__ = UserDict.__init__

    def __add__(self, other):
        if isinstance(other, ChromStats):
            for other_chrom, other_data in other.items():
                if other_chrom not in self.data.keys():
                    self.data[other_chrom] = other_data
                else:
                    if other_data["start"] < self.data[other_chrom]["start"]:
                        self.data[other_chrom]["start"] = other_data["start"]
                    if other_data["end"] > self.data[other_chrom]["end"]:
                        self.data[other_chrom]["end"] = other_data["end"]
        else:
            raise NotImplementedError
        return self

    @property
    def chroms_names(self):
        return list(self.data.keys())

    @property
    def chroms_lengths(self):
        return list(map(lambda data: data["end"] - data["start"], self.data.values()))

    @classmethod
    def from_batch(cls, data: pl.DataFrame):
        chrom_stats = {
            chrom: dict(start=start, end=end) for chrom, start, end in (
                data.group_by("chr", maintain_order=True)
                .agg(
                    start=pl.first("position"),
                    end=pl.last("position")
                )
                .iter_rows()
            )
        }

        return cls(chrom_stats)


@dataclass
class MethylationStats:
    context_stats: pl.DataFrame = dataclasses.field(
        default_factory=lambda: pl.DataFrame(schema=dict(
            context=pl.String,
            sum=pl.Float64,
            count=pl.UInt64
        ))
    )

    trinuc_stats: pl.DataFrame = dataclasses.field(
        default_factory=lambda: pl.DataFrame(schema=dict(
            trinuc=pl.String,
            sum=pl.Float64,
            count=pl.UInt64
        ))
    )

    def __add__(self, other):
        if isinstance(other, MethylationStats):
            self.context_stats = (
                pl.concat([self.context_stats, other.context_stats.cast(self.context_stats.schema)])
                .group_by("context")
                .agg(
                    sum=pl.sum("sum"),
                    count=pl.sum("count")
                )
            )

            self.trinuc_stats = (
                pl.concat([self.trinuc_stats, other.trinuc_stats.cast(self.trinuc_stats.schema)])
                .group_by("trinuc")
                .agg(
                    sum=pl.sum("sum"),
                    count=pl.sum("count")
                )
            )
        else:
            raise NotImplementedError()

        return self

    @classmethod
    def from_batch(cls, data: pl.DataFrame):
        trinuc_stats = (
            data
            .group_by("trinuc")
            .agg(
                context=pl.first("context"),
                sum=pl.sum("density"),
                count=pl.len()
            )
        )

        context_stats = (
            trinuc_stats
            .group_by("context")
            .agg(
                sum=pl.sum("sum"),
                count=pl.sum("count")
            )
        )

        return cls(context_stats, trinuc_stats.drop("context"))


class ReportStats:
    __slots__ = {
        "methylation_stats": MethylationStats,
        "chrom_stats": ChromStats
    }

    @classmethod
    def from_report(
            cls,
            file: str | Path,
            report_type: ReportTypes,
            use_threads: bool = True,
            **reader_kwargs
    ):
        out = cls()
        with UniversalReader(file, report_type, use_threads, **reader_kwargs) as reader:
            for batch in reader:
                data = batch.data.filter(pl.col("count_total") != 0)

                out.methylation_stats += MethylationStats.from_batch(data)
                out.chrom_stats += ChromStats.from_batch(data)

        return out


    def __init__(
            self,
            methylation_stats: MethylationStats | None = None,
            chrom_stats: ChromStats | None = None
    ):
        self.methylation_stats = methylation_stats if methylation_stats is not None else MethylationStats()
        self.chrom_stats = chrom_stats if chrom_stats is not None else ChromStats()