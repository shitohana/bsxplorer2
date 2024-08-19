from __future__ import annotations

import datetime
from collections import OrderedDict

import polars as pl
import pyarrow as pa
from progress.bar import Bar

polars2arrow = {
    pl.Int8: pa.int8(),
    pl.Int16: pa.int16(),
    pl.Int32: pa.int32(),
    pl.Int64: pa.int64(),
    pl.UInt8: pa.uint8(),
    pl.UInt16: pa.uint16(),
    pl.UInt32: pa.uint32(),
    pl.UInt64: pa.uint64(),
    pl.Float32: pa.float32(),
    pl.Float64: pa.float64(),
    pl.Boolean: pa.bool_(),
    pl.Binary: pa.binary(),
    pl.Utf8: pa.utf8(),
}
arrow2polars = {value: key for key, value in polars2arrow.items()}


def polars2arrow_convert(pl_schema: OrderedDict):
    if pl.Categorical in pl_schema.values():
        raise ValueError("You should write schema for Categorical manually")
    pa_schema = pa.schema([
        (key, polars2arrow[value]) for key, value in pl_schema.items()
    ])
    return pa_schema


def arrow2polars_convert(pa_schema: pa.Schema):
    pl_schema = OrderedDict()

    if not all(pa_type in arrow2polars.keys() for pa_type in pa_schema.types):
        raise KeyError("Not all field types have match in polars.")

    for name, pa_type in zip(pa_schema.names, pa_schema.types):
        pl_schema[name] = arrow2polars[pa_type]

    return pl_schema


class ReportBar(Bar):
    suffix = "%(progress2mb)d/%(max2mb)d Mb [%(elapsed_fmt)s | ETA: %(eta_fmt)s]"
    fill = "@"

    @property
    def progress2mb(self):
        return int(self.index) / (1024**2)

    @property
    def max2mb(self):
        return int(self.max) / (1024**2)

    @property
    def elapsed_fmt(self):
        return str(datetime.timedelta(seconds=self.elapsed))

    @property
    def eta_fmt(self):
        return str(datetime.timedelta(seconds=self.eta))
