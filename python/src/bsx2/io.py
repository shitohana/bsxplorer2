from __future__ import annotations

from typing import Any

from bsx2 import _bsx2 as x
from bsx2.guards import require_single_choice
from bsx2.validation import validate_context, validate_min_coverage, validate_strand

BsxFileReader = x.BsxFileReader
IpcCompression = x.IpcCompression
BsxFileWriter = x.BsxFileWriter
Compression = x.Compression
ReportWriter = x.ReportWriter
FilterOperation = x.FilterOperation
RegionReaderIterator = x.RegionReaderIterator
RegionReader = x.RegionReader

_CONTEXT_MAP = {
    "CG": x.Context.CG,
    "CHG": x.Context.CHG,
    "CHH": x.Context.CHH,
}
_STRAND_MAP = {
    "+": x.Strand.Forward,
    "-": x.Strand.Reverse,
}


def normalize_context(context: Any) -> x.Context:
    if isinstance(context, x.Context):
        return context
    normalized = validate_context(context)
    normalized = require_single_choice(
        normalized,
        name="context",
        allowed_hint="CG, CHG, or CHH",
    )
    return _CONTEXT_MAP[normalized]


def normalize_strand(strand: Any) -> x.Strand:
    if isinstance(strand, x.Strand):
        return strand
    normalized = validate_strand(strand, allow_both=False)
    return _STRAND_MAP[normalized]


def normalize_min_coverage(min_coverage: Any) -> int:
    return int(validate_min_coverage(min_coverage, integer=True))


def filter_context(target: Any, context: Any):
    return target.filter_context(normalize_context(context))


def filter_strand(target: Any, strand: Any):
    return target.filter_strand(normalize_strand(strand))


def filter_coverage_gt(target: Any, min_coverage: Any):
    return target.filter_coverage_gt(normalize_min_coverage(min_coverage))


__all__ = [
    "BsxFileReader",
    "IpcCompression",
    "BsxFileWriter",
    "Compression",
    "ReportWriter",
    "FilterOperation",
    "RegionReaderIterator",
    "RegionReader",
    "normalize_context",
    "normalize_strand",
    "normalize_min_coverage",
    "filter_context",
    "filter_strand",
    "filter_coverage_gt",
]
