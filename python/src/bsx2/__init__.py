"""BSX2 Python package exports.

When the Rust extension is available, this module exposes the full public API.
Source-tree documentation and smoke checks may run without the compiled
extension; in that case lightweight subpackages such as ``bsx2.analysis`` can
still be imported while the extension-backed names are simply unavailable.
"""

try:
    from ._bsx2 import (
        BsxFileReader,
        Compression,
        FilterOperation,
        RegionReader,
        RegionReaderIterator,
        ReportReader,
        ReportWriter,
    )
    from .types import (
        AggMethod,
        BatchIndex,
        BsxBatch,
        BsxColumns,
        Context,
        ContextData,
        Contig,
        GenomicPosition,
        GffEntry,
        GffEntryAttributes,
        HcAnnotStore,
        HcAnnotStoreIterator,
        LazyBsxBatch,
        ReportTypeSchema,
        Strand,
    )
except ModuleNotFoundError as exc:
    if exc.name not in {"bsx2._bsx2", "beartype"}:
        raise
    __all__: list[str] = []
else:
    __all__ = [
        "Strand",
        "Context",
        "BsxColumns",
        "AggMethod",
        "BsxBatch",
        "ContextData",
        "ReportTypeSchema",
        "LazyBsxBatch",
        "GenomicPosition",
        "Contig",
        "GffEntryAttributes",
        "GffEntry",
        "HcAnnotStore",
        "HcAnnotStoreIterator",
        "BatchIndex",
        "BsxFileReader",
        "Compression",
        "FilterOperation",
        "RegionReader",
        "RegionReaderIterator",
        "ReportReader",
        "ReportWriter",
    ]
