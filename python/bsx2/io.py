from bsx2 import _native # type: ignore

BsxFileReader = _native.io.BsxFileReader
BsxFileWriter = _native.io.BsxFileWriter
IpcCompression = _native.io.IpcCompression
Compression = _native.io.Compression
ReportReader = _native.io.ReportReader
ReportWriter = _native.io.ReportWriter
RegionReader = _native.io.RegionReader
FilterOperation = _native.io.FilterOperation
RegionReaderIterator = _native.io.RegionReaderIterator

__all__ = [
    "BsxFileReader",
    "BsxFileWriter",
    "IpcCompression",
    "Compression",
    "ReportReader",
    "ReportWriter",
    "RegionReader",
    "FilterOperation",
    "RegionReaderIterator"
]
