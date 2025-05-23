from bsx2 import _native # type: ignore

RegionReader = _native.io.RegionReader
Compression = _native.io.Compression
ReportReader = _native.io.ReportReader
ReportWriter = _native.io.ReportWriter
BsxFileReader = _native.io.BsxFileReader
IpcCompression = _native.io.IpcCompression
BsxFileWriter = _native.io.BsxFileWriter

__all__ = [
    "RegionReader",
    "Compression",
    "ReportReader",
    "ReportWriter",
    "BsxFileReader",
    "IpcCompression",
    "BsxFileWriter"
]
