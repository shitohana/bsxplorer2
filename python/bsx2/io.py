from bsx2 import _native # type: ignore

RegionReader = _native.io.RegionReader
Compression = _native.io.Compression
ReportReader = _native.io.ReportReader
ReportWriter = _native.io.ReportWriter
BsxFileReader = _native.io.BsxFileReader
IpcCompression = _native.io.IpcCompression
BsxIpcWriter = _native.io.BsxIpcWriter

__all__ = [
    "RegionReader",
    "Compression",
    "ReportReader",
    "ReportWriter",
    "BsxFileReader",
    "IpcCompression",
    "BsxIpcWriter"
]
