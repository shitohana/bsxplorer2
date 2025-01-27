import contextlib

with contextlib.suppress(ImportError):
    from .._lib import (
        RegionReader,
        BsxFileReader,
        IpcCompression,
        BsxIpcWriter
    )

__all__ = ["RegionReader", "BsxFileReader", "IpcCompression", "BsxIpcWriter"]