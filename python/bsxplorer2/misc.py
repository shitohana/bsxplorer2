import contextlib

with contextlib.suppress(ImportError):
    from ._lib import (
        RegionCoordinates,
        GenomicPosition,
        BsxBatch,
        EncodedBsxBatch,
        Context,
        Strand
    )
import io

# noinspection PyUnresolvedReferences
__all__ = [
    "RegionCoordinates",
    "GenomicPosition",
    "BsxBatch",
    "EncodedBsxBatch",
    "Context",
    "Strand"
]