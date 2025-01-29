from enum import Enum
from ..misc import BsxBatch, EncodedBsxBatch, RegionCoordinates


class IpcCompression(Enum):
    ZSTD = "ZSTD"
    LZ4 = "LZ4"

class BsxIpcWriter:
    def __init__(self, path: str, fai_path: str, compression: IpcCompression | None = None) -> 'BsxIpcWriter':
        """
        Initialize BSXplorer IPC file writer

        :param path: Path where the report will be saved
        :param fai_path: FASTA index path
        :param compression: Compression of the output file (default = None)
        """

    def write_batch(self, batch: BsxBatch): ...
    def write_encoded_batch(self, encoded_batch: EncodedBsxBatch): ...
    def close(self):
        """
        Write file footer and close the writer.
        """

class BsxFileReader:
    def __init__(self, path: str) -> 'BsxFileReader': ...

    def __iter__(self) -> 'BsxFileReader': ...
    def __next__(self) -> 'EncodedBsxBatch': ...

class RegionReader:
    def __init__(self, path: str, regions: list[RegionCoordinates]) -> 'RegionReader': ...
    
    def __iter__(self) -> 'RegionReader': ...
    def __next__(self) -> tuple['RegionCoordinates', 'EncodedBsxBatch']: ...