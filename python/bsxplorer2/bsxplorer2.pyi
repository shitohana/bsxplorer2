class RegionCoordinates:
    """
    A class, representing a genomic region [start, end].

    :param chr: name of the chromosome
    :param start: start of the region (included)
    :param end: end of the region (included)
    """

    def __init__(self, chr: str, start: int, end: int) -> 'RegionCoordinates': ...

    @property
    def chr(self) -> str: ...

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...

    @property
    def start_gpos(self) -> 'GenomicPosition':
        """
        :return: GenomicPosition instance
        """

    @property
    def end_gpos(self) -> 'GenomicPosition':
        """
        :return: GenomicPosition instance
        """
class GenomicPosition:
    def __init__(self, chr: str, position: int) -> 'GenomicPosition': ...
    def __add__(self, other) -> 'GenomicPosition': ...
    def __sub__(self, other) -> 'GenomicPosition': ...
    def __rshift__(self, other) -> 'RegionCoordinates':
        """
        Create a RegionCoordinates instance from two genomic positions.
        Fails if self and other chromosomes differ.

        :param other: GenomicPosition instance
        :return: 'RegionCoordinates'
        """
    @property
    def chr(self) -> str: ...

    @property
    def position(self) -> int: ...