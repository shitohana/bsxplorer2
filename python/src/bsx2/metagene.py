import heapq

import collections
import itertools
import math
from typing import Annotated, Callable

from beartype.typing import Hashable, Tuple, Sequence, Union, Iterator
from beartype import beartype
from beartype.vale import Is

from bsx2.io import RegionReader
from bsx2.types import Contig, BsxBatch
from src.bsx2 import _bsx2

# TODO: Move to specialized file
Fraction = Annotated[float, Is[lambda v: 0 <= v <= 1 or math.isnan(v)]]
Positions = Annotated[
    Sequence[Annotated[float, Is[lambda v: 0 <= v <= 1]]],
    Is[lambda seq: all(prev < current for prev, current in zip(seq[:-1], seq[1:]))]
]


class Metagene:
    def __init__(self):
        self.entries: dict[Hashable, Tuple[Positions, Sequence[Fraction]]] = dict()

    @classmethod
    @beartype
    def from_reader(
            cls,
            reader: RegionReader,
            contigs: Sequence[Contig],
            preprocess_fn: Callable[[BsxBatch], Tuple[Positions, Sequence[Fraction]]] = lambda batch: batch.normalized(),
    ) -> 'Metagene':
        contigs = reader.index().sort(list(contigs))
        new = cls()
        for batch, contig in zip(reader.iter_contigs(contigs), contigs):
            positions, density = preprocess_fn(batch)
            new.insert(str(contig), positions, density)

        return new

    @beartype
    def insert(self, name: Hashable, positions: Positions, density: Sequence[Fraction]):
        self.entries[name] = (positions, density)

    @beartype
    def remove(self, name: Hashable) -> Union[Tuple[Positions, Sequence[Fraction]], None]:
        return self.entries.pop(name, None)

    @beartype
    def get(self, name: Hashable) -> Union[Tuple[Positions, Sequence[Fraction]], None]:
        return self.entries.get(name, None)

    @beartype
    def union(self, other: 'Metagene'):
        self.entries |= other.entries

    def densities(self) -> Iterator[float]:
        """
        Chain all density arrays together

        Returns:
            Iterator over density values
        """
        return itertools.chain(*(density for _, density in self.entries.values()))

    def keys(self) -> collections.abc.Set:
        return self.entries.keys()

    def gather(self) -> Tuple[Positions, Sequence[Fraction]]:
        """
        Gather all positions and densities into two arrays

        Returns:
            Tuple of positions and densities lists
        """

        positions, densities = zip(*self.entries.values())
        # noinspection PyUnresolvedReferences
        return _bsx2.merge_metagene_values(positions, densities)

    def __len__(self):
        return len(self.entries)

    def __getitem__(self, item) -> Tuple[Positions, Sequence[Fraction]]:
        if (res := self.get(item)) is None:
            raise KeyError
        else:
            return res

    def __setitem__(self, key, value):
        if not isinstance(value, tuple):
            raise ValueError("Can't set entry with not a tuple")
        if len(value) > 2:
            raise ValueError("Value should be a tuple of [positions, densities]")
        self.insert(key, value[0], value[1])

    @beartype
    def __ior__(self, other: 'Metagene'):
        self.union(other)

