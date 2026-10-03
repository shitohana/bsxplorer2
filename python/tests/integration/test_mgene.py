import pytest
from bsx2.io import RegionReader
from bsx2.metagene import Metagene
from bsx2.types import Contig, Strand


def test_from_reader_keeps_contig_names_when_regions_are_missing(bsx_path):
    reader = RegionReader(bsx_path)
    missing = Contig("chr1", 1, 5, Strand.Null)
    present = Contig("chr1", 10, 31, Strand.Null)
    metagene = Metagene.from_reader(reader, [missing, present])
    assert len(metagene) == 1
    assert metagene.get(str(missing)) is None
    assert metagene[str(present)][1] == pytest.approx([0.2, 0.4, 0.6])
    positions, densities = metagene.gather()
    assert positions == sorted(positions)
    assert densities == pytest.approx([0.2, 0.4, 0.6])
