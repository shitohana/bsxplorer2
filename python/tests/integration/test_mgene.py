import pytest

from bsx2.io import RegionReader
from bsx2.metagene import Metagene
from bsx2.types import HcAnnotStore, Contig
from shared import RUSTTESTS_DIR


@pytest.fixture(scope="session")
def annotation() -> HcAnnotStore:
    return HcAnnotStore.from_gff(RUSTTESTS_DIR / "data" / "annot.gff")


@pytest.fixture(scope="session")
def region_reader() -> RegionReader:
    return RegionReader(str(RUSTTESTS_DIR / "data" / "report.bsx"))


@pytest.fixture(scope="session")
def contigs(annotation: HcAnnotStore, region_reader: RegionReader) -> list[Contig]:
    contigs = [entry.contig for _, entry in annotation.iter()]
    contigs = region_reader.index().sort(contigs)
    return contigs


@pytest.fixture(scope="session")
def metagene(contigs, region_reader: RegionReader) -> Metagene:
    mgene = Metagene()
    for batch, contig in zip(region_reader.iter_contigs(contigs), contigs):
        positions, density = batch.normalized()
        mgene.insert(str(contig), positions, density)

    return mgene


def test_basic_creation(metagene: Metagene) -> None:
    pass

def test_gather(metagene: Metagene) -> None:
    metagene.gather()