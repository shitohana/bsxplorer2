from __future__ import annotations

from io import StringIO

from bsx2.viz.compute import input_compat as compat
from bsx2.viz.compute.input_compat import (
    build_seqname_compatibility_report,
    normalize_seqname,
    read_bsx_seqnames,
    read_gff_seqnames,
)


class _FakeReader:
    def __init__(self, path: str, seqnames: tuple[str, ...]):
        self.path = path
        self._seqnames = seqnames

    def chr_order(self) -> list[str]:
        return list(self._seqnames)


def _reader_factory(seqnames: tuple[str, ...]):
    def _factory(path: str) -> _FakeReader:
        return _FakeReader(path, seqnames)

    return _factory


def _gff_text(seqnames: list[str]) -> str:
    rows = []
    for idx, seqname in enumerate(seqnames, start=1):
        rows.append(
            "\t".join(
                [
                    seqname,
                    "test",
                    "gene",
                    str(idx * 100),
                    str(idx * 100 + 99),
                    ".",
                    "+",
                    ".",
                    f"ID=gene_{idx}",
                ]
            )
        )
    return "\n".join(rows) + "\n"


def test_normalize_seqname_is_soft() -> None:
    assert normalize_seqname("chr1") == "1"
    assert normalize_seqname("Chr2") == "2"
    assert normalize_seqname("MT") == "MT"
    assert normalize_seqname("mitochondrion") == "MT"
    assert normalize_seqname("NC_003070.9") == "NC_003070.9"


def test_read_bsx_seqnames_uses_region_reader_chr_order() -> None:
    assert read_bsx_seqnames("report.bsx", reader_factory=_reader_factory(("chr2", "chr1", "chr1"))) == (
        "chr1",
        "chr2",
    )


def test_build_seqname_compatibility_report_exact_match(monkeypatch) -> None:
    monkeypatch.setattr(
        compat,
        "read_gff_seqnames",
        lambda path: ("chr1", "chr2"),
    )

    report = build_seqname_compatibility_report(
        "report.bsx",
        "annot.gff",
        reader_factory=_reader_factory(("chr1", "chr2")),
    )

    assert report.status == "exact"
    assert report.exact_match is True
    assert report.normalized_match is True
    assert report.mapping_hints == ()


def test_build_seqname_compatibility_report_soft_normalized_match(monkeypatch) -> None:
    monkeypatch.setattr(
        compat,
        "read_gff_seqnames",
        lambda path: ("1", "2", "MT"),
    )

    report = build_seqname_compatibility_report(
        "report.bsx",
        "annot.gff",
        reader_factory=_reader_factory(("chr1", "chr2", "chrM")),
    )

    assert report.status == "normalized"
    assert report.exact_match is False
    assert report.normalized_match is True
    assert report.only_report_normalized == ()
    assert report.only_annot_normalized == ()
    assert ("chr1", "1") in report.mapping_hints
    assert ("chrM", "MT") in report.mapping_hints


def test_build_seqname_compatibility_report_detects_real_mismatch(monkeypatch) -> None:
    monkeypatch.setattr(
        compat,
        "read_gff_seqnames",
        lambda path: ("1", "2"),
    )

    report = build_seqname_compatibility_report(
        "report.bsx",
        "annot.gff",
        reader_factory=_reader_factory(("NC_000001.11", "NC_000002.12")),
    )

    assert report.status == "mismatch"
    assert report.normalized_match is False
    assert "NC_000001.11" in report.only_report_normalized
    assert "1" in report.only_annot_normalized


def test_read_gff_seqnames_supports_gzip(monkeypatch) -> None:
    opened: dict[str, str] = {}

    def fake_gzip_open(path, mode="rt", encoding="utf-8", errors="replace"):
        opened["path"] = str(path)
        return StringIO(_gff_text(["chr1", "chr2"]))

    monkeypatch.setattr(compat.gzip, "open", fake_gzip_open)
    assert read_gff_seqnames("annot.gff.gz") == ("chr1", "chr2")
    assert opened["path"].endswith("annot.gff.gz")
