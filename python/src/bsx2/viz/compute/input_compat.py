from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from bsx2 import RegionReader

if TYPE_CHECKING:
    from collections.abc import Iterable


_MITOCHONDRIAL_NAMES = {
    "M",
    "MT",
    "Mt",
    "m",
    "mt",
    "mitochondria",
    "mitochondrion",
    "Mitochondria",
    "Mitochondrion",
}


@dataclass(frozen=True)
class SeqnameCompatibilityReport:
    """Summarise seqname compatibility between a BSX report and an annotation."""

    report_raw_seqnames: tuple[str, ...]
    annot_raw_seqnames: tuple[str, ...]
    report_normalized_seqnames: tuple[str, ...]
    annot_normalized_seqnames: tuple[str, ...]
    only_report_raw: tuple[str, ...]
    only_annot_raw: tuple[str, ...]
    only_report_normalized: tuple[str, ...]
    only_annot_normalized: tuple[str, ...]
    mapping_hints: tuple[tuple[str, str], ...]

    @property
    def exact_match(self) -> bool:
        return not self.only_report_raw and not self.only_annot_raw

    @property
    def normalized_match(self) -> bool:
        return not self.only_report_normalized and not self.only_annot_normalized

    @property
    def status(self) -> str:
        if self.exact_match:
            return "exact"
        if self.normalized_match:
            return "normalized"
        return "mismatch"


def normalize_seqname(name: str) -> str:
    """Apply a soft seqname normalization suitable for compatibility checks."""

    normalized = str(name).strip()
    normalized = re.sub(r"^chr", "", normalized, flags=re.IGNORECASE)
    if normalized in _MITOCHONDRIAL_NAMES:
        return "MT"
    return normalized


def _open_maybe_gzip(path: str | Path):
    resolved = Path(path)
    if resolved.suffix.lower() == ".gz":
        return gzip.open(resolved, "rt", encoding="utf-8", errors="replace")
    return resolved.open("r", encoding="utf-8", errors="replace")


def read_bsx_seqnames(
    bsx_path: str | Path,
    *,
    reader_factory=RegionReader,
) -> tuple[str, ...]:
    """Read seqnames from a BSX container via :class:`bsx2.RegionReader`."""

    reader = reader_factory(str(bsx_path))
    return tuple(sorted({str(name) for name in reader.chr_order()}))


def read_gff_seqnames(gff_path: str | Path) -> tuple[str, ...]:
    """Read raw seqnames from the first column of a GFF/GFF3 file."""

    seqnames: set[str] = set()
    with _open_maybe_gzip(gff_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqnames.add(parts[0].strip())
    return tuple(sorted(seqnames))


def _normalized_groups(seqnames: Iterable[str]) -> dict[str, tuple[str, ...]]:
    groups: dict[str, set[str]] = {}
    for raw_name in seqnames:
        key = normalize_seqname(raw_name)
        groups.setdefault(key, set()).add(str(raw_name))
    return {key: tuple(sorted(values)) for key, values in groups.items()}


def build_seqname_compatibility_report(
    bsx_path: str | Path,
    gff_path: str | Path,
    *,
    reader_factory=RegionReader,
) -> SeqnameCompatibilityReport:
    """Compare seqnames between a BSX report and a GFF annotation."""

    report_raw_seqnames = read_bsx_seqnames(bsx_path, reader_factory=reader_factory)
    annot_raw_seqnames = read_gff_seqnames(gff_path)

    report_raw_set = set(report_raw_seqnames)
    annot_raw_set = set(annot_raw_seqnames)

    report_groups = _normalized_groups(report_raw_seqnames)
    annot_groups = _normalized_groups(annot_raw_seqnames)
    report_normalized_set = set(report_groups)
    annot_normalized_set = set(annot_groups)

    mapping_hints: set[tuple[str, str]] = set()
    for normalized_name in sorted(report_normalized_set & annot_normalized_set):
        for report_raw in report_groups[normalized_name]:
            for annot_raw in annot_groups[normalized_name]:
                if report_raw != annot_raw:
                    mapping_hints.add((report_raw, annot_raw))

    return SeqnameCompatibilityReport(
        report_raw_seqnames=report_raw_seqnames,
        annot_raw_seqnames=annot_raw_seqnames,
        report_normalized_seqnames=tuple(sorted(report_normalized_set)),
        annot_normalized_seqnames=tuple(sorted(annot_normalized_set)),
        only_report_raw=tuple(sorted(report_raw_set - annot_raw_set)),
        only_annot_raw=tuple(sorted(annot_raw_set - report_raw_set)),
        only_report_normalized=tuple(sorted(report_normalized_set - annot_normalized_set)),
        only_annot_normalized=tuple(sorted(annot_normalized_set - report_normalized_set)),
        mapping_hints=tuple(sorted(mapping_hints)),
    )
