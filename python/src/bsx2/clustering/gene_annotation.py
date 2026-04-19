from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from beartype import beartype

from bsx2 import HcAnnotStore

from .config import AnnotationFormat
from .models import GeneAnnotation


@beartype
def infer_annotation_format(path: Path) -> AnnotationFormat:
    suffix = path.suffix.lower()
    if suffix in {".gff", ".gff3", ".gtf"}:
        return AnnotationFormat.GFF if suffix != ".gtf" else AnnotationFormat.GTF
    if suffix == ".bed":
        return AnnotationFormat.BED
    raise ValueError(f"Unsupported annotation format for path: {path}")


@beartype
def load_annotation_store(
    path: Path,
    annotation_format: AnnotationFormat | None = None,
) -> HcAnnotStore:
    fmt = annotation_format or infer_annotation_format(path)
    if fmt in {AnnotationFormat.GFF, AnnotationFormat.GTF}:
        return HcAnnotStore.from_gff(str(path))
    if fmt is AnnotationFormat.BED:
        return HcAnnotStore.from_bed(str(path))
    raise NotImplementedError(f"Unsupported annotation format: {fmt.value}")


def _entry_contig(entry) -> object | None:
    contig = getattr(entry, "contig", None)
    if contig is not None:
        return contig
    getter = getattr(entry, "get_contig", None)
    return getter() if callable(getter) else None


def _entry_feature_type(entry) -> str | None:
    feature_type = getattr(entry, "feature_type", None)
    try:
        return feature_type() if callable(feature_type) else feature_type
    except Exception:
        return None


def _entry_id(entry) -> str | None:
    entry_id = getattr(entry, "id", None)
    try:
        value = entry_id() if callable(entry_id) else entry_id
    except Exception:
        value = None
    normalized = _normalize_attr_value(value)
    if normalized:
        return normalized

    attrs = getattr(entry, "attributes", None)
    if attrs is None:
        return None
    for name in ("id", "alias", "name"):
        accessor = getattr(attrs, name, None)
        if accessor is None:
            continue
        try:
            attr_value = accessor() if callable(accessor) else accessor
        except Exception:
            continue
        normalized = _normalize_attr_value(attr_value)
        if normalized:
            return normalized
    return None


def _entry_name(entry) -> str | None:
    attrs = getattr(entry, "attributes", None)
    if attrs is None:
        return None
    for name in ("name", "alias"):
        accessor = getattr(attrs, name, None)
        if accessor is None:
            continue
        try:
            attr_value = accessor() if callable(accessor) else accessor
        except Exception:
            continue
        normalized = _normalize_attr_value(attr_value)
        if normalized:
            return normalized
    return None


def _normalize_attr_value(value) -> str | None:
    if value is None:
        return None
    if isinstance(value, list | tuple):
        if not value:
            return None
        return str(value[0])
    return str(value)


def _contig_strand(contig) -> str | None:
    strand = getattr(contig, "strand_str", None)
    try:
        return strand() if callable(strand) else strand
    except Exception:
        pass
    strand = getattr(contig, "strand", None)
    try:
        return str(strand() if callable(strand) else strand)
    except Exception:
        return None


def _contig_coord(contig, name: str) -> int | None:
    accessor = getattr(contig, name, None)
    try:
        value = accessor() if callable(accessor) else accessor
    except Exception:
        return None
    return None if value is None else int(value)


@beartype
def load_gene_annotations(
    path: Path,
    *,
    annotation_format: AnnotationFormat | None = None,
    min_gene_length_bp: int = 0,
    limit: int | None = None,
) -> list[GeneAnnotation]:
    store = load_annotation_store(path, annotation_format)
    annotation_format = annotation_format or infer_annotation_format(path)
    genes: list[GeneAnnotation] = []
    seen_counts: dict[str, int] = defaultdict(int)

    for item in store.iter():
        entry = item[1] if isinstance(item, tuple) and len(item) == 2 else item
        feature_type = _entry_feature_type(entry)
        if annotation_format is not AnnotationFormat.BED and feature_type != "gene":
            continue

        contig = _entry_contig(entry)
        if contig is None:
            continue

        chrom = str(getattr(contig, "seqname"))
        start = _contig_coord(contig, "start")
        end = _contig_coord(contig, "end")
        strand = _contig_strand(contig)
        if start is None or end is None or strand is None:
            continue
        if start < 0 or end <= start:
            continue

        gene_id = _entry_id(entry) or f"{chrom}:{start}-{end}"
        seen_counts[gene_id] += 1
        if seen_counts[gene_id] > 1:
            gene_id = f"{gene_id}#{seen_counts[gene_id]}"

        gene = GeneAnnotation(
            gene_id=gene_id,
            gene_name=_entry_name(entry),
            chrom=chrom,
            start=start,
            end=end,
            strand=strand,
        )
        if gene.length_bp < min_gene_length_bp:
            continue
        genes.append(gene)
        if limit is not None and len(genes) >= limit:
            break

    if not genes:
        raise ValueError(f"No genes were loaded from annotation: {path}")
    return genes
