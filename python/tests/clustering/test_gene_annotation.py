from __future__ import annotations

import shutil
import uuid
from pathlib import Path

from bsx2.clustering.config import AnnotationFormat
from bsx2.clustering.gene_annotation import infer_annotation_format, load_gene_annotations


def test_infer_annotation_format() -> None:
    assert infer_annotation_format(Path("genes.gff")) is AnnotationFormat.GFF
    assert infer_annotation_format(Path("genes.gtf")) is AnnotationFormat.GTF
    assert infer_annotation_format(Path("genes.bed")) is AnnotationFormat.BED


def test_load_gene_annotations_keeps_only_gene_entries() -> None:
    tmpdir = Path("python/tests/clustering") / f"_tmp_annot_{uuid.uuid4().hex}"
    tmpdir.mkdir(parents=True, exist_ok=False)
    try:
        annot = tmpdir / "genes.gff"
        annot.write_text(
            "\n".join(
                [
                    "chr1\t.\tgene\t10\t40\t.\t+\t.\tID=geneA;Name=GeneA",
                    "chr1\t.\tmRNA\t10\t40\t.\t+\t.\tID=txA;Parent=geneA",
                    "chr2\t.\tgene\t50\t90\t.\t-\t.\tID=geneB;Name=GeneB",
                ]
            ),
            encoding="utf-8",
        )

        genes = load_gene_annotations(annot)

        assert [gene.gene_id for gene in genes] == ["geneA", "geneB"]
        assert genes[0].gene_name == "GeneA"
        assert genes[1].strand == "-"
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
