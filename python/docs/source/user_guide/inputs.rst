Input formats
=============

BSXplorer2 plotting workflows rely on two principal inputs:

``report.bsx``
    BSX methylation report consumed by :class:`bsx2.RegionReader`.

``annot.gff``
    GFF or GFF3 annotation consumed by :class:`bsx2.HcAnnotStore`.

The plotting studio and gallery examples assume that both inputs refer to the
same genome build and compatible contig naming.

When files are uploaded into the interactive plot studio, the app runs a
seqname compatibility check between ``report.bsx`` and ``annot.gff``. The
result is reported as an exact match, a soft-normalized match (for example
``chr1`` versus ``1``), or a remaining mismatch that likely indicates
different assemblies or incompatible naming schemes.
