Input formats
=============

BSXplorer2 plotting workflows rely on two principal inputs:

``report.bsx``
    BSX methylation report consumed by :class:`bsx2.RegionReader`.

``annot.gff``
    GFF or GFF3 annotation consumed by :class:`bsx2.HcAnnotStore`.

The plotting studio and gallery examples assume that both inputs refer to the
same genome build and compatible contig naming.
