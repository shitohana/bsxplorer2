First metagene
==============

This tutorial shows the canonical AW25 metagene workflow from
``report.bsx`` and ``annot.gff``.

Annotation-driven assembly
--------------------------

.. code-block:: python

   from bsx2 import Context, HcAnnotStore, RegionReader
   from bsx2.viz import (
       AnnotProfileLayout,
       AnnotProfilePart,
       compute_from_annot,
       heatmap,
       line_plot,
   )

   reader = RegionReader("report.bsx")
   reader.clear_filters()
   reader.filter_context(Context.CG)

   annot = HcAnnotStore.from_gff("annot.gff")
   layout = AnnotProfileLayout(
       (
           AnnotProfilePart("up", 25, source="flank5", flank_bp=2000),
           AnnotProfilePart("body", 50, source="gene"),
           AnnotProfilePart("down", 25, source="flank3", flank_bp=2000),
       )
   )

   drd = compute_from_annot(reader, annot, layout=layout)
   line = line_plot(drd, segments=list(layout.segments), smooth=None)
   hm = heatmap(drd, segments=list(layout.segments), rank_rows=200)

Manual composed assembly
------------------------

The same metagene can be reproduced manually by preparing per-part contigs,
normalizing each part independently, and composing the parts into one profile.

That lower-level route is useful when the layout is already known or when the
inputs are not coming directly from annotation builders.

Expected result
---------------

For equivalent inputs and layout definitions, the annotation-driven and manual
composed routes should produce the same metagene-compatible
``DiscreteRegionData`` semantics.
