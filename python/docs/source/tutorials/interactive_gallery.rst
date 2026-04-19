Interactive plot studio
=======================

The interactive plot studio lets you upload ``report.bsx`` and ``annot.gff``,
configure parameters per plot, and build figures on demand.

Run the app
-----------

From the ``python`` directory:

.. code-block:: bash

   streamlit run aw25_interactive_plot_studio.py

Available plots
---------------

- Metagene line plot
- Metagene heatmap
- Segment box plot
- Chromosome methylation map
- Gene PCA
- Gene dendrogram
- Cluster metagene

Caching
-------

Each plot is cached by the combination of:

.. code-block:: text

   report.bsx + annot.gff + plot_id + parameters + theme

If the same plot was already built, the cached Plotly JSON is reused.
