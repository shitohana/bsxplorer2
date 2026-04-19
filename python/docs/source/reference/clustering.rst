Clustering
==========

``bsx2.clustering`` contains the feature-matrix builders, backend policy, and
result models used by PCA, dendrogram, and cluster metagene workflows.

Recommended flow:

1. Build a gene profile matrix.
2. Run clustering with the configured backend.
3. Convert the result into plot-ready data.
4. Render embedding, dendrogram, or cluster metagene views through ``bsx2.viz``.

.. automodule:: bsx2.clustering
   :members:
   :undoc-members:
   :show-inheritance:
