Development
===========

Repository notes
----------------

The Python subproject already contains additional markdown notes in
``python/docs``:

- ``release_notes_aw25.md``
- ``performance_notes.md``
- ``plotting_support_matrix.md``
- ``upgrade_notes.md``

These notes remain useful as project records, while the Sphinx site is the
canonical user-facing documentation surface.

Local build
-----------

.. code-block:: bash

   pip install -r docs/requirements.txt
   pip install -e .
   sphinx-build -b html docs/source docs/build/html
