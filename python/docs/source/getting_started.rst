Getting started
===============

Installation
------------

For local development on platforms where the Rust extension can be built,
install the package from the ``python`` subproject root:

.. code-block:: bash

   pip install -e .

If you only want to build the documentation locally, installing the editable
package is not required. The Sphinx configuration already imports modules from
``python/src`` directly.

Install documentation dependencies:

.. code-block:: bash

   pip install -r docs/requirements.txt

Basic usage
-----------

.. code-block:: python

   from bsx2 import Context, RegionReader

   reader = RegionReader("report.bsx")
   reader.clear_filters()
   reader.filter_context(Context.CG)

Input files
-----------

BSXplorer2 plotting and clustering workflows use two main input files:

``report.bsx``
    Methylation report in BSX format.

``annot.gff``
    Genome annotation in GFF or GFF3 format.

Interactive plot studio
-----------------------

Run the local AW25 studio from the ``python`` directory:

.. code-block:: bash

   streamlit run aw25_interactive_plot_studio.py

Build the Sphinx site
---------------------

.. code-block:: bash

   PYTHONPATH=src sphinx-build -b html docs/source docs/build/html

Then open:

.. code-block:: text

   docs/build/html/index.html

Windows note
------------

Editable installation may currently fail on Windows because the Rust dependency
``bsxplorer2 = 0.2.3`` imports ``std::os::fd`` in upstream code. This affects
``pip install -e .`` but does not block local Sphinx builds as long as the
prebuilt extension binaries already present in ``python/src/bsx2`` are usable.
