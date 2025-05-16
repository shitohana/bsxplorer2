"""BSX2 - Python bindings for bsxplorer2 for DNA methylation analysis

This package provides Python bindings to the bsxplorer2 Rust library.
It offers high-performance data structures and algorithms for analyzing
bisulfite sequencing data.
"""

from . import types  #type: ignore
from . import io  #type: ignore

__version__ = "0.1.0"
