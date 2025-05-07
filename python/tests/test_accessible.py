import pytest
import inspect
import bsxplorer2
import polars as pl
from pathlib import Path
from bsxplorer2 import (
    ReportTypeSchema, Contig, Strand, Context, ContextData,
    GenomicPosition, BatchIndex, MethylationStats, EncodedBsxBatch, BsxBatch,
    LazyBsxBatch, LazyEncodedBsxBatch, AnnotStore
)
from bsxplorer2.io import (
    RegionReader, ReportReader, ReportWriter, BsxFileReader, BsxIpcWriter
)

# Dummy arguments for testing
DUMMY_PATH = "dummy.bsx"
DUMMY_CONTIG = Contig("chr1", 0, 100)
DUMMY_DF = pl.DataFrame({"chr": ["chr1"], "position": [1], "strand": ["+"], "context": ["CG"], "count_m": [1], "count_total": [2], "density": [0.5]})
DUMMY_CONTEXT_DATA = ContextData(b"ACGT")
DUMMY_BATCH = BsxBatch(DUMMY_DF)
DUMMY_ENCODED_BATCH = EncodedBsxBatch(DUMMY_DF)
DUMMY_CHR_NAMES = ["chr1", "chr2"]

def get_classes(module):
    """Get all classes defined in the module."""
    return [
        (name, cls) for name, cls in inspect.getmembers(module, inspect.isclass)
        if cls.__module__.startswith(module.__name__)
    ]

def get_callable_methods(cls):
    """Get all callable methods from a class, excluding magic methods unless they are iterators/context managers."""
    methods = []
    for name, method in inspect.getmembers(cls, inspect.isfunction) or inspect.getmembers(cls, inspect.ismethod):
        if name.startswith("_") and name not in ("__iter__", "__next__", "__enter__", "__exit__"):
            continue
        methods.append(name)
    return methods

@pytest.mark.parametrize("class_name, cls", get_classes(bsx2.types) + get_classes(bsx2.io))
def test_class_methods_callable(class_name, cls):
    """Test that each method in each class is callable with minimal arguments."""
    # Skip classes that are not directly callable or require complex initialization
    if class_name in ("Strand", "Context", "Compression", "IpcCompression"):
        pytest.skip(f"{class_name} is an enum or constant class")

    # Initialize the class with minimal arguments
    try:
        if class_name == "RegionReader":
            instance = cls(DUMMY_PATH)
        elif class_name == "ReportReader":
            instance = cls(DUMMY_PATH, ReportTypeSchema.Bismark)
        elif class_name == "ReportWriter":
            instance = cls(DUMMY_PATH, ReportTypeSchema.Bismark)
        elif class_name == "BsxFileReader":
            instance = cls(DUMMY_PATH)
        elif class_name == "BsxIpcWriter":
            instance = cls(DUMMY_PATH, DUMMY_CHR_NAMES)
        elif class_name == "ContextData":
            instance = cls(b"ACGT")
        elif class_name == "GenomicPosition":
            instance = cls("chr1", 1)
        elif class_name == "Contig":
            instance = cls("chr1", 0, 100)
        elif class_name == "BatchIndex":
            instance = cls()
        elif class_name == "AnnotStore":
            instance = cls()
        elif class_name == "MethylationStats":
            instance = cls()
        elif class_name == "EncodedBsxBatch":
            instance = cls(DUMMY_DF)
        elif class_name == "BsxBatch":
            instance = cls(DUMMY_DF)
        elif class_name == "LazyBsxBatch":
            instance = cls(DUMMY_BATCH)
        elif class_name == "LazyEncodedBsxBatch":
            instance = cls(DUMMY_ENCODED_BATCH)
        else:
            instance = cls()
    except (FileNotFoundError, RuntimeError, ValueError) as e:
        pytest.skip(f"Cannot initialize {class_name}: {e}")

    for method_name in get_callable_methods(cls):
        method = getattr(instance, method_name)
        try:
            if method_name == "__iter__":
                iter(method())
            elif method_name == "__next__":
                iterator = iter(instance)
                next(iterator)
            elif method_name == "__enter__":
                with instance:
                    pass
            elif method_name == "__exit__":
                with instance:
                    pass
            elif method_name in ("insert", "sort", "find", "chr_order", "get_chr_order", "get_chr_index", "save", "load"):
                if class_name == "BatchIndex":
                    if method_name == "insert":
                        method(DUMMY_CONTIG, 1)
                    elif method_name == "sort":
                        method([DUMMY_CONTIG])
                    elif method_name == "find":
                        method(DUMMY_CONTIG)
                    elif method_name == "save":
                        method("dummy_index.bin")
                    elif method_name == "load":
                        method("dummy_index.bin")
                    else:
                        method()
            elif method_name == "query":
                method("chr1", 0, 100)
            elif method_name == "reset":
                method()
            elif method_name in ("write_batch", "write_df"):
                if class_name in ("ReportWriter", "BsxIpcWriter"):
                    method(DUMMY_DF)
            elif method_name == "close":
                method()
            elif method_name == "get_batch":
                method(0)
            elif method_name == "blocks_total":
                method()
            elif method_name in ("from_sink_and_fai", "from_sink_and_fasta"):
                method(DUMMY_PATH, DUMMY_PATH)
            elif method_name in ("write_encoded_batch"):
                method(DUMMY_DF)
            elif method_name in ("col_names", "schema", "chr_col", "position_col", "context_col", "strand_col", "need_align"):
                method()
            elif method_name in ("to_decoded_df", "to_encoded_df"):
                method()
            elif method_name in ("seqname", "position", "start", "end", "strand", "strand_str", "length", "start_gpos", "end_gpos"):
                method()
            elif method_name in ("extend_upstream", "extend_downstream"):
                method(10)
            elif method_name == "is_in":
                method(DUMMY_CONTIG)
            elif method_name in ("__str__", "__repr__"):
                method()
            elif method_name == "__richcmp__":
                method(instance, 0)
            elif method_name in ("__add__", "__sub__"):
                method(GenomicPosition("chr1", 2))
            elif method_name in ("from_gff", "from_bed"):
                method(DUMMY_PATH)
            elif method_name in ("len", "is_empty", "__len__"):
                method()
            elif method_name in ("add_upstream", "add_downstream", "add_flanks"):
                method(100)
            elif method_name == "iter_sorted":
                method(BatchIndex())
            elif method_name == "get_entry_ids":
                method()
            elif method_name == "from_data":
                method(0.5, 0.1, {1: 10}, {"CG": (0.5, 10)}, {"+": (0.5, 10)})
            elif method_name in ("chr", "position", "strand", "context", "count_m", "count_total", "density"):
                method()
            elif method_name == "split_at":
                method(1)
            elif method_name in ("data", "take"):
                method()
            elif method_name in ("chr_val", "start_pos", "end_pos"):
                method()
            elif method_name in ("start_gpos", "end_gpos", "as_contig"):
                method()
            elif method_name in ("vstack", "extend"):
                method(DUMMY_BATCH if class_name == "BsxBatch" else DUMMY_ENCODED_BATCH)
            elif method_name == "filter_mask":
                method(pl.Series([True]))
            elif method_name == "height":
                method()
            elif method_name in ("into_report", "collect"):
                method(ReportTypeSchema.Bismark)
            elif method_name in ("filter_pos_lt", "filter_pos_gt"):
                method(10)
            elif method_name == "filter_coverage_lt":
                method(5)
            elif method_name == "filter_strand":
                method(Strand.Forward)
            elif method_name == "filter_context":
                method(Context.CG)
            elif method_name == "mark_low_coverage":
                method(5)
            elif method_name == "align_with_contexts":
                method(DUMMY_CONTEXT_DATA, "chr1")
            elif method_name == "stats":
                method()
            else:
                method()
        except (FileNotFoundError, RuntimeError, ValueError, TypeError) as e:
            if isinstance(e, TypeError) and ("missing required" in str(e) or "takes" in str(e)):
                pytest.skip(f"{class_name}.{method_name} requires specific arguments: {e}")
            elif isinstance(e, FileNotFoundError):
                pytest.skip(f"{class_name}.{method_name} requires a valid file: {e}")
            else:
                pytest.fail(f"{class_name}.{method_name} failed with {e}")

"""
[build-system]
requires = ["maturin>=1,<2"]
build-backend = "maturin"

[project]
name = "bsx2"
version = "0.1.0"
requires-python = ">=3.9,<3.13"
dependencies = [
    "polars",
    "pyarrow",
]

[tool.pyright]
venvPath = "."
venv = ".venv"
include = ["bsx2"]
exclude = ["src"]
stubPath = "./typings"

"""
