"""Check every exported native symbol and member against shipped typing metadata."""

import ast
import importlib
from pathlib import Path

import pytest
from bsx2 import _bsx2

STUB_PATH = Path(_bsx2.__file__).with_name("_bsx2.pyi")
STUB = ast.parse(STUB_PATH.read_text())
CLASSES = [n for n in STUB.body if isinstance(n, ast.ClassDef)]
FUNCTIONS = [n for n in STUB.body if isinstance(n, ast.FunctionDef)]


@pytest.mark.parametrize("node", CLASSES + FUNCTIONS, ids=lambda n: n.name)
def test_every_native_export_is_accessible(node):
    symbol = getattr(_bsx2, node.name)
    assert callable(symbol)
    assert symbol.__module__ == "bsx2._bsx2"


def test_export_inventory_is_complete():
    names = {n.name for n in CLASSES + FUNCTIONS}
    assert set(_bsx2.__all__) == names
    assert len(CLASSES) == 24
    assert len(FUNCTIONS) == 1


@pytest.mark.parametrize("node", CLASSES, ids=lambda n: n.name)
def test_every_native_member_has_a_stub(node):
    runtime_class = getattr(_bsx2, node.name)
    declared = {
        n.name if isinstance(n, (ast.FunctionDef, ast.ClassDef)) else n.target.id
        for n in node.body
        if isinstance(n, (ast.FunctionDef, ast.ClassDef, ast.AnnAssign))
    }
    runtime = {name for name in vars(runtime_class) if not name.startswith("_")}
    assert runtime <= declared
    for name in declared:
        assert hasattr(runtime_class, name), f"{node.name}.{name} is inaccessible"


@pytest.mark.parametrize("module_name", ["bsx2.io", "bsx2.types"])
def test_python_reexports_are_the_native_objects(module_name):
    module = importlib.import_module(module_name)
    for name in module.__all__:
        assert getattr(module, name) is getattr(_bsx2, name)
