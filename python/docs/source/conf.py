from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"

sys.path.insert(0, str(SRC))

project = "BSXplorer2"
author = "BSXplorer2 contributors"
copyright = f"{datetime.now().year}, {author}"

release = "0.1.0"
version = "0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "myst_parser",
    "sphinx_copybutton",
]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"
templates_path = ["_templates"]
exclude_patterns: list[str] = []

html_theme = "sphinx_rtd_theme"
html_title = "BSXplorer2 documentation"
html_short_title = "BSXplorer2"
html_static_path = ["_static"]

html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "includehidden": True,
    "titles_only": False,
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
