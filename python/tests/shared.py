from pathlib import Path

PYTESTS_DIR = Path(__file__).parent
PYSRC_DIR = PYTESTS_DIR.parent
RUSTSRC_DIR = PYSRC_DIR.parent / "bsxplorer2"
RUSTTESTS_DIR = RUSTSRC_DIR / "tests"
