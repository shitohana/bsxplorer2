# bsx2 Python package

The `bsx2` package combines Python helpers with the `bsx2._bsx2` PyO3 extension.
Builds from this repository use the local `bsxplorer2` Rust core.

## Development setup

Install Rust, uv, and just. CPython 3.10–3.12 is supported; the local recipes use
Python 3.12. From the repository root, run:

```sh
just py-setup
```

This installs the dependencies pinned in `python/uv.lock` into `python/.venv`
and rebuilds the extension with the `pydev` Cargo profile. Repeat after Rust
changes. Run Python commands with `python/.venv/bin/python`.

```python
from bsx2.types import BsxBatch, Strand, Contig
from bsx2.io import BsxFileReader, ReportReader

batch = BsxBatch.empty()
assert batch.is_empty()
```

## Distributions

```sh
just py-build
```

Wheels and a source distribution are written to `dist/`. The Python release
profile uses panic unwinding for PyO3. The source distribution includes the local
Rust core, so installation does not silently substitute a registry release.
`just maturin debug`, `just maturin release`, and `just maturin build` remain
available.

To update dependencies intentionally, run `uv lock --project python`, review the
lockfile change, and repeat installation and validation. Normal recipes require
an unchanged lockfile.
