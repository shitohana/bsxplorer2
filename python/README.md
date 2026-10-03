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

## Tests and typing

```sh
just py-test
just py-typecheck
just test-full
```

`py-test` rebuilds and tests the installed extension. The pytest suite checks
access to every native export and its public members, plus synthetic batch,
annotation, region, iterator, and I/O behavior. `py-typecheck` runs ty and checks
the shipped stubs against runtime signatures with mypy stubtest. The package
includes `.pyi` files and `py.typed` for type checkers.

`test-full` runs Python tests, typing checks, distribution validation, Rust
workspace tests, and core doctests. `rs-test-fast` retains the fast Rust-only
nextest profile.

The existing Rust Polars adapter requires Python Polars 1.14.0, which is pinned
in package metadata and the lockfile. Native enum objects are PyO3 classes;
they are not Python `enum.Enum` subclasses. Schema methods return dictionaries.
Readers backed by memory maps accept paths and real file descriptors. Writers
also accept seekable binary streams such as `io.BytesIO`. Use writer context
managers or `close()` to finalize output and report any finalization errors.

## Distributions

```sh
just py-build
```

Wheels and a source distribution are written to `dist/`. The Python release
profile uses panic unwinding for PyO3. The source distribution includes the local
Rust core and workspace `Cargo.lock`, so installation does not silently substitute
a registry release. Repository wheel/editable builds use Cargo `--locked`.
For source distributions, Cargo may prune unused CLI entries from the archived
workspace lockfile; the distribution check rejects changed dependency versions.
Wheel builds bundle required external libraries using maturin's repair option.
`just maturin debug`, `just maturin release`, and `just maturin build` remain
available.

```sh
just py-package-test
```

This builds both distributions, installs the wheel in a temporary environment,
and checks imports, stubs, typing, and pytest outside the checkout. It then
extracts the source distribution, rebuilds its wheel, and repeats the checks.
Dependencies come from the uv lockfile; the source rebuild reuses the Cargo
dependency cache. The first release build can take several minutes. Validation
in this session is limited to local macOS arm64 with Python 3.12.

To update dependencies intentionally, run `uv lock --project python`, review the
lockfile change, and repeat installation and validation. Normal recipes require
an unchanged lockfile.
