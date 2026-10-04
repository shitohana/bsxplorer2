# bsx2 Python package

The `bsx2` package combines Python helpers with the `bsx2._bsx2` PyO3 extension.
Builds from this repository use the local `bsxplorer2` Rust core.

## Development setup

Install Nix and devenv. CPython 3.10–3.12 is supported; the development
environment uses Python 3.12. From the repository root, run:

```sh
devenv shell
py-setup
```

This installs the dependencies pinned in `python/uv.lock` into `python/.venv`
and rebuilds the extension with the `pydev` Cargo profile. Repeat after Rust
changes. Run Python commands with `python` inside devenv.

```python
from bsx2.types import BsxBatch, Strand, Contig
from bsx2.io import BsxFileReader, ReportReader

batch = BsxBatch.empty()
assert batch.is_empty()
```

## Tests and typing

```sh
py-test
py-typecheck
test-full
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
py-build
```

Devenv writes wheels and a source distribution to `dist/`, replacing previous
`bsx2` distributions. The Python release profile uses panic unwinding for PyO3. The source distribution includes the local
Rust core and workspace `Cargo.lock`, so installation does not silently substitute
a registry release. Repository wheel/editable builds use Cargo `--locked`.
For source distributions, Cargo may prune unused CLI entries from the archived
workspace lockfile; the distribution check rejects changed dependency versions.
Wheel builds bundle required external libraries using maturin's repair option.
`py-develop debug`, `py-develop release`, and `py-build` are available.
See [the development guide](../DEVELOPMENT.md) for tasks, hooks, and direnv.

```sh
py-package-test
```

This builds both distributions, installs the wheel in a temporary environment,
and checks imports, shipped type files, runtime stubs, and pytest outside the checkout. It then
extracts the source distribution, rebuilds its wheel, and repeats the checks.
Dependencies come from the uv lockfile; archive builds use a dedicated
Cargo cache in `target/package-check/`. The first release build can take several minutes. Validation
in this session is limited to local macOS arm64 with Python 3.12.

To update dependencies intentionally, run `uv lock --project python`, review the
lockfile change, and repeat installation and validation. Normal commands require
an unchanged lockfile.
