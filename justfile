default:
  just --list

format:
    cargo +nightly fmt --all -- --config-path rustfmt.toml

rs-coverage: test-full
    cargo +nightly llvm-cov --html --package bsxplorer2

test-full: py-test py-typecheck py-package-test
    cargo test --workspace --all-targets
    cargo test --doc -p bsxplorer2

rs-test-fast:
    cargo nextest run --profile unit

# Install locked dependencies and rebuild the local extension.
py-setup:
    uv sync --project python --python 3.12 --locked --no-install-project
    env -u CONDA_PREFIX VIRTUAL_ENV="{{justfile_directory()}}/python/.venv" python/.venv/bin/maturin develop --locked --manifest-path python/Cargo.toml --profile pydev

py-build:
    uv sync --project python --python 3.12 --locked --no-install-project --inexact
    python/.venv/bin/maturin build --locked --manifest-path python/Cargo.toml --profile pyrelease --auditwheel repair --interpreter python/.venv/bin/python --out dist
    python/.venv/bin/maturin sdist --manifest-path python/Cargo.toml --out dist

# Tests use the installed local extension, rebuilt before collection.
py-test: py-setup
    python/.venv/bin/python -m pytest python/tests

py-typecheck: py-setup
    python/.venv/bin/ty check python --python python/.venv/bin/python
    python/.venv/bin/python -m mypy.stubtest bsx2._bsx2 --concise
    python/.venv/bin/python -m mypy.stubtest bsx2.io bsx2.types --concise

# Verify wheel installation and a wheel rebuilt from the source distribution.
py-package-test: py-build
    python/.venv/bin/python python/scripts/check_package.py

[doc('debug, release, build')]
maturin mode="debug":
    #!/usr/bin/env sh
    set -eu
    case "{{mode}}" in
      debug) just py-setup ;;
      release)
        uv sync --project python --python 3.12 --locked --no-install-project
        env -u CONDA_PREFIX VIRTUAL_ENV="{{justfile_directory()}}/python/.venv" python/.venv/bin/maturin develop --locked --manifest-path python/Cargo.toml --profile pyrelease
        ;;
      build) just py-build ;;
      *) echo "Expected debug, release, or build" >&2; exit 2 ;;
    esac
