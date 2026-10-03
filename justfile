default:
  just --list

format:
    cargo +nightly fmt --all -- --config-path rustfmt.toml

rs-coverage: rs-test-full
    cargo +nightly llvm-cov --html --package bsxplorer2

rs-test-full:
    cargo nextest run --profile full

rs-test-fast:
    cargo nextest run --profile unit

# Install locked dependencies and rebuild the local extension.
py-setup:
    uv sync --project python --python 3.12 --locked --no-install-project
    env -u CONDA_PREFIX VIRTUAL_ENV="{{justfile_directory()}}/python/.venv" python/.venv/bin/maturin develop --manifest-path python/Cargo.toml --profile pydev

py-build:
    uv sync --project python --python 3.12 --locked --no-install-project
    python/.venv/bin/maturin build --manifest-path python/Cargo.toml --profile pyrelease --interpreter python/.venv/bin/python --out dist
    python/.venv/bin/maturin sdist --manifest-path python/Cargo.toml --out dist

[doc('debug, release, build')]
maturin mode="debug":
    #!/usr/bin/env sh
    set -eu
    case "{{mode}}" in
      debug) just py-setup ;;
      release)
        uv sync --project python --python 3.12 --locked --no-install-project
        env -u CONDA_PREFIX VIRTUAL_ENV="{{justfile_directory()}}/python/.venv" python/.venv/bin/maturin develop --manifest-path python/Cargo.toml --profile pyrelease
        ;;
      build) just py-build ;;
      *) echo "Expected debug, release, or build" >&2; exit 2 ;;
    esac
