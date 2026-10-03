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

[doc('debug, release, build')]
maturin mode="debug": rs-test-fast
    #!/usr/bin/env sh
    unset CONDA_PREFIX
    eval $(poetry env activate --project python)

    maturin \
        {{ if mode == "build" { "build" } else { "develop" } }} \
        --bindings pyo3 --manifest-path ./python/Cargo.toml \
        {{ if mode =~ "(release)|(build)" { "--release" } else { "" } }}
