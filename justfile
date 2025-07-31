coverage:
    cargo +nightly llvm-cov --html --package bsxplorer2

nextest:
    cargo nextest run --manifest-path bsxplorer2/Cargo.toml
