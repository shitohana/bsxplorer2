{ pkgs, lib, config, inputs, ... }:

let
  rustBin = inputs.rust-overlay.lib.mkRustBin { } pkgs;
  nightly = rustBin.nightly."2025-11-03".minimal.override {
    extensions = [ "rustfmt" "llvm-tools-preview" ];
  };
  root = lib.escapeShellArg config.devenv.root;
  python = "${pkgs.python312}/bin/python3.12";
  venv = "${config.devenv.root}/python/.venv";
  # Tasks and shell commands share implementations. Cargo and uv own packages.
  commands = {
    sync = ''
      uv sync --project python --python ${python} --locked --no-install-project --inexact
    '';
    develop = ''
      uv pip install --python "$VIRTUAL_ENV/bin/python" --no-deps --no-build-isolation --config-settings 'build-args=--locked --profile pydev' --editable python
    '';
    release = ''
      uv pip install --python "$VIRTUAL_ENV/bin/python" --no-deps --no-build-isolation --config-settings 'build-args=--locked --profile pyrelease' --editable python
    '';
    pythonTest = "python -m pytest python/tests";
    typecheck = ''
      ty check python --python "$VIRTUAL_ENV/bin/python"
      python -m mypy.stubtest bsx2._bsx2 --concise
      python -m mypy.stubtest bsx2.io bsx2.types --concise
    '';
    build = ''
      rm -f dist/bsx2-*.whl dist/bsx2-*.tar.gz
      maturin build --locked --manifest-path python/Cargo.toml --profile pyrelease --auditwheel repair --interpreter "$VIRTUAL_ENV/bin/python" --out dist
      maturin sdist --manifest-path python/Cargo.toml --out dist
    '';
    packageTest = "python python/scripts/check_package.py";
    lint = "ruff check python";
    pythonFormat = "ruff format python";
    pythonFormatCheck = "ruff format --check python";
    rustFormat = "RUSTFMT=${nightly}/bin/rustfmt ${nightly}/bin/cargo-fmt fmt --all -- --config-path rustfmt.toml";
    rustFormatCheck = "RUSTFMT=${nightly}/bin/rustfmt ${nightly}/bin/cargo-fmt fmt --all -- --check --config-path rustfmt.toml";
    rustFast = "cargo nextest run --locked --profile unit";
    rustTest = "cargo test --locked --workspace --all-targets";
    rustDocs = "cargo test --locked --doc -p bsxplorer2";
    coverage = ''
      export PATH="${nightly}/bin:$PATH"
      cargo llvm-cov --locked --html --package bsxplorer2
    '';
  };
  inProject = command: ''
    set -euo pipefail
    cd ${root}
    export VIRTUAL_ENV=${lib.escapeShellArg venv}
    export PATH="$VIRTUAL_ENV/bin:$PATH"
    ${command}
  '';
  task = command: after: { exec = inProject command; inherit after; };
  script = command: { exec = inProject command; };
  # Keep the full check ordered: build distributions before validating them.
  full = lib.concatStringsSep "\n" [
    commands.develop commands.lint
    commands.pythonTest commands.typecheck commands.build commands.packageTest
    commands.rustTest commands.rustDocs
  ];
  rustHook = stage: {
    enable = true;
    entry = toString (pkgs.writeShellScript "rust-test-${stage}" (inProject ''
      case ${lib.escapeShellArg stage} in
        pre-push) target_ref="''${PRE_COMMIT_REMOTE_BRANCH:-}" ;;
        pre-merge-commit) target_ref=$(git symbolic-ref --quiet HEAD) ;;
      esac
      if [ "$target_ref" != refs/heads/master ]; then
        echo "Skipping Rust tests: target is not master."
        exit 0
      fi
      exec rs-test-fast
    ''));
    language = "system";
    pass_filenames = false;
    always_run = true;
    require_serial = true;
    stages = [ stage ];
  };
in
{
  languages.rust = {
    enable = true;
    channel = "stable";
    version = "1.91.0";
    components = [ "rustc" "cargo" "clippy" "rustfmt" "rust-analyzer" "rust-src" ];
  };
  languages.python = {
    enable = true;
    package = pkgs.python312;
    lsp.enable = false;
    uv.enable = true;
  };
  packages = with pkgs; [
    git git-lfs direnv cargo-nextest cargo-llvm-cov pkg-config cmake xz
  ];
  env.UV_PROJECT_ENVIRONMENT = lib.mkForce venv;

  tasks = {
    "python:sync" = (task commands.sync [ ]) // { before = [ "devenv:enterShell" ]; };
    "python:setup" = task commands.develop [ "python:sync" ];
    "python:release" = task commands.release [ "python:sync" ];
    "python:test" = task commands.pythonTest [ "python:setup" ];
    "python:typecheck" = task commands.typecheck [ "python:setup" ];
    "python:lint" = task commands.lint [ "python:sync" ];
    "python:format" = task commands.pythonFormat [ "python:sync" ];
    "python:format-check" = task commands.pythonFormatCheck [ "python:sync" ];
    "python:build" = task commands.build [ "python:sync" ];
    "python:package-test" = task commands.packageTest [ "python:build" ];
    "rust:test-fast" = task commands.rustFast [ ];
    "rust:test" = task commands.rustTest [ ];
    "rust:doc-test" = task commands.rustDocs [ ];
    "rust:format" = task commands.rustFormat [ ];
    "rust:format-check" = task commands.rustFormatCheck [ ];
    "format:all" = task commands.rustFormat [ "python:format" ];
    "format:check" = task commands.rustFormatCheck [ "python:format-check" ];
    "check:full" = task full [ "python:sync" ];
    "rust:coverage" = task commands.coverage [ "check:full" ];
  };

  scripts = {
    py-setup = script commands.develop;
    py-test = script (commands.develop + commands.pythonTest + '' "$@"'');
    py-typecheck = script (commands.develop + commands.typecheck);
    py-build = script commands.build;
    py-package-test = script (commands.build + commands.packageTest);
    py-lint = script commands.lint;
    py-format = script commands.pythonFormat;
    py-format-check = script commands.pythonFormatCheck;
    rs-test-fast = script commands.rustFast;
    rs-format = script commands.rustFormat;
    rs-format-check = script commands.rustFormatCheck;
    test-full = script full;
    rs-coverage = script ''
      test-full
      ${commands.coverage}
    '';
    bsxplorer = script ''cargo run --locked -p bsxplorer-ci --bin bsxplorer -- "$@"'';
    py-develop = script ''
      case "''${1:-debug}" in
        debug) ${commands.develop} ;;
        release) ${commands.release} ;;
        *) echo "Expected debug or release" >&2; exit 2 ;;
      esac
    '';
    doctor = script ''
      for tool in cargo rustc rust-analyzer uv git prek; do
        command -v "$tool"
        "$tool" --version
      done
      ${nightly}/bin/rustfmt --version
      cargo nextest --version
      cargo llvm-cov --version
      python --version
      for tool in maturin ty ruff; do
        "$tool" --version
      done
      python - <<'PY'
      import importlib.metadata
      import importlib.util
      import sys
      print("Python:", sys.executable)
      for name in ("pytest", "mypy"):
          print(f"{name}: {importlib.metadata.version(name)}")
      spec = importlib.util.find_spec("bsx2")
      if spec is None:
          print("Extension: not installed; run py-setup")
      else:
          from bsx2 import _bsx2
          print("Extension:", _bsx2.__file__)
      PY
    '';
  };

  git-hooks = {
    package = pkgs.prek;
    configPath = ".devenv/git-hooks.yaml";
    default_stages = [ "pre-commit" ];
    hooks = {
      check-toml.enable = true;
      check-yaml.enable = true;
      check-added-large-files = {
        enable = true;
        stages = [ "pre-commit" ];
      };
      detect-private-keys.enable = true;
      typos.enable = true;
      rust-format = {
        enable = true;
        entry = "rs-format";
        language = "system";
        types = [ "rust" ];
        pass_filenames = false;
        require_serial = true;
        stages = [ "pre-commit" "pre-merge-commit" ];
      };
      rust-test = rustHook "pre-push";
      rust-test-merge = rustHook "pre-merge-commit";
    };
  };
  enterShell = ''
    export VIRTUAL_ENV=${lib.escapeShellArg venv}
    export PATH="$VIRTUAL_ENV/bin:$PATH"
    echo "BSXplorer2: doctor, py-setup, py-test, py-typecheck, test-full"
  '';
  enterTest = "test-full";
}
