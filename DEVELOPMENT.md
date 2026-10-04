# Development environment

All development uses devenv. Install Nix and devenv (2.1 or newer), then run
from the repository root:

```sh
devenv shell
doctor
py-setup
```

Devenv provides Rust 1.91.0, Python 3.12, uv, nextest, coverage tooling,
Git/Git LFS, and native build tools. Rust formatting and coverage use the pinned
2025-11-03 nightly directly, without requiring rustup. Cargo installs Rust
dependencies from `Cargo.lock`; uv installs Python dependencies and development
tools from `python/uv.lock`. Nix does not install application package dependencies.
`devenv.lock` pins the environment inputs. Review lock changes when updating tools.

Shell entry synchronizes `python/.venv` and installs Git hooks. It does not build
the extension. `py-setup` installs the local extension with the `pydev` profile;
Python test/type commands rebuild it before checking the installed package.
`doctor` prints tool versions, executable paths, and the extension location.

## Commands

The Justfile is replaced by commands available inside `devenv shell`:

- `py-setup`, `py-test [pytest arguments]`, `py-typecheck`
- `py-lint`, `py-format`, `py-format-check`
- `py-build`, `py-package-test`, `py-develop [debug|release]`
- `rs-test-fast`, `rs-format`, `rs-format-check`, `rs-coverage`
- `test-full`, `bsxplorer [CLI arguments]`, `doctor`

Tasks can also run directly without an interactive shell:

```sh
devenv tasks run python:test
devenv tasks run python:typecheck
devenv tasks run rust:test-fast
devenv tasks run format:check
devenv tasks run check:full
```

Other tasks include `python:setup`, `python:release`, `python:build`,
`python:package-test`, `python:lint`, `python:format`, `python:format-check`,
`rust:test`, `rust:doc-test`, `rust:format`, `rust:format-check`, and
`rust:coverage`. `format:all` formats Rust and Python; `format:check` modifies neither.
Ruff uses the uv development environment. Rust formatting uses `rustfmt.toml`.

`check:full`/`test-full` run Python lint, editable-package
tests/type checks, wheel/sdist validation, Rust workspace tests, and core doctests.
Distribution outputs use `dist/`. Each build removes previous `bsx2` wheels and
source archives so validation uses exactly the current build. Distribution checks
use `target/package-check/` for archive builds and never clean development
artifacts. The pinned Nix SDK currently produces macOS 14 arm64 wheels.
Coverage runs the full checks first. Run `devenv test` to also run the
configured commit hooks; the Rust formatting hook can modify Rust sources.

## Git hooks

Hook declarations live in `devenv.nix`. Devenv generates the ignored
`.devenv/git-hooks.yaml` and installs commit, push, and merge hooks using prek.
Existing Git LFS hooks are retained by the hook runner.

Commit checks validate TOML/YAML, large files, private keys, and spelling.
Rust formatting runs only when Rust files change, on commits and merges.
Nextest's `unit` profile runs only on pushes targeting `refs/heads/master` or
merge commits on `master`, including documentation-only changes. Other branches
skip the test gate. Fast-forward/conflicted-merge limitations of Git's merge hook
still apply. The existing coverage CI targets `main` and is unchanged.
Hooks call the same shell commands as manual checks and require an active
devenv environment.

## Optional automatic activation

The checked-in `.envrc` enables direnv integration. Add the direnv hook to your
shell (for zsh, `eval "$(direnv hook zsh)"`), make devenv and direnv available in
your normal PATH, and run `direnv allow` in this repository. This is optional;
`devenv shell` and `devenv tasks run` work without it.

## Intentional updates

Use `devenv update` for environment inputs, `uv lock --project python` for Python
dependencies, and Cargo's update commands for Rust dependencies. Normal setup
and build commands require existing lockfiles. Python 3.10–3.12 remains the
package support range; this environment selects 3.12. Development checks on
macOS do not establish Linux or distributable manylinux-wheel compatibility.
