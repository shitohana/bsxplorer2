# BSXplorer2 agent guidance

Applies throughout this repository. Follow the user's current task and preserve
unrelated worktree changes, staged changes, and stashes.

## Start here

- Read `ai/README.md`, then the relevant sections of `ai/project-overview.md`,
  `ai/documentation-audit.md`, and `ai/revival-plan.md` before substantive work.
- Treat their revision, branch, and test results as dated evidence. Verify the
  current checkout before using them as facts; roadmap items are not proof of
  working features. Work on the requested task rather than advancing the entire
  revival plan automatically.
- Inspect Git status before editing. Keep the task diff focused; do not restore,
  stage, or overwrite unrelated changes. Commit only when the user asks.
- Update relevant `ai/` notes as part of tasks that change behavior, verified
  status, or next steps. Date new verification evidence and preserve the
  distinction between historical findings and current results.

## Tools and exploration

- Follow `@RTK.md` (currently inherited from `~/RTK.md`). Prefix shell commands
  with `rtk`; use `rtk proxy <command>` when raw output is needed or a command is
  not supported directly.
- When `.codegraph/` exists, use the CodeGraph MCP tool or
  `rtk proxy codegraph explore "<symbol, file, or question>"` before text searches
  or reading source to understand or locate code. Reuse source already returned.
  If CodeGraph cannot answer the query, use targeted `rg` searches and reads.
- Without `.codegraph/`, skip CodeGraph. Do not create or rebuild its index
  unless requested. Prefer `rg` and `rg --files` for text and file searches.

## Architecture and scope

- `bsxplorer2/` owns the Rust domain model, BSX/report I/O, indexing, genomic
  queries, annotations, aggregation, and reusable analysis primitives.
- `console/` is package `bsxplorer-ci`, binary `bsxplorer`; keep command
  orchestration here and reusable computation in the core.
- `python/` is package `bsxplorer2-py`, exposing the `bsx2._bsx2` extension and
  Python analysis layer. Check that binding tests exercise the intended core
  version: its manifest currently declares a version dependency, not a path.
- Preserve existing architecture and public contracts. Propose substantial
  redesigns before implementing them; avoid unrelated refactors, dependency
  upgrades, and feature additions. Evaluate unmerged branches explicitly before
  integrating their work.

## Scientific and data correctness

- Establish expected behavior from initialization, transformations, callers,
  and tests. Resolve disagreements between code and prose explicitly.
- Do not silently choose unresolved scientific contracts: coordinate origins
  and interval bounds, null versus zero coverage, density units, overflow,
  replicate aggregation, missingness filters, statistical tests, or output
  semantics. Ask when the requested work depends on an undecided contract.
- Preserve schema types and encodings unless their change is part of the task.
  Check affected format, CLI, and Python boundaries when behavior changes.
- For algorithm repairs, use small synthetic fixtures with independently known
  results. Cover relevant missing-data, zero-coverage, boundary, and chromosome
  transition cases. Do not weaken assertions merely to make tests pass.
- Support scientific claims with implementation evidence and primary sources.
  Distinguish implemented, tested, scientifically validated, and released work.
- Profile before optimizing. Save reproducible benchmark commands, versions,
  dataset checksums, hardware, repetitions, and raw results. Separate index
  construction from reused-index queries and bound claims to measured conditions.

## Verification and reporting

- Run focused checks for the change first; broaden checks when its dependencies
  or affected interfaces warrant it. Plain Cargo commands default to the core
  crate, not the whole workspace.
- Rust integration checks: `rtk cargo test --workspace --all-targets`.
  Core documentation checks: `rtk cargo test --doc -p bsxplorer2`.
- Check Rust formatting with the pinned nightly and repository configuration:
  `rtk proxy devenv tasks run rust:format-check`.
  `devenv tasks run format:all` formats the Rust workspace and Python;
  nextest profiles select core tests. Git hooks format only when Rust files
  change and run `cargo nextest run --locked --profile unit` only
  for pushes to `master` or merge commits into `master`.
- For CLI changes, exercise affected help and success/failure workflows via
  `rtk cargo run -p bsxplorer-ci --bin bsxplorer -- <arguments>`.
- For binding or Python changes, build/install the intended extension and test
  the installed package. Source-layout imports and Rust compilation alone do
  not verify the Python API. Use the declared Python environment and report
  missing tooling rather than claiming unrun checks passed.
- Recheck relevant historical failures; distinguish pre-existing failures from
  regressions. Report commands run, outcomes, and checks left unrun with reasons.
  Documentation-only work does not require a Rust/Python test run.
- Keep explanations concise and concrete. Report what changed, why, validation,
  and unresolved limitations; do not present formatting as functional testing.

## Communication style

- Use the [Google Developer Documentation Style Guide](https://developers.google.com/style)
  for user-facing conversation, progress updates, and explanations, especially
  its [voice and tone guidance](https://developers.google.com/style/tone).
- Write clearly, directly, and concisely in a conversational, friendly, and
  respectful tone. Use active voice, consistent terminology, and simple words;
  explain necessary technical terms. Avoid slang, filler, and excessive formality.
- Adapt the guide's principles to the language the user uses; it does not
  require switching the conversation to English.
