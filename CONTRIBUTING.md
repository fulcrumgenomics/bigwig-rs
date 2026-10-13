# Contributing to bigwig

This document is for developers working on `bigwig` itself. User documentation is in [README.md](README.md) and the [API docs](https://docs.rs/bigwig).

## Toolchain and tools

`rust-toolchain.toml` pins the Rust release CI uses, with rustfmt and clippy; rustup installs it on the first build. The minimum supported Rust version is the `rust-version` in `Cargo.toml` (1.85), checked with that toolchain directly. Nightly is needed only for fuzzing.

```console
rustup toolchain install 1.85 nightly
cargo install cargo-nextest cargo-deny cargo-fuzz
```

`cargo-release` is needed only to cut a release (see [Releasing](#releasing)).

## Checks

CI runs the following on every pull request and push to `main` (`.github/workflows/check.yml`). The `ci-*` commands are aliases in `.cargo/config.toml`, so a local run is the same as CI's.

| Check | Command |
|---|---|
| Formatting (`cargo fmt --all` fixes it) | `cargo ci-fmt` |
| Clippy, pedantic, warnings as errors | `cargo ci-lint` and `cargo ci-lint-zlib` |
| Unit, integration and property tests | `cargo ci-test` and `cargo ci-test-zlib` |
| Doc tests (nextest does not run them) | `cargo ci-doctest` |
| API docs build without warnings | `cargo doc --all-features --no-deps --locked` |
| Advisories, licences, bans and sources | `cargo deny check` |
| Minimum supported Rust version | `cargo +1.85 check --all-features --locked` |

Lint and tests run twice because the two compression backends are exclusive: with both features on, libdeflate is used and the zlib-rs code is compiled out, so the `-zlib` aliases build with zlib-rs alone. `rustdocflags` in `.cargo/config.toml` makes `cargo doc` fail on any warning.

Every pull request also runs the [interop matrix](#interop-matrix) (`.github/workflows/interop.yml`), which must pass before merging.

## Tests

- **Unit tests** live in each module's `mod tests` under `src/`. Most write and read files in memory (`Cursor<Vec<u8>>`); `tempfile` is used only for the path-based API and the zoom spill file.
- **`tests/issues.rs`** holds regression tests derived from bugs reported against other wig and bigWig implementations, written through the public API.
- **`tests/properties.rs`** holds proptest round trips: random chromosome lists, data and writer options are written, read back and compared with a brute-force model. Each property runs 64 cases by default (fewer in the zlib-rs build); set `PROPTEST_CASES` for a longer soak, e.g. `PROPTEST_CASES=2000 cargo nextest run --all-features --test properties`. Failures are not persisted: a failing case is printed shrunk, and each bug it exposes gets a deterministic regression test at the end of the file.

## Fuzzing

`fuzz/` holds cargo-fuzz targets: `read_bigwig` opens arbitrary bytes as a bigWig and runs every read path; `parse_bedgraph` and `parse_wig` parse arbitrary text; `bedgraph_round_trip` checks that every bedGraph record that parses writes and re-parses unchanged; and `write_then_read` replays random push sequences, valid and invalid, against a model of the writer's validation rules and checks that the file reads back as the accepted items. `.github/workflows/fuzz.yml` runs each for a few minutes weekly and on demand; fuzzing is not part of the pull-request checks.

Run them from `fuzz/`. Every command names `+nightly`, because `rust-toolchain.toml` pins a stable toolchain that would otherwise win there. Seeds for `read_bigwig` are generated, never checked in:

```console
cd fuzz
cargo +nightly run --release --bin gen_seeds -- corpus/read_bigwig
cargo +nightly fuzz run read_bigwig -- -max_total_time=300 -rss_limit_mb=2048 -timeout=10 -max_len=262144
```

A crash leaves its input under `fuzz/artifacts/<target>/`. Turn it into a deterministic unit or regression test before fixing it.

## Interop matrix

`interop/` is a separate, unpublished crate and harness that checks the files `bigwig` writes against kent's utilities, pyBigWig (libBigWig), bigtools, `@gmod/bbi`, IGV desktop and igv.js, and checks that `bigwig` reads exactly what kent's utilities, pyBigWig and bigtools write. `.github/workflows/interop.yml` runs it on every pull request. [interop/README.md](interop/README.md) describes the fixtures, checks and tolerances, the known differences between implementations, and how to install the other tools and run it locally with `interop/run.sh`.

The interop crate has its own lint and build commands:

```console
cargo clippy --manifest-path interop/Cargo.toml --all-targets --locked -- -D warnings
cargo build --release --locked --manifest-path interop/Cargo.toml
```

`interop/corpus.sh` is a one-time validation against public bigWig files listed in `interop/corpus.tsv`, run by hand and never in CI because third-party downloads are too brittle to gate on. See [Public-file corpus](interop/README.md#public-file-corpus).

## Benchmarks

`benches/bigwig.rs` holds criterion benchmarks on synthetic, seeded data: writing per-base fixedStep, bedGraph runs and sparse variableStep, reading values over random regions, and summaries.

```console
cargo bench --bench bigwig
BIGWIG_BENCH_BASES=100000000 cargo bench --bench bigwig -- write_fixed_step
```

Write sizes default to 1 Mbp and 10 Mbp; `BIGWIG_BENCH_BASES` takes a comma-separated list of other sizes. Benchmark on a quiet machine, and compare an optimisation against the numbers before it: it must not regress any benchmark.

## Conventions

- **Generate test data in the test.** No data files are committed: tests build their inputs programmatically, so reviewers can see what each one exercises, and the public-file corpus is fetched from its manifest.
- **Many small tests, each named after the behaviour it asserts**, e.g. `zoom_records_never_span_gaps_longer_than_a_bin`, rather than table-driven tests. Test behaviour through the API rather than internals, so tests survive a refactor.
- **A bug found in another implementation gets a regression test** in `tests/issues.rs` whose doc comment cites the issue in the form `/// project #N (<url>): what went wrong`. Where the bug is a defective file from another writer, test both that ours never writes the defect and, where it makes sense, that our reader handles it.
- **The writer is strict; the reader is lenient.** The writer validates every push before changing any state, and never writes what some reader mishandles. The reader accepts what other writers produce wherever the meaning is unambiguous, returns `Error::Corrupt` where it is not, and must never panic on any input.
- **No `unsafe`.** The crate is `#![forbid(unsafe_code)]`.
- **Doc comments** on every public item, saying what it does, its inputs and outputs, and how it fails. Comments explain why, not what.
- **CHANGELOG.** A user-facing change (new API or option, behaviour change, a bug fix a user could hit, a notable performance change) adds a bullet under `## [Unreleased]` in `CHANGELOG.md`, in [Keep a Changelog](https://keepachangelog.com/) style. Internal refactors and test- or CI-only changes need none.

## Dependencies

Dependencies are kept deliberately few: libdeflate or zlib-rs for compression, `thiserror` and `log`. Prefer the standard library, and justify any new dependency. When adding or bumping one, check what it brings into each configuration's tree, including native build steps (today only libdeflate's C, in the default build):

```console
cargo tree --locked -e normal,build --prefix none --no-dedupe | sort -u
cargo tree --locked -e normal,build --prefix none --no-dedupe --no-default-features --features zlib-rs | sort -u
```

and run `cargo deny check`.

## Releasing

Releases use [cargo-release](https://github.com/crate-ci/cargo-release), configured in `release.toml` and under `[package.metadata.release]` in `Cargo.toml`. Before a release, make sure CI and the interop matrix pass on `main`.

Then cut the release from a clean, up-to-date `main`. `cargo release` is a dry run unless given `--execute`:

```console
cargo release X.Y.Z
cargo release X.Y.Z --execute
```

This moves the CHANGELOG's `[Unreleased]` notes under the new version, commits `Bump version to X.Y.Z`, tags `vX.Y.Z`, pushes, and publishes to crates.io, which needs a crates.io token (`cargo login`). A published version can't be replaced, only yanked.
