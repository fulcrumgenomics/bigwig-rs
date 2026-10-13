# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`BigWigWriter` and `WriterOptions`: a streaming bigWig writer with bounded memory.** It writes fixedStep, variableStep and bedGraph sections in kent's file layout, builds zoom levels in the same pass (spilling them to one temporary file beyond `zoom_memory_limit`), and validates every push before changing any state, so a rejected item leaves the writer usable. `finish()` is mandatory: an unfinished file keeps a zero magic number and is rejected by every reader. `create_with` requires an empty sink and rejects any other with an `InvalidInput` I/O error, so a reused sink can't leave stale bytes after the file.
- **`BigWigWriter::push_fixed_step_skip` writes dense values that mark missing data with NaN** (pyBigWig #132). Each NaN is a gap, and the file is exactly the one `push_fixed_step` writes when called once per run between gaps, so chunked calls write the same file as one call. Every other push still rejects NaN. The `BigWigWriter` docs give measured advice on choosing a push method for data with many gaps, where short fixedStep runs cost more than variableStep items.
- **`ChromIdOrder`** chooses chromosome IDs by first data arrival (the default, any push order) or by position in the caller's list (`Listed`, e.g. BAM header order), with data pushed in that order.
- **`ZoomPolicy`** chooses zoom levels: a data-derived ×4 ladder from a 256 bp first level by default (`Auto`), exact reductions (`Explicit`), or none.
- **`Compression`** chooses block compression, defaulting to deflate level 1. `Compression::MAX_LEVEL` reports the highest level the compiled backend supports.
- **`BigWigReader`: reads bigWig versions 1–4 in either byte order** from a path or any `Read + Seek` source. It provides info, chromosomes, streaming data items by region or for the whole file, per-base values into a caller buffer, exact or zoom-based summaries over a region or N bins (kent's level-choice rule), and raw zoom records. `reopen()` gives one handle per thread over shared metadata.
- **`bedgraph` and `wig` modules: streaming text readers and writers.** Wig positions are converted to 0-based on input, a `fixedStep` without `span=` uses `span = step` as kent does, and values are written in the shortest form that reads back to the same `f32`.
- Two compression backends behind features: `libdeflate` (default) and the pure-Rust `zlib-rs`.
- **`convert` module: streaming conversions between bigWig and text.** `bedgraph_to_bigwig` and `wig_to_bigwig` read sorted text once from any `BufRead` and write a bigWig, keeping wig blocks' section types and reporting what they did in a `Report` (items, skipped zero-length intervals, and the items `clip` truncated or dropped). An item the writer rejects becomes an `Error::AtLine` that gives the input line and wraps the writer's error, so callers can still match on the cause; text the readers cannot parse is an `Error::Parse`. A zero-length bedGraph line is skipped only if it would otherwise be valid; an unknown or finished chromosome, an overlap, a position past the end or a NaN is an error at its line. `bigwig_to_bedgraph` (optionally merging adjacent equal values) and `bigwig_to_wig` export a whole file, one chromosome or a region, chosen with `Selection`.
- **`ChromInfo::read_chrom_sizes`** parses a UCSC `chrom.sizes` file into a chromosome list, reporting a malformed line with its number.
