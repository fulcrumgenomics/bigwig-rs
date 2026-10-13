<h1 align="center">bigwig</h1>

<p align="center">
  <a href="https://github.com/fulcrumgenomics/bigwig-rs/actions/workflows/check.yml"><img src="https://github.com/fulcrumgenomics/bigwig-rs/actions/workflows/check.yml/badge.svg" alt="Check status"></a>
  <a href="https://github.com/fulcrumgenomics/bigwig-rs/blob/main/LICENSE"><img src="https://img.shields.io/github/license/fulcrumgenomics/bigwig-rs.svg" alt="Licence"></a>
  <a href="https://crates.io/crates/bigwig"><img src="https://img.shields.io/crates/v/bigwig.svg?colorB=319e8c" alt="Version info"></a>
  <a href="https://docs.rs/bigwig"><img src="https://docs.rs/bigwig/badge.svg" alt="Documentation"></a>
  <a href="https://www.fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics-badge.svg" alt="Fulcrum Genomics"></a>
  <br>
</p>

A general-purpose Rust library for reading and writing [bigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) files, and the [wig](https://genome.ucsc.edu/goldenPath/help/wiggle.html) and [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) text formats, with conversions between them. It writes files in the layout UCSC's own tools produce, which open in UCSC, IGV, igv.js, JBrowse 2, libBigWig/pyBigWig and bigtools, and reads the files those tools write.

## About Fulcrum Genomics

`bigwig` is developed and maintained by [Fulcrum Genomics](https://www.fulcrumgenomics.com). [Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with bigwig and beyond.

<p>
<a href="https://www.fulcrumgenomics.com"><picture><source media="(prefers-color-scheme: dark)" srcset="https://github.com/fulcrumgenomics/bigwig-rs/raw/HEAD/.github/logos/fulcrumgenomics-dark.svg?sanitize=true"><source media="(prefers-color-scheme: light)" srcset="https://github.com/fulcrumgenomics/bigwig-rs/raw/HEAD/.github/logos/fulcrumgenomics-light.svg?sanitize=true"><img alt="Fulcrum Genomics" src="https://github.com/fulcrumgenomics/bigwig-rs/raw/HEAD/.github/logos/fulcrumgenomics-light.svg?sanitize=true" height="36" align="middle"></picture></a>&nbsp;&nbsp;&nbsp;
<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white" alt="Email us" align="middle"></a>&nbsp;&nbsp;
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white" alt="Visit us" align="middle"></a>
</p>

## Contents

- [Features](#features)
- [Installing](#installing)
- [Writing](#writing)
- [Reading](#reading)
- [Conversions](#conversions)
- [Compatibility](#compatibility)
- [Defaults and differences from kent's tools](#defaults-and-differences-from-kents-tools)
- [Compression backends](#compression-backends)
- [Contributing](#contributing)
- [Licence](#licence)

The API reference is on [docs.rs](https://docs.rs/bigwig).

## Features

- A streaming writer with bounded memory. It writes all three section types (fixedStep, variableStep and bedGraph) and builds the zoom levels in the same pass, spilling them to a temporary file beyond a configurable limit.
- Every push is validated before it changes anything, so a rejected item can be skipped and writing carries on.
- A reader for bigWig versions 1 to 4 in either byte order, from a path or any `Read + Seek` source: data items, per-base values, exact and zoom-based summaries over a region or in N bins, and raw zoom records. `reopen()` gives one handle per thread over shared metadata.
- Streaming bedGraph and wig readers and writers, and conversions in both directions that read text once, so input can come from a pipe.
- Synchronous, with no async runtime, no thread pool, no `unsafe`, and a small dependency tree.

## Installing

```toml
[dependencies]
bigwig = "0.1"
```

The default compression backend compiles libdeflate's bundled C source, so it needs a C compiler; [Compression backends](#compression-backends) describes the pure-Rust alternative. The minimum supported Rust version is 1.85.

## Writing

```rust
use bigwig::{ChromInfo, WriterOptions};

fn main() -> bigwig::Result<()> {
    let chroms = [ChromInfo::new("chr1", 248_956_422), ChromInfo::new("chrM", 16_569)];
    let mut writer = WriterOptions::new().create("sample.bw", &chroms)?;

    // Per-base depth, pushed in chunks as it is produced.
    writer.push_fixed_step("chr1", 0, 1, 1, &[3.0; 1_000])?;
    writer.push_fixed_step("chr1", 1_000, 1, 1, &[4.0; 500])?;
    // A bedGraph-style interval: [start, end) with a value.
    writer.push_interval("chrM", 100, 200, 0.5)?;

    // Mandatory: writes the indexes, zoom levels and the real header.
    writer.finish()?;
    Ok(())
}
```

Each chromosome's data must be pushed sorted, without overlaps, and all together; chromosomes may come in any order. A rejected push changes nothing, so the caller may skip it and carry on. `finish()` is mandatory: an unfinished file keeps a zero magic number and every reader rejects it.

The push method chooses the section type: `push_fixed_step` for dense, evenly spaced values such as per-base depth or bins, `push_variable_step` for scattered values, and `push_interval`/`push_intervals` for piecewise-constant data. `WriterOptions` sets the section size, zoom levels, compression, chromosome ID order and clipping at chromosome ends; its docs list every option and its default.

### Missing values

bigWig has no missing-value marker: a base without data is simply not stored. The writer rejects NaN and infinite values, so missing data is written as gaps. For dense values that mark missing data with NaN, `push_fixed_step_skip` treats each NaN as a gap and writes each run between them as fixedStep:

```rust
// Bases 2 and 3 have no data, so this writes two runs: [0, 2) and [4, 5).
writer.push_fixed_step_skip("chr1", 0, 1, 1, &[1.0, 2.0, f32::NAN, f32::NAN, 3.0])?;
```

Each run is a section of its own, so where the runs between gaps are short, `push_variable_step` or bedGraph intervals give a smaller file. The [`BigWigWriter` docs](https://docs.rs/bigwig/latest/bigwig/struct.BigWigWriter.html#choosing-a-push-method) give measured guidance on choosing a push method.

## Reading

```rust
use bigwig::{BigWigReader, SummarySource};

fn main() -> bigwig::Result<()> {
    let mut reader = BigWigReader::open("sample.bw")?;
    for chrom in reader.chroms() {
        println!("{}\t{}", chrom.name, chrom.length);
    }

    // One value per base, NaN where there is no data.
    let values = reader.values("chr1", 0, 100)?;

    // The stored items overlapping a region.
    for item in reader.items("chr1", 0, 10_000)? {
        let item = item?;
        println!("{}\t{}\t{}", item.start, item.end, item.value);
    }

    // Mean signal in 100 bins, from a zoom level where one is fine enough.
    let bins = reader.summaries("chr1", 0, 1_000_000, 100, SummarySource::BestZoom)?;
    let means: Vec<Option<f64>> = bins.iter().map(|bin| bin.map(|summary| summary.mean())).collect();
    Ok(())
}
```

`SummarySource::Raw` summaries are exact. `SummarySource::BestZoom` uses the zoom level kent's `bigWigSummary` would choose and is approximate, `min` and `max` included. The reader is lenient where other writers deviate from kent's layout (chromosome tree after the data, unsorted chromosome keys, an undersized decompression buffer, oversized section ends) and strict where ambiguity would give wrong answers.

## Conversions

The `convert` module streams between text and bigWig. Text input is read once from any `BufRead`, so it can be a pipe, but must be sorted.

```rust
use std::fs::File;
use std::io::{BufReader, BufWriter};

use bigwig::convert::{Selection, bedgraph_to_bigwig, bigwig_to_bedgraph};
use bigwig::{BigWigReader, ChromInfo, WriterOptions};

fn main() -> bigwig::Result<()> {
    let chroms = ChromInfo::read_chrom_sizes(BufReader::new(File::open("hg38.chrom.sizes")?))?;
    let input = BufReader::new(File::open("signal.bedGraph")?);
    let sink = File::create("signal.bw")?;
    let report = bedgraph_to_bigwig(input, &chroms, &WriterOptions::new(), sink)?;
    println!("{} items", report.items);

    let mut reader = BigWigReader::open("signal.bw")?;
    let out = BufWriter::new(File::create("round_trip.bedGraph")?);
    bigwig_to_bedgraph(&mut reader, &Selection::All, true, out)?;
    Ok(())
}
```

`wig_to_bigwig` and `bigwig_to_wig` do the same for wig, keeping each block's section type. `Selection` exports everything, one chromosome, or a region. An item the writer rejects is reported with its input line.

### Replacing kent's tools

This crate ships no command-line tools. The operations kent's utilities perform are library functions that anyone can wrap in a CLI.

| kent tool | Library equivalent |
|---|---|
| `bigWigInfo` | `BigWigReader::info()`, `chroms()`, `zoom_levels()` |
| `bedGraphToBigWig` | `convert::bedgraph_to_bigwig` |
| `wigToBigWig` | `convert::wig_to_bigwig` |
| `bigWigToBedGraph` | `convert::bigwig_to_bedgraph` |
| `bigWigToWig` | `convert::bigwig_to_wig` |
| `bigWigSummary` | `BigWigReader::summaries()` |

`bigWigAverageOverBed` and `bigWigMerge` have no equivalent yet.

## Compatibility

The [interop matrix](interop/README.md) checks, in CI, that files this crate writes are read correctly by:

- kent's utilities (`bigWigInfo`, `bigWigToBedGraph`, `bigWigToWig`, `bigWigSummary`);
- libBigWig, through pyBigWig;
- bigtools;
- `@gmod/bbi`, the reader behind JBrowse 2;
- IGV desktop (2.19.8);
- igv.js (3.8.9).

It also checks the other direction: the reader must reproduce exactly what kent's `bedGraphToBigWig` and `wigToBigWig`, pyBigWig and bigtools write. Beyond the matrix, the reader has been run over a corpus of public bigWig files, including files attached to other projects' bug reports.

The writer avoids problems that specific readers have had with files from other writers:

- The chromosome tree comes before the data, as in kent's files. IGV 2.17.0–2.18.2, igv.js 3.0.5 and earlier, and `@gmod/bbi` before 7.0.3 failed to open files with the tree after the data ([bigtools #93](https://github.com/jackh726/bigtools/issues/93)).
- There is always at least one zoom level, unless `ZoomPolicy::None` asks for none: IGV 2.17.0–2.19.4 throw on a file with no zoom levels. No zoom reduction exceeds 2³¹ − 1, which IGV reads as negative ([bigtools #85](https://github.com/jackh726/bigtools/issues/85)).
- Zoom records are written only for bins with data, so a gap longer than a bin never gets a record claiming coverage, which IGV draws as signal ([bigtools PR #103](https://github.com/jackh726/bigtools/pull/103)). No record has zero covered bases.
- Up to 32,767 chromosomes, the chromosome tree is a single node of byte-sorted names, so a lookup never descends through internal nodes, where igv.js compares names with `localeCompare` rather than in byte order. Beyond that the tree has internal nodes, and igv.js can fail to find names that sort differently by bytes and by locale: of 20,000 `Scaffold_N` and 20,000 `contig_N` names, it found none of the `Scaffold_N`. kent's files have internal nodes beyond 256 chromosomes.
- The decompression buffer size in the header is the exact largest block, since `@gmod/bbi` decompresses into a buffer of exactly that size, and never 0 for a compressed file, which readers would take for an uncompressed one. Uncompressed files (`Compression::None`) are for testing and for pipelines that don't feed igv.js, which misreads them when a query spans more than one block.
- Creating a writer logs a warning for a chromosome longer than 2³¹ − 1 bases, since IGV, igv.js and JBrowse read coordinates as signed 32-bit integers.

Big-endian files can be written through `WriterOptions::byte_order`, for interoperability testing. Little-endian is the default, and the one to use for anything shared, because not every reader handles big-endian:

| Reads big-endian bigWig | |
|---|---|
| this crate | yes |
| kent | yes |
| bigtools | yes |
| pyBigWig | no |
| `@gmod/bbi` | no |
| IGV desktop | no |
| igv.js | no |

## Defaults and differences from kent's tools

The goal is files every reader handles, not byte-identical output. Three defaults differ from kent's, each measured on real data:

| Setting | Default | kent | Why |
|---|---|---|---|
| Deflate level | 1 | 6 | On dense per-base depth, about 11% larger for about 5× less CPU. `Compression::Deflate { level: 6 }` gives kent-sized output. |
| Items per section | 8,192 | 1,024 | On dense per-base depth, about 11% smaller and a third less compression CPU, while a single-base query still inflates only a 32 KiB fixedStep section. |
| First zoom level | 256 bp, or coarser for sparse data | about 35 bp for per-base data (`wigToBigWig`) | kent's 35 bp level alone is over half the size of the data; the whole ladder from 256 bp adds about 11%. |

The cost of the coarser first zoom level is that genome browsers read raw data over wider windows before switching to zoom data: up to roughly 128–512 kb rather than 18–70 kb, depending on the browser. A 500 kb window of per-base data is about 300 KB to fetch and about 1 ms to inflate. For a finer first level, e.g. 32 bp, use `ZoomPolicy::Auto { min_first_reduction: 32 }`. Each later level is 4× the previous, as in kent's files.

Other deliberate differences:

- Zoom records sit on a fixed grid (bin `k` of a level with reduction `r` is `[k·r, (k+1)·r)`) and cover the tight extent of the bases they hold, where kent's bins start wherever the data does.
- Every chromosome in the list passed to the writer goes in the chromosome tree, including those with no data. kent writes only chromosomes that have data.
- Text parsing is strict about numbers: a malformed line or a NaN value is an error naming the line. Zero-length bedGraph intervals are skipped and counted in the conversion `Report`, where kent writes them.
- With `WriterOptions::clip`, an item past a chromosome's end is truncated, so the bases that exist keep their data; `wigToBigWig -clip` drops it whole.
- Exports list chromosomes in `chromId` order, not name order, and print values in the shortest form that reads back to the same `f32`, not six significant digits. `bigwig_to_wig` clips items to the selection, writes one header per run rather than per stored section, and writes no empty headers.
- Zoom-based summaries keep the pro-rated number of covered bases, which may be fractional. kent rounds it up and scales the sums to match, so means agree but coverage and sums can be slightly lower than kent's.

## Compression backends

bigWig blocks are zlib streams. Exactly one backend is used:

| Feature | Backend | Notes |
|---|---|---|
| `libdeflate` (default) | [libdeflate](https://github.com/ebiggers/libdeflate) via `libdeflater` | Compiles bundled C with `cc`. |
| `zlib-rs` | [zlib-rs](https://github.com/trifectatechfoundation/zlib-rs) via `flate2` | Pure Rust; no C toolchain needed. |

At least one must be enabled. With both, libdeflate is used, so enabling `zlib-rs` elsewhere in a dependency graph never changes the output. To use the pure-Rust backend:

```toml
bigwig = { version = "0.1", default-features = false, features = ["zlib-rs"] }
```

libdeflate is the default because it measured better on real data. Nine datasets (ENCODE ChIP-seq, ATAC-seq, RNA-seq and WGBS tracks, two phastCons tracks, and per-base and binned depth) were written at several levels of each backend and three section sizes on aarch64 (Apple M4 Max). In 22 of the 27 dataset × section-size groups, some libdeflate level was at least as small and at least as cheap to write as every zlib-rs level, and reading was about 20% faster. x86_64 was not measured. At the default level 1, zlib-rs wrote files 1.3–1.85× the size of libdeflate's, so a zlib-rs build may want a higher level.

## Contributing

[CONTRIBUTING.md](CONTRIBUTING.md) covers building, the checks CI runs, the test suites (unit, property, fuzz and interop), benchmarks and conventions.

## Licence

MIT; see [LICENSE](LICENSE).
