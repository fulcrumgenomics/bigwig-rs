# bigWig interop matrix

Cross-implementation checks that the files `bigwig` writes are read correctly by kent's tools, pyBigWig, bigtools, `@gmod/bbi`, IGV desktop and igv.js, and that `bigwig` reads exactly what kent's tools, pyBigWig and bigtools write. This directory is a separate, unpublished crate with its own workspace, and it is not part of the published package.

## Pieces

- `src/`: the `bigwig-interop` binary.
  - `generate <seed> <outdir>` writes every fixture deterministically, using an inline splitmix64 generator.
  - `dump <info|chroms|items|values|summaries|zoom> <bw> ...` prints what our reader sees, as TSV. f32 values are widened to f64 and printed losslessly.
  - `dump <bedgraph|wig> <bw> [--chrom --start --end]` prints what `convert::bigwig_to_bedgraph` (with `merge_equal`) and `convert::bigwig_to_wig` write.
  - `scan` and `sample` serve the public-file corpus (see below).
- `compare.py`: runs every tool, compares what it reports against `dump`, and prints one PASS/FAIL line per check, followed by a summary. It runs under the Python that has pyBigWig.
- `bbi_dump.mjs`: reads a bigWig with `@gmod/bbi`, taking a batch of queries per file. It runs under node or bun.
- `IgvDump.java`: reads a bigWig with IGV desktop's own classes, headless, opening it as IGV opens a bigWig track, and takes a batch of requests per file. It needs IGV 2.19.x's jars and Java 21. A few IGV internals are reached reflectively; they are listed together in the file, for whoever moves it to a newer IGV.
- `igvjs_dump.mjs`: reads a bigWig with igv.js's own reader, from an igv.js source tree, taking a batch of queries per file. It runs under node or bun.
- `run.sh`: the driver. It builds the binary, fetches igv.js on its first run, generates the fixtures into `target/interop-run`, compiles `IgvDump.java`, and runs `compare.py`. It exits non-zero if any required check fails.
- `corpus.sh` and `corpus.py`: the one-time public-file validation, run by hand, and `corpus.tsv`, its manifest.
- `environment.yml` (conda) and `package.json` / `package-lock.json` (npm): the other implementations, pinned, apart from igv.js, which `run.sh` pins (see [igv.js](#igvjs)).

## Fixtures

Our bigWigs (`ours/`, listed in `manifest.tsv`):

| Fixture | What it covers |
|---|---|
| `sections` | fixedStep (including span < step), variableStep and bedGraph sections mixed on each chromosome, with abutting runs, abutting equal values and some full-precision f32 values; default options (Auto zoom) |
| `zoom_explicit` | the same data with Explicit zooms 50, 130, 1000, 4000 and 20000 (130 is not a multiple of 50), and 1,000 items per section |
| `uncompressed` | the same data, stored uncompressed |
| `big_endian` | the same data, written big-endian |
| `riker_per_base` | per-base depth as step-1 fixedStep pushed in chunks; BAM-order chromosomes, Listed IDs, and a zero-data chromosome |
| `riker_binned` | 100-base fixedStep bins over whole chromosomes |
| `riker_quantised` | abutting bedGraph runs of quantised values |
| `contigs_3000` | 3,000 contigs in BAM-header (non-lexicographic) order, mostly empty, with Listed IDs |
| `contigs_40000` | 40,000 contigs, which needs a multi-level B+ tree; data pushed in an order unrelated to the list |
| `locale_names` | contig names that `localeCompare` orders differently from bytes (`chrZ`/`chra`, `Chr3`/`chr1`, `chr_alt`/`chrEBV`, ...) |
| `cpg_sparse` | sparse one-base variableStep methylation fractions |
| `empty` | chromosomes and no data |

Text inputs for the other tools' writers (`text/`) are `input.bedGraph`, `fixed.wig`, `variable.wig` and `mixed.wig`. They share a `chrom.sizes` in BAM order, and each has an `.items.tsv` taken from the generated data, not from a parser.

## Checks

**Our files, read by them.** Every fixture is checked by:

- kent:
  - `bigWigInfo`: chromosomes with their IDs, zoom levels, byte order, compression, and the total summary;
  - `bigWigToBedGraph`: `convert::bigwig_to_bedgraph` with `merge_equal`, over the whole file and over a region of each sampled chromosome;
  - `bigWigToWig`: `convert::bigwig_to_wig` over the same selections. Every item must be in a block of the same kind (fixedStep with the same step), and our output must hold to its documented layout: no header without a value, runs merged under one header, and each `#bedGraph section` comment giving the extent of its lines;
  - `bigWigSummary`: every `-type`, over a raw-path query, a grid-aligned query for each zoom level, and an unaligned zoom query.
- pyBigWig: `chroms`, `header`, `intervals`, `values`, `stats(exact=True)`, and `stats(exact=False)` on the aligned queries.
- bigtools: `bigwiginfo` and `bigwigtobedgraph`.
- `@gmod/bbi`:
  - the header: chromosomes, zoom levels and the total summary;
  - `getFeatures` on raw data;
  - `getFeatures` at a scale that picks each zoom level, against our `zoom_records`.
- IGV desktop, through `IgvDump.java`:
  - the header: version, zoom levels, decompression buffer size and the total summary; every chromosome's name and ID, and IGV's lookup of each name;
  - `getRawData`, what IGV draws zoomed in, over the regions the other readers' raw queries cover;
  - every zoom level's records over each window, decoded as IGV decodes the level it picks, under its min, mean and max window functions, as are the records below;
  - what IGV draws at each of its zooms (chromosome length / (2^zoom × 700) bases per pixel), down to the first zoom where it uses no zoom level: the level must be the one IGV's rule picks, and the records ours;
  - the whole-genome view IGV first shows for a genome with several chromosomes: one zoom level across every long chromosome, in genome coordinates.
- igv.js, through `igvjs_dump.mjs`:
  - the header: type, byte order, version, decompression buffer size, chromosome count, zoom levels and the total summary;
  - the lookup of every chromosome name, and of a name the file does not hold;
  - `getFeatures` on raw data;
  - `getFeatures` at a resolution that picks each zoom level, under the mean, min and max window functions, over each window and over each whole sampled chromosome, where a query reads several blocks.

**Their files, read by us.** The text inputs are converted by:

- kent `bedGraphToBigWig`;
- kent `wigToBigWig`, from fixedStep, variableStep and mixed wig;
- pyBigWig `addEntries` in all three forms, plus a mixed file;
- bigtools `bedgraphtobigwig`.

On every converted file, our reader must reproduce the input items exactly. Every zoom level must also read back with the raw data's min and max.

**Tolerances.**

- Raw values, intervals, bases covered, min and max must be exact.
- Raw means and standard deviations must agree to a relative 1e-6.
- Zoom statistics must agree to f32 precision.

Where a tool prints fewer digits than that, the check also allows half a unit in the last printed digit: kent's `%g` has six significant digits, and its `%f` has six decimals. Every check also asserts that it compared something. A query on a chromosome with data must find data, and a query on an empty chromosome must find none.

## Known differences encoded in `compare.py`

Each is commented in the code where it is applied.

- **Printing.** kent prints values with `%g` and summaries with `%f`, so values are compared as printed text.
- **Chromosome order.** kent and bigtools list chromosomes in name order, while we list them in chromId order, so comparisons sort both.
- **kent's bigWigToWig layout.** kent writes one header per stored section, including empty ones, and does not clip items to a region; `convert::bigwig_to_wig` merges runs, writes no empty header, and clips. Items and their block kinds are compared, with kent's items clipped to the region, plus the `#bedGraph section` comments over the whole file.
- **kent's whole-number base counts.** kent's `normalizeCount` rounds each bin's pro-rated base count up to a whole number and rescales the sums.
  - This changes zoom coverage and std where records straddle bins.
  - It changes raw coverage by one base where the count's floating-point error lands just above a whole number. This happened locally: `1603.0000000000002` was counted as 1604.
  - The check repeats kent's arithmetic on our items and records. It also requires our count to equal kent's unrounded count.
- **kent's coverage formula.** kent's coverage divides by the mean bin width, not by the bin's own width.
- **kent's pixel replication.** kent replicates bases across pixels when there are more bins than bases. Queries avoid that case.
- **Zoom-level choice.** kent picks the largest reduction strictly below half the bin width, and libBigWig accepts one equal to it. Aligned queries use bin widths where both rules pick the same level, and our `BestZoom` must pick kent's level.
- **pyBigWig's header.** `header()` truncates floats to integers.
- **pyBigWig's zoom sum.** The zoom-based `sum` is not compared, because libBigWig truncates overlap fractions ([pyBigWig #146](https://github.com/deeptools/pyBigWig/issues/146)).
- **Std near zero.** A std computed from sums is accurate only to about 1e-6 of the RMS value once the values are nearly constant, so stds also match when their variances agree to within 1e-12 of the mean square.
  - Below that floor, kent prints `n/a`: it takes an unguarded sqrt of a variance rounded below zero.
  - libBigWig's zoom std is NaN for at most one covered base, or for a variance below −1e-8.
  - On Linux, the square root of that negative variance also sets `errno`, which libBigWig treats as an error, so the whole `stats()` call fails. The check accepts that only for a query with a bin whose variance is within the floor of zero.
  - libBigWig's exact std is two-pass, so it is exactly 0 for constant values.
- **Std of an empty file.** kent and bigtools print NaN, and our `Summary::std` is 0.
- **bigtools' input order.** bigtools `bedgraphtobigwig` needs `--sorted start` for BAM-order chromosomes; its default requires byte order.
- **pyBigWig's fixedStep end.** pyBigWig takes a fixedStep block to end at `start + count·step` and rejects a following block that starts sooner, so the text fixtures leave that room. Its files' fixedStep sections do carry the overshooting `end` ([pyBigWig #166](https://github.com/deeptools/pyBigWig/issues/166)). With pyBigWig 0.3.26 it was 6 steps plus `step − span` past the last item, and our reader ignores it as designed.
- **Zoom data in their files.** Only min and max are compared: kent's `wigToBigWig` over-counts coverage (`bbiAddToSummary` in `lib/bbiWrite.c`), libBigWig zeroes some sums, and bigtools writes phantom records over gaps ([bigtools PR #103](https://github.com/jackh726/bigtools/pull/103)).
- **Items ending at the query start.** IGV and igv.js return raw items and zoom records whose end equals the query start: their decoders skip only those ending before it. The checks drop them and say how many.
- **IGV's zoom records on chromId 0.** `BBFile.decodeZoomData` filters records by region only `if (chrIdx > 0)`, so a zoom query on the chromosome with ID 0 returns every record of the blocks it reads, other chromosomes' included, since zoom blocks span chromosomes (in kent's files as in ours). The checks drop the records outside the query on chromId 0, and say how many; anywhere else such a record fails the check.
- **IGV's zoom mean.** IGV computes a record's mean in single precision, as `sumData / validCount`, so the expected mean is that f32 quotient.
- **igv.js's chromosome lookups.** Decoding features caches the name of every chromosome in the tree leaves it reads, and a cached name resolves even where igv.js's tree search would miss it, so `igvjs_dump.mjs` looks every name up before any feature query. igv.js returns no features both for a name it cannot find and for a region without data, so the empty chromosome's check also requires igv.js to find its name.

## Big-endian results

These were recorded locally, on macOS arm64. kent and our reader are required to read the file. The other results are recorded, not required.

| Reader | Reads our big-endian file |
|---|---|
| `bigwig` (this crate) | yes |
| kent utilities (built from commit `233187ca`) | yes (`isSwapped: 1`) |
| bigtools 0.5.8 | yes (`bigwiginfo`, `bigwigtobedgraph`) |
| pyBigWig 0.3.26 | no: `pyBigWig.open` fails (libBigWig is native-endian only) |
| `@gmod/bbi` 11.3.2 | no: "not a BigWig/BigBed file" (little-endian only) |
| IGV 2.19.8 | no: "Bad magic number". `BBFile.loadBinaryBuffer` ignores the byte order it is given, so the big-endian retry reads the header little-endian again. |
| igv.js 3.8.9 | no: the big-endian magic number is compared with a mistyped constant, `0x26FC8F66` ([igv.js #786](https://github.com/igvteam/igv.js/issues/786)), and the decoders ignore the byte order. |

## Uncompressed results

igv.js is not required to read the `uncompressed` fixture, and its results there are recorded. When a query reads more than one block, igv.js decodes each from the start of the buffer holding them all, which is right only for compressed blocks, so it returns duplicated items and misses others. With igv.js 3.8.9, the raw queries over the two whole chromosomes stored in several blocks failed, and every other check passed. Every other reader is required to read the fixture, and does.

## Running locally

The matrix needs the six other implementations. The quickest way to get most of them is the conda environment CI uses, created with micromamba, mamba or conda (`conda env create -f interop/environment.yml`), which has kent's tools, pyBigWig, bigtools, IGV's jars and Java 21, plus `@gmod/bbi` from npm. igv.js comes from GitHub: `run.sh` fetches it on its first run, which needs network access, and reuses it after that.

```bash
micromamba create -f interop/environment.yml
micromamba activate bigwig-interop
(cd interop && npm ci)
interop/run.sh
```

Every tool defaults to its name on PATH, and IGV's jars to the environment's `lib/igv`, so inside the environment nothing else is needed. bioconda has every tool for x86_64 Linux, where CI runs; on other platforms some may be missing. Tools built or installed elsewhere can be named with environment variables:

```bash
KENT_DIR=/path/to/kent/bin \
BIGTOOLS_DIR=$HOME/.cargo/bin \
PYTHON=/path/to/venv/bin/python \
JS_RUNTIME=bun \
JAVA=/path/to/jdk-21/bin/java \
IGV_LIB=/path/to/IGV_2.19.8/lib \
interop/run.sh
```

- `KENT_DIR`: the directory holding kent's `bigWigInfo`, `bigWigToBedGraph`, `bigWigToWig`, `bigWigSummary`, `bedGraphToBigWig` and `wigToBigWig`.
- `BIGTOOLS_DIR`: the directory holding the bigtools binaries, e.g. from `cargo install bigtools`.
- `PYTHON`: a Python that can `import pyBigWig`.
- `JS_RUNTIME`: `node` (the default) or `bun`. Where only bun is installed, `bunx --bun npm ci` in this directory installs `@gmod/bbi`.
- `JAVA`: a Java 21 `java` from a JDK, with `javac` in the same directory.
- `IGV_LIB`: the directory holding IGV 2.19.x's jars: by default `$CONDA_PREFIX/lib/igv`, where bioconda's `igv` package puts them, or the `lib` directory of IGV's own release ([IGV_2.19.8.zip](https://data.broadinstitute.org/igv/projects/downloads/2.19/IGV_2.19.8.zip)), which gave the same results.
- `IGVJS_DIR`: an igv.js source tree to use instead of the pinned release, such as a patched copy, laid out as `igvjs_dump.mjs` describes.

`run.sh` also accepts `SEED`, `OUT_DIR`, `INTEROP_BIN` and per-tool overrides, which are listed at its top. Arguments restrict the run to the named fixtures; `theirs` selects the other tools' writers. For example, `interop/run.sh sections theirs` runs those two.

### igv.js

igv.js is not on conda, and its npm package can't be used here: it ships only a bundle, which needs a browser's `document`. Instead, `run.sh` fetches igv.js's source at the commit tagged v3.8.9, and the source of igv-utils at the commit igv.js's `package.json` names, from GitHub. It keeps only igv.js's `js/` and `package.json`, and igv-utils' `src/` at `node_modules/igv-utils`, where igv.js imports it from, and installs nothing. The kept files must match a sha256 of their paths and contents, rather than of the tarballs, whose bytes GitHub does not promise to keep stable. Only then are they moved into `target/igvjs-3.8.9`, which later runs reuse; delete it to fetch again. The commits, version and hash are at the top of `run.sh`.

The file is handed to igv.js as a `File`, as a browser's file picker does, so no web server is needed.

## CI

`.github/workflows/interop.yml` runs `run.sh` on Linux x86_64 for every pull request, and on demand. It sets up the conda environment with `setup-micromamba`, runs `npm ci`, builds the binary, and runs the matrix, which fetches igv.js from GitHub.

## Public-file corpus

`corpus.sh` validates our reader and writer against public bigWig files. It is run by hand, never in CI, because third-party downloads are too brittle to gate on. The files are never committed: `corpus.tsv` lists each one with what is needed to fetch it again and to know why it is there.

### The manifest

`corpus.tsv` has one row per downloadable file, either cited by an issue in another implementation's tracker or chosen directly (curated) to cover a use case or a writer. Files that are unreachable, too large or behind a file-sharing service are listed too, so it is a complete inventory. The scripts read it by column name:

| Column | What it holds |
|---|---|
| `name` | A short, unique name. The file is fetched to `<download-dir>/<name>/`. |
| `url` | The download URL. |
| `sha256`, `size` | The file's digest and size in bytes, blank if it has never been fetched. A download must match `sha256`, and `size` decides whether it fits the download caps. |
| `archive_member` | For a zip or gzip download, the bigWig inside. Archives are unpacked and every bigWig in them is checked, so this is for reference. |
| `provider` | Who serves the file: ENCODE, UCSC, Ensembl, modENCODE, a repo fixture, an issue attachment, a lab server, an S3 bucket, Dropbox, Google Drive or ownCloud. |
| `source_issue` | The issue that cites the file, as `owner/repo#N`; blank for curated files. |
| `source_url` | Where the file was found: the issue, or for a curated file its ENCODE file page, UCSC download directory or repository page. |
| `writer` | The software that wrote the file, as claimed or presumed; blank if unknown. |
| `format` | What the file's header shows: bigWig version, byte order, compression (with the uncompress buffer size), whether the chrom tree comes before or after the data, and libBigWig's 10 reserved zoom slots or a missing total summary. Blank for files never fetched. |
| `expectation` | What the validation should find. A row whose expectation starts with `skip:` is never fetched: it is behind a file-sharing service, or is bedGraph text. |
| `notes` | What the file exercises, and whether it could be fetched on the last run. |

Results are not committed. A file that fails becomes a regression test with generated data, in `tests/issues.rs` or the unit tests, so the permanent suite never depends on the downloads.

### Fetching the corpus

```bash
interop/corpus.sh --fetch-only interop/corpus.tsv interop/target/corpus interop/target/corpus-filled.tsv
```

This fetches every row not marked `skip:` to `<download-dir>/<name>/<the URL's file name>` and checks it against `sha256`. It needs only curl and Python: no kent tools and no Rust build.

- A file already at that path is kept. A file under `CORPUS_LOCAL_DIRS` (colon-separated directories) whose name ends with the URL's file name and whose sha256 matches is hard-linked instead of downloaded, or copied across file systems. An earlier download directory, or copies under other names such as `07-test.bw` or `pyBigWig_test.bw` for `test.bw`, can be reused this way.
- Each download is capped at `CORPUS_MAX_FILE_BYTES` (1 GiB by default), and all downloads together at `CORPUS_BUDGET_BYTES` (4 GB by default). A row whose `size` is over a cap is skipped. For a full rebuild, raise either to a number of bytes or lift it with `none`: the 14 rows over the default cap come to about 46 GB, the largest being 9.6 GB.
- The filled manifest is the manifest with `sha256` and `size` filled in for everything fetched. Every other column is carried through. A new row therefore needs only a `name`, `url` and `expectation`, and fetching it records its digest and size.
- The exit status is 1 only if a file does not match its `sha256`. Unreachable and skipped rows are reported, not failures.

### Running the validation

```bash
KENT_DIR=/path/to/kent/bin \
interop/corpus.sh interop/corpus.tsv interop/target/corpus interop/target/corpus-filled.tsv interop/target/corpus.md
```

This fetches the corpus as above, builds the `bigwig-interop` binary unless `INTEROP_BIN` is set, and needs kent's `bigWigToBedGraph`, `bigWigToWig` and `bigWigSummary`.

- Zip, gzip and tar downloads are unpacked.
- For each bigWig, `scan` opens it, reads every item, recomputes the total summary, and rewrites the items with our writer, which must read back identically.
- On four windows around randomly sampled items, kent `bigWigToBedGraph` and `bigWigToWig` are compared with `convert`. Our `values` and raw summaries are compared with kent's bedGraph, and kent `bigWigSummary` with our summaries on a raw-path and a zoom-path query, using the tolerances and known differences above.
- A panic in our binary is reported as such.
- The results file is a Markdown table of every file and check, with the claimed writer from the `writer` column, plus the entries that were not fetched and why. `--annotations <tsv>` (columns `url` and `writer`) overrides the claimed writers.
