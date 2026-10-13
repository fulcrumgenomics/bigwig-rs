//! [`BigWigWriter`] and [`WriterOptions`]: streaming bigWig output with validation, zoom levels
//! and a bounded-memory zoom spill.
//!
//! The file is laid out as kent's `bedGraphToBigWig` lays it out: the fixed header, reserved
//! zoom-header slots, the total summary, the chromosome B+ tree, the data sections and their
//! R-tree, then each zoom level's blocks and R-tree, and the trailing magic. `create` writes
//! everything up to the data, with a zeroed header (magic included) that `finish` patches, so a
//! file that is never finished fails every reader's magic check.
//!
//! A fixedStep section's values reach the total summary and the zoom accumulators when the
//! section is flushed; bedGraph and variableStep items reach them as they are pushed. Either way
//! the order of every f64 sum depends only on the data and `items_per_section`, never on how the
//! caller sliced its pushes. Zoom records are packed into blocks and compressed at once; the
//! compressed blocks stay in memory up to a limit, beyond which they spill to a temporary file
//! that `finish` copies back.

use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Read, Seek, SeekFrom, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::bbi::{
    BIGWIG_MAGIC, BIGWIG_VERSION, ByteOrder, ChromInfo, FIXED_HEADER_SIZE, FixedHeader,
    TOTAL_SUMMARY_SIZE, ZOOM_HEADER_SIZE, ZoomLevel, encode_total_summary,
};
use crate::bptree::{self, ChromEntry, EncodedChromTree};
use crate::compress::{BlockCompressor, Compression, CompressionJob, MAX_INFLATED_SIZE};
use crate::error::{Error, Result};
use crate::rtree::{self, LeafEntry};
use crate::section::{SectionBuilder, SectionKind};
use crate::zoom::{
    AUTO_MAX_LEVELS, PitchSampler, RawSummaryAccumulator, Summary, ZoomBuilder, ZoomPolicy,
    ZoomRecord, kept_level_count,
};

/// Default maximum items per data section.
const DEFAULT_ITEMS_PER_SECTION: u16 = 8_192;
/// Default zoom records per zoom block, as kent writes.
const DEFAULT_ZOOM_ITEMS_PER_BLOCK: u32 = 1_024;
/// Default fan-out of the R-tree indexes, as kent writes.
const DEFAULT_INDEX_BLOCK_SIZE: u32 = 256;
/// Default bytes of compressed zoom blocks held in memory before spilling to a temporary file.
const DEFAULT_ZOOM_MEMORY_LIMIT: usize = 64 << 20;
/// Capacity of the buffer between the writer and its sink.
const SINK_BUFFER_SIZE: usize = 256 << 10;
/// Longest chromosome that readers using signed 32-bit coordinates (IGV, igv.js, JBrowse) handle.
const LONGEST_SIGNED_CHROM_LENGTH: u32 = i32::MAX.unsigned_abs();
/// Attempts at finding an unused spill-file name before giving up.
const SPILL_NAME_ATTEMPTS: u32 = 1_000;
/// Stands in for the output file name in spill-file names when the sink has no path.
const UNNAMED_OUTPUT: &str = "bigwig";

/// Numbers the spill files of every writer in the process, so their names never collide.
static NEXT_SPILL_FILE: AtomicU64 = AtomicU64::new(0);

/// Options for creating a [`BigWigWriter`], in the style of [`std::fs::OpenOptions`].
///
/// Start from [`WriterOptions::new`] (every default), change what you need with the chainable
/// setters, then call [`create`](Self::create) for a file path or
/// [`create_with`](Self::create_with) for any seekable sink. None of the options is needed for a
/// correct file. The value is `Clone`, so one can be kept in a configuration and reused for several
/// outputs. Every option is validated by `create`/`create_with`, which return
/// [`Error::InvalidOption`] naming the option.
///
/// | Option | Default |
/// |---|---|
/// | [`items_per_section`](Self::items_per_section) | 8,192 |
/// | [`zoom`](Self::zoom) | `ZoomPolicy::Auto { min_first_reduction: 256 }` |
/// | [`zoom_items_per_block`](Self::zoom_items_per_block) | 1,024 |
/// | [`index_block_size`](Self::index_block_size) | 256 |
/// | [`compression`](Self::compression) | `Compression::Deflate { level: 1 }` |
/// | [`zoom_memory_limit`](Self::zoom_memory_limit) | 64 MiB |
/// | [`temp_dir`](Self::temp_dir) | the output file's directory, or `std::env::temp_dir()` |
/// | [`chrom_id_order`](Self::chrom_id_order) | [`ChromIdOrder::Arrival`] |
/// | [`clip`](Self::clip) | false |
/// | [`byte_order`](Self::byte_order) | [`ByteOrder::Little`] |
///
/// # Example
///
/// ```
/// use std::io::Cursor;
///
/// use bigwig::{ChromInfo, Compression, WriterOptions};
///
/// let chroms = [ChromInfo::new("chr1", 1_000), ChromInfo::new("chr2", 500)];
/// let options = WriterOptions::new()
///     .items_per_section(1_024)
///     .compression(Compression::Deflate { level: 6 });
///
/// let mut writer = options.create_with(Cursor::new(Vec::new()), &chroms)?;
/// writer.push_fixed_step("chr1", 100, 1, 1, &[1.0, 2.0, 3.0])?;
/// writer.push_intervals("chr2", &[(0, 50, 0.5), (60, 70, 1.5)])?;
/// let bytes = writer.finish()?.into_inner();
///
/// // A finished file starts and ends with the bigWig magic number.
/// assert_eq!(bytes[..4], 0x888F_FC26_u32.to_le_bytes());
/// assert_eq!(bytes[bytes.len() - 4..], 0x888F_FC26_u32.to_le_bytes());
/// # Ok::<(), bigwig::Error>(())
/// ```
#[derive(Clone, Debug)]
pub struct WriterOptions {
    items_per_section: u16,
    zoom: ZoomPolicy,
    zoom_items_per_block: u32,
    index_block_size: u32,
    compression: Compression,
    zoom_memory_limit: usize,
    temp_dir: Option<PathBuf>,
    chrom_id_order: ChromIdOrder,
    clip: bool,
    byte_order: ByteOrder,
}

impl WriterOptions {
    /// Options with every default (see the table above).
    pub fn new() -> Self {
        Self {
            items_per_section: DEFAULT_ITEMS_PER_SECTION,
            zoom: ZoomPolicy::default(),
            zoom_items_per_block: DEFAULT_ZOOM_ITEMS_PER_BLOCK,
            index_block_size: DEFAULT_INDEX_BLOCK_SIZE,
            compression: Compression::default(),
            zoom_memory_limit: DEFAULT_ZOOM_MEMORY_LIMIT,
            temp_dir: None,
            chrom_id_order: ChromIdOrder::default(),
            clip: false,
            byte_order: ByteOrder::default(),
        }
    }

    /// Maximum items per data section: fixedStep values, variableStep items or bedGraph
    /// intervals. Must be at least 1. Larger sections compress better; smaller ones make
    /// single-base queries decompress less.
    #[must_use]
    pub fn items_per_section(mut self, items: u16) -> Self {
        self.items_per_section = items;
        self
    }

    /// How zoom levels are chosen (see [`ZoomPolicy`]).
    #[must_use]
    pub fn zoom(mut self, policy: ZoomPolicy) -> Self {
        self.zoom = policy;
        self
    }

    /// Zoom records per compressed zoom block: 1 to 2,097,152 (a block of at most 64 MiB).
    #[must_use]
    pub fn zoom_items_per_block(mut self, records: u32) -> Self {
        self.zoom_items_per_block = records;
        self
    }

    /// Fan-out of the R-tree indexes over the data and every zoom level: 2 to 65,535. The
    /// chromosome B+ tree is not affected; it is a single node for up to 32,767 chromosomes.
    #[must_use]
    pub fn index_block_size(mut self, children: u32) -> Self {
        self.index_block_size = children;
        self
    }

    /// Block compression. The level must be supported by the compiled backend
    /// ([`Compression::MAX_LEVEL`]).
    #[must_use]
    pub fn compression(mut self, compression: Compression) -> Self {
        self.compression = compression;
        self
    }

    /// Bytes of compressed zoom blocks kept in memory before they spill to a temporary file. Any
    /// value is valid; 0 spills every block as soon as it is compressed.
    #[must_use]
    pub fn zoom_memory_limit(mut self, bytes: usize) -> Self {
        self.zoom_memory_limit = bytes;
        self
    }

    /// Directory for the zoom spill file. Defaults to the output file's directory, or to
    /// [`std::env::temp_dir`] when writing to a sink with no path. It is used only if spilling
    /// starts, and must then exist and be writable.
    #[must_use]
    pub fn temp_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.temp_dir = Some(dir.into());
        self
    }

    /// How `chromId`s are assigned: by first data arrival (the default) or by position in the
    /// chromosome list passed to `create` (see [`ChromIdOrder`]).
    #[must_use]
    pub fn chrom_id_order(mut self, order: ChromIdOrder) -> Self {
        self.chrom_id_order = order;
        self
    }

    /// With `true`, items extending past their chromosome's end are truncated at the end, and
    /// items starting at or past it are dropped, instead of failing with
    /// [`Error::OutOfBounds`]. A truncated fixedStep or variableStep item is written in a
    /// one-item section of its own, because its span differs. Clipping changes only what is
    /// written: items are still checked for order against the earlier items' full extents,
    /// dropped items included, so splitting a push differently never changes what is accepted.
    #[must_use]
    pub fn clip(mut self, clip: bool) -> Self {
        self.clip = clip;
        self
    }

    /// Byte order of the written file. Leave it little-endian unless testing interoperability:
    /// several readers (libBigWig on little-endian hosts, `@gmod/bbi`, IGV) can't read big-endian
    /// files.
    #[must_use]
    pub fn byte_order(mut self, order: ByteOrder) -> Self {
        self.byte_order = order;
        self
    }

    /// Validates the options and chromosome list, then creates (or truncates) the file at `path`
    /// and writes the placeholder header, the reserved zoom-header slots, the total-summary
    /// placeholder and the chromosome tree.
    ///
    /// `chroms` may be in any order. Nothing is created when validation fails. Fails with
    /// [`Error::InvalidOption`], [`Error::InvalidChromList`] or [`Error::Io`]; an error creating
    /// the file names `path` and keeps the underlying error's kind.
    pub fn create(
        &self,
        path: impl AsRef<Path>,
        chroms: &[ChromInfo],
    ) -> Result<BigWigWriter<File>> {
        let path = path.as_ref();
        self.validate()?;
        let chrom_list = ChromList::new(chroms, self.chrom_id_order)?;
        let file = File::create(path).map_err(|error| {
            Error::Io(io::Error::new(
                error.kind(),
                format!("cannot create {}: {error}", path.display()),
            ))
        })?;
        let spill_dir = self.temp_dir.clone().unwrap_or_else(|| match path.parent() {
            Some(parent) if !parent.as_os_str().is_empty() => parent.to_path_buf(),
            _ => PathBuf::from("."),
        });
        let output_name = path
            .file_name()
            .map_or_else(|| UNNAMED_OUTPUT.to_owned(), |name| name.to_string_lossy().into_owned());
        let spill_location = SpillLocation { dir: spill_dir, output_name };
        BigWigWriter::new(file, chrom_list, self, spill_location, Some(path.to_path_buf()))
    }

    /// As [`create`](Self::create), writing to any seekable sink, which must be empty. The file
    /// is written from offset 0 of the sink, because bigWig offsets are absolute, and the sink is
    /// never truncated, so a longer sink would keep stale bytes after the file's trailing magic,
    /// which readers such as kent's reject. Fails with [`Error::Io`] of kind
    /// [`InvalidInput`](io::ErrorKind::InvalidInput) for a sink that isn't empty, after validating
    /// the options and chromosome list.
    pub fn create_with<W: Write + Seek>(
        &self,
        sink: W,
        chroms: &[ChromInfo],
    ) -> Result<BigWigWriter<W>> {
        self.validate()?;
        let chrom_list = ChromList::new(chroms, self.chrom_id_order)?;
        let spill_dir = self.temp_dir.clone().unwrap_or_else(std::env::temp_dir);
        let spill_location =
            SpillLocation { dir: spill_dir, output_name: UNNAMED_OUTPUT.to_owned() };
        BigWigWriter::new(sink, chrom_list, self, spill_location, None)
    }

    /// Checks every option, returning [`Error::InvalidOption`] naming the first invalid one.
    fn validate(&self) -> Result<()> {
        if self.items_per_section == 0 {
            return Err(Error::InvalidOption("items_per_section must be at least 1".to_owned()));
        }
        self.zoom.validate()?;
        if let ZoomPolicy::Explicit(reductions) = &self.zoom {
            if reductions.len() > usize::from(u16::MAX) {
                return Err(Error::InvalidOption(format!(
                    "zoom: the header holds at most {} levels, got {} Explicit reductions",
                    u16::MAX,
                    reductions.len()
                )));
            }
        }
        // A larger zoom block would inflate past the size every reader of this crate accepts.
        let max_records_per_block = (MAX_INFLATED_SIZE / ZoomRecord::ENCODED_LEN) as u64;
        if self.zoom_items_per_block == 0
            || u64::from(self.zoom_items_per_block) > max_records_per_block
        {
            return Err(Error::InvalidOption(format!(
                "zoom_items_per_block must be from 1 to {max_records_per_block}, got {}",
                self.zoom_items_per_block
            )));
        }
        // The R-tree's node `count` is a u16, and a fan-out below 2 never reduces to a root.
        if !(2..=u32::from(u16::MAX)).contains(&self.index_block_size) {
            return Err(Error::InvalidOption(format!(
                "index_block_size must be from 2 to {}, got {}",
                u16::MAX,
                self.index_block_size
            )));
        }
        self.compression.validate()
    }
}

impl Default for WriterOptions {
    fn default() -> Self {
        Self::new()
    }
}

/// How a [`BigWigWriter`] assigns `chromId`s, the indices data and zoom records refer to.
///
/// Either way the IDs are dense (`0..n`), every listed chromosome is in the file's chromosome
/// tree, and the file's data, its `chromId`s and its index are all in the same order, which some
/// readers rely on.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum ChromIdOrder {
    /// IDs follow the order in which chromosomes first receive data; chromosomes that never
    /// receive data take the remaining IDs at `finish()`, in chromosome-list order. Data may be
    /// pushed in any chromosome order.
    #[default]
    Arrival,
    /// IDs are positions in the chromosome list passed to `create`, so the list states the order
    /// (e.g. BAM header order, where `chromId` equals the BAM `tid`). Data must be pushed in list
    /// order; chromosomes may be skipped. A push to a chromosome earlier in the list than the
    /// current one fails with [`Error::ChromOutOfOrder`].
    Listed,
}

/// Writes a bigWig file from streamed data. **[`finish`](Self::finish) must be called** to write
/// the indexes, the zoom levels and the real header: a writer dropped without it leaves a file
/// whose magic number is 0, which every reader rejects, and logs a warning.
///
/// Create one with [`WriterOptions::create`] or [`WriterOptions::create_with`]. Push each
/// chromosome's data sorted by position and without overlaps, one chromosome at a time, with the
/// methods under [Choosing a push method](#choosing-a-push-method). Section types may be mixed on
/// one chromosome. Pushing to a different chromosome finishes the current one, which can then
/// receive no more data. Chromosomes with no data need no calls.
///
/// Every push is validated before it changes anything: a push rejected with a data-validation
/// error ([`Error::ChromOutOfOrder`], [`Error::ChromAlreadyFinished`],
/// [`Error::UnsortedOrOverlapping`], [`Error::OutOfBounds`], [`Error::InvalidItem`],
/// [`Error::NonFiniteValue`], or [`Error::UnknownChrom`]) writes nothing, and the writer stays
/// usable. An I/O or compression error is returned by the call that hit it; every later call,
/// `finish` included, returns [`Error::WriterFailed`].
///
/// Memory is bounded: one open section, one open zoom block per level, compressed zoom blocks up
/// to [`WriterOptions::zoom_memory_limit`], and the index entries (32 bytes per data section or
/// zoom block). `BigWigWriter<W>` is `Send` when `W` is, so it can live on a writer thread.
///
/// # Choosing a push method
///
/// - [`push_fixed_step`](Self::push_fixed_step) writes fixedStep sections, the most compact form
///   for dense, evenly spaced values such as per-base depth or binned signal.
/// - [`push_fixed_step_skip`](Self::push_fixed_step_skip) takes dense values in which NaN marks
///   missing data, and writes each run between gaps as fixedStep.
/// - [`push_variable_step`](Self::push_variable_step) writes variableStep sections, for scattered
///   values.
/// - [`push_interval`](Self::push_interval) and [`push_intervals`](Self::push_intervals) write
///   bedGraph sections, for piecewise-constant data such as depth given as runs of equal values.
///
/// For data with many gaps the choice matters, because each fixedStep run is a section of its
/// own, costing roughly 65 bytes (its header, an index entry and a compressed block) besides its
/// values. With default options, variableStep wrote smaller files than `push_fixed_step_skip`
/// when the runs between gaps were shorter than roughly 20 to 50 values, the cut-off rising the
/// better the values compress (from noisy floating-point signal to integer depth). For integer
/// depth with a gap at least every few hundred values, one bedGraph interval per run of equal
/// values was smaller than either; without gaps, fixedStep was as small or smaller.
///
/// # Example
///
/// ```
/// use std::io::Cursor;
///
/// use bigwig::{ChromInfo, WriterOptions};
///
/// let chroms = [ChromInfo::new("chr2", 300), ChromInfo::new("chr1", 1_000)];
/// let mut writer = WriterOptions::new().create_with(Cursor::new(Vec::new()), &chroms)?;
///
/// // Per-base depth for chr1, pushed in chunks as it is produced: consecutive chunks extend
/// // the same fixedStep run.
/// writer.push_fixed_step("chr1", 0, 1, 1, &[3.0; 400])?;
/// writer.push_fixed_step("chr1", 400, 1, 1, &[4.0; 600])?;
///
/// // A rejected push changes nothing, so the caller may skip it and carry on.
/// assert!(writer.push_interval("chr2", 50, 20, 1.0).is_err());
/// writer.push_interval("chr2", 20, 50, 1.0)?;
///
/// let bytes = writer.finish()?.into_inner();
/// assert_eq!(bytes[..4], 0x888F_FC26_u32.to_le_bytes());
/// # Ok::<(), bigwig::Error>(())
/// ```
pub struct BigWigWriter<W: Write + Seek> {
    /// The sink, with the position of the next byte written.
    out: OutputSink<W>,
    byte_order: ByteOrder,
    items_per_section: u16,
    index_block_size: u32,
    clip: bool,
    chrom_id_order: ChromIdOrder,
    compression: Compression,
    chroms: ChromList,
    /// The chromosome receiving data; `None` before the first item is written.
    current: Option<CurrentChrom>,
    /// Where `create` put each region.
    layout: Layout,
    /// With [`ChromIdOrder::Arrival`], the chromosome tree as `create` wrote it, with placeholder
    /// IDs that `finish` replaces.
    arrival_chrom_tree: Option<EncodedChromTree>,
    /// The open data section, if `section_open`.
    section: SectionBuilder,
    /// Whether `section` holds a section not yet flushed.
    section_open: bool,
    /// The open section's values when it is fixedStep, fed to the summaries when it is flushed.
    section_values: Vec<f32>,
    /// Compresses data sections, through the batch seam a parallel executor will use.
    compressor: BlockCompressor,
    /// The batch of one section passed to the compressor, its buffers reused across sections.
    compression_job: CompressionJob,
    /// One leaf entry per data section, in file order.
    data_index: Vec<LeafEntry>,
    /// Largest decompressed data section so far.
    largest_data_section: usize,
    /// Items that clipping truncated at their chromosome's end and wrote.
    clipped_items: u64,
    /// Items that clipping dropped because they start at or past their chromosome's end.
    dropped_items: u64,
    total_summary: RawSummaryAccumulator,
    zoom: ZoomWriter,
    /// Set by an I/O or compression error; every later call then fails with `WriterFailed`.
    failed: bool,
    /// Logs a warning if the writer is dropped before `finish`.
    unfinished: UnfinishedWarning,
}

impl<W: Write + Seek> BigWigWriter<W> {
    /// Writes the placeholder header, reserved zoom-header slots and total summary, the
    /// chromosome tree and the `sectionCount` placeholder. The options and the chromosome list
    /// must already be validated.
    fn new(
        sink: W,
        chroms: ChromList,
        options: &WriterOptions,
        spill_location: SpillLocation,
        output_path: Option<PathBuf>,
    ) -> Result<Self> {
        let byte_order = options.byte_order;
        let zoom_store = ZoomBlockStore::new(options, spill_location)?;
        let zoom = ZoomWriter::new(&options.zoom, chroms.longest_length(), zoom_store);
        let reserved_zoom_slots = zoom.reserved_header_slots();

        let mut out = OutputSink::new(sink)?;
        let header_len = FIXED_HEADER_SIZE + reserved_zoom_slots * ZOOM_HEADER_SIZE;
        let total_summary_offset = header_len as u64;
        let chrom_tree_offset = total_summary_offset + TOTAL_SUMMARY_SIZE as u64;
        // The magic stays 0 until `finish`, so an unfinished file is rejected by every reader.
        out.write_all(&vec![0; header_len + TOTAL_SUMMARY_SIZE])?;

        // With `Arrival` the IDs written here are placeholders, replaced at `finish`.
        let entries: Vec<ChromEntry<'_>> = chroms
            .chroms
            .iter()
            .enumerate()
            .map(|(position, chrom)| ChromEntry {
                name: &chrom.name,
                chrom_id: list_position_id(position),
                length: chrom.length,
            })
            .collect();
        let chrom_tree = bptree::encode(&entries, chrom_tree_offset, byte_order);
        out.write_all(&chrom_tree.bytes)?;
        let full_data_offset = out.position();
        out.write_all(&[0; 8])?;

        let arrival_chrom_tree =
            (options.chrom_id_order == ChromIdOrder::Arrival).then_some(chrom_tree);
        Ok(Self {
            out,
            byte_order,
            items_per_section: options.items_per_section,
            index_block_size: options.index_block_size,
            clip: options.clip,
            chrom_id_order: options.chrom_id_order,
            compression: options.compression,
            chroms,
            current: None,
            layout: Layout {
                reserved_zoom_slots,
                total_summary_offset,
                chrom_tree_offset,
                full_data_offset,
            },
            arrival_chrom_tree,
            section: SectionBuilder::new(byte_order),
            section_open: false,
            section_values: Vec::new(),
            compressor: BlockCompressor::new(options.compression)?,
            compression_job: CompressionJob::default(),
            data_index: Vec::new(),
            largest_data_section: 0,
            clipped_items: 0,
            dropped_items: 0,
            total_summary: RawSummaryAccumulator::unbounded(),
            zoom,
            failed: false,
            unfinished: UnfinishedWarning { output: output_path, armed: true },
        })
    }

    /// Pushes a fixedStep run: `values[i]` covers `[start + i·step, start + i·step + span)`.
    ///
    /// Requires `1 ≤ span ≤ step` and a non-empty slice. A call on the same chromosome with the
    /// same step and span whose `start` continues exactly where the previous call's run ended
    /// extends the same run, and the same section until it is full, so data may be pushed in
    /// chunks of any size with identical output.
    pub fn push_fixed_step(
        &mut self,
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Result<()> {
        self.push_fixed_step_indexed(chrom, start, step, span, values)
            .map_err(|rejected| rejected.error)
    }

    /// Pushes fixedStep runs from dense values in which NaN marks missing data: `values[i]`
    /// covers `[start + i·step, start + i·step + span)`, except that a NaN writes nothing there.
    ///
    /// Writes exactly what [`push_fixed_step`](Self::push_fixed_step) writes when called once per
    /// maximal run of values without NaN, in order, so runs continue across calls as they do
    /// there and chunked calls write the same file as one call. Requires what `push_fixed_step`
    /// requires; ±infinity is still rejected, and the first run, not `start`, must not begin
    /// before the previous item's end. The call is validated whole before anything is written.
    /// A call whose values are all NaN is checked only for its chromosome name and arguments,
    /// and then writes nothing and changes nothing.
    ///
    /// Each run is a section of its own, so where runs between gaps are short (below roughly 20
    /// to 50 values) [`push_variable_step`](Self::push_variable_step) writes a smaller file; see
    /// [Choosing a push method](Self#choosing-a-push-method).
    pub fn push_fixed_step_skip(
        &mut self,
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Result<()> {
        self.check_usable()?;
        if values.iter().all(|value| value.is_nan()) {
            // Gaps alone are no data, so they neither enter nor finish a chromosome.
            self.chroms.lookup(chrom).ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))?;
            return Self::check_fixed_step_arguments(chrom, start, step, span, values);
        }
        let target = self.resolve_chrom(chrom)?;
        let checked = self
            .check_fixed_step(&target, chrom, start, step, span, values, true)
            .map_err(|rejected| rejected.error)?;
        let result =
            self.write_fixed_step_runs(&target, start, step, span, &values[..checked.kept]);
        self.fail_on_error(result)?;
        self.record_push(&target, checked);
        Ok(())
    }

    /// Pushes variableStep items: each `(start, value)` covers `[start, start + span)`.
    ///
    /// Requires `span ≥ 1` and a non-empty slice, with items sorted and not overlapping.
    pub fn push_variable_step(
        &mut self,
        chrom: &str,
        span: u32,
        items: &[(u32, f32)],
    ) -> Result<()> {
        self.push_variable_step_indexed(chrom, span, items).map_err(|rejected| rejected.error)
    }

    /// Pushes one bedGraph interval `[start, end)` with `value`. Requires `start < end`.
    pub fn push_interval(&mut self, chrom: &str, start: u32, end: u32, value: f32) -> Result<()> {
        self.push_intervals(chrom, &[(start, end, value)])
    }

    /// Pushes bedGraph intervals, each `(start, end, value)` covering `[start, end)`.
    ///
    /// Requires a non-empty slice of intervals with `start < end`, sorted and not overlapping.
    pub fn push_intervals(&mut self, chrom: &str, intervals: &[(u32, u32, f32)]) -> Result<()> {
        self.check_usable()?;
        let target = self.resolve_chrom(chrom)?;
        let checked = self.check_intervals(&target, chrom, intervals, false)?;
        let result = self.write_intervals(&target, &intervals[..checked.kept]);
        self.fail_on_error(result)?;
        self.record_push(&target, checked);
        Ok(())
    }

    /// Finishes the file and returns the sink, positioned at the end of the file.
    ///
    /// Flushes the open section and zoom records, assigns the remaining `chromId`s, writes the
    /// data index, each zoom level and the trailing magic, then patches the chromosome IDs, the
    /// header, the zoom headers, the total summary and `sectionCount`. Zoom levels that don't
    /// have fewer records than the level before are dropped, and their blocks are never copied
    /// into the file. Fails with [`Error::WriterFailed`] if an earlier call failed with an I/O or
    /// compression error.
    pub fn finish(mut self) -> Result<W> {
        // Whatever happens now, the caller has called `finish` and will see any error.
        self.unfinished.armed = false;
        self.check_usable()?;
        self.write_trailer()?;
        self.out.into_inner()
    }

    /// As [`push_fixed_step`](Self::push_fixed_step), with the index of the value an error
    /// concerns (see [`PushError`]).
    pub(crate) fn push_fixed_step_indexed(
        &mut self,
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> std::result::Result<(), PushError> {
        self.check_usable()?;
        let target = self.resolve_chrom(chrom)?;
        let checked = self.check_fixed_step(&target, chrom, start, step, span, values, false)?;
        let result = self.write_fixed_step(&target, start, step, span, &values[..checked.kept]);
        self.fail_on_error(result)?;
        self.record_push(&target, checked);
        Ok(())
    }

    /// As [`push_variable_step`](Self::push_variable_step), with the index of the item an error
    /// concerns (see [`PushError`]).
    pub(crate) fn push_variable_step_indexed(
        &mut self,
        chrom: &str,
        span: u32,
        items: &[(u32, f32)],
    ) -> std::result::Result<(), PushError> {
        self.check_usable()?;
        let target = self.resolve_chrom(chrom)?;
        let checked = self.check_variable_step(&target, chrom, span, items)?;
        let result = self.write_variable_step(&target, span, &items[..checked.kept]);
        self.fail_on_error(result)?;
        self.record_push(&target, checked);
        Ok(())
    }

    /// Validates the zero-length bedGraph interval `[at, at)` as [`push_interval`] validates any
    /// other (its chromosome, order, bounds and value), then skips it: nothing is written, but
    /// later items on its chromosome may not start before `at`, as kent requires. Fails as
    /// `push_interval` does.
    ///
    /// [`push_interval`]: Self::push_interval
    pub(crate) fn skip_empty_interval(&mut self, chrom: &str, at: u32, value: f32) -> Result<()> {
        self.check_usable()?;
        let target = self.resolve_chrom(chrom)?;
        let checked = self.check_intervals(&target, chrom, &[(at, at, value)], true)?;
        self.chroms.last_end[target.index] = checked.end;
        Ok(())
    }

    /// Items that clipping has truncated at their chromosome's end, and written.
    pub(crate) fn clipped_items(&self) -> u64 {
        self.clipped_items
    }

    /// Items that clipping has dropped because they start at or past their chromosome's end.
    pub(crate) fn dropped_items(&self) -> u64 {
        self.dropped_items
    }

    /// Fails with [`Error::WriterFailed`] if an earlier I/O or compression error left the writer
    /// unusable.
    fn check_usable(&self) -> Result<()> {
        if self.failed { Err(Error::WriterFailed) } else { Ok(()) }
    }

    /// Finds the chromosome named `chrom` and checks that it may receive data now, without
    /// changing anything.
    ///
    /// Fails with [`Error::UnknownChrom`]; with [`Error::ChromOutOfOrder`] under
    /// [`ChromIdOrder::Listed`] for a chromosome earlier in the list than the current one; and
    /// with [`Error::ChromAlreadyFinished`] for one that received data before another did.
    fn resolve_chrom(&self, chrom: &str) -> Result<ChromTarget> {
        let current = self.current.map(|current| current.index);
        let index = match current {
            Some(index) if self.chroms.chroms[index].name == chrom => index,
            _ => self.chroms.lookup(chrom).ok_or_else(|| Error::UnknownChrom(chrom.to_owned()))?,
        };
        if current != Some(index) {
            if let Some(current) = current
                .filter(|&current| self.chrom_id_order == ChromIdOrder::Listed && index < current)
            {
                return Err(Error::ChromOutOfOrder {
                    chrom: chrom.to_owned(),
                    current: self.chroms.chroms[current].name.clone(),
                });
            }
            if self.chroms.has_data[index] {
                return Err(Error::ChromAlreadyFinished(chrom.to_owned()));
            }
        }
        Ok(ChromTarget {
            index,
            length: self.chroms.chroms[index].length,
            previous_end: self.chroms.last_end[index],
        })
    }

    /// Validates a fixedStep push and returns how many leading values to write (all of them, or
    /// with clipping those that start before the chromosome end) and where its last item ends.
    /// With `nan_is_gap`, a NaN value is a gap rather than an item: only the other values are
    /// checked for order and bounds, and only they count as dropped by clipping. An error names
    /// the first value rejected, as if each had been pushed alone.
    #[expect(clippy::too_many_arguments, reason = "a fixedStep run, its target and the NaN rule")]
    fn check_fixed_step(
        &self,
        target: &ChromTarget,
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
        nan_is_gap: bool,
    ) -> std::result::Result<CheckedPush, PushError> {
        Self::check_fixed_step_arguments(chrom, start, step, span, values)?;
        // In u64 nothing below can overflow, except where saturating: every other term is below
        // 2^33. Saturating at u64::MAX still orders every later item, all of which start below
        // 2^32.
        let (start, step, span) = (u64::from(start), u64::from(step), u64::from(span));
        let item_start = |index: usize| start.saturating_add((index as u64).saturating_mul(step));
        let is_item = |value: &f32| !(nan_is_gap && value.is_nan());

        let first_item = values.iter().position(is_item);
        if let Some(index) = first_item.filter(|&index| item_start(index) < target.previous_end) {
            return Err(PushError {
                index,
                error: Error::UnsortedOrOverlapping {
                    chrom: chrom.to_owned(),
                    start: item_start(index),
                    previous_end: target.previous_end,
                },
            });
        }

        let length = u64::from(target.length);
        let count = values.len() as u64;
        let first_non_finite = values.iter().position(|value| is_item(value) && !value.is_finite());
        let kept = if self.clip {
            // Values starting before the chromosome end form a prefix of the slice.
            if start >= length { 0 } else { ((length - 1 - start) / step + 1).min(count) }
        } else {
            // Values ending within the chromosome form a prefix; the first item past it is the
            // error.
            let fitting = if start + span > length {
                0
            } else {
                ((length - start - span) / step + 1).min(count)
            };
            // `fitting` is at most `count`, a slice length.
            let fitting = usize::try_from(fitting).unwrap_or(values.len());
            let first_past_end =
                values[fitting..].iter().position(is_item).map(|offset| fitting + offset);
            // An item past the end is rejected unless a non-finite value comes before it.
            if let Some(index) = first_past_end
                .filter(|&index| first_non_finite.is_none_or(|non_finite| non_finite >= index))
            {
                return Err(PushError {
                    index,
                    error: Error::OutOfBounds {
                        chrom: chrom.to_owned(),
                        start: item_start(index),
                        end: item_start(index).saturating_add(span),
                        length: target.length,
                    },
                });
            }
            count
        };
        if let Some(index) = first_non_finite {
            return Err(PushError {
                index,
                error: Error::NonFiniteValue {
                    chrom: chrom.to_owned(),
                    position: item_start(index),
                    value: values[index],
                },
            });
        }
        // `kept` is at most `count`, a slice length.
        let kept = usize::try_from(kept).unwrap_or(values.len());
        let last_item = values.iter().rposition(is_item);
        Ok(CheckedPush {
            kept,
            dropped: values[kept..].iter().filter(|value| is_item(value)).count(),
            end: last_item
                .map_or(target.previous_end, |index| item_start(index).saturating_add(span)),
        })
    }

    /// Checks a fixedStep push's arguments on their own: a non-empty slice and `1 ≤ span ≤ step`.
    fn check_fixed_step_arguments(
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Result<()> {
        let invalid = |reason| Error::InvalidItem {
            chrom: chrom.to_owned(),
            start: u64::from(start),
            reason,
        };
        if values.is_empty() {
            return Err(invalid("empty slice of values"));
        }
        if step == 0 {
            return Err(invalid("fixedStep step must be at least 1"));
        }
        if span == 0 {
            return Err(invalid("span must be at least 1"));
        }
        if span > step {
            return Err(invalid("fixedStep span must not exceed step, or its items would overlap"));
        }
        Ok(())
    }

    /// Validates a variableStep push and returns how many leading items to write (all of them,
    /// or with clipping those that start before the chromosome end) and where the last one ends.
    /// An error names the first item rejected.
    fn check_variable_step(
        &self,
        target: &ChromTarget,
        chrom: &str,
        span: u32,
        items: &[(u32, f32)],
    ) -> std::result::Result<CheckedPush, PushError> {
        let Some(&(first_start, _)) = items.first() else {
            return Err(PushError::from(Error::InvalidItem {
                chrom: chrom.to_owned(),
                start: 0,
                reason: "empty slice of items",
            }));
        };
        if span == 0 {
            return Err(PushError::from(Error::InvalidItem {
                chrom: chrom.to_owned(),
                start: u64::from(first_start),
                reason: "span must be at least 1",
            }));
        }
        let length = u64::from(target.length);
        let mut previous_end = target.previous_end;
        let mut kept = 0;
        for (index, &(start, value)) in items.iter().enumerate() {
            let end = u64::from(start) + u64::from(span);
            let rejected = |error| PushError { index, error };
            if u64::from(start) < previous_end {
                return Err(rejected(Error::UnsortedOrOverlapping {
                    chrom: chrom.to_owned(),
                    start: u64::from(start),
                    previous_end,
                }));
            }
            if end > length && !self.clip {
                return Err(rejected(Error::OutOfBounds {
                    chrom: chrom.to_owned(),
                    start: u64::from(start),
                    end,
                    length: target.length,
                }));
            }
            if !value.is_finite() {
                return Err(rejected(Error::NonFiniteValue {
                    chrom: chrom.to_owned(),
                    position: u64::from(start),
                    value,
                }));
            }
            // Items are sorted, so those starting before the chromosome end form a prefix.
            if u64::from(start) < length {
                kept += 1;
            }
            previous_end = end;
        }
        Ok(CheckedPush { kept, dropped: items.len() - kept, end: previous_end })
    }

    /// Validates a bedGraph push and returns how many leading intervals to write (all of them,
    /// or with clipping those that start before the chromosome end) and where the last one ends.
    /// Zero-length intervals are rejected unless `allow_empty`.
    fn check_intervals(
        &self,
        target: &ChromTarget,
        chrom: &str,
        intervals: &[(u32, u32, f32)],
        allow_empty: bool,
    ) -> Result<CheckedPush> {
        if intervals.is_empty() {
            return Err(Error::InvalidItem {
                chrom: chrom.to_owned(),
                start: 0,
                reason: "empty slice of intervals",
            });
        }
        let mut previous_end = target.previous_end;
        let mut kept = 0;
        for &(start, end, value) in intervals {
            if start > end || (start == end && !allow_empty) {
                return Err(Error::InvalidItem {
                    chrom: chrom.to_owned(),
                    start: u64::from(start),
                    reason: "interval start must be below its end",
                });
            }
            if u64::from(start) < previous_end {
                return Err(Error::UnsortedOrOverlapping {
                    chrom: chrom.to_owned(),
                    start: u64::from(start),
                    previous_end,
                });
            }
            if end > target.length && !self.clip {
                return Err(Error::OutOfBounds {
                    chrom: chrom.to_owned(),
                    start: u64::from(start),
                    end: u64::from(end),
                    length: target.length,
                });
            }
            if !value.is_finite() {
                return Err(Error::NonFiniteValue {
                    chrom: chrom.to_owned(),
                    position: u64::from(start),
                    value,
                });
            }
            // Intervals are sorted, so those starting before the chromosome end form a prefix.
            if start < target.length {
                kept += 1;
            }
            previous_end = u64::from(end);
        }
        Ok(CheckedPush { kept, dropped: intervals.len() - kept, end: previous_end })
    }

    /// Records an accepted push: where its last item ends, which later items on the chromosome
    /// may not start before, and how many items clipping dropped.
    fn record_push(&mut self, target: &ChromTarget, checked: CheckedPush) {
        self.chroms.last_end[target.index] = checked.end;
        self.dropped_items += checked.dropped as u64;
    }

    /// Writes validated fixedStep values with NaN values as gaps, each run between gaps as
    /// [`write_fixed_step`](Self::write_fixed_step) writes a push of it.
    fn write_fixed_step_runs(
        &mut self,
        target: &ChromTarget,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Result<()> {
        let mut run_index = 0;
        for run in values.split(|value| value.is_nan()) {
            if !run.is_empty() {
                let run_start = u64::from(start) + run_index as u64 * u64::from(step);
                self.write_fixed_step(target, within_chrom(run_start), step, span, run)?;
            }
            // The NaN after the run.
            run_index += run.len() + 1;
        }
        Ok(())
    }

    /// Writes validated fixedStep values, each starting before the chromosome end. Only the last
    /// can extend past it (with clipping); it is truncated into a one-item section of its own.
    fn write_fixed_step(
        &mut self,
        target: &ChromTarget,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Result<()> {
        let Some((&last_value, leading)) = values.split_last() else { return Ok(()) };
        let chrom_id = self.enter_chrom(target)?;
        let last_start = u64::from(start) + leading.len() as u64 * u64::from(step);
        let length = u64::from(target.length);
        if last_start + u64::from(span) > length {
            self.clipped_items += 1;
            self.write_fixed_step_run(chrom_id, start, step, span, leading)?;
            let clipped_span = within_chrom(length - last_start);
            self.write_fixed_step_run(
                chrom_id,
                within_chrom(last_start),
                step,
                clipped_span,
                &[last_value],
            )
        } else {
            self.write_fixed_step_run(chrom_id, start, step, span, values)
        }
    }

    /// Appends a fixedStep run to the open section if it continues it, and otherwise to new
    /// sections, flushing each section that fills.
    fn write_fixed_step_run(
        &mut self,
        chrom_id: u32,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Result<()> {
        let kind = SectionKind::FixedStep { step, span };
        let mut next_start = u64::from(start);
        let mut remaining = values;
        while !remaining.is_empty() {
            self.open_section(chrom_id, kind, within_chrom(next_start))?;
            let room = usize::from(self.items_per_section - self.section.item_count());
            let (now, later) = remaining.split_at(room.min(remaining.len()));
            self.section.push_fixed_step_values(now);
            self.section_values.extend_from_slice(now);
            next_start += now.len() as u64 * u64::from(step);
            remaining = later;
            if self.section.item_count() == self.items_per_section {
                self.flush_section()?;
            }
        }
        Ok(())
    }

    /// Writes validated variableStep items, each starting before the chromosome end. An item
    /// extending past it (with clipping) is truncated into a section with a shorter span.
    fn write_variable_step(
        &mut self,
        target: &ChromTarget,
        span: u32,
        items: &[(u32, f32)],
    ) -> Result<()> {
        if items.is_empty() {
            return Ok(());
        }
        let chrom_id = self.enter_chrom(target)?;
        for &(start, value) in items {
            let item_span = span.min(target.length - start);
            if item_span < span {
                self.clipped_items += 1;
            }
            self.open_section(chrom_id, SectionKind::VariableStep { span: item_span }, start)?;
            self.section.push_variable_step(start, value);
            self.add_item(chrom_id, start, start + item_span, value)?;
        }
        Ok(())
    }

    /// Writes validated bedGraph intervals, each starting before the chromosome end, truncating
    /// one that extends past it (with clipping).
    fn write_intervals(
        &mut self,
        target: &ChromTarget,
        intervals: &[(u32, u32, f32)],
    ) -> Result<()> {
        if intervals.is_empty() {
            return Ok(());
        }
        let chrom_id = self.enter_chrom(target)?;
        for &(start, end, value) in intervals {
            if end > target.length {
                self.clipped_items += 1;
            }
            let end = end.min(target.length);
            self.open_section(chrom_id, SectionKind::BedGraph, start)?;
            self.section.push_bed_graph(start, end, value);
            self.add_item(chrom_id, start, end, value)?;
        }
        Ok(())
    }

    /// Makes `target` the current chromosome, first finishing the previous one by flushing its
    /// open section. Returns the chromosome's `chromId`, assigning it under
    /// [`ChromIdOrder::Arrival`] if this is its first data.
    fn enter_chrom(&mut self, target: &ChromTarget) -> Result<u32> {
        if let Some(current) = self.current.filter(|current| current.index == target.index) {
            return Ok(current.chrom_id);
        }
        self.flush_section()?;
        let chrom_id = self.chroms.start_data(target.index);
        self.current = Some(CurrentChrom { index: target.index, chrom_id });
        Ok(chrom_id)
    }

    /// Makes the open section one that can take an item of `kind` on `chrom_id` starting at
    /// `start`, flushing the open section and beginning a new one unless it is of the same kind
    /// on the same chromosome and, for fixedStep, `start` continues its run.
    fn open_section(&mut self, chrom_id: u32, kind: SectionKind, start: u32) -> Result<()> {
        if self.section_open
            && self.section.chrom_id() == chrom_id
            && self.section.kind() == kind
            && self.section.item_count() < self.items_per_section
        {
            let continues = match kind {
                SectionKind::FixedStep { step, .. } => {
                    u64::from(self.section.start())
                        + u64::from(self.section.item_count()) * u64::from(step)
                        == u64::from(start)
                }
                SectionKind::BedGraph | SectionKind::VariableStep { .. } => true,
            };
            if continues {
                return Ok(());
            }
        }
        self.flush_section()?;
        match kind {
            SectionKind::BedGraph => self.section.begin_bed_graph(chrom_id),
            SectionKind::VariableStep { span } => self.section.begin_variable_step(chrom_id, span),
            SectionKind::FixedStep { step, span } => {
                self.section.begin_fixed_step(chrom_id, start, step, span);
            }
        }
        self.section_open = true;
        Ok(())
    }

    /// Feeds one written bedGraph or variableStep item to the total summary and the zoom levels,
    /// then flushes the open section if it is full.
    fn add_item(&mut self, chrom_id: u32, start: u32, end: u32, value: f32) -> Result<()> {
        self.total_summary.add_item(start, end, value);
        self.zoom.add_item(chrom_id, start, end, value)?;
        if self.section.item_count() == self.items_per_section {
            self.flush_section()?;
        }
        Ok(())
    }

    /// Encodes, compresses and writes the open section, if any, and records its index entry. A
    /// fixedStep section's values are fed to the total summary and the zoom levels here.
    fn flush_section(&mut self) -> Result<()> {
        if !self.section_open {
            return Ok(());
        }
        self.section_open = false;
        let (chrom_id, start, end, kind) = (
            self.section.chrom_id(),
            self.section.start(),
            self.section.end(),
            self.section.kind(),
        );

        let raw = self.section.finish();
        self.largest_data_section = self.largest_data_section.max(raw.len());
        let job = &mut self.compression_job;
        job.raw.clear();
        job.raw.extend_from_slice(raw);
        self.compressor.compress_batch(std::slice::from_mut(job))?;
        let data_offset = self.out.position();
        self.out.write_all(&job.compressed)?;
        self.data_index.push(LeafEntry {
            start_chrom_ix: chrom_id,
            start_base: start,
            end_chrom_ix: chrom_id,
            end_base: end,
            data_offset,
            data_size: job.compressed.len() as u64,
        });

        if let SectionKind::FixedStep { step, span } = kind {
            self.total_summary.add_fixed_step(start, step, span, &self.section_values);
            self.zoom.add_fixed_step(chrom_id, start, step, span, &self.section_values)?;
            self.section_values.clear();
        }
        Ok(())
    }

    /// Writes everything after the data sections and patches what `create` left as placeholders.
    fn write_trailer(&mut self) -> Result<()> {
        self.flush_section()?;
        self.zoom.finish_records()?;
        self.chroms.assign_remaining_ids();

        let full_index_offset = self.out.position();
        let data_index = rtree::encode(
            &self.data_index,
            full_index_offset,
            self.index_block_size,
            1,
            self.byte_order,
        );
        self.out.write_all(&data_index)?;
        let zoom_levels = self.zoom.write_kept_levels(&mut self.out, self.index_block_size)?;
        self.out.write_all(&self.byte_order.u32_bytes(BIGWIG_MAGIC))?;

        self.patch_chrom_ids()?;
        self.patch_header(full_index_offset, &zoom_levels)
    }

    /// With [`ChromIdOrder::Arrival`], rewrites the chromosome tree with every chromosome's final
    /// `chromId` in place of the placeholder. The tree's size doesn't depend on the IDs, so it
    /// is rewritten in one piece.
    fn patch_chrom_ids(&mut self) -> Result<()> {
        let Some(tree) = &mut self.arrival_chrom_tree else { return Ok(()) };
        let tree_offset = self.layout.chrom_tree_offset;
        for (position, &field_offset) in tree.chrom_id_offsets.iter().enumerate() {
            let at = usize::try_from(field_offset - tree_offset)
                .expect("the tree was encoded in memory, so its offsets fit in usize");
            let chrom_id = self.chroms.chrom_ids[position].expect("every chromId is assigned");
            tree.bytes[at..at + 4].copy_from_slice(&self.byte_order.u32_bytes(chrom_id));
        }
        self.out.overwrite(tree_offset, &tree.bytes)
    }

    /// Writes the real fixed header, the kept levels' zoom headers (the other reserved slots
    /// stay zero), the total summary and `sectionCount`.
    fn patch_header(&mut self, full_index_offset: u64, zoom_levels: &[ZoomLevel]) -> Result<()> {
        let uncompress_buf_size = match self.compression {
            Compression::None => 0,
            Compression::Deflate { .. } => {
                // Kent's floor for a file without blocks: readers treat 0 as "uncompressed".
                let floor = self.zoom.store.records_per_block as usize * ZoomRecord::ENCODED_LEN;
                let largest = self.largest_data_section.max(self.zoom.largest_kept_block());
                let largest = if largest == 0 { floor } else { largest };
                u32::try_from(largest)
                    .expect("sections and zoom blocks are validated to stay below 64 MiB")
            }
        };
        let header = FixedHeader {
            version: BIGWIG_VERSION,
            zoom_levels: u16::try_from(zoom_levels.len())
                .expect("the zoom policy is validated to fit the header"),
            chrom_tree_offset: self.layout.chrom_tree_offset,
            full_data_offset: self.layout.full_data_offset,
            full_index_offset,
            field_count: 0,
            defined_field_count: 0,
            auto_sql_offset: 0,
            total_summary_offset: self.layout.total_summary_offset,
            uncompress_buf_size,
            extension_offset: 0,
        };
        // A file with no data gets an all-zero summary, as kent's `bedGraphToBigWig` writes.
        let total_summary = self.total_summary.summary().unwrap_or(Summary {
            covered_bases: 0.0,
            min: 0.0,
            max: 0.0,
            sum: 0.0,
            sum_squares: 0.0,
        });

        let mut prefix = Vec::with_capacity(
            FIXED_HEADER_SIZE
                + self.layout.reserved_zoom_slots * ZOOM_HEADER_SIZE
                + TOTAL_SUMMARY_SIZE,
        );
        header.encode(self.byte_order, &mut prefix);
        for level in zoom_levels {
            level.encode(self.byte_order, &mut prefix);
        }
        prefix.resize(FIXED_HEADER_SIZE + self.layout.reserved_zoom_slots * ZOOM_HEADER_SIZE, 0);
        encode_total_summary(&total_summary, self.byte_order, &mut prefix);
        debug_assert_eq!(prefix.len() as u64, self.layout.chrom_tree_offset);
        self.out.overwrite(0, &prefix)?;

        let section_count = self.data_index.len() as u64;
        self.out.overwrite(self.layout.full_data_offset, &self.byte_order.u64_bytes(section_count))
    }

    /// Marks the writer unusable if `result` is an error, and returns it.
    fn fail_on_error(&mut self, result: Result<()>) -> Result<()> {
        if result.is_err() {
            self.failed = true;
        }
        result
    }
}

/// Where `create` put the regions before the data, for `finish` to fill in.
#[derive(Clone, Copy, Debug)]
struct Layout {
    /// Zoom-header slots reserved after the fixed header; unused ones stay zero.
    reserved_zoom_slots: usize,
    /// Offset of the total summary, right after the reserved slots.
    total_summary_offset: u64,
    /// Offset of the chromosome B+ tree, right after the total summary.
    chrom_tree_offset: u64,
    /// Offset of `sectionCount`, right after the tree; the data sections follow it.
    full_data_offset: u64,
}

/// The validated chromosome list, with each chromosome's `chromId` and whether it has data.
#[derive(Debug)]
struct ChromList {
    /// The chromosomes in the caller's order.
    chroms: Vec<ChromInfo>,
    /// List position of each name.
    index_by_name: HashMap<String, usize>,
    /// `chromId` of each chromosome by list position: all assigned at creation with
    /// [`ChromIdOrder::Listed`], at first data (or at `finish`) with [`ChromIdOrder::Arrival`].
    chrom_ids: Vec<Option<u32>>,
    /// Whether each chromosome has received data, by list position.
    has_data: Vec<bool>,
    /// End of the last item accepted for each chromosome before clipping, by list position, or 0
    /// for none. Items that clipping dropped count, so where a push is split doesn't change
    /// which pushes are accepted.
    last_end: Vec<u64>,
    /// The next `chromId` to assign under [`ChromIdOrder::Arrival`].
    next_arrival_id: u32,
}

impl ChromList {
    /// Validates `chroms` and indexes them by name.
    ///
    /// Fails with [`Error::InvalidChromList`] for an empty list, a duplicate name, a name that is
    /// empty or has a byte outside printable ASCII other than space (0x21 to 0x7E), or a length
    /// of 0. Logs a warning for each chromosome longer than 2³¹ − 1, which readers using signed
    /// 32-bit coordinates mishandle.
    fn new(chroms: &[ChromInfo], order: ChromIdOrder) -> Result<Self> {
        if chroms.is_empty() {
            return Err(Error::InvalidChromList("the list is empty".to_owned()));
        }
        if u32::try_from(chroms.len()).is_err() {
            return Err(Error::InvalidChromList(format!(
                "{} chromosomes is more than the {} chromIds can number",
                chroms.len(),
                u32::MAX
            )));
        }
        let mut index_by_name = HashMap::with_capacity(chroms.len());
        for (position, chrom) in chroms.iter().enumerate() {
            if chrom.name.is_empty() {
                return Err(Error::InvalidChromList(format!(
                    "the name at position {position} is empty"
                )));
            }
            if let Some(byte) = chrom.name.bytes().find(|byte| !(0x21..=0x7E).contains(byte)) {
                return Err(Error::InvalidChromList(format!(
                    "name {:?} contains byte {byte:#04x}; names must be printable ASCII other \
                     than space",
                    chrom.name
                )));
            }
            if chrom.length == 0 {
                return Err(Error::InvalidChromList(format!("{} has length 0", chrom.name)));
            }
            if index_by_name.insert(chrom.name.clone(), position).is_some() {
                return Err(Error::InvalidChromList(format!(
                    "{} is listed more than once",
                    chrom.name
                )));
            }
        }
        for chrom in chroms.iter().filter(|chrom| chrom.length > LONGEST_SIGNED_CHROM_LENGTH) {
            log::warn!(
                "chromosome {} is {} bases long, longer than {LONGEST_SIGNED_CHROM_LENGTH}; \
                 IGV, igv.js and JBrowse read coordinates as signed 32-bit integers and may \
                 mishandle data beyond that",
                chrom.name,
                chrom.length
            );
        }
        let chrom_ids = match order {
            ChromIdOrder::Arrival => vec![None; chroms.len()],
            ChromIdOrder::Listed => (0..chroms.len()).map(|p| Some(list_position_id(p))).collect(),
        };
        Ok(Self {
            chroms: chroms.to_vec(),
            index_by_name,
            chrom_ids,
            has_data: vec![false; chroms.len()],
            last_end: vec![0; chroms.len()],
            next_arrival_id: 0,
        })
    }

    /// List position of the chromosome named `name`.
    fn lookup(&self, name: &str) -> Option<usize> {
        self.index_by_name.get(name).copied()
    }

    /// Length of the longest chromosome.
    fn longest_length(&self) -> u32 {
        self.chroms.iter().map(|chrom| chrom.length).max().unwrap_or(0)
    }

    /// Marks the chromosome at list position `index` as having data and returns its `chromId`,
    /// assigning the next arrival ID if it has none yet.
    fn start_data(&mut self, index: usize) -> u32 {
        self.has_data[index] = true;
        *self.chrom_ids[index].get_or_insert_with(|| {
            let chrom_id = self.next_arrival_id;
            self.next_arrival_id += 1;
            chrom_id
        })
    }

    /// Assigns arrival IDs to the chromosomes that received no data, in list order.
    fn assign_remaining_ids(&mut self) {
        for chrom_id in self.chrom_ids.iter_mut().filter(|chrom_id| chrom_id.is_none()) {
            *chrom_id = Some(self.next_arrival_id);
            self.next_arrival_id += 1;
        }
    }
}

/// A chromosome a push is for, as resolved before anything is written.
#[derive(Clone, Copy, Debug)]
struct ChromTarget {
    /// Position in the chromosome list.
    index: usize,
    /// Length in bases.
    length: u32,
    /// End of the last item accepted for it before clipping, or 0 if it has none; new items may
    /// not start before.
    previous_end: u64,
}

/// A validated push: how many leading items to write, how many clipping dropped, and where its
/// last item ends.
#[derive(Clone, Copy, Debug)]
struct CheckedPush {
    /// Leading items (or fixedStep values, gaps included) to write: all of them, or with
    /// clipping those starting before the chromosome end.
    kept: usize,
    /// Items after the kept ones, which clipping dropped.
    dropped: usize,
    /// End of the last item before clipping, which later items on the chromosome may not start
    /// before; for a push without items, the previous end.
    end: u64,
}

/// An error from a push, with the index in the pushed slice of the item it concerns, so a caller
/// can say where in its input the error is.
#[derive(Debug)]
pub(crate) struct PushError {
    /// The item rejected, for a data-validation error about one item; 0 for any other error,
    /// such as one about the chromosome or the whole push.
    pub(crate) index: usize,
    pub(crate) error: Error,
}

impl From<Error> for PushError {
    fn from(error: Error) -> Self {
        Self { index: 0, error }
    }
}

/// The chromosome currently receiving data.
#[derive(Clone, Copy, Debug)]
struct CurrentChrom {
    /// Position in the chromosome list.
    index: usize,
    chrom_id: u32,
}

/// The writer's sink behind a buffer, with the absolute offset of the next byte written.
///
/// The position is tracked here because seeking a `BufWriter` (even to ask its position) flushes
/// it.
struct OutputSink<W: Write + Seek> {
    inner: BufWriter<W>,
    /// Offset just past the last byte appended; overwrites before it don't move it.
    position: u64,
}

impl<W: Write + Seek> OutputSink<W> {
    /// Buffers `sink`, which must be empty: the file is written from offset 0 and the sink never
    /// truncated, so bytes beyond the file's end would remain. Fails with [`Error::Io`] of kind
    /// [`InvalidInput`](io::ErrorKind::InvalidInput) otherwise.
    fn new(mut sink: W) -> Result<Self> {
        // Seeking to the end of an empty sink also positions it at offset 0.
        let len = sink.seek(SeekFrom::End(0))?;
        if len != 0 {
            return Err(Error::Io(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("the sink must be empty, but it holds {len} bytes"),
            )));
        }
        Ok(Self { inner: BufWriter::with_capacity(SINK_BUFFER_SIZE, sink), position: 0 })
    }

    /// Offset of the next byte appended.
    fn position(&self) -> u64 {
        self.position
    }

    /// Appends `bytes`.
    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        self.inner.write_all(bytes)?;
        self.position += bytes.len() as u64;
        Ok(())
    }

    /// Appends exactly `len` bytes read from `source`.
    fn copy_from(&mut self, source: &mut impl Read, len: u64) -> Result<()> {
        let copied = io::copy(&mut source.take(len), &mut self.inner)?;
        if copied != len {
            return Err(Error::Io(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("zoom spill file ended after {copied} of {len} bytes"),
            )));
        }
        self.position += len;
        Ok(())
    }

    /// Overwrites the bytes at `offset`, which must lie before the current position.
    fn overwrite(&mut self, offset: u64, bytes: &[u8]) -> Result<()> {
        debug_assert!(offset + bytes.len() as u64 <= self.position);
        self.inner.seek(SeekFrom::Start(offset))?;
        self.inner.write_all(bytes)?;
        Ok(())
    }

    /// Flushes the buffer and returns the sink, positioned at the end of what was written.
    fn into_inner(mut self) -> Result<W> {
        self.inner.seek(SeekFrom::Start(self.position))?;
        self.inner.into_inner().map_err(|error| Error::Io(error.into_error()))
    }
}

/// Turns the data into zoom records and stores them in blocks: the zoom policy's ladder, the
/// pitch sample it may need first, and the block store.
struct ZoomWriter {
    policy: ZoomPolicy,
    /// Length of the longest listed chromosome, which bounds the ladder.
    longest_chrom_length: u32,
    stage: ZoomStage,
    /// Kept apart from `stage` so the builder can emit records into it.
    store: ZoomBlockStore,
}

impl ZoomWriter {
    /// Prepares zoom output for `policy`. A ladder that doesn't depend on the data is fixed now.
    fn new(policy: &ZoomPolicy, longest_chrom_length: u32, mut store: ZoomBlockStore) -> Self {
        let stage = if policy.samples_pitch() {
            ZoomStage::Sampling(PitchSample::default())
        } else {
            let reductions = policy.reductions(&PitchSampler::new(), longest_chrom_length);
            if reductions.is_empty() {
                ZoomStage::Disabled
            } else {
                store.add_levels(&reductions);
                ZoomStage::Accumulating(ZoomBuilder::new(&reductions))
            }
        };
        Self { policy: policy.clone(), longest_chrom_length, stage, store }
    }

    /// Zoom-header slots to reserve after the fixed header: none without zoom levels, otherwise
    /// at least the ten kent's `bedGraphToBigWig` reserves, and one per level of a longer
    /// `Explicit` ladder.
    fn reserved_header_slots(&self) -> usize {
        match self.stage {
            ZoomStage::Disabled => 0,
            ZoomStage::Sampling(_) | ZoomStage::Accumulating(_) => {
                self.store.levels.len().max(AUTO_MAX_LEVELS)
            }
        }
    }

    /// Adds a written fixedStep section's values.
    fn add_fixed_step(
        &mut self,
        chrom_id: u32,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Result<()> {
        let (mut start, mut values) = (start, values);
        if let ZoomStage::Sampling(sample) = &mut self.stage {
            let taken = sample.add_fixed_step(chrom_id, start, step, span, values);
            if !sample.sampler.is_full() {
                return Ok(());
            }
            self.fix_ladder()?;
            if taken == values.len() {
                return Ok(());
            }
            start = within_chrom(u64::from(start) + taken as u64 * u64::from(step));
            values = &values[taken..];
        }
        let Self { stage, store, .. } = self;
        if let ZoomStage::Accumulating(builder) = stage {
            builder.push_fixed_step(
                chrom_id,
                start,
                step,
                span,
                values,
                &mut |level, record| store.add_record(level, record),
            )?;
        }
        Ok(())
    }

    /// Adds one written bedGraph or variableStep item.
    fn add_item(&mut self, chrom_id: u32, start: u32, end: u32, value: f32) -> Result<()> {
        let Self { stage, store, .. } = self;
        match stage {
            ZoomStage::Disabled => Ok(()),
            ZoomStage::Sampling(sample) => {
                sample.add_item(chrom_id, start, end, value);
                if sample.sampler.is_full() { self.fix_ladder() } else { Ok(()) }
            }
            ZoomStage::Accumulating(builder) => {
                builder.push_item(chrom_id, start, end, value, &mut |level, record| {
                    store.add_record(level, record)
                })
            }
        }
    }

    /// Completes every record and block: fixes the ladder from whatever was sampled if that
    /// hasn't happened yet, closes the open records and flushes the open blocks.
    fn finish_records(&mut self) -> Result<()> {
        self.fix_ladder()?;
        let Self { stage, store, .. } = self;
        if let ZoomStage::Accumulating(builder) = stage {
            builder.close(&mut |level, record| store.add_record(level, record))?;
        }
        self.store.flush_open_blocks()
    }

    /// Ends sampling: fixes the ladder from the sampled pitch and replays the sample through a
    /// new builder. Does nothing unless sampling.
    fn fix_ladder(&mut self) -> Result<()> {
        if !matches!(self.stage, ZoomStage::Sampling(_)) {
            return Ok(());
        }
        let ZoomStage::Sampling(sample) = std::mem::replace(&mut self.stage, ZoomStage::Disabled)
        else {
            unreachable!("the stage was just checked to be sampling");
        };
        let reductions = self.policy.reductions(&sample.sampler, self.longest_chrom_length);
        self.store.add_levels(&reductions);
        let mut builder = ZoomBuilder::new(&reductions);
        let store = &mut self.store;
        let mut emit = |level: usize, record: ZoomRecord| store.add_record(level, record);
        for unit in &sample.units {
            match unit {
                SampledUnit::FixedStep { chrom_id, start, step, span, values } => {
                    let values = &sample.values[values.clone()];
                    builder.push_fixed_step(*chrom_id, *start, *step, *span, values, &mut emit)?;
                }
                SampledUnit::Items(items) => {
                    for item in &sample.items[items.clone()] {
                        builder.push_item(
                            item.chrom_id,
                            item.start,
                            item.end,
                            item.value,
                            &mut emit,
                        )?;
                    }
                }
            }
        }
        self.stage = ZoomStage::Accumulating(builder);
        Ok(())
    }

    /// Writes each kept level's `zoomCount`, blocks and index, returning the levels' headers.
    fn write_kept_levels<W: Write + Seek>(
        &mut self,
        out: &mut OutputSink<W>,
        index_block_size: u32,
    ) -> Result<Vec<ZoomLevel>> {
        (0..self.store.kept_level_count())
            .map(|level| self.store.write_level(level, out, index_block_size))
            .collect()
    }

    /// Largest decompressed block among the kept levels.
    fn largest_kept_block(&self) -> usize {
        self.store.levels[..self.store.kept_level_count()]
            .iter()
            .map(|level| level.largest_block)
            .max()
            .unwrap_or(0)
    }
}

/// Where zoom accumulation stands.
#[derive(Debug)]
enum ZoomStage {
    /// [`ZoomPolicy::None`]: no zoom levels.
    Disabled,
    /// The ladder depends on the data's pitch, which is still being measured.
    Sampling(PitchSample),
    /// The ladder is fixed and records are being accumulated.
    Accumulating(ZoomBuilder),
}

/// The leading data, kept while the `Auto` policy measures its pitch, to replay into the zoom
/// builder once the ladder is fixed. It holds at most [`crate::zoom::PITCH_SAMPLE_ITEMS`] items.
#[derive(Debug, Default)]
struct PitchSample {
    sampler: PitchSampler,
    /// The sampled data in the order it arrived.
    units: Vec<SampledUnit>,
    /// Values of the sampled fixedStep runs.
    values: Vec<f32>,
    /// The sampled bedGraph and variableStep items.
    items: Vec<SampledItem>,
}

impl PitchSample {
    /// Samples the leading values of a fixedStep section that fit the sample, returning how many.
    fn add_fixed_step(
        &mut self,
        chrom_id: u32,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> usize {
        let taken = self.sampler.offer_fixed_step(chrom_id, start, step, span, values.len());
        if taken > 0 {
            let first = self.values.len();
            self.values.extend_from_slice(&values[..taken]);
            let values = first..self.values.len();
            self.units.push(SampledUnit::FixedStep { chrom_id, start, step, span, values });
        }
        taken
    }

    /// Samples one item. The sample must not be full.
    fn add_item(&mut self, chrom_id: u32, start: u32, end: u32, value: f32) {
        let sampled = self.sampler.offer_item(chrom_id, start, end);
        debug_assert!(sampled, "the ladder is fixed as soon as the sample is full");
        let index = self.items.len();
        self.items.push(SampledItem { chrom_id, start, end, value });
        match self.units.last_mut() {
            Some(SampledUnit::Items(items)) => items.end = index + 1,
            _ => self.units.push(SampledUnit::Items(index..index + 1)),
        }
    }
}

/// A run of sampled data, in arrival order.
#[derive(Clone, Debug)]
enum SampledUnit {
    /// A fixedStep section's sampled values: a range of [`PitchSample::values`].
    FixedStep { chrom_id: u32, start: u32, step: u32, span: u32, values: Range<usize> },
    /// Consecutive items: a range of [`PitchSample::items`].
    Items(Range<usize>),
}

/// A sampled bedGraph or variableStep item.
#[derive(Clone, Copy, Debug)]
struct SampledItem {
    chrom_id: u32,
    start: u32,
    end: u32,
    value: f32,
}

/// Packs each level's zoom records into compressed blocks and holds them until `finish`: in
/// memory up to a limit, then in one spill file shared by every level.
struct ZoomBlockStore {
    byte_order: ByteOrder,
    /// Records per block (`zoom_items_per_block`).
    records_per_block: u32,
    /// Compressed bytes held in memory, over all levels, beyond which they spill.
    memory_limit: usize,
    /// One entry per level of the ladder, in increasing reduction order; empty until it is fixed.
    levels: Vec<ZoomLevelBlocks>,
    compressor: BlockCompressor,
    /// Reused output buffer for compressing a block.
    compressed: Vec<u8>,
    /// Compressed bytes currently held in memory, over all levels.
    in_memory: usize,
    /// The spill file, created once memory first exceeds the limit.
    spill: Option<SpillFile>,
    spill_location: SpillLocation,
}

impl ZoomBlockStore {
    /// Creates an empty store with the options' block size, compression and memory limit.
    fn new(options: &WriterOptions, spill_location: SpillLocation) -> Result<Self> {
        Ok(Self {
            byte_order: options.byte_order,
            records_per_block: options.zoom_items_per_block,
            memory_limit: options.zoom_memory_limit,
            levels: Vec::new(),
            compressor: BlockCompressor::new(options.compression)?,
            compressed: Vec::new(),
            in_memory: 0,
            spill: None,
            spill_location,
        })
    }

    /// Adds an empty level for each of `reductions`.
    fn add_levels(&mut self, reductions: &[u32]) {
        self.levels.extend(reductions.iter().map(|&reduction| ZoomLevelBlocks::new(reduction)));
    }

    /// Appends a completed record to its level's open block, starting a new block on a new
    /// chromosome so that no block spans two, and flushing the block once it is full.
    fn add_record(&mut self, level_index: usize, record: ZoomRecord) -> Result<()> {
        let level = &self.levels[level_index];
        if level.open_records > 0 && level.open_chrom_id != record.chrom_id {
            self.flush_block(level_index)?;
        }
        self.levels[level_index].push_record(&record, self.byte_order);
        if self.levels[level_index].open_records == self.records_per_block {
            self.flush_block(level_index)?;
        }
        Ok(())
    }

    /// Flushes every level's open block.
    fn flush_open_blocks(&mut self) -> Result<()> {
        for level_index in 0..self.levels.len() {
            self.flush_block(level_index)?;
        }
        Ok(())
    }

    /// Compresses a level's open block, if it has records, into the level's stored blocks and
    /// records its index entry; spills every level once memory exceeds the limit.
    fn flush_block(&mut self, level_index: usize) -> Result<()> {
        let level = &mut self.levels[level_index];
        if level.open_records == 0 {
            return Ok(());
        }
        self.compressor.compress(&level.open_block, &mut self.compressed)?;
        level.entries.push(LeafEntry {
            start_chrom_ix: level.open_chrom_id,
            start_base: level.open_start,
            end_chrom_ix: level.open_chrom_id,
            end_base: level.open_end,
            data_offset: level.stored_len,
            data_size: self.compressed.len() as u64,
        });
        level.stored_len += self.compressed.len() as u64;
        level.largest_block = level.largest_block.max(level.open_block.len());
        // An exactly sized copy per block. Keeping the compression buffer itself would hold its
        // capacity, about the raw block's size and often several times the bytes counted; and
        // appending to one growing buffer would let its capacity overshoot the limit by double.
        level.in_memory.push(self.compressed.clone());
        level.open_block.clear();
        level.open_records = 0;
        self.in_memory += self.compressed.len();
        if self.in_memory > self.memory_limit {
            self.spill_all()?;
        }
        Ok(())
    }

    /// Appends every level's in-memory blocks to the spill file, creating it if needed.
    fn spill_all(&mut self) -> Result<()> {
        let spill = match &mut self.spill {
            Some(spill) => spill,
            None => self.spill.insert(SpillFile::create(&self.spill_location)?),
        };
        for level in self.levels.iter_mut().filter(|level| !level.in_memory.is_empty()) {
            level.spilled.push(spill.append(&level.in_memory)?);
            level.in_memory.clear();
        }
        self.in_memory = 0;
        Ok(())
    }

    /// Number of leading levels to write: the first level is always kept, and the first level
    /// that doesn't have fewer records than the one before is dropped with every level after it.
    fn kept_level_count(&self) -> usize {
        let record_counts: Vec<u64> = self.levels.iter().map(|level| level.record_count).collect();
        kept_level_count(&record_counts)
    }

    /// Appends level `level_index`'s `zoomCount`, its blocks (from the spill file, then from
    /// memory) and its R-tree to `out`, returning the level's header. Every block must be flushed.
    fn write_level<W: Write + Seek>(
        &mut self,
        level_index: usize,
        out: &mut OutputSink<W>,
        index_block_size: u32,
    ) -> Result<ZoomLevel> {
        let level = &self.levels[level_index];
        debug_assert_eq!(level.open_records, 0);
        let data_offset = out.position();
        // Readers reach blocks through the index and ignore `zoomCount`, so a count beyond u32
        // (possible only past 2^32 bases at reduction 1) saturates rather than failing.
        let zoom_count = u32::try_from(level.record_count).unwrap_or(u32::MAX);
        out.write_all(&self.byte_order.u32_bytes(zoom_count))?;
        let blocks_offset = out.position();
        for &run in &level.spilled {
            let spill = self.spill.as_mut().expect("a level with spilled runs has a spill file");
            spill.copy_run(run, out)?;
        }
        for block in &level.in_memory {
            out.write_all(block)?;
        }

        let index_offset = out.position();
        debug_assert_eq!(index_offset, blocks_offset + level.stored_len);
        let entries: Vec<LeafEntry> = level
            .entries
            .iter()
            .map(|entry| LeafEntry { data_offset: blocks_offset + entry.data_offset, ..*entry })
            .collect();
        let index = rtree::encode(
            &entries,
            index_offset,
            index_block_size,
            self.records_per_block,
            self.byte_order,
        );
        out.write_all(&index)?;
        Ok(ZoomLevel::new(level.reduction, data_offset, index_offset))
    }
}

/// One zoom level's blocks: the open one being filled and the compressed ones stored so far.
#[derive(Debug)]
struct ZoomLevelBlocks {
    reduction: u32,
    /// Records emitted so far, the level's `zoomCount`.
    record_count: u64,
    /// Encoded records of the open block.
    open_block: Vec<u8>,
    /// Records in the open block.
    open_records: u32,
    /// Chromosome of the open block's records.
    open_chrom_id: u32,
    /// Start of the open block's first record.
    open_start: u32,
    /// End of the open block's last record.
    open_end: u32,
    /// One leaf entry per stored block, with offsets relative to the level's first block.
    entries: Vec<LeafEntry>,
    /// Total size of the stored blocks, which is the relative offset of the next.
    stored_len: u64,
    /// Stored blocks not yet spilled, in order, following those in `spilled`.
    in_memory: Vec<Vec<u8>>,
    /// Runs of stored blocks in the spill file, in order.
    spilled: Vec<SpillRun>,
    /// Largest decompressed block so far.
    largest_block: usize,
}

impl ZoomLevelBlocks {
    /// Creates a level with no records.
    fn new(reduction: u32) -> Self {
        Self {
            reduction,
            record_count: 0,
            open_block: Vec::new(),
            open_records: 0,
            open_chrom_id: 0,
            open_start: 0,
            open_end: 0,
            entries: Vec::new(),
            stored_len: 0,
            in_memory: Vec::new(),
            spilled: Vec::new(),
            largest_block: 0,
        }
    }

    /// Appends `record` to the open block, which must be empty or on the record's chromosome.
    fn push_record(&mut self, record: &ZoomRecord, byte_order: ByteOrder) {
        if self.open_records == 0 {
            self.open_chrom_id = record.chrom_id;
            self.open_start = record.start;
        }
        record.encode(byte_order, &mut self.open_block);
        self.open_end = record.end;
        self.open_records += 1;
        self.record_count += 1;
    }
}

/// The temporary file zoom blocks spill to, shared by every level.
///
/// On Unix it is unlinked as soon as it is opened, so nothing is left behind even if the
/// process dies; elsewhere it is removed when dropped.
#[derive(Debug)]
struct SpillFile {
    file: File,
    /// Bytes appended so far.
    len: u64,
    /// The path to remove on drop, or `None` once the file is unlinked.
    path: Option<PathBuf>,
}

impl SpillFile {
    /// Creates a new, uniquely named file `.<output>.zoom-<pid>-<n>.tmp` in the spill directory.
    fn create(location: &SpillLocation) -> Result<Self> {
        for _ in 0..SPILL_NAME_ATTEMPTS {
            let number = NEXT_SPILL_FILE.fetch_add(1, Ordering::Relaxed);
            let name =
                format!(".{}.zoom-{}-{number}.tmp", location.output_name, std::process::id());
            let path = location.dir.join(name);
            match OpenOptions::new().read(true).write(true).create_new(true).open(&path) {
                Ok(file) => return Ok(Self::opened(file, path)),
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {}
                Err(error) => return Err(Error::Io(error)),
            }
        }
        Err(Error::Io(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "found no unused zoom spill-file name in {} after {SPILL_NAME_ATTEMPTS} attempts",
                location.dir.display()
            ),
        )))
    }

    /// Wraps a newly created spill file, unlinking it at once where the platform allows.
    fn opened(file: File, path: PathBuf) -> Self {
        // An unlinked file stays usable through its handle on Unix. If unlinking fails, the
        // file is removed on drop instead.
        #[cfg(unix)]
        let path = fs::remove_file(&path).is_err().then_some(path);
        #[cfg(not(unix))]
        let path = Some(path);
        Self { file, len: 0, path }
    }

    /// Appends `blocks` as one run.
    fn append(&mut self, blocks: &[Vec<u8>]) -> Result<SpillRun> {
        let offset = self.len;
        self.file.seek(SeekFrom::Start(offset))?;
        let mut writer = BufWriter::with_capacity(SINK_BUFFER_SIZE, &mut self.file);
        for block in blocks {
            writer.write_all(block)?;
        }
        writer.flush()?;
        let len = blocks.iter().map(|block| block.len() as u64).sum();
        self.len += len;
        Ok(SpillRun { offset, len })
    }

    /// Copies a run of blocks to `out`.
    fn copy_run<W: Write + Seek>(&mut self, run: SpillRun, out: &mut OutputSink<W>) -> Result<()> {
        self.file.seek(SeekFrom::Start(run.offset))?;
        out.copy_from(&mut self.file, run.len)
    }
}

impl Drop for SpillFile {
    fn drop(&mut self) {
        if let Some(path) = self.path.take() {
            if let Err(error) = fs::remove_file(&path) {
                log::warn!("could not remove zoom spill file {}: {error}", path.display());
            }
        }
    }
}

/// A run of one level's blocks in the spill file.
#[derive(Clone, Copy, Debug)]
struct SpillRun {
    /// Offset in the spill file.
    offset: u64,
    /// Length in bytes.
    len: u64,
}

/// Where a spill file goes and what it is named after.
#[derive(Clone, Debug)]
struct SpillLocation {
    /// Directory to create it in.
    dir: PathBuf,
    /// The output file's name, or a stand-in for a sink with no path.
    output_name: String,
}

/// Logs a warning when dropped while armed: the writer was dropped without `finish`.
#[derive(Debug)]
struct UnfinishedWarning {
    /// The output file, when the writer was created from a path.
    output: Option<PathBuf>,
    /// Cleared by `finish`.
    armed: bool,
}

impl Drop for UnfinishedWarning {
    fn drop(&mut self) {
        if !self.armed {
            return;
        }
        match &self.output {
            Some(path) => log::warn!(
                "bigWig writer for {} dropped without finish(): the file is incomplete and its \
                 magic number is 0, so readers will reject it",
                path.display()
            ),
            None => log::warn!(
                "bigWig writer dropped without finish(): the output is incomplete and its magic \
                 number is 0, so readers will reject it"
            ),
        }
    }
}

/// The `chromId` equal to list position `position`; the list was validated to be numbered by u32.
fn list_position_id(position: usize) -> u32 {
    u32::try_from(position)
        .expect("the chromosome list is validated to have at most u32::MAX entries")
}

/// Narrows a position that validation has shown lies within a chromosome, so within u32.
fn within_chrom(position: u64) -> u32 {
    u32::try_from(position).expect("validated items lie within their chromosome")
}

#[cfg(test)]
mod tests {
    use std::cell::{Cell, RefCell};
    use std::io::Cursor;
    use std::rc::Rc;

    use super::*;
    use crate::bbi::{BPTREE_MAGIC, ByteReader, RTREE_MAGIC, decode_total_summary};
    use crate::compress::BlockDecompressor;
    use crate::rtree::RTree;
    use crate::section::{DataItem, Section};

    /// A writer over an in-memory file.
    type TestWriter = BigWigWriter<Cursor<Vec<u8>>>;

    // ---- Writing --------------------------------------------------------------------------

    /// Chromosomes for most tests, deliberately listed out of byte order.
    fn test_chroms() -> Vec<ChromInfo> {
        vec![
            ChromInfo::new("chr2", 50_000),
            ChromInfo::new("chr10", 20_000),
            ChromInfo::new("chr1", 100_000),
            ChromInfo::new("chrM", 16_569),
        ]
    }

    /// A single chromosome of 100 bases, short enough to reach its end.
    fn short_chrom() -> Vec<ChromInfo> {
        vec![ChromInfo::new("chr1", 100)]
    }

    /// A writer with `options` over `chroms`, writing to memory.
    fn new_writer(options: &WriterOptions, chroms: &[ChromInfo]) -> TestWriter {
        options.create_with(Cursor::new(Vec::new()), chroms).unwrap()
    }

    /// Runs `pushes` on a new writer over [`test_chroms`] and returns the finished file.
    fn write_file(options: &WriterOptions, pushes: impl FnOnce(&mut TestWriter)) -> Vec<u8> {
        write_file_over(options, &test_chroms(), pushes)
    }

    /// Runs `pushes` on a new writer over `chroms` and returns the finished file.
    fn write_file_over(
        options: &WriterOptions,
        chroms: &[ChromInfo],
        pushes: impl FnOnce(&mut TestWriter),
    ) -> Vec<u8> {
        let mut writer = new_writer(options, chroms);
        pushes(&mut writer);
        writer.finish().unwrap().into_inner()
    }

    /// `count` varied values, all multiples of 0.25, so every sum in a test is exact.
    fn values(count: usize, seed: usize) -> Vec<f32> {
        (0..count)
            .map(|i| f32::from(u16::try_from((i * 37 + seed * 11) % 400).unwrap()) * 0.25)
            .collect()
    }

    /// `count` values whose sums round in f64, so summation order shows in the bytes.
    fn inexact_values(count: usize, seed: usize) -> Vec<f32> {
        (0..count)
            .map(|i| f32::from(u16::try_from((i * 7_919 + seed) % 10_007).unwrap()) / 3.0)
            .collect()
    }

    /// The items a fixedStep run on `chrom_id` stands for.
    fn fixed_step_items(
        chrom_id: u32,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Vec<DataItem> {
        (0..)
            .zip(values)
            .map(|(i, &value)| DataItem {
                chrom_id,
                start: start + i * step,
                end: start + i * step + span,
                value,
            })
            .collect()
    }

    /// Asserts that `rejected` fails and leaves no trace: running `before`, the rejected call,
    /// then `after` writes exactly the file `before` then `after` write. Returns the error.
    fn rejection_without_trace(
        options: &WriterOptions,
        chroms: &[ChromInfo],
        before: impl Fn(&mut TestWriter),
        rejected: impl FnOnce(&mut TestWriter) -> Result<()>,
        after: impl Fn(&mut TestWriter),
    ) -> Error {
        let mut writer = new_writer(options, chroms);
        before(&mut writer);
        let error = rejected(&mut writer).expect_err("the push should be rejected");
        after(&mut writer);
        let with_rejection = writer.finish().unwrap().into_inner();
        let without_rejection = write_file_over(options, chroms, |writer| {
            before(writer);
            after(writer);
        });
        assert_eq!(with_rejection, without_rejection);
        error
    }

    /// As [`rejection_without_trace`] with default options over [`test_chroms`], where `before`
    /// pushes 100 per-base values to chr1 at 0 and `after` 100 more at 100.
    fn rejection_after_chr1_data(rejected: impl FnOnce(&mut TestWriter) -> Result<()>) -> Error {
        rejection_without_trace(
            &WriterOptions::new(),
            &test_chroms(),
            |writer| writer.push_fixed_step("chr1", 0, 1, 1, &values(100, 1)).unwrap(),
            rejected,
            |writer| writer.push_fixed_step("chr1", 100, 1, 1, &values(100, 2)).unwrap(),
        )
    }

    /// A sink shared with the test, so its bytes can be inspected after the writer is dropped.
    #[derive(Clone, Default)]
    struct SharedSink(Rc<RefCell<Cursor<Vec<u8>>>>);

    impl Write for SharedSink {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.borrow_mut().write(buf)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    impl Seek for SharedSink {
        fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
            self.0.borrow_mut().seek(pos)
        }
    }

    /// A sink whose writes fail while the shared flag is set.
    struct FailingSink {
        inner: Cursor<Vec<u8>>,
        failing: Rc<Cell<bool>>,
    }

    impl Write for FailingSink {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.failing.get() {
                return Err(io::Error::other("disk full"));
            }
            self.inner.write(buf)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    impl Seek for FailingSink {
        fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
            self.inner.seek(pos)
        }
    }

    // ---- Decoding -------------------------------------------------------------------------

    /// A finished file decoded with the crate's own codecs.
    struct DecodedFile {
        byte_order: ByteOrder,
        header: FixedHeader,
        zoom_levels: Vec<ZoomLevel>,
        total_summary: Summary,
        /// Chromosomes ordered by `chromId`.
        chroms: Vec<ChromInfo>,
        section_count: u64,
        data_index: IndexHeader,
        /// Data sections in index order.
        sections: Vec<DecodedSection>,
        zoom: Vec<DecodedZoomLevel>,
    }

    impl DecodedFile {
        /// Every item, in index order.
        fn items(&self) -> Vec<DataItem> {
            self.sections.iter().flat_map(|section| section.items.iter().copied()).collect()
        }

        /// The `chromId` of the chromosome named `name`.
        fn chrom_id(&self, name: &str) -> u32 {
            let position = self.chroms.iter().position(|chrom| chrom.name == name).unwrap();
            u32::try_from(position).unwrap()
        }

        /// The zoom levels' reductions, in header order.
        fn reductions(&self) -> Vec<u32> {
            self.zoom_levels.iter().map(ZoomLevel::reduction).collect()
        }

        /// Item count of each section, in index order.
        fn section_sizes(&self) -> Vec<usize> {
            self.sections.iter().map(|section| section.items.len()).collect()
        }

        /// Largest decompressed data section or zoom block.
        fn largest_block(&self) -> usize {
            let sections = self.sections.iter().map(|section| section.raw_len);
            let blocks =
                self.zoom.iter().flat_map(|level| level.blocks.iter().map(|block| block.raw_len));
            sections.chain(blocks).max().unwrap_or(0)
        }
    }

    /// One decoded data section.
    struct DecodedSection {
        entry: LeafEntry,
        raw_len: usize,
        chrom_id: u32,
        start: u32,
        end: u32,
        kind: SectionKind,
        items: Vec<DataItem>,
    }

    /// One decoded zoom level.
    struct DecodedZoomLevel {
        zoom_count: u32,
        index: IndexHeader,
        blocks: Vec<DecodedZoomBlock>,
    }

    impl DecodedZoomLevel {
        /// Every record, in index order.
        fn records(&self) -> Vec<ZoomRecord> {
            self.blocks.iter().flat_map(|block| block.records.iter().copied()).collect()
        }
    }

    /// One decoded zoom block.
    struct DecodedZoomBlock {
        entry: LeafEntry,
        raw_len: usize,
        records: Vec<ZoomRecord>,
    }

    /// An R-tree's header fields and its root node's header.
    struct IndexHeader {
        block_size: u32,
        item_count: u64,
        end_file_offset: u64,
        items_per_slot: u32,
        root_is_leaf: u8,
        root_count: u16,
    }

    /// The chromosome tree's header fields and root node.
    struct ChromTreeLayout {
        block_size: u32,
        key_size: u32,
        value_size: u32,
        item_count: u64,
        root_is_leaf: u8,
        /// Keys of the root's used slots, NULs trimmed, in slot order.
        root_keys: Vec<String>,
        /// Size of a single-leaf tree in bytes.
        len: u64,
    }

    /// A reader over `bytes` from absolute `offset`.
    fn reader_at(bytes: &[u8], offset: u64, byte_order: ByteOrder) -> ByteReader<'_> {
        ByteReader::new(&bytes[usize::try_from(offset).unwrap()..], byte_order, offset)
    }

    /// The stored bytes of the block an index entry points at.
    fn block_of(bytes: &[u8], entry: LeafEntry) -> &[u8] {
        let start = usize::try_from(entry.data_offset).unwrap();
        &bytes[start..start + usize::try_from(entry.data_size).unwrap()]
    }

    /// Every leaf entry of the R-tree at `offset`.
    fn index_leaves(bytes: &[u8], offset: u64) -> Vec<LeafEntry> {
        let mut source = Cursor::new(bytes);
        let tree = RTree::read(&mut source, offset, bytes.len() as u64).unwrap();
        tree.leaves(&mut source).unwrap()
    }

    /// Reads the header and root node header of the R-tree at `offset`.
    fn index_header(bytes: &[u8], offset: u64, byte_order: ByteOrder) -> IndexHeader {
        let mut reader = reader_at(bytes, offset, byte_order);
        assert_eq!(reader.u32().unwrap(), RTREE_MAGIC);
        let block_size = reader.u32().unwrap();
        let item_count = reader.u64().unwrap();
        reader.skip(16).unwrap();
        let end_file_offset = reader.u64().unwrap();
        let items_per_slot = reader.u32().unwrap();
        reader.skip(4).unwrap();
        let root_is_leaf = reader.u8().unwrap();
        reader.skip(1).unwrap();
        let root_count = reader.u16().unwrap();
        IndexHeader {
            block_size,
            item_count,
            end_file_offset,
            items_per_slot,
            root_is_leaf,
            root_count,
        }
    }

    /// Reads the chromosome tree's header and root node directly.
    fn chrom_tree_layout(bytes: &[u8]) -> ChromTreeLayout {
        let (header, byte_order) = FixedHeader::decode(bytes).unwrap();
        let mut reader = reader_at(bytes, header.chrom_tree_offset, byte_order);
        assert_eq!(reader.u32().unwrap(), BPTREE_MAGIC);
        let block_size = reader.u32().unwrap();
        let key_size = reader.u32().unwrap();
        let value_size = reader.u32().unwrap();
        let item_count = reader.u64().unwrap();
        reader.skip(8).unwrap();
        let root_is_leaf = reader.u8().unwrap();
        reader.skip(1).unwrap();
        let root_count = reader.u16().unwrap();
        let key_len = usize::try_from(key_size).unwrap();
        let root_keys = (0..root_count)
            .map(|_| {
                let key = reader.bytes(key_len).unwrap();
                reader.skip(8).unwrap();
                let name = key.iter().copied().take_while(|&byte| byte != 0).collect();
                String::from_utf8(name).unwrap()
            })
            .collect();
        let len = 32 + 4 + u64::from(block_size) * (u64::from(key_size) + 8);
        ChromTreeLayout {
            block_size,
            key_size,
            value_size,
            item_count,
            root_is_leaf,
            root_keys,
            len,
        }
    }

    /// Decodes a finished file.
    fn decode(bytes: &[u8]) -> DecodedFile {
        let (header, byte_order) = FixedHeader::decode(bytes).unwrap();
        let mut zoom_headers = reader_at(bytes, FIXED_HEADER_SIZE as u64, byte_order);
        let zoom_levels: Vec<ZoomLevel> = (0..header.zoom_levels)
            .map(|_| ZoomLevel::decode(&mut zoom_headers).unwrap())
            .collect();
        let total_summary =
            decode_total_summary(&mut reader_at(bytes, header.total_summary_offset, byte_order))
                .unwrap();
        let file_len = bytes.len() as u64;
        let chroms =
            bptree::read(&mut Cursor::new(bytes), header.chrom_tree_offset, file_len).unwrap();
        let section_count = reader_at(bytes, header.full_data_offset, byte_order).u64().unwrap();

        let mut decompressor = BlockDecompressor::new(header.uncompress_buf_size);
        let sections = index_leaves(bytes, header.full_index_offset)
            .into_iter()
            .map(|entry| {
                let block = block_of(bytes, entry);
                let raw = decompressor.decompress(block, entry.data_offset).unwrap();
                let section = Section::parse(raw, byte_order, entry.data_offset).unwrap();
                DecodedSection {
                    entry,
                    raw_len: raw.len(),
                    chrom_id: section.chrom_id(),
                    start: section.start(),
                    end: section.end_hint(),
                    kind: section.kind(),
                    items: section.items().map(Result::unwrap).collect(),
                }
            })
            .collect();
        let zoom = zoom_levels
            .iter()
            .map(|level| {
                let zoom_count = reader_at(bytes, level.data_offset(), byte_order).u32().unwrap();
                let blocks = index_leaves(bytes, level.index_offset())
                    .into_iter()
                    .map(|entry| {
                        let block = block_of(bytes, entry);
                        let raw = decompressor.decompress(block, entry.data_offset).unwrap();
                        assert_eq!(raw.len() % ZoomRecord::ENCODED_LEN, 0);
                        let mut reader = ByteReader::new(raw, byte_order, entry.data_offset);
                        let records = (0..raw.len() / ZoomRecord::ENCODED_LEN)
                            .map(|_| ZoomRecord::decode(&mut reader).unwrap())
                            .collect();
                        DecodedZoomBlock { entry, raw_len: raw.len(), records }
                    })
                    .collect();
                let index = index_header(bytes, level.index_offset(), byte_order);
                DecodedZoomLevel { zoom_count, index, blocks }
            })
            .collect();
        DecodedFile {
            byte_order,
            header,
            zoom_levels,
            total_summary,
            chroms,
            section_count,
            data_index: index_header(bytes, header.full_index_offset, byte_order),
            sections,
            zoom,
        }
    }

    /// Asserts the structural rules every file the writer produces follows: one chromosome per
    /// section and zoom block, index bounds matching the blocks, file and `chromId` order
    /// agreeing, items and records sorted without overlaps, no record without data, `zoomCount`
    /// equal to the records, and each level smaller than the last.
    fn assert_structure(file: &DecodedFile) {
        assert_eq!(file.section_count, file.sections.len() as u64);
        assert_eq!(file.data_index.item_count, file.sections.len() as u64);
        for pair in file.sections.windows(2) {
            assert!(pair[0].entry.data_offset < pair[1].entry.data_offset);
            assert!(pair[0].chrom_id <= pair[1].chrom_id, "file order is chromId order");
        }
        for section in &file.sections {
            let entry = section.entry;
            assert_eq!(
                (entry.start_chrom_ix, entry.end_chrom_ix),
                (section.chrom_id, section.chrom_id)
            );
            assert_eq!((entry.start_base, entry.end_base), (section.start, section.end));
            assert!(section.items.iter().all(|item| item.chrom_id == section.chrom_id));
            assert!(section.items.iter().all(|item| item.value.is_finite()));
        }
        for pair in file.items().windows(2) {
            if pair[0].chrom_id == pair[1].chrom_id {
                assert!(pair[0].end <= pair[1].start, "items are sorted without overlaps");
            }
        }

        for (level, header) in file.zoom.iter().zip(&file.zoom_levels) {
            let records = level.records();
            assert_eq!(u64::from(level.zoom_count), records.len() as u64);
            assert_eq!(level.index.item_count, level.blocks.len() as u64);
            for block in &level.blocks {
                let chrom_id = block.records[0].chrom_id;
                assert!(block.records.iter().all(|record| record.chrom_id == chrom_id));
                let entry = block.entry;
                assert_eq!((entry.start_chrom_ix, entry.end_chrom_ix), (chrom_id, chrom_id));
                assert_eq!(entry.start_base, block.records[0].start);
                assert_eq!(entry.end_base, block.records[block.records.len() - 1].end);
            }
            for record in &records {
                assert!(record.summary.covered_bases >= 1.0, "no record without data");
                assert!(record.start < record.end);
                assert!(record.summary.covered_bases <= f64::from(record.end - record.start));
                assert!(record.end - record.start <= header.reduction());
            }
            for pair in records.windows(2) {
                assert!((pair[0].chrom_id, pair[0].end) <= (pair[1].chrom_id, pair[1].start));
            }
        }
        let reductions = file.reductions();
        assert!(reductions.windows(2).all(|pair| pair[0] < pair[1]), "{reductions:?}");
        let counts: Vec<u32> = file.zoom.iter().map(|level| level.zoom_count).collect();
        assert!(counts.windows(2).all(|pair| pair[1] < pair[0]), "{counts:?}");
    }

    /// A [`Summary`] of nothing, to accumulate into.
    fn empty_summary() -> Summary {
        Summary {
            covered_bases: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            sum: 0.0,
            sum_squares: 0.0,
        }
    }

    /// The zoom records a level with `reduction` must hold, computed base by base from `items`
    /// over chromosomes of `lengths` (by `chromId`), with values rounded to f32 as stored.
    fn expected_records(items: &[DataItem], lengths: &[u32], reduction: u32) -> Vec<ZoomRecord> {
        let bin_len = usize::try_from(reduction).unwrap();
        let mut records = Vec::new();
        for (chrom_id, &length) in (0..).zip(lengths) {
            let mut bases = vec![None; usize::try_from(length).unwrap()];
            for item in items.iter().filter(|item| item.chrom_id == chrom_id) {
                for base in item.start..item.end {
                    bases[usize::try_from(base).unwrap()] = Some(f64::from(item.value));
                }
            }
            for (bin_start, bin) in (0..).step_by(bin_len).zip(bases.chunks(bin_len)) {
                let covered: Vec<(u32, f64)> = (bin_start..)
                    .zip(bin)
                    .filter_map(|(base, value)| value.map(|value| (base, value)))
                    .collect();
                let (Some(first), Some(last)) = (covered.first(), covered.last()) else {
                    continue;
                };
                let mut summary = empty_summary();
                for &(_, value) in &covered {
                    summary.covered_bases += 1.0;
                    summary.min = summary.min.min(value);
                    summary.max = summary.max.max(value);
                    summary.sum += value;
                    summary.sum_squares += value * value;
                }
                let record = ZoomRecord { chrom_id, start: first.0, end: last.0 + 1, summary };
                let mut encoded = Vec::new();
                record.encode(ByteOrder::Little, &mut encoded);
                let mut reader = ByteReader::new(&encoded, ByteOrder::Little, 0);
                records.push(ZoomRecord::decode(&mut reader).unwrap());
            }
        }
        records
    }

    /// The exact summary of `items`.
    fn summary_of(items: &[DataItem]) -> Summary {
        let mut summary = empty_summary();
        for item in items {
            let (bases, value) = (f64::from(item.end - item.start), f64::from(item.value));
            summary.covered_bases += bases;
            summary.min = summary.min.min(value);
            summary.max = summary.max.max(value);
            summary.sum += bases * value;
            summary.sum_squares += bases * value * value;
        }
        summary
    }

    /// Pushes a mix of all three section types to three chromosomes.
    fn push_mixed_data<W: Write + Seek>(writer: &mut BigWigWriter<W>) {
        writer.push_fixed_step("chr1", 1_000, 1, 1, &values(20_000, 1)).unwrap();
        writer.push_fixed_step("chr1", 30_000, 10, 5, &values(500, 2)).unwrap();
        writer.push_intervals("chr1", &[(40_000, 40_100, 1.5), (40_200, 40_250, -2.0)]).unwrap();
        writer.push_variable_step("chrM", 25, &[(100, 3.0), (200, 4.0), (16_000, 5.0)]).unwrap();
        let interval_values = values(2_000, 3);
        let intervals: Vec<(u32, u32, f32)> =
            (0..).zip(interval_values).map(|(i, value)| (i * 20, i * 20 + 15, value)).collect();
        writer.push_intervals("chr2", &intervals).unwrap();
    }

    // ---- Options --------------------------------------------------------------------------

    #[test]
    fn default_options_match_the_documented_defaults() {
        let options = WriterOptions::default();
        assert_eq!(options.items_per_section, 8_192);
        assert_eq!(options.zoom, ZoomPolicy::Auto { min_first_reduction: 256 });
        assert_eq!(options.zoom_items_per_block, 1_024);
        assert_eq!(options.index_block_size, 256);
        assert_eq!(options.compression, Compression::Deflate { level: 1 });
        assert_eq!(options.zoom_memory_limit, 64 << 20);
        assert_eq!(options.temp_dir, None);
        assert_eq!(options.chrom_id_order, ChromIdOrder::Arrival);
        assert!(!options.clip);
        assert_eq!(options.byte_order, ByteOrder::Little);
    }

    /// The message of the `InvalidOption` error from creating a writer with `options`.
    fn invalid_option_message(options: &WriterOptions) -> String {
        match options.create_with(Cursor::new(Vec::new()), &test_chroms()) {
            Err(Error::InvalidOption(message)) => message,
            Err(other) => panic!("expected InvalidOption, got {other:?}"),
            Ok(_) => panic!("expected InvalidOption, got a writer"),
        }
    }

    #[test]
    fn zero_items_per_section_is_an_invalid_option() {
        let message = invalid_option_message(&WriterOptions::new().items_per_section(0));
        assert!(message.contains("items_per_section"), "{message}");
    }

    #[test]
    fn zero_zoom_items_per_block_is_an_invalid_option() {
        let message = invalid_option_message(&WriterOptions::new().zoom_items_per_block(0));
        assert!(message.contains("zoom_items_per_block"), "{message}");
    }

    #[test]
    fn zoom_blocks_larger_than_the_inflate_cap_are_an_invalid_option() {
        let message = invalid_option_message(&WriterOptions::new().zoom_items_per_block(2_097_153));
        assert!(message.contains("zoom_items_per_block"), "{message}");
        let at_the_cap = WriterOptions::new().zoom_items_per_block(2_097_152);
        assert!(at_the_cap.create_with(Cursor::new(Vec::new()), &test_chroms()).is_ok());
    }

    #[test]
    fn index_block_size_below_two_is_an_invalid_option() {
        let message = invalid_option_message(&WriterOptions::new().index_block_size(1));
        assert!(message.contains("index_block_size"), "{message}");
    }

    #[test]
    fn index_block_size_above_65535_is_an_invalid_option() {
        let message = invalid_option_message(&WriterOptions::new().index_block_size(65_536));
        assert!(message.contains("index_block_size"), "{message}");
    }

    #[test]
    fn index_block_size_of_65535_is_accepted() {
        let bytes = write_file(&WriterOptions::new().index_block_size(65_535), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        });
        assert_eq!(decode(&bytes).data_index.block_size, 65_535);
    }

    #[test]
    fn unsupported_compression_level_is_an_invalid_option() {
        let level = Compression::MAX_LEVEL + 1;
        let options = WriterOptions::new().compression(Compression::Deflate { level });
        assert!(invalid_option_message(&options).contains("compression"));
    }

    #[test]
    fn invalid_zoom_policy_is_an_invalid_option() {
        let options = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![]));
        assert!(invalid_option_message(&options).contains("zoom"));
    }

    #[test]
    fn invalid_options_create_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bw");
        let result = WriterOptions::new().index_block_size(0).create(&path, &test_chroms());
        assert!(matches!(result, Err(Error::InvalidOption(_))));
        assert!(!path.exists());
    }

    // ---- Chromosome list ------------------------------------------------------------------

    /// The message of the `InvalidChromList` error from creating a writer over `chroms`.
    fn chrom_list_message(chroms: &[ChromInfo]) -> String {
        match WriterOptions::new().create_with(Cursor::new(Vec::new()), chroms) {
            Err(Error::InvalidChromList(message)) => message,
            Err(other) => panic!("expected InvalidChromList, got {other:?}"),
            Ok(_) => panic!("expected InvalidChromList, got a writer"),
        }
    }

    #[test]
    fn empty_chrom_list_is_invalid() {
        assert!(chrom_list_message(&[]).contains("empty"));
    }

    #[test]
    fn empty_chrom_name_is_invalid() {
        let message = chrom_list_message(&[ChromInfo::new("chr1", 10), ChromInfo::new("", 10)]);
        assert!(message.contains("position 1"), "{message}");
    }

    #[test]
    fn chrom_name_with_a_space_is_invalid() {
        assert!(chrom_list_message(&[ChromInfo::new("chr 1", 10)]).contains("0x20"));
    }

    #[test]
    fn chrom_name_with_a_nul_is_invalid() {
        assert!(chrom_list_message(&[ChromInfo::new("chr\0x", 10)]).contains("0x00"));
    }

    #[test]
    fn chrom_name_beyond_ascii_is_invalid() {
        assert!(chrom_list_message(&[ChromInfo::new("chré", 10)]).contains("printable ASCII"));
    }

    #[test]
    fn duplicate_chrom_name_is_invalid() {
        let chroms =
            [ChromInfo::new("chr1", 10), ChromInfo::new("chr2", 10), ChromInfo::new("chr1", 5)];
        assert!(chrom_list_message(&chroms).contains("more than once"));
    }

    #[test]
    fn zero_length_chrom_is_invalid() {
        assert!(chrom_list_message(&[ChromInfo::new("chr1", 0)]).contains("length 0"));
    }

    #[test]
    fn invalid_chrom_list_creates_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bw");
        let result = WriterOptions::new().create(&path, &[]);
        assert!(matches!(result, Err(Error::InvalidChromList(_))));
        assert!(!path.exists());
    }

    #[test]
    fn chrom_longer_than_i32_max_is_accepted() {
        let chroms = [ChromInfo::new("huge", u32::MAX)];
        let bytes = write_file_over(&WriterOptions::new(), &chroms, |writer| {
            writer.push_interval("huge", u32::MAX - 10, u32::MAX, 2.0).unwrap();
        });
        let file = decode(&bytes);
        assert_eq!(file.chroms, chroms);
        let item = DataItem { chrom_id: 0, start: u32::MAX - 10, end: u32::MAX, value: 2.0 };
        assert_eq!(file.items(), [item]);
    }

    // ---- Layout: region order, header slots, total summary ---------------------------------

    #[test]
    fn regions_are_in_kents_order() {
        let bytes = write_file(&WriterOptions::new(), push_mixed_data);
        let file = decode(&bytes);
        let header = file.header;
        assert!(FIXED_HEADER_SIZE as u64 <= header.total_summary_offset);
        assert!(header.total_summary_offset < header.chrom_tree_offset);
        assert!(header.chrom_tree_offset < header.full_data_offset);
        assert!(header.full_data_offset < file.sections[0].entry.data_offset);
        let last_section = file.sections.last().unwrap().entry;
        assert!(last_section.data_offset + last_section.data_size <= header.full_index_offset);
        assert!(!file.zoom_levels.is_empty());
        let mut previous_end = header.full_index_offset;
        for level in &file.zoom_levels {
            assert!(previous_end < level.data_offset());
            assert!(level.data_offset() < level.index_offset());
            previous_end = level.index_offset();
        }
        assert_eq!(bytes[bytes.len() - 4..], BIGWIG_MAGIC.to_le_bytes());
        assert_structure(&file);
    }

    #[test]
    fn total_summary_immediately_follows_ten_reserved_zoom_header_slots() {
        let bytes = write_file(&WriterOptions::new(), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        });
        let header = decode(&bytes).header;
        assert_eq!(header.total_summary_offset, 64 + 10 * 24);
        assert_eq!(header.chrom_tree_offset, 64 + 10 * 24 + 40);
    }

    #[test]
    fn explicit_ladder_longer_than_ten_levels_reserves_a_slot_per_level() {
        let options = WriterOptions::new().zoom(ZoomPolicy::Explicit((1..=12).collect()));
        let bytes = write_file(&options, |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        });
        assert_eq!(decode(&bytes).header.total_summary_offset, 64 + 12 * 24);
    }

    #[test]
    fn zoom_policy_none_writes_no_zoom_levels_and_reserves_no_slots() {
        let options = WriterOptions::new().zoom(ZoomPolicy::None);
        let file = decode(&write_file(&options, push_mixed_data));
        assert_eq!(file.header.zoom_levels, 0);
        assert_eq!(file.header.total_summary_offset, 64);
        assert_eq!(file.total_summary, summary_of(&file.items()));
        assert_structure(&file);
    }

    #[test]
    fn unused_zoom_header_slots_are_zero() {
        let bytes = write_file(&WriterOptions::new(), push_mixed_data);
        let header = decode(&bytes).header;
        let used = FIXED_HEADER_SIZE + usize::from(header.zoom_levels) * ZOOM_HEADER_SIZE;
        assert!(used < 64 + 10 * 24);
        let total_summary_offset = usize::try_from(header.total_summary_offset).unwrap();
        assert!(bytes[used..total_summary_offset].iter().all(|&byte| byte == 0));
    }

    #[test]
    fn total_summary_is_exact_over_the_raw_items() {
        let file = decode(&write_file(&WriterOptions::new(), push_mixed_data));
        assert_ne!(file.header.total_summary_offset, 0);
        assert_eq!(file.total_summary, summary_of(&file.items()));
    }

    #[test]
    fn total_summary_max_of_a_single_value_is_that_value() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_interval("chr1", 0, 10, 5.0).unwrap();
        }));
        assert_eq!((file.total_summary.min, file.total_summary.max), (5.0, 5.0));
    }

    #[test]
    fn total_summary_max_of_negative_values_is_negative() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &[-1.0, -2.0]).unwrap();
        }));
        assert_eq!((file.total_summary.min, file.total_summary.max), (-2.0, -1.0));
    }

    // ---- Zoom levels ----------------------------------------------------------------------

    #[test]
    fn auto_ladder_for_per_base_data_grows_by_four_from_256() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(100_000, 1)).unwrap();
        }));
        assert_eq!(file.reductions(), [256, 1_024, 4_096, 16_384, 65_536]);
        assert_structure(&file);
    }

    #[test]
    fn auto_ladder_for_binned_data_starts_at_ten_bins() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 100, 100, &values(1_000, 1)).unwrap();
        }));
        assert_eq!(file.reductions(), [1_000, 4_000, 16_000, 64_000]);
    }

    #[test]
    fn auto_ladder_measures_pitch_across_chromosomes() {
        // Ten 1 kb items every 2 kb on each of two chromosomes: each chromosome's extent is
        // 19 kb, and the space between chromosomes doesn't count, so the pitch is 1.9 kb and
        // the first level 19 kb.
        let items: Vec<(u32, u32, f32)> =
            (0..10).map(|i| (i * 2_000, i * 2_000 + 1_000, 1.0)).collect();
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &items).unwrap();
            writer.push_intervals("chr2", &items).unwrap();
        }));
        assert_eq!(file.reductions()[0], 19_000);
    }

    #[test]
    fn zoom_levels_longer_than_the_longest_chromosome_are_not_written() {
        let options = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![100, 1_000, 200_000]));
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(100_000, 1)).unwrap();
        }));
        assert_eq!(file.reductions(), [100, 1_000]);
    }

    /// Single-base items every 2 kb on chr1: every level up to 1 kb has a record per item.
    fn push_sparse_sites(writer: &mut TestWriter) {
        let sites: Vec<(u32, u32, f32)> =
            (0..50).map(|i| (i * 2_000, i * 2_000 + 1, 1.0)).collect();
        writer.push_intervals("chr1", &sites).unwrap();
    }

    #[test]
    fn zoom_levels_that_do_not_shrink_are_dropped() {
        let options = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10, 100, 1_000]));
        let file = decode(&write_file(&options, push_sparse_sites));
        assert_eq!(file.reductions(), [10]);
        assert_eq!(file.zoom[0].zoom_count, 50);
    }

    #[test]
    fn dropped_zoom_levels_leave_no_bytes_in_the_file() {
        let ladder = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10, 100, 1_000]));
        let single = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10]));
        assert_eq!(write_file(&ladder, push_sparse_sites), write_file(&single, push_sparse_sites));
    }

    #[test]
    fn zoom_records_match_a_base_by_base_model() {
        // 10 is not a multiple of 4, so it is accumulated from the raw items; 16 and 64 are
        // merged from finer levels.
        let options = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![4, 10, 16, 64]));
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr10", 3, 1, 1, &values(700, 1)).unwrap();
            writer.push_fixed_step("chr10", 1_000, 7, 3, &values(300, 2)).unwrap();
            writer.push_intervals("chr10", &[(4_000, 4_050, 2.5), (4_053, 4_054, -1.0)]).unwrap();
            writer.push_variable_step("chrM", 9, &[(5, 1.0), (17, 2.0), (40, 3.0)]).unwrap();
        }));
        assert_eq!(file.reductions(), [4, 10, 16, 64]);
        let lengths: Vec<u32> = file.chroms.iter().map(|chrom| chrom.length).collect();
        for (level, reduction) in file.zoom.iter().zip(file.reductions()) {
            assert_eq!(level.records(), expected_records(&file.items(), &lengths, reduction));
        }
        assert_structure(&file);
    }

    #[test]
    fn each_zoom_block_holds_one_chromosome() {
        let options =
            WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10])).zoom_items_per_block(100_000);
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(500, 1)).unwrap();
            writer.push_fixed_step("chr2", 0, 1, 1, &values(500, 2)).unwrap();
            writer.push_fixed_step("chrM", 0, 1, 1, &values(500, 3)).unwrap();
        }));
        assert_eq!(file.zoom[0].blocks.len(), 3);
        assert_structure(&file);
    }

    #[test]
    fn zoom_blocks_hold_zoom_items_per_block_records() {
        let options =
            WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10])).zoom_items_per_block(7);
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(500, 1)).unwrap();
        }));
        let sizes: Vec<usize> =
            file.zoom[0].blocks.iter().map(|block| block.records.len()).collect();
        assert_eq!(sizes, [7, 7, 7, 7, 7, 7, 7, 1]);
        assert_eq!(file.zoom[0].index.items_per_slot, 7);
    }

    #[test]
    fn index_items_per_slot_is_one_for_data_and_the_block_size_for_zoom() {
        let file = decode(&write_file(&WriterOptions::new(), push_mixed_data));
        assert_eq!(file.data_index.items_per_slot, 1);
        assert_eq!(file.zoom[0].index.items_per_slot, 1_024);
    }

    #[test]
    fn index_block_size_sets_the_fan_out_of_every_index() {
        let options = WriterOptions::new()
            .items_per_section(10)
            .zoom(ZoomPolicy::Explicit(vec![10]))
            .zoom_items_per_block(10)
            .index_block_size(4);
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(1_000, 1)).unwrap();
        }));
        assert_eq!((file.data_index.block_size, file.data_index.root_is_leaf), (4, 0));
        assert!(file.data_index.root_count <= 4);
        assert_eq!(file.sections.len(), 100);
        assert_eq!((file.zoom[0].index.block_size, file.zoom[0].index.root_is_leaf), (4, 0));
        assert_eq!(file.zoom[0].blocks.len(), 10);
        assert_structure(&file);
    }

    #[test]
    fn index_end_file_offset_is_the_index_offset() {
        let file = decode(&write_file(&WriterOptions::new(), push_mixed_data));
        assert_eq!(file.data_index.end_file_offset, file.header.full_index_offset);
        for (level, header) in file.zoom.iter().zip(&file.zoom_levels) {
            assert_eq!(level.index.end_file_offset, header.index_offset());
        }
    }

    // ---- Chromosome tree and IDs -----------------------------------------------------------

    #[test]
    fn chrom_tree_keys_are_byte_sorted_in_a_single_leaf() {
        let layout = chrom_tree_layout(&write_file(&WriterOptions::new(), push_mixed_data));
        assert_eq!(layout.root_is_leaf, 1);
        assert_eq!((layout.block_size, layout.item_count), (4, 4));
        assert_eq!(layout.root_keys, ["chr1", "chr10", "chr2", "chrM"]);
    }

    #[test]
    fn chrom_tree_values_are_eight_bytes_and_keys_fit_the_longest_name() {
        let layout = chrom_tree_layout(&write_file(&WriterOptions::new(), push_mixed_data));
        assert_eq!((layout.value_size, layout.key_size), (8, 5));
    }

    #[test]
    fn three_thousand_chromosomes_fit_a_single_leaf() {
        let chroms: Vec<ChromInfo> =
            (1..=3_000).map(|n| ChromInfo::new(format!("contig{n}"), 1_000)).collect();
        let bytes = write_file_over(&WriterOptions::new(), &chroms, |writer| {
            writer.push_interval("contig2999", 0, 10, 1.0).unwrap();
        });
        let layout = chrom_tree_layout(&bytes);
        assert_eq!((layout.root_is_leaf, layout.block_size), (1, 3_000));
        let file = decode(&bytes);
        assert_eq!(file.chroms.len(), 3_000);
        assert_eq!(file.chroms[0].name, "contig2999");
    }

    #[test]
    fn every_listed_chromosome_is_in_the_tree() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        }));
        let mut written = file.chroms.clone();
        let mut listed = test_chroms();
        written.sort_by(|a, b| a.name.cmp(&b.name));
        listed.sort_by(|a, b| a.name.cmp(&b.name));
        assert_eq!(written, listed);
    }

    #[test]
    fn arrival_ids_follow_first_data_then_list_order() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
            writer.push_interval("chrM", 0, 10, 1.0).unwrap();
        }));
        let names: Vec<&str> = file.chroms.iter().map(|chrom| chrom.name.as_str()).collect();
        // chr2 and chr10 have no data, so they take the remaining IDs in list order.
        assert_eq!(names, ["chr1", "chrM", "chr2", "chr10"]);
        let section_chroms: Vec<u32> =
            file.sections.iter().map(|section| section.chrom_id).collect();
        assert_eq!(section_chroms, [0, 1]);
    }

    #[test]
    fn chromosomes_may_be_written_in_any_order() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chrM", 10, 1, 1, &values(5, 1)).unwrap();
            writer.push_intervals("chr10", &[(5, 15, 2.0)]).unwrap();
            writer.push_fixed_step("chr1", 0, 1, 1, &values(5, 2)).unwrap();
        }));
        let mut expected = fixed_step_items(file.chrom_id("chrM"), 10, 1, 1, &values(5, 1));
        expected.push(DataItem { chrom_id: file.chrom_id("chr10"), start: 5, end: 15, value: 2.0 });
        expected.extend(fixed_step_items(file.chrom_id("chr1"), 0, 1, 1, &values(5, 2)));
        assert_eq!(file.items(), expected);
        assert_structure(&file);
    }

    #[test]
    fn listed_ids_are_list_positions_around_chromosomes_without_data() {
        let options = WriterOptions::new().chrom_id_order(ChromIdOrder::Listed);
        let file = decode(&write_file(&options, |writer| {
            writer.push_interval("chr10", 0, 10, 1.0).unwrap();
            writer.push_interval("chrM", 0, 10, 1.0).unwrap();
        }));
        assert_eq!(file.chroms, test_chroms());
        let section_chroms: Vec<u32> =
            file.sections.iter().map(|section| section.chrom_id).collect();
        assert_eq!(section_chroms, [1, 3]);
        assert_structure(&file);
    }

    #[test]
    fn listed_order_rejects_a_chromosome_earlier_in_the_list_without_trace() {
        let options = WriterOptions::new().chrom_id_order(ChromIdOrder::Listed);
        let error = rejection_without_trace(
            &options,
            &test_chroms(),
            |writer| writer.push_interval("chr1", 0, 10, 1.0).unwrap(),
            |writer| writer.push_interval("chr10", 0, 10, 1.0),
            |writer| writer.push_interval("chrM", 0, 10, 1.0).unwrap(),
        );
        assert!(matches!(
            error,
            Error::ChromOutOfOrder { chrom, current } if chrom == "chr10" && current == "chr1"
        ));
    }

    #[test]
    fn chromosomes_without_data_appear_only_in_the_tree() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr2", 0, 1, 1, &values(10, 1)).unwrap();
        }));
        assert_eq!(file.chroms.len(), 4);
        let chr2 = file.chrom_id("chr2");
        assert!(file.sections.iter().all(|section| section.chrom_id == chr2));
        assert!(file.zoom[0].records().iter().all(|record| record.chrom_id == chr2));
    }

    // ---- Sections -------------------------------------------------------------------------

    #[test]
    fn items_of_every_section_type_round_trip() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 100, 10, 4, &values(30, 1)).unwrap();
            writer.push_variable_step("chr1", 3, &[(500, 1.0), (510, -1.5)]).unwrap();
            writer.push_intervals("chr1", &[(600, 650, 0.25), (650, 700, 8.0)]).unwrap();
        }));
        let chrom_id = file.chrom_id("chr1");
        let mut expected = fixed_step_items(chrom_id, 100, 10, 4, &values(30, 1));
        expected.extend([
            DataItem { chrom_id, start: 500, end: 503, value: 1.0 },
            DataItem { chrom_id, start: 510, end: 513, value: -1.5 },
            DataItem { chrom_id, start: 600, end: 650, value: 0.25 },
            DataItem { chrom_id, start: 650, end: 700, value: 8.0 },
        ]);
        assert_eq!(file.items(), expected);
        let kinds: Vec<SectionKind> = file.sections.iter().map(|section| section.kind).collect();
        let expected_kinds = [
            SectionKind::FixedStep { step: 10, span: 4 },
            SectionKind::VariableStep { span: 3 },
            SectionKind::BedGraph,
        ];
        assert_eq!(kinds, expected_kinds);
    }

    #[test]
    fn section_ends_follow_kents_formulas() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 100, 10, 5, &[1.0, 2.0, 3.0]).unwrap();
            writer.push_variable_step("chr1", 7, &[(200, 1.0), (300, 2.0)]).unwrap();
            writer.push_intervals("chr1", &[(400, 410, 1.0), (420, 490, 2.0)]).unwrap();
        }));
        let extents: Vec<(u32, u32)> =
            file.sections.iter().map(|section| (section.start, section.end)).collect();
        assert_eq!(extents, [(100, 125), (200, 307), (400, 490)]);
    }

    #[test]
    fn fixed_step_sections_split_at_items_per_section() {
        let options = WriterOptions::new().items_per_section(4);
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(10, 1)).unwrap();
        }));
        assert_eq!(file.section_sizes(), [4, 4, 2]);
        let starts: Vec<u32> = file.sections.iter().map(|section| section.start).collect();
        assert_eq!(starts, [0, 4, 8]);
        assert_eq!(file.items(), fixed_step_items(file.chrom_id("chr1"), 0, 1, 1, &values(10, 1)));
    }

    #[test]
    fn exactly_items_per_section_values_fill_one_section() {
        let options = WriterOptions::new().items_per_section(4);
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(4, 1)).unwrap();
        }));
        assert_eq!(file.section_sizes(), [4]);
    }

    #[test]
    fn items_per_section_of_one_writes_a_section_per_item() {
        let options = WriterOptions::new().items_per_section(1);
        let file = decode(&write_file(&options, |writer| {
            writer.push_intervals("chr1", &[(0, 5, 1.0), (5, 10, 2.0), (20, 30, 3.0)]).unwrap();
        }));
        assert_eq!(file.section_sizes(), [1, 1, 1]);
    }

    #[test]
    fn items_per_section_of_65535_holds_65535_values() {
        let options = WriterOptions::new().items_per_section(u16::MAX);
        let file = decode(&write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(65_536, 1)).unwrap();
        }));
        assert_eq!(file.section_sizes(), [65_535, 1]);
    }

    #[test]
    fn bed_graph_sections_split_at_items_per_section() {
        let options = WriterOptions::new().items_per_section(3);
        let intervals: Vec<(u32, u32, f32)> = (0..7).map(|i| (i * 10, i * 10 + 5, 1.0)).collect();
        let file = decode(&write_file(&options, |writer| {
            writer.push_intervals("chr1", &intervals).unwrap();
        }));
        assert_eq!(file.section_sizes(), [3, 3, 1]);
    }

    #[test]
    fn variable_step_sections_split_at_items_per_section() {
        let options = WriterOptions::new().items_per_section(2);
        let items: Vec<(u32, f32)> = (0..5).map(|i| (i * 10, 1.0)).collect();
        let file = decode(&write_file(&options, |writer| {
            writer.push_variable_step("chr1", 5, &items).unwrap();
        }));
        assert_eq!(file.section_sizes(), [2, 2, 1]);
    }

    #[test]
    fn continuing_fixed_step_pushes_extend_the_open_section() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(100, 1)).unwrap();
            writer.push_fixed_step("chr1", 100, 1, 1, &values(100, 2)).unwrap();
        }));
        assert_eq!(file.section_sizes(), [200]);
    }

    #[test]
    fn continuing_fixed_step_pushes_fill_sections_as_one_push_would() {
        let options = WriterOptions::new().items_per_section(150);
        let chunked = write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 2, 1, &values(100, 1)).unwrap();
            writer.push_fixed_step("chr1", 200, 2, 1, &values(100, 2)).unwrap();
        });
        let mut all = values(100, 1);
        all.extend(values(100, 2));
        let whole = write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 2, 1, &all).unwrap();
        });
        assert_eq!(chunked, whole);
        assert_eq!(decode(&whole).section_sizes(), [150, 50]);
    }

    #[test]
    fn fixed_step_gap_starts_a_new_section() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(100, 1)).unwrap();
            writer.push_fixed_step("chr1", 101, 1, 1, &values(100, 2)).unwrap();
        }));
        assert_eq!(file.section_sizes(), [100, 100]);
    }

    #[test]
    fn changed_step_or_span_starts_a_new_section() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(10, 1)).unwrap();
            writer.push_fixed_step("chr1", 10, 2, 1, &values(10, 2)).unwrap();
            writer.push_fixed_step("chr1", 30, 2, 2, &values(10, 3)).unwrap();
        }));
        assert_eq!(file.section_sizes(), [10, 10, 10]);
    }

    #[test]
    fn changed_section_type_starts_a_new_section() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &[(0, 5, 1.0)]).unwrap();
            writer.push_variable_step("chr1", 5, &[(10, 1.0)]).unwrap();
            writer.push_intervals("chr1", &[(20, 25, 1.0)]).unwrap();
        }));
        assert_eq!(file.section_sizes(), [1, 1, 1]);
    }

    #[test]
    fn chromosome_change_starts_a_new_section() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &[(0, 5, 1.0)]).unwrap();
            writer.push_intervals("chr2", &[(5, 10, 1.0)]).unwrap();
        }));
        assert_eq!(file.section_sizes(), [1, 1]);
        assert_structure(&file);
    }

    #[test]
    fn single_interval_pushes_share_a_section() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            for i in 0..100 {
                writer.push_interval("chr1", i * 10, i * 10 + 5, 1.0).unwrap();
            }
        }));
        assert_eq!(file.section_sizes(), [100]);
    }

    // ---- Compression and uncompressBufSize --------------------------------------------------

    #[test]
    fn uncompress_buf_size_is_the_largest_decompressed_block() {
        let file = decode(&write_file(&WriterOptions::new(), push_mixed_data));
        assert_eq!(file.header.uncompress_buf_size as usize, file.largest_block());
    }

    #[test]
    fn uncompress_buf_size_can_come_from_a_zoom_block() {
        let options = WriterOptions::new().items_per_section(1).zoom(ZoomPolicy::Explicit(vec![1]));
        let intervals: Vec<(u32, u32, f32)> = (0..100).map(|i| (i * 2, i * 2 + 1, 1.0)).collect();
        let file = decode(&write_file(&options, |writer| {
            writer.push_intervals("chr1", &intervals).unwrap();
        }));
        assert_eq!(file.header.uncompress_buf_size, 100 * 32);
        assert_eq!(file.header.uncompress_buf_size as usize, file.largest_block());
    }

    #[test]
    fn uncompress_buf_size_ignores_dropped_zoom_levels() {
        // Bases 14 and 15 share a 10-base bin but not a 15-base one, so the 15-base level has
        // more records than the 10-base level and is dropped, with its 64-byte block.
        let options =
            WriterOptions::new().items_per_section(1).zoom(ZoomPolicy::Explicit(vec![10, 15]));
        let file = decode(&write_file(&options, |writer| {
            writer.push_intervals("chr1", &[(14, 15, 1.0), (15, 16, 2.0)]).unwrap();
        }));
        assert_eq!(file.reductions(), [10]);
        assert_eq!(file.header.uncompress_buf_size, 24 + 12);
    }

    #[test]
    fn blocks_are_zlib_streams() {
        let bytes = write_file(&WriterOptions::new(), push_mixed_data);
        let file = decode(&bytes);
        let section_entries = file.sections.iter().map(|section| section.entry);
        let zoom_entries =
            file.zoom.iter().flat_map(|level| level.blocks.iter().map(|block| block.entry));
        for entry in section_entries.chain(zoom_entries) {
            let block = block_of(&bytes, entry);
            assert_eq!(block[0], 0x78, "deflate with a 32 KiB window");
            assert_eq!((u16::from(block[0]) * 256 + u16::from(block[1])) % 31, 0);
        }
    }

    #[test]
    fn compression_none_stores_raw_blocks_with_zero_uncompress_buf_size() {
        let options = WriterOptions::new().compression(Compression::None);
        let file = decode(&write_file(&options, push_mixed_data));
        assert_eq!(file.header.uncompress_buf_size, 0);
        for section in &file.sections {
            assert_eq!(section.entry.data_size, section.raw_len as u64);
        }
        assert_eq!(file.total_summary, summary_of(&file.items()));
    }

    // ---- Byte order -----------------------------------------------------------------------

    #[test]
    fn files_are_little_endian_by_default() {
        let bytes = write_file(&WriterOptions::new(), push_mixed_data);
        assert_eq!(bytes[..4], BIGWIG_MAGIC.to_le_bytes());
        assert_eq!(decode(&bytes).byte_order, ByteOrder::Little);
    }

    #[test]
    fn big_endian_file_holds_the_same_content_as_little_endian() {
        let little = decode(&write_file(&WriterOptions::new(), push_mixed_data));
        let bytes = write_file(&WriterOptions::new().byte_order(ByteOrder::Big), push_mixed_data);
        assert_eq!(bytes[..4], BIGWIG_MAGIC.to_be_bytes());
        assert_eq!(bytes[bytes.len() - 4..], BIGWIG_MAGIC.to_be_bytes());
        let big = decode(&bytes);
        assert_eq!(big.byte_order, ByteOrder::Big);
        assert_eq!(big.chroms, little.chroms);
        assert_eq!(big.items(), little.items());
        assert_eq!(big.total_summary, little.total_summary);
        assert_eq!(big.reductions(), little.reductions());
        for (big_level, little_level) in big.zoom.iter().zip(&little.zoom) {
            assert_eq!(big_level.records(), little_level.records());
        }
        assert_structure(&big);
    }

    // ---- Validation -----------------------------------------------------------------------

    #[test]
    fn unknown_chromosome_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| writer.push_interval("chrX", 0, 10, 1.0));
        assert!(matches!(error, Error::UnknownChrom(name) if name == "chrX"));
    }

    #[test]
    fn data_for_a_finished_chromosome_is_rejected_without_trace() {
        let error = rejection_without_trace(
            &WriterOptions::new(),
            &test_chroms(),
            |writer| {
                writer.push_interval("chr1", 0, 10, 1.0).unwrap();
                writer.push_interval("chr2", 0, 10, 1.0).unwrap();
            },
            |writer| writer.push_interval("chr1", 100, 110, 1.0),
            |writer| writer.push_interval("chr2", 100, 110, 1.0).unwrap(),
        );
        assert!(matches!(error, Error::ChromAlreadyFinished(name) if name == "chr1"));
    }

    #[test]
    fn rejected_push_to_another_chromosome_does_not_finish_the_current_one() {
        // The pushes after the rejected one add chr1 data, which fails if chr1 was finished.
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step("chr2", 0, 1, 1, &[1.0, f32::NAN])
        });
        assert!(matches!(error, Error::NonFiniteValue { .. }));
    }

    #[test]
    fn item_starting_before_the_previous_end_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| writer.push_interval("chr1", 99, 120, 1.0));
        assert!(matches!(error, Error::UnsortedOrOverlapping { start: 99, previous_end: 100, .. }));
    }

    #[test]
    fn unsorted_intervals_within_one_push_are_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_intervals("chr1", &[(300, 400, 1.0), (350, 500, 2.0)])
        });
        assert!(matches!(
            error,
            Error::UnsortedOrOverlapping { start: 350, previous_end: 400, .. }
        ));
    }

    #[test]
    fn overlapping_variable_step_items_are_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_variable_step("chr1", 10, &[(300, 1.0), (305, 2.0)])
        });
        assert!(matches!(
            error,
            Error::UnsortedOrOverlapping { start: 305, previous_end: 310, .. }
        ));
    }

    #[test]
    fn overlap_across_section_types_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_variable_step("chr1", 1, &[(99, 1.0)]));
        assert!(matches!(error, Error::UnsortedOrOverlapping { start: 99, .. }));
    }

    #[test]
    fn fixed_step_starting_before_the_previous_end_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_fixed_step("chr1", 50, 1, 1, &[1.0]));
        assert!(matches!(error, Error::UnsortedOrOverlapping { start: 50, .. }));
    }

    #[test]
    fn empty_interval_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| writer.push_interval("chr1", 500, 500, 1.0));
        assert!(matches!(error, Error::InvalidItem { start: 500, .. }));
    }

    #[test]
    fn inverted_interval_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| writer.push_interval("chr1", 500, 400, 1.0));
        assert!(matches!(error, Error::InvalidItem { start: 500, .. }));
    }

    #[test]
    fn zero_step_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_fixed_step("chr1", 500, 0, 0, &[1.0]));
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("step")));
    }

    #[test]
    fn zero_fixed_step_span_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_fixed_step("chr1", 500, 1, 0, &[1.0]));
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("span")));
    }

    #[test]
    fn zero_variable_step_span_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_variable_step("chr1", 0, &[(500, 1.0)]));
        assert!(matches!(
            error,
            Error::InvalidItem { start: 500, reason, .. } if reason.contains("span")
        ));
    }

    #[test]
    fn span_above_step_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_fixed_step("chr1", 500, 2, 3, &[1.0]));
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("exceed")));
    }

    #[test]
    fn empty_fixed_step_slice_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_fixed_step("chr1", 500, 1, 1, &[]));
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("empty")));
    }

    #[test]
    fn empty_variable_step_slice_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| writer.push_variable_step("chr1", 1, &[]));
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("empty")));
    }

    #[test]
    fn empty_interval_slice_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| writer.push_intervals("chr1", &[]));
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("empty")));
    }

    #[test]
    fn fixed_step_past_the_chromosome_end_is_rejected_at_its_first_item_past_it() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step("chrM", 16_560, 1, 1, &values(20, 1))
        });
        assert!(matches!(
            error,
            Error::OutOfBounds { chrom, start: 16_569, end: 16_570, length: 16_569 }
                if chrom == "chrM"
        ));
    }

    #[test]
    fn variable_step_past_the_chromosome_end_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_variable_step("chrM", 10, &[(100, 1.0), (16_565, 2.0)])
        });
        assert!(matches!(error, Error::OutOfBounds { start: 16_565, end: 16_575, .. }));
    }

    #[test]
    fn interval_past_the_chromosome_end_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_interval("chrM", 16_000, 16_570, 1.0));
        assert!(matches!(error, Error::OutOfBounds { start: 16_000, end: 16_570, .. }));
    }

    #[test]
    fn fixed_step_run_past_u32_max_is_out_of_bounds() {
        let chroms = [ChromInfo::new("huge", u32::MAX)];
        let error = rejection_without_trace(
            &WriterOptions::new(),
            &chroms,
            |_| {},
            |writer| writer.push_fixed_step("huge", u32::MAX - 5, 1_000, 1, &[1.0, 2.0, 3.0]),
            |writer| writer.push_interval("huge", 0, 10, 1.0).unwrap(),
        );
        let past_u32 = u64::from(u32::MAX) - 5 + 1_000;
        assert!(matches!(
            error,
            Error::OutOfBounds { start, end, .. } if start == past_u32 && end == past_u32 + 1
        ));
    }

    #[test]
    fn nan_is_rejected_without_trace_at_its_position() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step("chr1", 500, 10, 1, &[1.0, 2.0, f32::NAN])
        });
        assert!(
            matches!(error, Error::NonFiniteValue { position: 520, value, .. } if value.is_nan())
        );
    }

    #[test]
    fn infinite_interval_value_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_interval("chr1", 500, 510, f32::INFINITY)
        });
        assert!(matches!(error, Error::NonFiniteValue { position: 500, .. }));
    }

    #[test]
    fn infinite_variable_step_value_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_variable_step("chr1", 1, &[(500, 1.0), (600, f32::NEG_INFINITY)])
        });
        assert!(matches!(error, Error::NonFiniteValue { position: 600, .. }));
    }

    #[test]
    fn fixed_step_push_reports_a_non_finite_value_before_the_first_value_past_the_end() {
        let mut writer = new_writer(&WriterOptions::new(), &short_chrom());
        let mut values = vec![1.0; 20];
        values[3] = f32::NAN;
        let error = writer.push_fixed_step("chr1", 90, 1, 1, &values).unwrap_err();
        assert!(matches!(error, Error::NonFiniteValue { position: 93, .. }));
    }

    // ---- Indexed pushes -------------------------------------------------------------------

    #[test]
    fn indexed_fixed_step_push_gives_the_index_of_the_first_value_past_the_end() {
        let mut writer = new_writer(&WriterOptions::new(), &short_chrom());
        let rejected = writer.push_fixed_step_indexed("chr1", 90, 1, 1, &[1.0; 20]).unwrap_err();
        assert_eq!(rejected.index, 10);
        assert!(matches!(rejected.error, Error::OutOfBounds { start: 100, .. }));
    }

    #[test]
    fn indexed_fixed_step_push_gives_the_index_of_a_non_finite_value() {
        let mut writer = new_writer(&WriterOptions::new(), &short_chrom());
        let rejected =
            writer.push_fixed_step_indexed("chr1", 0, 10, 1, &[1.0, 2.0, f32::NAN]).unwrap_err();
        assert_eq!(rejected.index, 2);
        assert!(matches!(rejected.error, Error::NonFiniteValue { position: 20, .. }));
    }

    #[test]
    fn indexed_variable_step_push_gives_the_index_of_the_rejected_item() {
        let mut writer = new_writer(&WriterOptions::new(), &short_chrom());
        let items = [(0, 1.0), (10, 2.0), (15, 3.0)];
        let rejected = writer.push_variable_step_indexed("chr1", 10, &items).unwrap_err();
        assert_eq!(rejected.index, 2);
        assert!(matches!(rejected.error, Error::UnsortedOrOverlapping { start: 15, .. }));
    }

    #[test]
    fn indexed_push_gives_index_zero_for_an_error_about_the_chromosome() {
        let mut writer = new_writer(&WriterOptions::new(), &short_chrom());
        let rejected = writer.push_variable_step_indexed("chrX", 1, &[(0, 1.0), (5, 2.0)]);
        let rejected = rejected.unwrap_err();
        assert_eq!(rejected.index, 0);
        assert!(matches!(rejected.error, Error::UnknownChrom(_)));
    }

    #[test]
    fn indexed_push_rejection_leaves_no_trace() {
        let options = WriterOptions::new();
        rejection_without_trace(
            &options,
            &short_chrom(),
            |writer| writer.push_interval("chr1", 0, 10, 1.0).unwrap(),
            |writer| {
                let items = [(20, 1.0), (30, 2.0), (35, f32::NAN)];
                writer.push_variable_step_indexed("chr1", 1, &items).map_err(|e| e.error)
            },
            |writer| writer.push_interval("chr1", 20, 30, 2.0).unwrap(),
        );
    }

    // ---- fixedStep pushes with NaN gaps ---------------------------------------------------

    /// `values` with NaN over each of `gaps`.
    fn with_gaps(mut values: Vec<f32>, gaps: &[Range<usize>]) -> Vec<f32> {
        for gap in gaps {
            values[gap.clone()].fill(f32::NAN);
        }
        values
    }

    /// Pushes each maximal run of `values` without NaN in a `push_fixed_step` call of its own.
    fn push_each_run(
        writer: &mut TestWriter,
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) {
        let mut index = 0;
        while index < values.len() {
            let run_len = values[index..].iter().take_while(|value| !value.is_nan()).count();
            if run_len > 0 {
                let run_start = start + u32::try_from(index).unwrap() * step;
                let run = &values[index..index + run_len];
                writer.push_fixed_step(chrom, run_start, step, span, run).unwrap();
            }
            index += run_len + 1;
        }
    }

    /// Asserts that pushing `values` to chr1 from 1,000 with step 2 and span 1, with
    /// `push_fixed_step_skip` in chunks ending at each of `splits`, writes the same file as one
    /// push. Sections hold four values, so runs fill several.
    fn assert_split_gapped_push_writes_the_same_file(values: &[f32], splits: &[usize]) {
        let options = WriterOptions::new().items_per_section(4);
        let (start, step, span) = (1_000, 2, 1);
        let whole = write_file(&options, |writer| {
            writer.push_fixed_step_skip("chr1", start, step, span, values).unwrap();
        });
        let split = write_file(&options, |writer| {
            let mut from = 0;
            for &to in splits.iter().chain([&values.len()]) {
                let chunk_start = start + u32::try_from(from).unwrap() * step;
                let chunk = &values[from..to];
                writer.push_fixed_step_skip("chr1", chunk_start, step, span, chunk).unwrap();
                from = to;
            }
        });
        assert_eq!(split, whole);
    }

    #[test]
    fn nan_values_of_a_gapped_push_are_gaps_between_runs() {
        let nan = f32::NAN;
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            let values = [1.0, nan, nan, 2.0, 3.0, nan, 4.0];
            writer.push_fixed_step_skip("chr1", 100, 10, 3, &values).unwrap();
        }));
        let chrom_id = file.chrom_id("chr1");
        let mut expected = fixed_step_items(chrom_id, 100, 10, 3, &[1.0]);
        expected.extend(fixed_step_items(chrom_id, 130, 10, 3, &[2.0, 3.0]));
        expected.extend(fixed_step_items(chrom_id, 160, 10, 3, &[4.0]));
        assert_eq!(file.items(), expected);
        assert_eq!(file.section_sizes(), [1, 2, 1]);
    }

    #[test]
    fn gapped_push_writes_what_one_fixed_step_push_per_run_writes() {
        let options = WriterOptions::new().items_per_section(4);
        let gapped = with_gaps(values(50, 1), &[0..2, 10..14, 20..21, 48..50]);
        let one_push = write_file(&options, |writer| {
            writer.push_fixed_step_skip("chr1", 1_000, 5, 3, &gapped).unwrap();
        });
        let per_run = write_file(&options, |writer| {
            push_each_run(writer, "chr1", 1_000, 5, 3, &gapped);
        });
        assert_eq!(one_push, per_run);
    }

    #[test]
    fn gapped_push_split_inside_a_gap_writes_the_same_file() {
        let gapped = with_gaps(values(40, 1), &[10..14, 25..26]);
        assert_split_gapped_push_writes_the_same_file(&gapped, &[12]);
    }

    #[test]
    fn gapped_push_split_at_the_edges_of_gaps_writes_the_same_file() {
        let gapped = with_gaps(values(40, 1), &[10..14, 25..26]);
        assert_split_gapped_push_writes_the_same_file(&gapped, &[10, 26]);
    }

    #[test]
    fn gapped_push_split_inside_a_run_writes_the_same_file() {
        let gapped = with_gaps(values(40, 1), &[10..14, 25..26]);
        assert_split_gapped_push_writes_the_same_file(&gapped, &[3, 20]);
    }

    #[test]
    fn gapped_push_split_around_an_all_nan_chunk_writes_the_same_file() {
        let gapped = with_gaps(values(40, 1), &[10..14, 25..26]);
        assert_split_gapped_push_writes_the_same_file(&gapped, &[10, 14]);
    }

    #[test]
    fn leading_gaps_of_a_gapped_push_may_overlap_the_previous_item() {
        let file = decode(&write_file(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &values(100, 1)).unwrap();
            let mut gapped = vec![1.0; 20];
            gapped[..10].fill(f32::NAN);
            writer.push_fixed_step_skip("chr1", 90, 1, 1, &gapped).unwrap();
        }));
        // The first value continues the run that ended at 100.
        assert_eq!(file.section_sizes(), [110]);
    }

    #[test]
    fn trailing_gaps_of_a_gapped_push_do_not_order_later_items() {
        let mut writer = new_writer(&WriterOptions::new(), &test_chroms());
        writer.push_fixed_step_skip("chr1", 0, 1, 1, &[1.0, 2.0, f32::NAN, f32::NAN]).unwrap();
        writer.push_interval("chr1", 2, 3, 3.0).unwrap();
    }

    #[test]
    fn gapped_push_whose_first_run_overlaps_the_previous_item_is_rejected_without_trace() {
        let mut gapped = vec![1.0; 10];
        gapped[..5].fill(f32::NAN);
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step_skip("chr1", 90, 1, 1, &gapped)
        });
        assert!(matches!(error, Error::UnsortedOrOverlapping { start: 95, previous_end: 100, .. }));
    }

    #[test]
    fn gapped_push_reports_a_first_item_past_u32_max_exactly() {
        let chroms = [ChromInfo::new("huge", u32::MAX)];
        let mut writer = new_writer(&WriterOptions::new().clip(true), &chroms);
        writer.push_variable_step("huge", 100, &[(u32::MAX - 10, 1.0)]).unwrap();
        // The leading gap puts the first item at u32::MAX + 5, before the unclipped end above.
        let error =
            writer.push_fixed_step_skip("huge", u32::MAX - 5, 10, 1, &[f32::NAN, 2.0]).unwrap_err();
        assert!(matches!(
            error,
            Error::UnsortedOrOverlapping { start, previous_end, .. }
                if start == u64::from(u32::MAX) + 5 && previous_end == u64::from(u32::MAX) + 90
        ));
    }

    #[test]
    fn gapped_push_whose_last_run_is_past_the_end_is_rejected_without_trace() {
        // Runs cover 16,500-16,510 and 16,550-16,580 of chrM, which ends at 16,569.
        let mut gapped = vec![1.0; 80];
        gapped[10..50].fill(f32::NAN);
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step_skip("chrM", 16_500, 1, 1, &gapped)
        });
        assert!(matches!(error, Error::OutOfBounds { start: 16_569, end: 16_570, .. }));
    }

    #[test]
    fn infinity_in_a_gapped_push_is_rejected_without_trace_at_its_position() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step_skip("chr1", 500, 10, 1, &[1.0, f32::NAN, f32::INFINITY])
        });
        assert!(matches!(
            error,
            Error::NonFiniteValue { position: 520, value, .. } if value == f32::INFINITY
        ));
    }

    #[test]
    fn infinity_in_a_gapped_push_is_rejected_where_clipping_would_drop_it() {
        let error = rejection_without_trace(
            &WriterOptions::new().clip(true),
            &test_chroms(),
            |writer| writer.push_fixed_step("chrM", 0, 1, 1, &values(100, 1)).unwrap(),
            |writer| {
                let values = [1.0, f32::NAN, f32::NEG_INFINITY];
                writer.push_fixed_step_skip("chrM", 16_560, 10, 1, &values)
            },
            |writer| writer.push_fixed_step("chrM", 100, 1, 1, &values(100, 2)).unwrap(),
        );
        assert!(matches!(error, Error::NonFiniteValue { position: 16_580, .. }));
    }

    #[test]
    fn clip_of_a_gapped_push_counts_and_writes_as_one_push_per_run() {
        let options = WriterOptions::new().clip(true);
        // Values start every 4 bases from 16,540 and chrM ends at 16,569: the value at 16,568 is
        // truncated, those at 16,572 and 16,584 are dropped, and the NaN ones are gaps.
        let nan = f32::NAN;
        let gapped = [1.0, 2.0, nan, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, nan, nan, 9.0];
        let mut one_push = new_writer(&options, &test_chroms());
        one_push.push_fixed_step_skip("chrM", 16_540, 4, 3, &gapped).unwrap();
        let mut per_run = new_writer(&options, &test_chroms());
        push_each_run(&mut per_run, "chrM", 16_540, 4, 3, &gapped);
        let counts = |writer: &TestWriter| (writer.clipped_items(), writer.dropped_items());
        assert_eq!(counts(&one_push), (1, 2));
        assert_eq!(counts(&one_push), counts(&per_run));
        assert_eq!(one_push.finish().unwrap().into_inner(), per_run.finish().unwrap().into_inner());
    }

    #[test]
    fn all_nan_gapped_push_writes_nothing_and_leaves_the_current_chromosome_open() {
        let push_chr1 = |writer: &mut TestWriter, start: u32, seed: usize| {
            writer.push_fixed_step("chr1", start, 1, 1, &values(100, seed)).unwrap();
        };
        let with_nan_push = write_file(&WriterOptions::new(), |writer| {
            push_chr1(writer, 0, 1);
            writer.push_fixed_step_skip("chr2", 0, 1, 1, &[f32::NAN; 10]).unwrap();
            push_chr1(writer, 100, 2);
        });
        let without_nan_push = write_file(&WriterOptions::new(), |writer| {
            push_chr1(writer, 0, 1);
            push_chr1(writer, 100, 2);
        });
        assert_eq!(with_nan_push, without_nan_push);
    }

    #[test]
    fn all_nan_gapped_push_to_a_finished_chromosome_is_accepted() {
        let mut writer = new_writer(&WriterOptions::new(), &test_chroms());
        writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        writer.push_interval("chr2", 0, 10, 1.0).unwrap();
        writer.push_fixed_step_skip("chr1", 100, 1, 1, &[f32::NAN; 3]).unwrap();
        writer.push_interval("chr2", 10, 20, 1.0).unwrap();
    }

    #[test]
    fn all_nan_gapped_push_to_an_earlier_listed_chromosome_is_accepted() {
        let options = WriterOptions::new().chrom_id_order(ChromIdOrder::Listed);
        let mut writer = new_writer(&options, &test_chroms());
        writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        // chr2 is listed before chr1.
        writer.push_fixed_step_skip("chr2", 0, 1, 1, &[f32::NAN; 3]).unwrap();
        writer.push_interval("chr1", 10, 20, 1.0).unwrap();
    }

    #[test]
    fn all_nan_gapped_push_does_not_order_later_items() {
        let mut writer = new_writer(&WriterOptions::new(), &test_chroms());
        writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        writer.push_fixed_step_skip("chr1", 500, 1, 1, &[f32::NAN; 100]).unwrap();
        writer.push_interval("chr1", 10, 20, 1.0).unwrap();
    }

    #[test]
    fn all_nan_gapped_push_to_an_unknown_chromosome_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step_skip("chrX", 0, 1, 1, &[f32::NAN])
        });
        assert!(matches!(error, Error::UnknownChrom(name) if name == "chrX"));
    }

    #[test]
    fn all_nan_gapped_push_with_span_above_step_is_rejected_without_trace() {
        let error = rejection_after_chr1_data(|writer| {
            writer.push_fixed_step_skip("chr1", 500, 2, 3, &[f32::NAN])
        });
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("exceed")));
    }

    #[test]
    fn empty_gapped_push_is_rejected_without_trace() {
        let error =
            rejection_after_chr1_data(|writer| writer.push_fixed_step_skip("chr1", 500, 1, 1, &[]));
        assert!(matches!(error, Error::InvalidItem { reason, .. } if reason.contains("empty")));
    }

    // ---- Skipped empty intervals -----------------------------------------------------------

    #[test]
    fn skipping_an_empty_interval_writes_nothing() {
        let with_skip = write_file_over(&WriterOptions::new(), &short_chrom(), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
            writer.skip_empty_interval("chr1", 20, 2.0).unwrap();
            writer.push_interval("chr1", 30, 40, 3.0).unwrap();
        });
        let without_skip = write_file_over(&WriterOptions::new(), &short_chrom(), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
            writer.push_interval("chr1", 30, 40, 3.0).unwrap();
        });
        assert_eq!(with_skip, without_skip);
    }

    #[test]
    fn skipping_an_empty_interval_checks_it_as_a_push() {
        let error = rejection_without_trace(
            &WriterOptions::new(),
            &short_chrom(),
            |writer| writer.push_interval("chr1", 0, 10, 1.0).unwrap(),
            |writer| writer.skip_empty_interval("chr1", 5, 2.0),
            |writer| writer.push_interval("chr1", 30, 40, 3.0).unwrap(),
        );
        assert!(matches!(error, Error::UnsortedOrOverlapping { start: 5, previous_end: 10, .. }));
    }

    #[test]
    fn items_may_not_start_before_a_skipped_empty_interval() {
        let mut writer = new_writer(&WriterOptions::new(), &short_chrom());
        writer.skip_empty_interval("chr1", 25, 1.0).unwrap();
        let error = writer.push_interval("chr1", 22, 30, 2.0).unwrap_err();
        assert!(matches!(error, Error::UnsortedOrOverlapping { start: 22, previous_end: 25, .. }));
    }

    // ---- Clipping -------------------------------------------------------------------------

    #[test]
    fn clip_truncates_an_interval_at_the_chromosome_end() {
        let file = decode(&write_file(&WriterOptions::new().clip(true), |writer| {
            writer.push_intervals("chrM", &[(16_000, 16_100, 1.0), (16_500, 16_600, 2.0)]).unwrap();
        }));
        let chrom_id = file.chrom_id("chrM");
        let expected = [
            DataItem { chrom_id, start: 16_000, end: 16_100, value: 1.0 },
            DataItem { chrom_id, start: 16_500, end: 16_569, value: 2.0 },
        ];
        assert_eq!(file.items(), expected);
    }

    #[test]
    fn clip_drops_intervals_starting_at_or_past_the_end() {
        let file = decode(&write_file(&WriterOptions::new().clip(true), |writer| {
            let intervals = [(16_000, 16_100, 1.0), (16_569, 16_600, 2.0), (17_000, 17_100, 3.0)];
            writer.push_intervals("chrM", &intervals).unwrap();
        }));
        let chrom_id = file.chrom_id("chrM");
        assert_eq!(file.items(), [DataItem { chrom_id, start: 16_000, end: 16_100, value: 1.0 }]);
    }

    #[test]
    fn clip_truncates_a_fixed_step_item_into_a_section_of_its_own() {
        let file = decode(&write_file(&WriterOptions::new().clip(true), |writer| {
            writer.push_fixed_step("chrM", 16_560, 4, 3, &[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        }));
        let chrom_id = file.chrom_id("chrM");
        let expected = [
            DataItem { chrom_id, start: 16_560, end: 16_563, value: 1.0 },
            DataItem { chrom_id, start: 16_564, end: 16_567, value: 2.0 },
            DataItem { chrom_id, start: 16_568, end: 16_569, value: 3.0 },
        ];
        assert_eq!(file.items(), expected);
        let kinds: Vec<SectionKind> = file.sections.iter().map(|section| section.kind).collect();
        let expected_kinds = [
            SectionKind::FixedStep { step: 4, span: 3 },
            SectionKind::FixedStep { step: 4, span: 1 },
        ];
        assert_eq!(kinds, expected_kinds);
    }

    #[test]
    fn clip_truncates_a_variable_step_item_into_a_section_of_its_own() {
        let file = decode(&write_file(&WriterOptions::new().clip(true), |writer| {
            let items = [(16_550, 1.0), (16_565, 2.0), (16_580, 3.0)];
            writer.push_variable_step("chrM", 10, &items).unwrap();
        }));
        let chrom_id = file.chrom_id("chrM");
        let expected = [
            DataItem { chrom_id, start: 16_550, end: 16_560, value: 1.0 },
            DataItem { chrom_id, start: 16_565, end: 16_569, value: 2.0 },
        ];
        assert_eq!(file.items(), expected);
        let kinds: Vec<SectionKind> = file.sections.iter().map(|section| section.kind).collect();
        let expected_kinds =
            [SectionKind::VariableStep { span: 10 }, SectionKind::VariableStep { span: 4 }];
        assert_eq!(kinds, expected_kinds);
    }

    #[test]
    fn clip_of_every_item_writes_nothing_and_leaves_the_chromosome_unfinished() {
        let file = decode(&write_file(&WriterOptions::new().clip(true), |writer| {
            writer.push_intervals("chrM", &[(17_000, 17_100, 1.0)]).unwrap();
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
            writer.push_interval("chrM", 17_100, 17_200, 2.0).unwrap();
        }));
        assert_eq!(file.items().len(), 1);
        assert_eq!(file.chroms[0].name, "chr1");
    }

    #[test]
    fn clip_writes_the_same_file_for_a_split_push_as_for_the_single_push() {
        let options = WriterOptions::new().clip(true);
        let intervals = [(80, 90, 1.0), (90, 110, 2.0), (110, 120, 3.0)];
        let single = write_file_over(&options, &short_chrom(), |writer| {
            writer.push_intervals("chr1", &intervals).unwrap();
        });
        let split = write_file_over(&options, &short_chrom(), |writer| {
            for interval in &intervals {
                writer.push_intervals("chr1", std::slice::from_ref(interval)).unwrap();
            }
        });
        assert_eq!(split, single);
    }

    #[test]
    fn clip_rejects_overlap_with_a_truncated_interval_whether_or_not_the_push_is_split() {
        let options = WriterOptions::new().clip(true);
        let intervals = [(90, 110, 1.0), (105, 120, 2.0)];
        let single_error =
            new_writer(&options, &short_chrom()).push_intervals("chr1", &intervals).unwrap_err();
        let mut split = new_writer(&options, &short_chrom());
        split.push_intervals("chr1", &intervals[..1]).unwrap();
        let split_error = split.push_intervals("chr1", &intervals[1..]).unwrap_err();
        assert!(matches!(
            split_error,
            Error::UnsortedOrOverlapping { start: 105, previous_end: 110, .. }
        ));
        assert_eq!(split_error.to_string(), single_error.to_string());
    }

    #[test]
    fn clip_rejects_overlap_with_a_truncated_variable_step_item_whether_or_not_the_push_is_split() {
        let options = WriterOptions::new().clip(true);
        let items = [(95, 1.0), (101, 2.0)];
        let single_error = new_writer(&options, &short_chrom())
            .push_variable_step("chr1", 10, &items)
            .unwrap_err();
        let mut split = new_writer(&options, &short_chrom());
        split.push_variable_step("chr1", 10, &items[..1]).unwrap();
        let split_error = split.push_variable_step("chr1", 10, &items[1..]).unwrap_err();
        assert!(matches!(
            split_error,
            Error::UnsortedOrOverlapping { start: 101, previous_end: 105, .. }
        ));
        assert_eq!(split_error.to_string(), single_error.to_string());
    }

    #[test]
    fn clip_rejects_an_interval_starting_before_one_it_dropped() {
        let mut writer = new_writer(&WriterOptions::new().clip(true), &short_chrom());
        writer.push_interval("chr1", 80, 90, 1.0).unwrap();
        writer.push_interval("chr1", 150, 160, 2.0).unwrap();
        let error = writer.push_interval("chr1", 95, 99, 3.0).unwrap_err();
        assert!(matches!(error, Error::UnsortedOrOverlapping { start: 95, previous_end: 160, .. }));
    }

    #[test]
    fn clip_rejects_an_item_starting_inside_a_dropped_fixed_step_item() {
        let mut writer = new_writer(&WriterOptions::new().clip(true), &short_chrom());
        // Items cover 98-103 (truncated), 108-113 and 118-123 (dropped).
        writer.push_fixed_step("chr1", 98, 10, 5, &[1.0, 2.0, 3.0]).unwrap();
        let error = writer.push_variable_step("chr1", 1, &[(120, 4.0)]).unwrap_err();
        assert!(matches!(
            error,
            Error::UnsortedOrOverlapping { start: 120, previous_end: 123, .. }
        ));
    }

    #[test]
    fn clip_rejects_an_item_starting_before_a_wholly_dropped_first_push() {
        let mut writer = new_writer(&WriterOptions::new().clip(true), &test_chroms());
        writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        writer.push_interval("chrM", 17_000, 17_100, 2.0).unwrap();
        let error = writer.push_interval("chrM", 16_000, 16_100, 3.0).unwrap_err();
        assert!(matches!(
            error,
            Error::UnsortedOrOverlapping { start: 16_000, previous_end: 17_100, .. }
        ));
    }

    #[test]
    fn clip_reports_a_previous_end_past_u32_max_exactly() {
        let chroms = [ChromInfo::new("huge", u32::MAX)];
        let mut writer = new_writer(&WriterOptions::new().clip(true), &chroms);
        writer.push_variable_step("huge", 100, &[(u32::MAX - 10, 1.0)]).unwrap();
        let error = writer.push_variable_step("huge", 1, &[(u32::MAX, 2.0)]).unwrap_err();
        let previous_end = u64::from(u32::MAX) + 90;
        assert!(matches!(
            error,
            Error::UnsortedOrOverlapping { start, previous_end: end, .. }
                if start == u64::from(u32::MAX) && end == previous_end
        ));
    }

    #[test]
    fn clip_counts_a_truncated_fixed_step_value_and_the_dropped_values_after_it() {
        let mut writer = new_writer(&WriterOptions::new().clip(true), &test_chroms());
        // Values start at 16_560, 16_564, 16_568 (truncated), 16_572 and 16_576 (dropped).
        writer.push_fixed_step("chrM", 16_560, 4, 3, &[1.0; 5]).unwrap();
        assert_eq!((writer.clipped_items(), writer.dropped_items()), (1, 2));
    }

    #[test]
    fn clip_counts_a_truncated_variable_step_item_and_the_dropped_items_after_it() {
        let mut writer = new_writer(&WriterOptions::new().clip(true), &test_chroms());
        let items = [(16_550, 1.0), (16_565, 2.0), (16_580, 3.0)];
        writer.push_variable_step("chrM", 10, &items).unwrap();
        assert_eq!((writer.clipped_items(), writer.dropped_items()), (1, 1));
    }

    #[test]
    fn clip_counts_a_truncated_interval_and_the_dropped_intervals_after_it() {
        let mut writer = new_writer(&WriterOptions::new().clip(true), &test_chroms());
        let intervals = [(16_000, 16_100, 1.0), (16_500, 16_600, 2.0), (17_000, 17_100, 3.0)];
        writer.push_intervals("chrM", &intervals).unwrap();
        assert_eq!((writer.clipped_items(), writer.dropped_items()), (1, 1));
    }

    #[test]
    fn clip_counts_nothing_for_a_rejected_push() {
        let mut writer = new_writer(&WriterOptions::new().clip(true), &test_chroms());
        let intervals = [(16_500, 16_600, 1.0), (17_000, 17_100, f32::NAN)];
        writer.push_intervals("chrM", &intervals).unwrap_err();
        assert_eq!((writer.clipped_items(), writer.dropped_items()), (0, 0));
    }

    #[test]
    fn clipped_items_reach_the_summary_and_zoom_levels_truncated() {
        let file = decode(&write_file(&WriterOptions::new().clip(true), |writer| {
            writer.push_fixed_step("chrM", 16_560, 4, 3, &[1.0, 2.0, 3.0, 4.0]).unwrap();
        }));
        assert_eq!(file.total_summary, summary_of(&file.items()));
        assert_eq!(file.total_summary.covered_bases, 7.0);
        assert_eq!(file.zoom[0].records().last().unwrap().end, 16_569);
    }

    // ---- Empty files, finishing and failure ------------------------------------------------

    #[test]
    fn file_without_data_finishes_valid() {
        let file = decode(&write_file(&WriterOptions::new(), |_| {}));
        assert_eq!(file.section_count, 0);
        assert!(file.sections.is_empty());
        assert_eq!(file.data_index.item_count, 0);
        assert_eq!((file.data_index.root_is_leaf, file.data_index.root_count), (1, 0));
        assert_eq!(file.reductions(), [256]);
        assert_eq!(file.zoom[0].zoom_count, 0);
        assert_eq!((file.zoom[0].index.root_is_leaf, file.zoom[0].index.root_count), (1, 0));
        assert_ne!(file.header.total_summary_offset, 0);
        let zero = Summary { covered_bases: 0.0, min: 0.0, max: 0.0, sum: 0.0, sum_squares: 0.0 };
        assert_eq!(file.total_summary, zero);
        assert_eq!(file.chroms.len(), 4);
    }

    #[test]
    fn empty_compressed_file_has_a_nonzero_uncompress_buf_size() {
        let file = decode(&write_file(&WriterOptions::new(), |_| {}));
        assert_eq!(file.header.uncompress_buf_size, 1_024 * 32);
    }

    #[test]
    fn empty_uncompressed_file_has_a_zero_uncompress_buf_size() {
        let options = WriterOptions::new().compression(Compression::None);
        let file = decode(&write_file(&options, |_| {}));
        assert_eq!(file.header.uncompress_buf_size, 0);
    }

    #[test]
    fn file_without_data_keeps_one_explicit_level() {
        let options = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10, 100]));
        let file = decode(&write_file(&options, |_| {}));
        assert_eq!(file.reductions(), [10]);
        assert_eq!(file.zoom[0].zoom_count, 0);
    }

    #[test]
    fn finish_patches_every_header_field() {
        let bytes = write_file(&WriterOptions::new(), push_mixed_data);
        let file = decode(&bytes);
        let header = file.header;
        assert_eq!(header.version, 4);
        assert!(header.zoom_levels >= 1);
        assert_eq!(usize::from(header.zoom_levels), file.zoom.len());
        assert_eq!(header.total_summary_offset, 64 + 10 * 24);
        assert_eq!(header.chrom_tree_offset, header.total_summary_offset + 40);
        let tree_len = chrom_tree_layout(&bytes).len;
        assert_eq!(header.full_data_offset, header.chrom_tree_offset + tree_len);
        assert_eq!(file.sections[0].entry.data_offset, header.full_data_offset + 8);
        assert_eq!(file.section_count, file.sections.len() as u64);
        assert_eq!(file.data_index.item_count, file.section_count);
        assert_eq!((header.field_count, header.defined_field_count), (0, 0));
        assert_eq!((header.auto_sql_offset, header.extension_offset), (0, 0));
        assert_eq!(header.uncompress_buf_size as usize, file.largest_block());
        assert_eq!(file.total_summary, summary_of(&file.items()));
        // The placeholder chromIds were replaced: chr1 received data first.
        assert_eq!(file.chroms[0].name, "chr1");
    }

    #[test]
    fn finish_returns_the_sink_positioned_at_the_end() {
        let mut writer = new_writer(&WriterOptions::new(), &test_chroms());
        push_mixed_data(&mut writer);
        let sink = writer.finish().unwrap();
        assert_eq!(sink.position(), sink.get_ref().len() as u64);
    }

    #[test]
    fn create_with_rejects_a_sink_that_is_not_empty() {
        let result = WriterOptions::new().create_with(Cursor::new(vec![0xAA; 10]), &test_chroms());
        let Err(Error::Io(error)) = result else { panic!("expected an I/O error") };
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(error.to_string(), "the sink must be empty, but it holds 10 bytes");
    }

    #[test]
    fn create_with_writes_from_offset_zero_of_an_empty_sink_whatever_its_position() {
        let mut sink = Cursor::new(Vec::new());
        sink.set_position(10);
        let mut writer = WriterOptions::new().create_with(sink, &test_chroms()).unwrap();
        push_mixed_data(&mut writer);
        let bytes = writer.finish().unwrap().into_inner();
        assert_eq!(bytes, write_file(&WriterOptions::new(), push_mixed_data));
    }

    #[test]
    fn dropping_an_unfinished_writer_leaves_a_zero_magic() {
        let sink = SharedSink::default();
        let mut writer = WriterOptions::new().create_with(sink.clone(), &test_chroms()).unwrap();
        writer.push_fixed_step("chr1", 0, 1, 1, &values(100_000, 1)).unwrap();
        drop(writer);
        let bytes = sink.0.borrow().get_ref().clone();
        assert!(bytes.len() > 10_000, "the data written so far reached the sink");
        assert_eq!(bytes[..FIXED_HEADER_SIZE], [0; FIXED_HEADER_SIZE]);
        assert!(matches!(FixedHeader::decode(&bytes), Err(Error::NotBigWig { found: 0 })));
    }

    #[test]
    fn io_error_fails_every_later_call() {
        let failing = Rc::new(Cell::new(false));
        let sink = FailingSink { inner: Cursor::new(Vec::new()), failing: Rc::clone(&failing) };
        let options = WriterOptions::new().compression(Compression::None);
        let mut writer = options.create_with(sink, &test_chroms()).unwrap();
        failing.set(true);
        // 400 KB of raw sections overflow the writer's 256 KiB buffer, so the sink is written.
        let error = writer.push_fixed_step("chr1", 0, 1, 1, &values(100_000, 1)).unwrap_err();
        assert!(matches!(error, Error::Io(_)));
        failing.set(false);
        assert!(matches!(writer.push_interval("chr2", 0, 10, 1.0), Err(Error::WriterFailed)));
        assert!(matches!(writer.push_interval("chrX", 0, 10, 1.0), Err(Error::WriterFailed)));
        assert!(matches!(writer.finish(), Err(Error::WriterFailed)));
    }

    #[test]
    fn finish_returns_an_io_error_it_hits() {
        let failing = Rc::new(Cell::new(false));
        let sink = FailingSink { inner: Cursor::new(Vec::new()), failing: Rc::clone(&failing) };
        let mut writer = WriterOptions::new().create_with(sink, &test_chroms()).unwrap();
        writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        failing.set(true);
        assert!(matches!(writer.finish(), Err(Error::Io(_))));
    }

    // ---- Zoom spill -----------------------------------------------------------------------

    /// Options that make many small zoom blocks, so spilling has plenty to do.
    fn small_zoom_blocks() -> WriterOptions {
        WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10, 40, 160])).zoom_items_per_block(16)
    }

    /// The names of the entries in `dir`, sorted.
    fn dir_entries(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn spilled_zoom_blocks_give_the_same_file_as_blocks_kept_in_memory() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bw");
        let options = small_zoom_blocks().zoom_memory_limit(0);
        let mut writer = options.create(&path, &test_chroms()).unwrap();
        push_mixed_data(&mut writer);
        assert!(writer.zoom.store.spill.is_some(), "the zoom blocks spilled");
        assert!(writer.zoom.store.levels.iter().all(|level| !level.spilled.is_empty()));
        writer.finish().unwrap();

        let in_memory = write_file(&small_zoom_blocks(), push_mixed_data);
        assert_eq!(fs::read(&path).unwrap(), in_memory);
        assert_structure(&decode(&in_memory));
    }

    #[test]
    fn spilling_starts_only_beyond_the_memory_limit() {
        let options = small_zoom_blocks().zoom_memory_limit(usize::MAX);
        let mut writer = new_writer(&options, &test_chroms());
        push_mixed_data(&mut writer);
        assert!(writer.zoom.store.spill.is_none());
    }

    #[test]
    fn spill_file_is_removed_by_finish() {
        let dir = tempfile::tempdir().unwrap();
        let options = small_zoom_blocks().zoom_memory_limit(0);
        let mut writer = options.create(dir.path().join("out.bw"), &test_chroms()).unwrap();
        push_mixed_data(&mut writer);
        writer.finish().unwrap();
        assert_eq!(dir_entries(dir.path()), ["out.bw"]);
    }

    #[test]
    fn spill_file_is_removed_when_an_unfinished_writer_is_dropped() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bw");
        let options = small_zoom_blocks().zoom_memory_limit(0);
        let mut writer = options.create(&path, &test_chroms()).unwrap();
        push_mixed_data(&mut writer);
        assert!(writer.zoom.store.spill.is_some());
        drop(writer);
        assert_eq!(dir_entries(dir.path()), ["out.bw"]);
        assert_eq!(fs::read(&path).unwrap()[..4], [0; 4]);
    }

    #[cfg(unix)]
    #[test]
    fn spill_file_is_unlinked_as_soon_as_it_is_created() {
        let dir = tempfile::tempdir().unwrap();
        let options = small_zoom_blocks().zoom_memory_limit(0);
        let mut writer = options.create(dir.path().join("out.bw"), &test_chroms()).unwrap();
        push_mixed_data(&mut writer);
        assert!(writer.zoom.store.spill.as_ref().is_some_and(|spill| spill.path.is_none()));
        assert_eq!(dir_entries(dir.path()), ["out.bw"]);
        writer.finish().unwrap();
    }

    #[test]
    fn spill_file_goes_next_to_the_output_named_after_it() {
        let dir = tempfile::tempdir().unwrap();
        let writer =
            WriterOptions::new().create(dir.path().join("out.bw"), &test_chroms()).unwrap();
        let location = &writer.zoom.store.spill_location;
        assert_eq!(location.dir, dir.path());
        assert_eq!(location.output_name, "out.bw");
    }

    #[test]
    fn spill_file_of_a_sink_without_a_path_goes_to_the_system_temp_dir() {
        let writer = new_writer(&WriterOptions::new(), &test_chroms());
        let location = &writer.zoom.store.spill_location;
        assert_eq!(location.dir, std::env::temp_dir());
        assert_eq!(location.output_name, UNNAMED_OUTPUT);
    }

    #[test]
    fn spill_file_goes_to_the_temp_dir_option() {
        let spill_dir = tempfile::tempdir().unwrap();
        let output_dir = tempfile::tempdir().unwrap();
        let options = small_zoom_blocks().zoom_memory_limit(0).temp_dir(spill_dir.path());
        let mut writer = options.create(output_dir.path().join("out.bw"), &test_chroms()).unwrap();
        assert_eq!(writer.zoom.store.spill_location.dir, spill_dir.path());
        push_mixed_data(&mut writer);
        assert!(writer.zoom.store.spill.is_some());
        writer.finish().unwrap();
        assert!(dir_entries(spill_dir.path()).is_empty());
    }

    #[test]
    fn spill_into_a_missing_directory_fails_the_writer() {
        let dir = tempfile::tempdir().unwrap();
        let options = WriterOptions::new()
            .zoom(ZoomPolicy::Explicit(vec![10]))
            .zoom_items_per_block(1)
            .zoom_memory_limit(0)
            .temp_dir(dir.path().join("missing"));
        let mut writer = new_writer(&options, &test_chroms());
        // The first section, and with it the first zoom block, is flushed when chr1 is finished.
        writer.push_fixed_step("chr1", 0, 1, 1, &values(100, 1)).unwrap();
        let error = writer.push_interval("chr2", 0, 10, 1.0).unwrap_err();
        assert!(matches!(error, Error::Io(_)));
        assert!(matches!(writer.push_interval("chr2", 0, 10, 1.0), Err(Error::WriterFailed)));
    }

    // ---- Paths ----------------------------------------------------------------------------

    #[test]
    fn create_writes_the_file_at_the_path() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bw");
        let mut writer = WriterOptions::new().create(&path, &test_chroms()).unwrap();
        push_mixed_data(&mut writer);
        writer.finish().unwrap();
        assert_eq!(fs::read(&path).unwrap(), write_file(&WriterOptions::new(), push_mixed_data));
    }

    #[test]
    fn create_error_names_the_path_and_keeps_the_error_kind() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("missing").join("out.bw");
        let result = WriterOptions::new().create(&path, &test_chroms());
        let Err(Error::Io(error)) = result else { panic!("expected an I/O error") };
        assert_eq!(error.kind(), io::ErrorKind::NotFound);
        let message = error.to_string();
        assert!(message.starts_with(&format!("cannot create {}: ", path.display())), "{message}");
    }

    #[test]
    fn create_truncates_an_existing_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bw");
        fs::write(&path, vec![0xAA; 1 << 20]).unwrap();
        let mut writer = WriterOptions::new().create(&path, &test_chroms()).unwrap();
        writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        writer.finish().unwrap();
        let expected = write_file(&WriterOptions::new(), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0).unwrap();
        });
        assert_eq!(fs::read(&path).unwrap(), expected);
    }

    // ---- Determinism ----------------------------------------------------------------------

    /// Pushes `values` per base from 0, as fixedStep calls of varying sizes.
    fn push_in_varying_chunks(writer: &mut TestWriter, chrom: &str, values: &[f32]) {
        let sizes = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1_597];
        let mut start = 0;
        for size in sizes.iter().cycle() {
            if start == values.len() {
                break;
            }
            let end = (start + size).min(values.len());
            let position = u32::try_from(start).unwrap();
            writer.push_fixed_step(chrom, position, 1, 1, &values[start..end]).unwrap();
            start = end;
        }
    }

    /// The same per-base data pushed as one slice per chromosome or as many small slices gives
    /// byte-identical files: the zoom ladder (which depended on chunking in bigtools #63), the
    /// records and every f64 sum depend on the data alone. The values make f64 sums round, so
    /// a change in summation order would show.
    #[test]
    fn zoom_ladder_is_independent_of_push_granularity() {
        let chr1 = inexact_values(100_000, 1);
        let chr2 = inexact_values(50_000, 2);
        let options = WriterOptions::new().items_per_section(1_000);
        let whole = write_file(&options, |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &chr1).unwrap();
            writer.push_fixed_step("chr2", 0, 1, 1, &chr2).unwrap();
        });
        let chunked = write_file(&options, |writer| {
            push_in_varying_chunks(writer, "chr1", &chr1);
            push_in_varying_chunks(writer, "chr2", &chr2);
        });
        assert_eq!(whole, chunked);
        assert_eq!(decode(&whole).reductions(), [256, 1_024, 4_096, 16_384, 65_536]);
    }

    #[test]
    fn same_input_gives_identical_bytes() {
        let options = WriterOptions::new().byte_order(ByteOrder::Big);
        assert_eq!(write_file(&options, push_mixed_data), write_file(&options, push_mixed_data));
    }

    // ---- Threads --------------------------------------------------------------------------

    #[test]
    fn writer_is_send_when_its_sink_is() {
        fn assert_send<T: Send>() {}
        assert_send::<BigWigWriter<File>>();
        assert_send::<BigWigWriter<Cursor<Vec<u8>>>>();
    }
}
