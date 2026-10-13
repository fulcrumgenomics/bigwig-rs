//! `BigWigReader`: opening a bigWig and querying its items, values, summaries and zoom levels.

use std::cmp::Reverse;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fmt;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::iter::FusedIterator;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::bbi::{
    ByteOrder, ByteReader, ChromInfo, FIXED_HEADER_SIZE, FixedHeader, TOTAL_SUMMARY_SIZE,
    ZOOM_HEADER_SIZE, ZoomLevel, decode_total_summary, read_exact_at, read_exact_at_into,
};
use crate::bptree;
use crate::compress::BlockDecompressor;
use crate::error::{Error, Result};
use crate::rtree::{LeafEntry, RTree, WalkBuffers};
use crate::section::{DataItem, Section, SectionKind};
use crate::zoom::{
    RawSummaryAccumulator, Summary, ZoomRecord, ZoomSummaryAccumulator, best_zoom, bin_bounds,
};

/// Blocks separated by at most this many bytes are fetched with one read, as kent, igv.js and
/// bbi-js merge them: a few wasted bytes cost less than another round trip on a remote source.
const MAX_MERGE_GAP: u64 = 8 * 1024;
/// A merged read stops growing at this size, so iterating a whole file holds a bounded amount of
/// compressed data rather than a chromosome's worth. A single larger block is still read whole.
const MAX_MERGED_READ: u64 = 1 << 20;
/// [`BigWigReader::summaries`] allows more bins than bases up to this many, for display; past it
/// every extra bin is empty, and the allocation could not otherwise be bounded.
const MIN_BIN_LIMIT: u32 = 1 << 20;

/// Decodes one decompressed block into `out` (replacing its contents), keeping what lies in the
/// region when one is given. Decoded values are left in `out` even when a later one fails.
type DecodeFn<T> = fn(&Metadata, &[u8], u64, Option<Region>, &mut Vec<T>) -> Result<()>;

/// A bigWig file opened for reading.
///
/// Opening parses the header, the zoom headers, the total summary and the whole chromosome tree;
/// queries then read only the index nodes and data blocks they need. Coordinates are 0-based and
/// half-open, and chromosomes are named as the file names them.
///
/// Queries take `&mut self` because they move the source's position. For parallel queries, open
/// one handle per thread with [`reopen`](Self::reopen) or [`with_source`](Self::with_source),
/// which share the parsed metadata. `Clone` is deliberately not implemented: a cloned
/// [`File`] shares its OS cursor, so cloned handles would race on seeks. A failed query leaves the
/// handle usable.
///
/// # Examples
///
/// ```no_run
/// use bigwig::{BigWigReader, SummarySource};
///
/// let mut reader = BigWigReader::open("signal.bw")?;
/// for chrom in reader.chroms() {
///     println!("{}\t{}", chrom.name, chrom.length);
/// }
///
/// // The stored items overlapping chr1:0-10,000, unclipped.
/// for item in reader.items("chr1", 0, 10_000)? {
///     let item = item?;
///     println!("{}\t{}\t{}", item.start, item.end, item.value);
/// }
///
/// // One value per base, NaN where there is no data.
/// let values = reader.values("chr1", 1_000, 1_100)?;
/// assert_eq!(values.len(), 100);
///
/// // Mean signal in 100 bins, from a zoom level where one is fine enough.
/// let means: Vec<Option<f64>> = reader
///     .summaries("chr1", 0, 1_000_000, 100, SummarySource::BestZoom)?
///     .iter()
///     .map(|bin| bin.map(|summary| summary.mean()))
///     .collect();
/// # Ok::<(), bigwig::Error>(())
/// ```
pub struct BigWigReader<R> {
    /// The file being read.
    source: R,
    /// Everything parsed at open, shared by every handle on the same file.
    metadata: Arc<Metadata>,
    /// Scratch space reused from one query to the next.
    buffers: Buffers,
}

impl BigWigReader<BufReader<File>> {
    /// Opens the bigWig file at `path`.
    ///
    /// Fails with [`Error::Io`], naming the path, if the file can't be opened, and otherwise as
    /// [`new`](BigWigReader::new) does.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = BufReader::new(open_file(path)?);
        Self::from_source(file, Some(path.to_path_buf()))
    }

    /// Opens another independent handle on the same file, sharing the already-parsed metadata,
    /// so nothing is parsed again. This is the way to query one file from several threads.
    ///
    /// Fails with [`Error::Io`] if this handle was not created by [`open`](Self::open) (so has no
    /// path), if the file can't be opened, or if its length has changed since it was opened.
    pub fn reopen(&self) -> Result<Self> {
        let Some(path) = &self.metadata.path else {
            return Err(Error::Io(std::io::Error::other(
                "reopen needs a handle created by BigWigReader::open",
            )));
        };
        let file = open_file(path)?;
        let file_len = file.metadata()?.len();
        if file_len != self.metadata.file_len {
            return Err(Error::Io(std::io::Error::other(format!(
                "{} is now {file_len} bytes but was {} when opened",
                path.display(),
                self.metadata.file_len
            ))));
        }
        Ok(self.with_source(BufReader::new(file)))
    }
}

impl<R: Read + Seek> BigWigReader<R> {
    /// Opens a bigWig from any seekable source: a file, a `Cursor<Vec<u8>>`, a remote range
    /// reader.
    ///
    /// Reads the header, zoom headers, total summary, chromosome tree and index headers, and
    /// checks that every offset they hold lies inside the source. Fails with
    /// [`Error::NotBigWig`] when the magic matches in neither byte order, [`Error::IsBigBed`] for
    /// a bigBed, [`Error::UnsupportedVersion`] for a version outside 1 to 4, [`Error::Corrupt`]
    /// for an inconsistent structure (including two chromosomes with the same name, which would
    /// make lookups by name ambiguous), and [`Error::Io`] when the source fails.
    pub fn new(source: R) -> Result<Self> {
        Self::from_source(source, None)
    }

    /// Creates a handle on `source`, which must hold the same file as this one, sharing this
    /// handle's parsed metadata. Nothing is read, and the source is not checked.
    pub fn with_source<S: Read + Seek>(&self, source: S) -> BigWigReader<S> {
        BigWigReader {
            source,
            metadata: Arc::clone(&self.metadata),
            buffers: Buffers::new(self.metadata.info.uncompress_buf_size),
        }
    }

    /// Returns the underlying source.
    pub fn into_inner(self) -> R {
        self.source
    }

    /// The header fields, zoom levels and total summary.
    pub fn info(&self) -> &Info {
        &self.metadata.info
    }

    /// The chromosomes, ordered by `chromId`, so `chroms()[item.chrom_id]` is an item's
    /// chromosome.
    ///
    /// Writers differ: kent's tools number chromosomes in name order, libBigWig and bigtools in
    /// input order. The order of the keys in the file's chromosome tree plays no part.
    pub fn chroms(&self) -> &[ChromInfo] {
        &self.metadata.chroms
    }

    /// The chromosome called `name`, or `None` if the file has none.
    pub fn chrom(&self, name: &str) -> Option<&ChromInfo> {
        let &chrom_id = self.metadata.chrom_ids.get(name)?;
        self.metadata.chrom(chrom_id)
    }

    /// The name of the chromosome with `chrom_id` (e.g. `reader.chrom_name(item.chrom_id)`), or
    /// `None` if there is no such chromosome.
    pub fn chrom_name(&self, chrom_id: u32) -> Option<&str> {
        self.metadata.chrom(chrom_id).map(|chrom| chrom.name.as_str())
    }

    /// Recomputes the total summary from the raw data in one full pass, for files whose header
    /// summary is known to be wrong (libBigWig's `maxVal`, for example). With no data, every
    /// field is 0, as kent's `bedGraphToBigWig` writes it.
    ///
    /// Bases covered by two overlapping items, which only lax writers produce, count twice.
    pub fn compute_total_summary(&mut self) -> Result<Summary> {
        let mut total = RawSummaryAccumulator::unbounded();
        for item in self.all_items()? {
            let item = item?;
            total.add_item(item.start, item.end, item.value);
        }
        Ok(total.summary().unwrap_or(Summary {
            covered_bases: 0.0,
            min: 0.0,
            max: 0.0,
            sum: 0.0,
            sum_squares: 0.0,
        }))
    }

    /// The items overlapping `[start, end)` on `chrom`, in file order and not clipped to the
    /// region (as pyBigWig's `intervals()`). An empty region overlaps nothing.
    ///
    /// fixedStep and variableStep items come back with the extent their step and span give them.
    /// A zero-length item (`start == end`) holds no bases and is skipped. An item that runs past
    /// the end of its chromosome is clipped to it, and one that starts at or past the end is
    /// dropped. Skipping and clipping each log a warning once per file.
    ///
    /// Fails with [`Error::UnknownChrom`] for an unknown chromosome and [`Error::InvalidRegion`]
    /// when `start > end` or `end` is past the chromosome's length. Reading may fail with
    /// [`Error::Corrupt`], [`Error::Decompress`] or [`Error::Io`]; the iterator then yields that
    /// error and stops.
    pub fn items(&mut self, chrom: &str, start: u32, end: u32) -> Result<Items<'_, R>> {
        let region = self.resolve_region(chrom, start, end)?;
        let index = self.metadata.data_index;
        let blocks = self.blocks_in_region(&index, region)?;
        Ok(Items(BlockStream::new(self, blocks, Some(region), Metadata::decode_items)))
    }

    /// Every item of the chromosome `chrom`, in file order, as [`items`](Self::items) over the
    /// whole chromosome.
    pub fn chrom_items(&mut self, chrom: &str) -> Result<Items<'_, R>> {
        let length =
            self.chrom(chrom).ok_or_else(|| Error::UnknownChrom(chrom.to_string()))?.length;
        self.items(chrom, 0, length)
    }

    /// Every item in the file, chromosome by chromosome in `chromId` order and by start within
    /// each, whatever order the data sections are stored in.
    ///
    /// Items are checked, clipped and dropped as [`items`](Self::items) describes, and an item whose
    /// `chromId` is not in the chromosome tree is [`Error::Corrupt`].
    pub fn all_items(&mut self) -> Result<Items<'_, R>> {
        let blocks = self.all_data_blocks()?;
        Ok(Items(BlockStream::new(self, blocks, None, Metadata::decode_items)))
    }

    /// Fills `out[i]` with the value of base `start + i` of `chrom`, or with `missing` where no
    /// item covers it, and returns the number of bases covered. Where items overlap, which only
    /// lax writers produce, the later item in file order wins.
    ///
    /// Once the reader's buffers have grown to fit the query, nothing is allocated. Fails as
    /// [`items`](Self::items) does, and with [`Error::InvalidRegion`] when `out.len()` is not
    /// `end - start`.
    pub fn values_into(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        out: &mut [f32],
        missing: f32,
    ) -> Result<u64> {
        let region = self.resolve_region(chrom, start, end)?;
        if usize::try_from(end - start).ok() != Some(out.len()) {
            return Err(invalid_region(
                chrom,
                start,
                end,
                "the output buffer's length differs from the region's",
            ));
        }
        out.fill(missing);
        if start == end {
            return Ok(0);
        }

        let Self { source, metadata, buffers } = self;
        let Buffers { walk, blocks, read, decompressor, items, covered } = buffers;
        metadata.data_index.overlapping_into(source, region.chrom_id, start, end, walk, blocks)?;
        sort_and_dedupe_blocks(blocks);
        // A bitmap of covered bases, since a stored value may equal `missing`.
        covered.clear();
        covered.resize(out.len().div_ceil(64), 0);
        let mut cursor = BlockCursor::default();
        while let Some(block) =
            cursor.next_block(blocks, source, metadata.file_len, read, decompressor)
        {
            let (bytes, offset) = block?;
            metadata.decode_items(bytes, offset, Some(region), items)?;
            for item in items.iter() {
                let (from, to) = region.relative_overlap(item.start, item.end);
                out[from..to].fill(item.value);
                mark_bits(covered, from, to);
            }
        }
        Ok(covered.iter().map(|word| u64::from(word.count_ones())).sum())
    }

    /// The value of every base of `[start, end)` on `chrom`, with NaN where there is no data (as
    /// pyBigWig's `values()`). Fails as [`values_into`](Self::values_into) does.
    pub fn values(&mut self, chrom: &str, start: u32, end: u32) -> Result<Vec<f32>> {
        self.resolve_region(chrom, start, end)?;
        let length = usize::try_from(end - start)
            .map_err(|_| invalid_region(chrom, start, end, "the region is too long to hold"))?;
        let mut out = vec![f32::NAN; length];
        self.values_into(chrom, start, end, &mut out, f32::NAN)?;
        Ok(out)
    }

    /// Summary of the data in `[start, end)` on `chrom`, or `None` when no base in the region has
    /// data. As [`summaries`](Self::summaries) with one bin.
    pub fn summary(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        source: SummarySource,
    ) -> Result<Option<Summary>> {
        Ok(self.summaries(chrom, start, end, 1, source)?.pop().flatten())
    }

    /// Summaries of `n_bins` equal bins of `[start, end)` on `chrom` (kent's `bigWigSummary`,
    /// pyBigWig's `stats` with `nBins`), each `None` when no base in the bin has data.
    ///
    /// Bin *i* is `[start + ⌊L·i/n⌋, start + ⌊L·(i+1)/n⌋)` with `L = end − start`, as kent and
    /// libBigWig compute bins; with more bins than bases, some bins are empty and `None`. kent's
    /// widening of empty bins and its replication of bases across pixels are display conveniences
    /// that are not reproduced.
    ///
    /// [`SummarySource::Raw`] is exact. [`SummarySource::BestZoom`] uses the zoom level kent's
    /// `bigWigSummary` would (see [`best_zoom_for`](Self::best_zoom_for)), with each zoom record
    /// contributing the fraction of it that overlaps a bin, and falls back to raw data when no
    /// level is fine enough. Zoom-based results are approximate, `min` and `max` included.
    /// Records with no covered bases contribute nothing. Zoom-based `covered_bases` is the
    /// pro-rated count, which may be fractional: kent rounds it up to a whole number and scales
    /// `sum` and `sum_squares` to match, so means agree with kent's but coverage and sums can be
    /// slightly lower.
    ///
    /// Bases covered by two overlapping items, which only lax writers produce, count twice in raw
    /// summaries, as they do in kent's and libBigWig's.
    ///
    /// Fails as [`items`](Self::items) does, and with [`Error::InvalidRegion`] when `n_bins` is
    /// 0, when it exceeds the larger of the region's length and 2²⁰ (1,048,576), or when the
    /// allocator refuses memory for that many bins.
    pub fn summaries(
        &mut self,
        chrom: &str,
        start: u32,
        end: u32,
        n_bins: u32,
        source: SummarySource,
    ) -> Result<Vec<Option<Summary>>> {
        let region = self.resolve_region(chrom, start, end)?;
        if n_bins == 0 {
            return Err(invalid_region(chrom, start, end, "zero bins requested"));
        }
        if n_bins > (end - start).max(MIN_BIN_LIMIT) {
            return Err(invalid_region(
                chrom,
                start,
                end,
                "more bins requested than the larger of the region's length and 2^20",
            ));
        }
        let mut bins = self.bin_storage(region, usize::try_from(n_bins).unwrap_or(usize::MAX))?;
        bins.extend((0..n_bins).map(|index| bin_bounds(start, end, n_bins, index)));
        let level = match source {
            SummarySource::Raw => None,
            SummarySource::BestZoom => self.best_zoom_index((end - start) / n_bins),
        };
        match level {
            None => self.raw_summaries(region, &bins),
            Some(level_index) => self.zoom_summaries(level_index, region, &bins),
        }
    }

    /// The file's zoom levels, in header order (which need not be increasing).
    pub fn zoom_levels(&self) -> &[ZoomLevel] {
        &self.metadata.info.zoom_levels
    }

    /// The records of zoom level `level` overlapping `[start, end)` on `chrom`, in file order and
    /// not clipped to the region. Overlap is half-open, so a record ending at `start` is not
    /// included.
    ///
    /// Fails with [`Error::ForeignZoomLevel`] when `level` is not one of this file's levels,
    /// otherwise as [`items`](Self::items) does; a record whose `chromId` is not in the
    /// chromosome tree is [`Error::Corrupt`].
    pub fn zoom_records(
        &mut self,
        level: ZoomLevel,
        chrom: &str,
        start: u32,
        end: u32,
    ) -> Result<ZoomRecords<'_, R>> {
        let level_index = self.metadata.level_index(level)?;
        let region = self.resolve_region(chrom, start, end)?;
        Ok(ZoomRecords(self.level_records(level_index, region)?))
    }

    /// The zoom level [`SummarySource::BestZoom`] uses for bins `bin_width` bases wide, or `None`
    /// when it would use raw data.
    ///
    /// This is kent's rule: the level with the largest reduction strictly less than
    /// `⌊bin_width / 2⌋`, and none when no level qualifies or `⌊bin_width / 2⌋ ≤ 1`. Levels are
    /// not assumed to be sorted. libBigWig accepts a reduction equal to half the bin width, so the
    /// two can choose different levels.
    pub fn best_zoom_for(&self, bin_width: u32) -> Option<ZoomLevel> {
        self.best_zoom_index(bin_width).map(|index| self.metadata.info.zoom_levels[index])
    }

    /// Calls `visit` with the kind and items of each data section that may hold data in
    /// `region` (a chromosome and a range in it), or of every section when `region` is `None`.
    /// Sections come in the order [`items`](Self::items) and [`all_items`](Self::all_items) read
    /// them, with their items clipped and filtered as those methods describe, so a section's
    /// slice may be empty. The public iterators drop section types, which wig output keeps.
    ///
    /// Fails as `items` does, or with the first error `visit` returns.
    pub(crate) fn for_each_section(
        &mut self,
        region: Option<(&str, u32, u32)>,
        mut visit: impl FnMut(SectionKind, &[DataItem]) -> Result<()>,
    ) -> Result<()> {
        let (blocks, region) = match region {
            Some((chrom, start, end)) => {
                let region = self.resolve_region(chrom, start, end)?;
                let index = self.metadata.data_index;
                (self.blocks_in_region(&index, region)?, Some(region))
            }
            None => (self.all_data_blocks()?, None),
        };
        let Self { source, metadata, buffers } = self;
        let Buffers { read, decompressor, items, .. } = buffers;
        let mut cursor = BlockCursor::default();
        while let Some(block) =
            cursor.next_block(&blocks, source, metadata.file_len, read, decompressor)
        {
            let (bytes, offset) = block?;
            let kind = metadata.decode_section(bytes, offset, region, items)?;
            visit(kind, items)?;
        }
        Ok(())
    }

    /// Opens a handle on `source`, remembering `path` for [`BigWigReader::reopen`].
    fn from_source(mut source: R, path: Option<PathBuf>) -> Result<Self> {
        let metadata = Metadata::read(&mut source, path)?;
        let buffers = Buffers::new(metadata.info.uncompress_buf_size);
        Ok(Self { source, metadata: Arc::new(metadata), buffers })
    }

    /// Looks up `chrom` and checks that `[start, end)` lies within it.
    fn resolve_region(&self, chrom: &str, start: u32, end: u32) -> Result<Region> {
        let unknown = || Error::UnknownChrom(chrom.to_string());
        let &chrom_id = self.metadata.chrom_ids.get(chrom).ok_or_else(unknown)?;
        let length = self.metadata.chrom(chrom_id).ok_or_else(unknown)?.length;
        if start > end {
            return Err(invalid_region(chrom, start, end, "start is after end"));
        }
        if end > length {
            return Err(invalid_region(chrom, start, end, "end is past the chromosome's length"));
        }
        Ok(Region { chrom_id, start, end })
    }

    /// Index into the zoom levels of the level kent's rule picks for bins `bin_width` wide.
    fn best_zoom_index(&self, bin_width: u32) -> Option<usize> {
        best_zoom(self.metadata.info.zoom_levels.iter().map(ZoomLevel::reduction), bin_width)
    }

    /// Exact summaries of `bins`, which partition `region` in order.
    fn raw_summaries(
        &mut self,
        region: Region,
        bins: &[(u32, u32)],
    ) -> Result<Vec<Option<Summary>>> {
        let mut bin_summaries = self.bin_storage(region, bins.len())?;
        bin_summaries
            .extend(bins.iter().map(|&(start, end)| RawSummaryAccumulator::new(start, end)));
        let index = self.metadata.data_index;
        let blocks = self.blocks_in_region(&index, region)?;
        for item in BlockStream::new(self, blocks, Some(region), Metadata::decode_items) {
            let item = item?;
            for bin in &mut bin_summaries[overlapping_bins(bins, item.start, item.end)] {
                bin.add_item(item.start, item.end, item.value);
            }
        }
        let mut summaries = self.bin_storage(region, bins.len())?;
        summaries.extend(bin_summaries.iter().map(RawSummaryAccumulator::summary));
        Ok(summaries)
    }

    /// Summaries of `bins`, which partition `region` in order, from zoom level `level_index`.
    fn zoom_summaries(
        &mut self,
        level_index: usize,
        region: Region,
        bins: &[(u32, u32)],
    ) -> Result<Vec<Option<Summary>>> {
        let mut bin_summaries = self.bin_storage(region, bins.len())?;
        bin_summaries
            .extend(bins.iter().map(|&(start, end)| ZoomSummaryAccumulator::new(start, end)));
        for record in self.level_records(level_index, region)? {
            let record = record?;
            for bin in &mut bin_summaries[overlapping_bins(bins, record.start, record.end)] {
                bin.add_record(&record);
            }
        }
        let mut summaries = self.bin_storage(region, bins.len())?;
        summaries.extend(bin_summaries.iter().map(ZoomSummaryAccumulator::summary));
        Ok(summaries)
    }

    /// An empty vector with room for `n_bins` entries, one per bin of `region`. The count comes
    /// from the caller and may be more than memory holds, so a refused allocation is
    /// [`Error::InvalidRegion`] rather than an abort.
    fn bin_storage<T>(&self, region: Region, n_bins: usize) -> Result<Vec<T>> {
        let mut storage = Vec::new();
        storage.try_reserve_exact(n_bins).map_err(|_| {
            let chrom =
                self.metadata.chrom(region.chrom_id).map_or("", |chrom| chrom.name.as_str());
            invalid_region(chrom, region.start, region.end, "too many bins to allocate")
        })?;
        Ok(storage)
    }

    /// Streams the records of zoom level `level_index` that overlap `region`.
    fn level_records(
        &mut self,
        level_index: usize,
        region: Region,
    ) -> Result<BlockStream<'_, R, ZoomRecord>> {
        let index = self.metadata.zoom_indexes[level_index];
        let blocks = self.blocks_in_region(&index, region)?;
        Ok(BlockStream::new(self, blocks, Some(region), Metadata::decode_records))
    }

    /// The blocks of `index` that may hold data in `region`, in file order, each once (see
    /// [`sort_and_dedupe_blocks`]). An empty region has none.
    fn blocks_in_region(&mut self, index: &RTree, region: Region) -> Result<Vec<LeafEntry>> {
        if region.start == region.end {
            return Ok(Vec::new());
        }
        let Region { chrom_id, start, end } = region;
        let mut blocks = index.overlapping(&mut self.source, chrom_id, start, end)?;
        sort_and_dedupe_blocks(&mut blocks);
        Ok(blocks)
    }

    /// Every data block, each once (see [`sort_and_dedupe_blocks`]), ordered by chromosome, then
    /// start, then file offset.
    fn all_data_blocks(&mut self) -> Result<Vec<LeafEntry>> {
        let mut blocks = self.metadata.data_index.leaves(&mut self.source)?;
        sort_and_dedupe_blocks(&mut blocks);
        blocks.sort_unstable_by_key(|block| {
            (block.start_chrom_ix, block.start_base, block.data_offset)
        });
        Ok(blocks)
    }
}

impl<R> fmt::Debug for BigWigReader<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BigWigReader")
            .field("info", &self.metadata.info)
            .field("chrom_count", &self.metadata.chroms.len())
            .finish_non_exhaustive()
    }
}

/// The header fields, zoom levels and total summary of an open file.
///
/// Fields may be added in later versions, so an `Info` comes only from a reader:
///
/// ```compile_fail,E0639
/// let info = bigwig::Info {
///     version: 4,
///     byte_order: bigwig::ByteOrder::Little,
///     zoom_levels: Vec::new(),
///     total_summary: None,
///     uncompress_buf_size: 0,
///     section_count: 0,
/// };
/// ```
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Info {
    /// Format version, 1 to 4.
    pub version: u16,
    /// Byte order, detected from the magic number.
    pub byte_order: ByteOrder,
    /// The zoom levels, in header order.
    pub zoom_levels: Vec<ZoomLevel>,
    /// The total summary as stored, or `None` when the file has none (`totalSummaryOffset` is
    /// 0). Values are reported as the writer computed them, even when wrong; see
    /// [`BigWigReader::compute_total_summary`].
    pub total_summary: Option<Summary>,
    /// Size of the largest decompressed block as the writer declared it, or 0 when blocks are
    /// stored uncompressed.
    pub uncompress_buf_size: u32,
    /// Number of data sections (not items), from the count before the data. Version 1 files
    /// don't store it consistently, so for them it may be meaningless.
    pub section_count: u64,
}

/// Where [`BigWigReader::summary`] and [`BigWigReader::summaries`] take their data from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SummarySource {
    /// The raw data: exact (pyBigWig's `exact=True`).
    Raw,
    /// A zoom level when one is fine enough for the bin width, otherwise the raw data (pyBigWig's
    /// default and kent's `bigWigSummary`). Faster, but approximate.
    BestZoom,
}

/// Iterator over data items, from [`BigWigReader::items`], [`BigWigReader::chrom_items`] and
/// [`BigWigReader::all_items`].
///
/// It holds one merged read of compressed blocks and one decoded block at a time. After an error
/// it yields nothing more.
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct Items<'a, R>(BlockStream<'a, R, DataItem>);

impl<R: Read + Seek> Iterator for Items<'_, R> {
    type Item = Result<DataItem>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}

impl<R: Read + Seek> FusedIterator for Items<'_, R> {}

/// Iterator over the records of one zoom level, from [`BigWigReader::zoom_records`].
///
/// It holds one merged read of compressed blocks and one decoded block at a time. After an error
/// it yields nothing more.
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct ZoomRecords<'a, R>(BlockStream<'a, R, ZoomRecord>);

impl<R: Read + Seek> Iterator for ZoomRecords<'_, R> {
    type Item = Result<ZoomRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}

impl<R: Read + Seek> FusedIterator for ZoomRecords<'_, R> {}

/// Everything parsed when a file is opened, shared by every handle on it.
#[derive(Debug)]
struct Metadata {
    /// The header fields reported to callers.
    info: Info,
    /// Chromosomes indexed by `chromId`.
    chroms: Vec<ChromInfo>,
    /// `chromId` of every chromosome name.
    chrom_ids: HashMap<String, u32>,
    /// The index over the data sections.
    data_index: RTree,
    /// The index of each zoom level, parallel to `info.zoom_levels`.
    zoom_indexes: Vec<RTree>,
    /// Length of the file; every block read must lie inside it.
    file_len: u64,
    /// Path the file was opened from, for [`BigWigReader::reopen`].
    path: Option<PathBuf>,
    /// Whether the warning about zero-length items has been logged.
    warned_about_empty_items: AtomicBool,
    /// Whether the warning about items running past their chromosome's end has been logged.
    warned_about_clipping: AtomicBool,
}

impl Metadata {
    /// Reads and checks everything needed to query `source`, which was opened from `path`.
    fn read<R: Read + Seek>(source: &mut R, path: Option<PathBuf>) -> Result<Self> {
        let file_len = source.seek(SeekFrom::End(0))?;
        // Files shorter than the header still reach the magic check, so short input that isn't a
        // bigWig, empty input included, is `NotBigWig` rather than `Corrupt`.
        let header_len =
            usize::try_from(file_len).map_or(FIXED_HEADER_SIZE, |len| len.min(FIXED_HEADER_SIZE));
        let (header, byte_order) =
            FixedHeader::decode(&read_exact_at(source, 0, header_len, file_len)?)?;

        let zoom_levels = read_zoom_headers(source, header.zoom_levels, byte_order, file_len)?;
        let total_summary = match header.total_summary_offset {
            0 => None,
            offset => {
                let bytes = read_exact_at(source, offset, TOTAL_SUMMARY_SIZE, file_len)?;
                Some(decode_total_summary(&mut ByteReader::new(&bytes, byte_order, offset))?)
            }
        };
        let section_count = match header.full_data_offset {
            0 => 0,
            offset => {
                let bytes = read_exact_at(source, offset, 8, file_len)?;
                ByteReader::new(&bytes, byte_order, offset).u64()?
            }
        };
        let chroms = bptree::read(source, header.chrom_tree_offset, file_len)?;
        let chrom_ids = index_chrom_names(&chroms, header.chrom_tree_offset)?;
        let data_index = RTree::read(source, header.full_index_offset, file_len)?;
        let zoom_indexes = zoom_levels
            .iter()
            .map(|level| RTree::read(source, level.index_offset(), file_len))
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            info: Info {
                version: header.version,
                byte_order,
                zoom_levels,
                total_summary,
                uncompress_buf_size: header.uncompress_buf_size,
                section_count,
            },
            chroms,
            chrom_ids,
            data_index,
            zoom_indexes,
            file_len,
            path,
            warned_about_empty_items: AtomicBool::new(false),
            warned_about_clipping: AtomicBool::new(false),
        })
    }

    /// The chromosome with `chrom_id`, if there is one.
    fn chrom(&self, chrom_id: u32) -> Option<&ChromInfo> {
        self.chroms.get(usize::try_from(chrom_id).ok()?)
    }

    /// Index into the zoom levels of `level`, or [`Error::ForeignZoomLevel`] if it is not one of
    /// this file's.
    fn level_index(&self, level: ZoomLevel) -> Result<usize> {
        self.info
            .zoom_levels
            .iter()
            .position(|&own| own == level)
            .ok_or(Error::ForeignZoomLevel { reduction: level.reduction() })
    }

    /// A [`DecodeFn`] for data sections: decodes the section in `bytes`, read from `offset`, into
    /// its items, without zero-length ones and clipped to their chromosome (see
    /// [`Metadata::clip_or_drop`]).
    ///
    /// Fails with [`Error::Corrupt`] for a malformed section or a `chromId` that is not in the
    /// chromosome tree.
    fn decode_items(
        &self,
        bytes: &[u8],
        offset: u64,
        region: Option<Region>,
        out: &mut Vec<DataItem>,
    ) -> Result<()> {
        self.decode_section(bytes, offset, region, out).map(|_kind| ())
    }

    /// As [`Metadata::decode_items`], also returning the section's kind.
    fn decode_section(
        &self,
        bytes: &[u8],
        offset: u64,
        region: Option<Region>,
        out: &mut Vec<DataItem>,
    ) -> Result<SectionKind> {
        out.clear();
        let section = Section::parse(bytes, self.info.byte_order, offset)?;
        let chrom_id = section.chrom_id();
        let Some(chrom) = self.chrom(chrom_id) else {
            return Err(Error::corrupt(
                offset,
                format!("section chromId {chrom_id} is not in the chromosome tree"),
            ));
        };
        // A section holds one chromosome, and the index may hand back a neighbour's.
        if region.is_some_and(|region| region.chrom_id != chrom_id) {
            return Ok(section.kind());
        }
        for item in section.items() {
            let Some(item) = self.clip_or_drop(item?, chrom, offset) else { continue };
            if region.is_none_or(|region| region.overlaps(item.start, item.end)) {
                out.push(item);
            }
        }
        Ok(section.kind())
    }

    /// Returns `item`, from the block at `offset`, as a query sees it, or `None` when it holds no
    /// data. Two kinds of item that other writers store are tolerated, each logging a warning once
    /// per file, so that the rest of the file stays readable:
    ///
    /// - A zero-length item (`start == end`) covers no bases and is dropped. kent's
    ///   `bedGraphToBigWig` stores them, rejecting only `start > end`, and kent and bigtools read
    ///   past them.
    /// - An item running past the end of its chromosome `chrom` is shortened to it, and one
    ///   starting at or past the end has nothing left and is dropped. libBigWig writes such items
    ///   unchecked, and bigtools' merge has produced them (bigtools #78).
    fn clip_or_drop(&self, item: DataItem, chrom: &ChromInfo, offset: u64) -> Option<DataItem> {
        if item.start == item.end {
            if !self.warned_about_empty_items.swap(true, Ordering::Relaxed) {
                log::warn!(
                    "item {}:{}-{} in the block at offset {offset} is zero-length; skipping it \
                     and any others (reported once per file)",
                    chrom.name,
                    item.start,
                    item.end
                );
            }
            return None;
        }
        if item.end <= chrom.length {
            return Some(item);
        }
        if !self.warned_about_clipping.swap(true, Ordering::Relaxed) {
            log::warn!(
                "item {}:{}-{} in the block at offset {offset} extends past the chromosome's \
                 length {}; clipping it and any others to the chromosome, dropping those that \
                 start past its end (reported once per file)",
                chrom.name,
                item.start,
                item.end,
                chrom.length
            );
        }
        (item.start < chrom.length).then_some(DataItem { end: chrom.length, ..item })
    }

    /// A [`DecodeFn`] for zoom blocks: decodes the records in `bytes`, read from `offset`.
    /// Bytes after the last whole record are ignored, as kent ignores them (it divides the block
    /// length by the record size) and as trailing bytes after a section's items are.
    ///
    /// Fails with [`Error::Corrupt`] when a record's `chromId` is not in the chromosome tree.
    fn decode_records(
        &self,
        bytes: &[u8],
        offset: u64,
        region: Option<Region>,
        out: &mut Vec<ZoomRecord>,
    ) -> Result<()> {
        out.clear();
        let mut reader = ByteReader::new(bytes, self.info.byte_order, offset);
        while reader.remaining() >= ZoomRecord::ENCODED_LEN {
            let record_offset = reader.file_offset();
            let record = ZoomRecord::decode(&mut reader)?;
            if self.chrom(record.chrom_id).is_none() {
                return Err(Error::corrupt(
                    record_offset,
                    format!(
                        "zoom record chromId {} is not in the chromosome tree",
                        record.chrom_id
                    ),
                ));
            }
            let keep = region.is_none_or(|region| {
                region.chrom_id == record.chrom_id && region.overlaps(record.start, record.end)
            });
            if keep {
                out.push(record);
            }
        }
        Ok(())
    }
}

/// Scratch space a handle reuses from one query to the next.
struct Buffers {
    /// For walking the R-trees.
    walk: WalkBuffers,
    /// The blocks a query reads.
    blocks: Vec<LeafEntry>,
    /// One merged read of (compressed) blocks.
    read: Vec<u8>,
    /// Inflates blocks into a reusable buffer.
    decompressor: BlockDecompressor,
    /// The items of one block.
    items: Vec<DataItem>,
    /// One bit per base of a [`BigWigReader::values_into`] region: whether it has data.
    covered: Vec<u64>,
}

impl Buffers {
    /// Creates empty buffers, apart from an inflate buffer of `uncompress_buf_size` bytes.
    fn new(uncompress_buf_size: u32) -> Self {
        Self {
            walk: WalkBuffers::default(),
            blocks: Vec::new(),
            read: Vec::new(),
            decompressor: BlockDecompressor::new(uncompress_buf_size),
            items: Vec::new(),
            covered: Vec::new(),
        }
    }
}

/// A streaming iterator over what a list of blocks decodes to: data items or zoom records. It
/// borrows its handle's source and buffers, and holds one decoded block of its own.
struct BlockStream<'a, R, T> {
    /// The handle whose source and buffers are used.
    reader: &'a mut BigWigReader<R>,
    /// The blocks to read, in order; emptied after an error so nothing more is read.
    blocks: Vec<LeafEntry>,
    /// Position in `blocks`.
    cursor: BlockCursor,
    /// The region to keep, or `None` for everything.
    region: Option<Region>,
    /// Turns one decompressed block into values.
    decode: DecodeFn<T>,
    /// The current block's values.
    decoded: Vec<T>,
    /// Index of the next value of `decoded` to yield.
    position: usize,
    /// An error to yield once the values decoded before it are exhausted.
    pending_error: Option<Error>,
}

impl<'a, R: Read + Seek, T: Copy> BlockStream<'a, R, T> {
    /// Creates a stream over `blocks`, keeping only what `decode` finds inside `region`.
    fn new(
        reader: &'a mut BigWigReader<R>,
        blocks: Vec<LeafEntry>,
        region: Option<Region>,
        decode: DecodeFn<T>,
    ) -> Self {
        Self {
            reader,
            blocks,
            cursor: BlockCursor::default(),
            region,
            decode,
            decoded: Vec::new(),
            position: 0,
            pending_error: None,
        }
    }
}

impl<R: Read + Seek, T: Copy> Iterator for BlockStream<'_, R, T> {
    type Item = Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(&value) = self.decoded.get(self.position) {
                self.position += 1;
                return Some(Ok(value));
            }
            if let Some(error) = self.pending_error.take() {
                return Some(Err(error));
            }

            let Self { reader, blocks, cursor, region, decode, decoded, position, pending_error } =
                self;
            let BigWigReader { source, metadata, buffers } = &mut **reader;
            let Buffers { read, decompressor, .. } = buffers;
            let outcome =
                cursor.next_block(blocks, source, metadata.file_len, read, decompressor)?;
            decoded.clear();
            *position = 0;
            let decoded_block = outcome
                .and_then(|(bytes, offset)| decode(metadata, bytes, offset, *region, decoded));
            if let Err(error) = decoded_block {
                *pending_error = Some(error);
                blocks.clear();
            }
        }
    }
}

/// Position in a list of blocks being read in order, where runs of nearby blocks are fetched with
/// one read.
#[derive(Debug, Default)]
struct BlockCursor {
    /// Index of the next block.
    next: usize,
    /// Blocks before this index are in the current read.
    read_end: usize,
    /// File offset of the first byte of the current read.
    read_offset: u64,
}

impl BlockCursor {
    /// The next block's decompressed bytes and file offset, or `None` after the last block.
    ///
    /// When the next block is not in the current read, it and the blocks after it that are
    /// contiguous or separated by at most [`MAX_MERGE_GAP`] bytes are read into `read_buf`, up to
    /// [`MAX_MERGED_READ`] bytes. Every block must lie inside the file, as
    /// [`RTree`] leaf entries are checked to.
    fn next_block<'b, R: Read + Seek>(
        &mut self,
        blocks: &[LeafEntry],
        source: &mut R,
        file_len: u64,
        read_buf: &'b mut Vec<u8>,
        decompressor: &'b mut BlockDecompressor,
    ) -> Option<Result<(&'b [u8], u64)>> {
        let block = *blocks.get(self.next)?;
        if self.next >= self.read_end {
            if let Err(error) = self.read_run(blocks, source, file_len, read_buf) {
                self.next = blocks.len();
                return Some(Err(error));
            }
        }
        self.next += 1;
        // The run holds the whole block, and its length fits in usize, so these can't fail.
        let block_bytes = usize::try_from(block.data_offset - self.read_offset)
            .ok()
            .zip(usize::try_from(block.data_size).ok())
            .and_then(|(start, len)| read_buf.get(start..start.checked_add(len)?));
        let Some(block_bytes) = block_bytes else {
            return Some(Err(Error::corrupt(block.data_offset, "block lies outside its read")));
        };
        Some(
            decompressor
                .decompress(block_bytes, block.data_offset)
                .map(|bytes| (bytes, block.data_offset)),
        )
    }

    /// Reads the run of blocks starting at `blocks[self.next]` into `read_buf`.
    fn read_run<R: Read + Seek>(
        &mut self,
        blocks: &[LeafEntry],
        source: &mut R,
        file_len: u64,
        read_buf: &mut Vec<u8>,
    ) -> Result<()> {
        let first = blocks[self.next];
        // Leaf entries are checked to lie inside the file, so these sums can't overflow.
        let mut run_end = first.data_offset + first.data_size;
        let mut run_blocks = 1;
        for block in &blocks[self.next + 1..] {
            let block_end = block.data_offset + block.data_size;
            let mergeable = block.data_offset >= run_end
                && block.data_offset - run_end <= MAX_MERGE_GAP
                && block_end - first.data_offset <= MAX_MERGED_READ;
            if !mergeable {
                break;
            }
            run_end = block_end;
            run_blocks += 1;
        }
        let run_len = usize::try_from(run_end - first.data_offset).map_err(|_| {
            Error::corrupt(first.data_offset, "block is too large to read on this platform")
        })?;
        read_exact_at_into(source, first.data_offset, run_len, file_len, read_buf)?;
        self.read_offset = first.data_offset;
        self.read_end = self.next + run_blocks;
        Ok(())
    }
}

/// A query region on one chromosome, already checked against the chromosome's length.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Region {
    /// Chromosome of the region.
    chrom_id: u32,
    /// First base of the region.
    start: u32,
    /// One past the last base of the region.
    end: u32,
}

impl Region {
    /// Whether `[start, end)` shares at least one base with the region.
    fn overlaps(self, start: u32, end: u32) -> bool {
        start < self.end && end > self.start
    }

    /// The part of `[start, end)` inside the region, as offsets from the region's start. The two
    /// must overlap.
    fn relative_overlap(self, start: u32, end: u32) -> (usize, usize) {
        let offset = |position: u32| {
            // The region's length fits in usize (it is the length of the caller's buffer).
            usize::try_from(position - self.start).unwrap_or(usize::MAX)
        };
        (offset(start.max(self.start)), offset(end.min(self.end)))
    }
}

/// Opens the file at `path` for reading. A failure is [`Error::Io`] with the error's kind kept and
/// the path added to its message, since a bare "No such file or directory" doesn't say which.
fn open_file(path: &Path) -> Result<File> {
    File::open(path).map_err(|error| {
        Error::Io(std::io::Error::new(error.kind(), format!("{}: {error}", path.display())))
    })
}

/// Reads `count` zoom headers, which follow the fixed header, and checks that each level's data
/// offset lies inside the file (its index is read separately, which checks its offset).
fn read_zoom_headers<R: Read + Seek>(
    source: &mut R,
    count: u16,
    byte_order: ByteOrder,
    file_len: u64,
) -> Result<Vec<ZoomLevel>> {
    let offset = FIXED_HEADER_SIZE as u64;
    let bytes = read_exact_at(source, offset, usize::from(count) * ZOOM_HEADER_SIZE, file_len)?;
    let mut reader = ByteReader::new(&bytes, byte_order, offset);
    let mut levels = Vec::with_capacity(usize::from(count));
    for _ in 0..count {
        let header_offset = reader.file_offset();
        let level = ZoomLevel::decode(&mut reader)?;
        if level.data_offset() > file_len {
            return Err(Error::corrupt(
                header_offset + 8,
                format!(
                    "zoom data offset {} is past the end of the {file_len}-byte file",
                    level.data_offset()
                ),
            ));
        }
        levels.push(level);
    }
    Ok(levels)
}

/// Maps each chromosome name to its `chromId` (its index in `chroms`). A name used twice is
/// [`Error::Corrupt`] at the tree's offset, since a lookup by that name would see only part of the
/// data.
fn index_chrom_names(chroms: &[ChromInfo], tree_offset: u64) -> Result<HashMap<String, u32>> {
    let mut chrom_ids = HashMap::with_capacity(chroms.len());
    for (chrom_id, chrom) in (0u32..).zip(chroms) {
        match chrom_ids.entry(chrom.name.clone()) {
            Entry::Vacant(entry) => {
                entry.insert(chrom_id);
            }
            Entry::Occupied(entry) => {
                return Err(Error::corrupt(
                    tree_offset,
                    format!(
                        "chromosome name {:?} is used by chromIds {} and {chrom_id}",
                        chrom.name,
                        entry.get()
                    ),
                ));
            }
        }
    }
    Ok(chrom_ids)
}

/// Sorts `blocks` by file offset and drops each block whose bytes lie wholly inside a block kept
/// before it, so a block named by several leaf entries is kept once. Writers never store blocks
/// that overlap, but a hostile index can name one large block thousands of times, and each copy
/// would otherwise be read and inflated again.
fn sort_and_dedupe_blocks(blocks: &mut Vec<LeafEntry>) {
    // The largest block at each offset comes first, so it is the one kept.
    blocks.sort_unstable_by_key(|block| (block.data_offset, Reverse(block.data_size)));
    let mut kept_end: Option<u64> = None;
    blocks.retain(|block| {
        // Leaf entries are checked to lie inside the file, so this can't overflow.
        let end = block.data_offset + block.data_size;
        // Every kept block starts at or before this one, so ending within the furthest kept end
        // means lying inside the block that reaches it.
        let inside_a_kept_block = kept_end.is_some_and(|kept_end| end <= kept_end);
        if !inside_a_kept_block {
            kept_end = Some(end);
        }
        !inside_a_kept_block
    });
}

/// The range of `bins` (in order, partitioning a region) that `[start, end)` overlaps. Empty bins
/// inside the range overlap nothing, which the accumulators handle.
fn overlapping_bins(bins: &[(u32, u32)], start: u32, end: u32) -> Range<usize> {
    let first = bins.partition_point(|&(_, bin_end)| bin_end <= start);
    let last = bins.partition_point(|&(bin_start, _)| bin_start < end);
    first..last.max(first)
}

/// Sets bits `[from, to)` of the bitmap `bits` (bit *i* is bit `i % 64` of word `i / 64`).
fn mark_bits(bits: &mut [u64], from: usize, to: usize) {
    if from >= to {
        return;
    }
    let (first_word, last_word) = (from / 64, (to - 1) / 64);
    let first_mask = u64::MAX << (from % 64);
    let last_mask = u64::MAX >> (63 - (to - 1) % 64);
    if first_word == last_word {
        bits[first_word] |= first_mask & last_mask;
    } else {
        bits[first_word] |= first_mask;
        bits[first_word + 1..last_word].fill(u64::MAX);
        bits[last_word] |= last_mask;
    }
}

/// Builds an [`Error::InvalidRegion`].
fn invalid_region(chrom: &str, start: u32, end: u32, reason: &'static str) -> Error {
    Error::InvalidRegion {
        chrom: chrom.to_string(),
        start: u64::from(start),
        end: u64::from(end),
        reason,
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::bbi::{BIGWIG_MAGIC, BPTREE_MAGIC, encode_total_summary};
    use crate::bptree::ChromEntry;
    use crate::compress::{BlockCompressor, Compression};
    use crate::rtree;
    use crate::section::SectionBuilder;
    use crate::zoom::ZoomBuilder;

    /// Size of an R-tree index header.
    const RTREE_HEADER_SIZE: usize = 48;
    /// A summary with every field 0.
    const ZERO_SUMMARY: Summary =
        Summary { covered_bases: 0.0, min: 0.0, max: 0.0, sum: 0.0, sum_squares: 0.0 };
    /// The sample's chromosomes: `chr3` has no data.
    const SAMPLE_CHROMS: [(&str, u32); 3] = [("chr1", 1000), ("chr2", 500), ("chr3", 100)];

    // ---- Building test files -----------------------------------------------------------------

    /// One data section of a test file.
    #[derive(Clone, Debug)]
    enum TestSection {
        /// Items `(start, end, value)`.
        BedGraph { chrom_id: u32, items: Vec<(u32, u32, f32)> },
        /// Items `(start, value)`, each covering `span` bases.
        VariableStep { chrom_id: u32, span: u32, items: Vec<(u32, f32)> },
        /// Value *i* covers `[start + i·step, start + i·step + span)`.
        FixedStep { chrom_id: u32, start: u32, step: u32, span: u32, values: Vec<f32> },
    }

    impl TestSection {
        /// Encodes the section, returning its bytes and its extent `(chromId, start, end)`.
        fn encode(&self, order: ByteOrder) -> (Vec<u8>, u32, u32, u32) {
            let mut builder = SectionBuilder::new(order);
            match self {
                Self::BedGraph { chrom_id, items } => {
                    builder.begin_bed_graph(*chrom_id);
                    for &(start, end, value) in items {
                        builder.push_bed_graph(start, end, value);
                    }
                }
                Self::VariableStep { chrom_id, span, items } => {
                    builder.begin_variable_step(*chrom_id, *span);
                    for &(start, value) in items {
                        builder.push_variable_step(start, value);
                    }
                }
                Self::FixedStep { chrom_id, start, step, span, values } => {
                    builder.begin_fixed_step(*chrom_id, *start, *step, *span);
                    for &value in values {
                        builder.push_fixed_step(value);
                    }
                }
            }
            let bytes = builder.finish().to_vec();
            (bytes, builder.chrom_id(), builder.start(), builder.end())
        }

        /// The items the section holds, with stepped extents expanded.
        fn items(&self) -> Vec<DataItem> {
            match self {
                Self::BedGraph { chrom_id, items } => items
                    .iter()
                    .map(|&(start, end, value)| item(*chrom_id, start, end, value))
                    .collect(),
                Self::VariableStep { chrom_id, span, items } => items
                    .iter()
                    .map(|&(start, value)| item(*chrom_id, start, start + span, value))
                    .collect(),
                Self::FixedStep { chrom_id, start, step, span, values } => (0u32..)
                    .zip(values)
                    .map(|(index, &value)| {
                        let item_start = start + index * step;
                        item(*chrom_id, item_start, item_start + span, value)
                    })
                    .collect(),
            }
        }
    }

    /// A bedGraph section on `chrom_id`.
    fn bed_graph(chrom_id: u32, items: &[(u32, u32, f32)]) -> TestSection {
        TestSection::BedGraph { chrom_id, items: items.to_vec() }
    }

    /// A variableStep section on `chrom_id`.
    fn variable_step(chrom_id: u32, span: u32, items: &[(u32, f32)]) -> TestSection {
        TestSection::VariableStep { chrom_id, span, items: items.to_vec() }
    }

    /// A fixedStep section on `chrom_id`.
    fn fixed_step(chrom_id: u32, start: u32, step: u32, span: u32, values: &[f32]) -> TestSection {
        TestSection::FixedStep { chrom_id, start, step, span, values: values.to_vec() }
    }

    /// How the chromosome tree is laid out.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum ChromTreeStyle {
        /// Keys byte-sorted, as kent and this crate write them.
        Sorted,
        /// One leaf with keys in `chromId` order, unsorted, as libBigWig writes it.
        Unsorted,
    }

    /// How the data index is written.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum DataIndexStyle {
        /// A complete R-tree.
        Full,
        /// No index: `fullIndexOffset` is 0, as libBigWig writes empty files.
        OffsetZero,
        /// The index header alone, with no root node, as kent writes empty files.
        HeaderOnly,
    }

    /// What the total summary holds.
    #[derive(Clone, Copy, Debug, PartialEq)]
    enum TotalSummaryStyle {
        /// Computed from the sections; all zeros without data.
        FromData,
        /// The given values.
        Given(Summary),
        /// None at all: `totalSummaryOffset` is 0.
        Absent,
    }

    /// Where the chromosome tree goes and how zoom data starts.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Layout {
        /// The chromosome tree before the data, and each zoom level's record count before its
        /// blocks.
        Kent,
        /// The chromosome tree after the data, and zoom blocks with no count before them.
        Bigtools,
    }

    /// The pieces of a bigWig file for a test to vary, assembled by [`TestFile::build`] with the
    /// crate's own encoders.
    #[derive(Clone, Debug)]
    struct TestFile {
        /// Chromosomes as (name, length), indexed by `chromId`.
        chroms: Vec<(String, u32)>,
        /// Data sections, in file order.
        sections: Vec<TestSection>,
        /// Reductions (strictly increasing) of zoom levels computed from the sections.
        zoom_reductions: Vec<u32>,
        /// Zoom levels given record by record, after the computed ones.
        zoom_levels: Vec<(u32, Vec<ZoomRecord>)>,
        byte_order: ByteOrder,
        version: u16,
        /// Blocks are zlib streams rather than raw bytes.
        compressed: bool,
        /// `uncompressBufSize` to declare instead of the largest decompressed block.
        uncompress_buf_size: Option<u32>,
        /// Children per R-tree node.
        rtree_block_size: u32,
        /// Records per zoom block.
        zoom_block_records: usize,
        chrom_tree: ChromTreeStyle,
        data_index: DataIndexStyle,
        total_summary: TotalSummaryStyle,
        layout: Layout,
        /// Zoom-header slots to reserve; `bedGraphToBigWig` reserves 10 whatever it writes.
        zoom_header_slots: usize,
        /// Added to every section's header `end` and R-tree bound (libBigWig's overstated end).
        section_end_excess: u32,
        /// Zero bytes written after each data section.
        section_gap: usize,
        trailing_magic: bool,
        /// The data starts at least this far into the file; the gap reads as zeros.
        data_start: u64,
        /// How many times each block's leaf entry is repeated in its R-tree; 1 in a sound file.
        leaf_entry_copies: usize,
    }

    impl TestFile {
        /// A little-endian, compressed, version 4 file in kent's layout, with `sections` on
        /// `chroms` and no zoom levels.
        fn new(chroms: &[(&str, u32)], sections: Vec<TestSection>) -> Self {
            Self {
                chroms: chroms.iter().map(|&(name, length)| (name.to_string(), length)).collect(),
                sections,
                zoom_reductions: Vec::new(),
                zoom_levels: Vec::new(),
                byte_order: ByteOrder::Little,
                version: 4,
                compressed: true,
                uncompress_buf_size: None,
                rtree_block_size: 256,
                zoom_block_records: 4,
                chrom_tree: ChromTreeStyle::Sorted,
                data_index: DataIndexStyle::Full,
                total_summary: TotalSummaryStyle::FromData,
                layout: Layout::Kent,
                zoom_header_slots: 0,
                section_end_excess: 0,
                section_gap: 0,
                trailing_magic: true,
                data_start: 0,
                leaf_entry_copies: 1,
            }
        }

        /// The whole file.
        fn build(&self) -> Vec<u8> {
            let (mut file, tail_offset, tail) = self.build_parts();
            file.resize(usize::try_from(tail_offset).unwrap(), 0);
            file.extend(tail);
            file
        }

        /// A handle on the whole file.
        fn open(&self) -> BigWigReader<Cursor<Vec<u8>>> {
            BigWigReader::new(Cursor::new(self.build())).unwrap()
        }

        /// Every item of the sections, in file order.
        fn items(&self) -> Vec<DataItem> {
            self.sections.iter().flat_map(TestSection::items).collect()
        }

        /// The file as the bytes before the data, the offset of the data region, and the bytes
        /// from there to the end. The two parts are adjacent unless `data_start` leaves a gap.
        fn build_parts(&self) -> (Vec<u8>, u64, Vec<u8>) {
            let order = self.byte_order;
            let zoom_levels = self.all_zoom_levels();
            let zoom_header_slots = self.zoom_header_slots.max(zoom_levels.len());
            let mut head = vec![0; FIXED_HEADER_SIZE + zoom_header_slots * ZOOM_HEADER_SIZE];
            let total_summary_offset = match self.total_summary {
                TotalSummaryStyle::Absent => 0,
                TotalSummaryStyle::FromData => {
                    let offset = head.len() as u64;
                    encode_total_summary(&self.summary_of_data(), order, &mut head);
                    offset
                }
                TotalSummaryStyle::Given(summary) => {
                    let offset = head.len() as u64;
                    encode_total_summary(&summary, order, &mut head);
                    offset
                }
            };
            let mut chrom_tree_offset = head.len() as u64;
            if self.layout == Layout::Kent {
                head.extend(self.chrom_tree(chrom_tree_offset));
            }

            let compression =
                if self.compressed { Compression::Deflate { level: 6 } } else { Compression::None };
            let mut tail = Tail {
                offset: (head.len() as u64).max(self.data_start),
                bytes: Vec::new(),
                compressor: BlockCompressor::new(compression).unwrap(),
                byte_order: order,
                largest_block: 0,
            };
            let full_data_offset = tail.position();
            let entries = self.put_sections(&mut tail);
            if self.layout == Layout::Bigtools {
                chrom_tree_offset = tail.position();
                let tree = self.chrom_tree(chrom_tree_offset);
                tail.bytes.extend(tree);
            }
            let full_index_offset = match self.data_index {
                DataIndexStyle::OffsetZero => 0,
                DataIndexStyle::Full => {
                    tail.put_index(&self.leaf_entries(&entries), self.rtree_block_size, 1)
                }
                DataIndexStyle::HeaderOnly => {
                    let offset = tail.put_index(&[], self.rtree_block_size, 1);
                    let header_end =
                        usize::try_from(offset - tail.offset).unwrap() + RTREE_HEADER_SIZE;
                    tail.bytes.truncate(header_end);
                    offset
                }
            };
            let levels: Vec<ZoomLevel> = zoom_levels
                .iter()
                .map(|(reduction, records)| self.put_zoom_level(&mut tail, *reduction, records))
                .collect();
            if self.trailing_magic {
                order.put_u32(&mut tail.bytes, BIGWIG_MAGIC);
            }

            let header = FixedHeader {
                version: self.version,
                zoom_levels: u16::try_from(levels.len()).unwrap(),
                chrom_tree_offset,
                full_data_offset,
                full_index_offset,
                field_count: 0,
                defined_field_count: 0,
                auto_sql_offset: 0,
                total_summary_offset,
                uncompress_buf_size: self.uncompress_buf_size.unwrap_or(if self.compressed {
                    u32::try_from(tail.largest_block).unwrap()
                } else {
                    0
                }),
                extension_offset: 0,
            };
            let mut fixed = Vec::new();
            header.encode(order, &mut fixed);
            for level in &levels {
                level.encode(order, &mut fixed);
            }
            head[..fixed.len()].copy_from_slice(&fixed);
            (head, tail.offset, tail.bytes)
        }

        /// Appends the section count and the sections to `tail`, returning their leaf entries.
        fn put_sections(&self, tail: &mut Tail) -> Vec<LeafEntry> {
            let order = self.byte_order;
            order.put_u64(&mut tail.bytes, self.sections.len() as u64);
            let mut entries = Vec::new();
            for section in &self.sections {
                let (mut bytes, chrom_id, start, end) = section.encode(order);
                let end = end + self.section_end_excess;
                bytes[8..12].copy_from_slice(&order.u32_bytes(end));
                entries.push(tail.put_block(&bytes, (chrom_id, start), (chrom_id, end)));
                tail.bytes.resize(tail.bytes.len() + self.section_gap, 0);
            }
            entries
        }

        /// Appends a zoom level's data and index to `tail` and returns its header.
        fn put_zoom_level(
            &self,
            tail: &mut Tail,
            reduction: u32,
            records: &[ZoomRecord],
        ) -> ZoomLevel {
            let data_offset = tail.position();
            if self.layout == Layout::Kent {
                self.byte_order.put_u32(&mut tail.bytes, u32::try_from(records.len()).unwrap());
            }
            let mut entries = Vec::new();
            for block in records.chunks(self.zoom_block_records) {
                let mut raw = Vec::new();
                for record in block {
                    record.encode(self.byte_order, &mut raw);
                }
                let (first, last) = (block[0], block[block.len() - 1]);
                entries.push(tail.put_block(
                    &raw,
                    (first.chrom_id, first.start),
                    (last.chrom_id, last.end),
                ));
            }
            let records_per_block = u32::try_from(self.zoom_block_records).unwrap();
            let entries = self.leaf_entries(&entries);
            let index_offset = tail.put_index(&entries, self.rtree_block_size, records_per_block);
            ZoomLevel::new(reduction, data_offset, index_offset)
        }

        /// The leaf entries an R-tree over `blocks` holds: each repeated `leaf_entry_copies`
        /// times.
        fn leaf_entries(&self, blocks: &[LeafEntry]) -> Vec<LeafEntry> {
            blocks
                .iter()
                .flat_map(|&block| std::iter::repeat_n(block, self.leaf_entry_copies))
                .collect()
        }

        /// The computed zoom levels followed by the given ones.
        fn all_zoom_levels(&self) -> Vec<(u32, Vec<ZoomRecord>)> {
            let mut levels: Vec<(u32, Vec<ZoomRecord>)> =
                self.zoom_reductions.iter().map(|&reduction| (reduction, Vec::new())).collect();
            if !levels.is_empty() {
                let mut builder = ZoomBuilder::new(&self.zoom_reductions);
                let mut emit = |level: usize, record: ZoomRecord| {
                    levels[level].1.push(record);
                    Ok(())
                };
                for item in self.items() {
                    let DataItem { chrom_id, start, end, value } = item;
                    builder.push_item(chrom_id, start, end, value, &mut emit).unwrap();
                }
                builder.close(&mut emit).unwrap();
            }
            levels.extend(self.zoom_levels.iter().cloned());
            levels
        }

        /// The exact summary of every item, or all zeros without data.
        fn summary_of_data(&self) -> Summary {
            let mut total = RawSummaryAccumulator::unbounded();
            for item in self.items() {
                total.add_item(item.start, item.end, item.value);
            }
            total.summary().unwrap_or(ZERO_SUMMARY)
        }

        /// The chromosome tree, to be written at `offset`.
        fn chrom_tree(&self, offset: u64) -> Vec<u8> {
            let entries: Vec<ChromEntry<'_>> = self
                .chroms
                .iter()
                .zip(0u32..)
                .map(|((name, length), chrom_id)| ChromEntry { name, chrom_id, length: *length })
                .collect();
            match self.chrom_tree {
                ChromTreeStyle::Sorted => bptree::encode(&entries, offset, self.byte_order).bytes,
                ChromTreeStyle::Unsorted => unsorted_chrom_tree(&entries, self.byte_order),
            }
        }
    }

    /// The part of a test file from the data on, appended block by block.
    struct Tail {
        /// File offset of the first byte.
        offset: u64,
        bytes: Vec<u8>,
        compressor: BlockCompressor,
        byte_order: ByteOrder,
        /// Size of the largest block before compression, for `uncompressBufSize`.
        largest_block: usize,
    }

    impl Tail {
        /// File offset of the next byte.
        fn position(&self) -> u64 {
            self.offset + self.bytes.len() as u64
        }

        /// Appends `raw` as one stored block, returning its leaf entry with bounds `start` to
        /// `end`, each a (chromId, base) pair.
        fn put_block(&mut self, raw: &[u8], start: (u32, u32), end: (u32, u32)) -> LeafEntry {
            self.largest_block = self.largest_block.max(raw.len());
            let mut stored = Vec::new();
            self.compressor.compress(raw, &mut stored).unwrap();
            let entry = LeafEntry {
                start_chrom_ix: start.0,
                start_base: start.1,
                end_chrom_ix: end.0,
                end_base: end.1,
                data_offset: self.position(),
                data_size: stored.len() as u64,
            };
            self.bytes.extend(stored);
            entry
        }

        /// Appends an R-tree over `entries` and returns its offset.
        fn put_index(
            &mut self,
            entries: &[LeafEntry],
            block_size: u32,
            items_per_slot: u32,
        ) -> u64 {
            let offset = self.position();
            let order = self.byte_order;
            self.bytes.extend(rtree::encode(entries, offset, block_size, items_per_slot, order));
            offset
        }
    }

    /// libBigWig's chromosome tree: one leaf holding `entries` in the order given, unsorted.
    fn unsorted_chrom_tree(entries: &[ChromEntry<'_>], order: ByteOrder) -> Vec<u8> {
        let key_size = entries.iter().map(|entry| entry.name.len()).max().unwrap_or(0);
        let count = u32::try_from(entries.len()).unwrap();
        let mut tree = Vec::new();
        order.put_u32(&mut tree, BPTREE_MAGIC);
        order.put_u32(&mut tree, count);
        order.put_u32(&mut tree, u32::try_from(key_size).unwrap());
        order.put_u32(&mut tree, 8);
        order.put_u64(&mut tree, u64::from(count));
        order.put_u64(&mut tree, 0);
        if entries.is_empty() {
            return tree;
        }
        tree.extend([1, 0]);
        order.put_u16(&mut tree, u16::try_from(count).unwrap());
        for entry in entries {
            tree.extend(entry.name.as_bytes());
            tree.resize(tree.len() + key_size - entry.name.len(), 0);
            order.put_u32(&mut tree, entry.chrom_id);
            order.put_u32(&mut tree, entry.length);
        }
        tree
    }

    /// A `Read + Seek` source of `head`, then zeros up to `tail_offset`, then `tail`, which never
    /// materialises the zeros, so a test can place data past 4 GiB.
    struct SparseSource {
        head: Vec<u8>,
        tail_offset: u64,
        tail: Vec<u8>,
        position: u64,
    }

    impl SparseSource {
        /// Total length of the source.
        fn len(&self) -> u64 {
            self.tail_offset + self.tail.len() as u64
        }
    }

    impl Read for SparseSource {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let copy = |from: &[u8], at: u64, buf: &mut [u8]| {
                let at = usize::try_from(at).unwrap();
                let len = buf.len().min(from.len() - at);
                buf[..len].copy_from_slice(&from[at..at + len]);
                len
            };
            let len = if self.position >= self.len() {
                0
            } else if self.position < self.head.len() as u64 {
                copy(&self.head, self.position, buf)
            } else if self.position < self.tail_offset {
                let gap = usize::try_from(self.tail_offset - self.position).unwrap_or(usize::MAX);
                let len = buf.len().min(gap);
                buf[..len].fill(0);
                len
            } else {
                copy(&self.tail, self.position - self.tail_offset, buf)
            };
            self.position += len as u64;
            Ok(len)
        }
    }

    impl Seek for SparseSource {
        fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
            let target = match position {
                SeekFrom::Start(offset) => Some(offset),
                SeekFrom::End(delta) => self.len().checked_add_signed(delta),
                SeekFrom::Current(delta) => self.position.checked_add_signed(delta),
            };
            self.position = target.ok_or_else(|| std::io::Error::other("seek out of range"))?;
            Ok(self.position)
        }
    }

    /// A `Read + Seek` source that logs the file range of every read from it.
    struct ReadLog<R> {
        inner: R,
        position: u64,
        reads: Vec<Range<u64>>,
    }

    impl<R> ReadLog<R> {
        /// A source reading from `inner`, with an empty log.
        fn new(inner: R) -> Self {
            Self { inner, position: 0, reads: Vec::new() }
        }
    }

    impl<R: Read> Read for ReadLog<R> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let len = self.inner.read(buf)?;
            self.reads.push(self.position..self.position + len as u64);
            self.position += len as u64;
            Ok(len)
        }
    }

    impl<R: Seek> Seek for ReadLog<R> {
        fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
            self.position = self.inner.seek(position)?;
            Ok(self.position)
        }
    }

    // ---- Expectations ------------------------------------------------------------------------

    /// Sections of every type on `chr1` (chromId 0), and a bedGraph section on `chr2`. Values
    /// are multiples of 0.5, so every sum is exact whatever the order of addition.
    fn sample_sections() -> Vec<TestSection> {
        vec![
            bed_graph(0, &[(10, 20, 1.0), (20, 25, 2.0), (40, 60, -3.5)]),
            variable_step(0, 5, &[(100, 4.0), (110, 5.0), (200, 6.5)]),
            fixed_step(0, 300, 10, 5, &[7.0, 8.0, 9.0, 10.0]),
            bed_graph(1, &[(0, 100, 0.5), (400, 500, 1.5)]),
        ]
    }

    /// The sample as a file.
    fn sample_file() -> TestFile {
        TestFile::new(&SAMPLE_CHROMS, sample_sections())
    }

    /// The item `[start, end)` on `chrom_id` with `value`.
    fn item(chrom_id: u32, start: u32, end: u32, value: f32) -> DataItem {
        DataItem { chrom_id, start, end, value }
    }

    /// Collects an iterator of results, failing the test on the first error.
    fn collect<T>(iter: impl Iterator<Item = Result<T>>) -> Vec<T> {
        iter.collect::<Result<Vec<T>>>().unwrap()
    }

    /// The items of `file` on `chrom_id` that overlap `[start, end)`, in file order.
    fn expected_items(file: &TestFile, chrom_id: u32, start: u32, end: u32) -> Vec<DataItem> {
        file.items()
            .into_iter()
            .filter(|item| item.chrom_id == chrom_id && item.start < end && item.end > start)
            .collect()
    }

    /// The value of every base of `[start, end)` on `chrom_id`, NaN where no item covers it, with
    /// later items overwriting earlier ones.
    fn dense_values(file: &TestFile, chrom_id: u32, start: u32, end: u32) -> Vec<f32> {
        let mut values = vec![f32::NAN; usize::try_from(end - start).unwrap()];
        for item in file.items().into_iter().filter(|item| item.chrom_id == chrom_id) {
            for position in item.start.max(start)..item.end.min(end) {
                values[usize::try_from(position - start).unwrap()] = item.value;
            }
        }
        values
    }

    /// The exact summary of `values`, skipping NaN, or `None` if every value is NaN.
    #[allow(clippy::cast_precision_loss, reason = "test regions are far below 2^53 bases")]
    fn dense_summary(values: &[f32]) -> Option<Summary> {
        let covered: Vec<f64> =
            values.iter().filter(|value| !value.is_nan()).map(|&value| f64::from(value)).collect();
        (!covered.is_empty()).then(|| Summary {
            covered_bases: covered.len() as f64,
            min: covered.iter().copied().fold(f64::INFINITY, f64::min),
            max: covered.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            sum: covered.iter().sum(),
            sum_squares: covered.iter().map(|value| value * value).sum(),
        })
    }

    /// Asserts that `left` and `right` are equal element by element, treating NaN as equal to NaN.
    fn assert_same_values(left: &[f32], right: &[f32]) {
        assert_eq!(left.len(), right.len());
        for (index, (a, b)) in left.iter().zip(right).enumerate() {
            assert!(a == b || (a.is_nan() && b.is_nan()), "base {index}: {a} != {b}");
        }
    }

    /// The data index's leaf entries in the little-endian `file`.
    fn data_blocks(file: &[u8]) -> Vec<LeafEntry> {
        let full_index_offset = u64::from_le_bytes(file[24..32].try_into().unwrap());
        let mut source = Cursor::new(file);
        let index = RTree::read(&mut source, full_index_offset, file.len() as u64).unwrap();
        index.leaves(&mut source).unwrap()
    }

    /// A handle on `file` whose source logs its reads, with the reads made while opening it
    /// already cleared.
    fn open_logged(file: &TestFile) -> BigWigReader<ReadLog<Cursor<Vec<u8>>>> {
        let mut reader = BigWigReader::new(ReadLog::new(Cursor::new(file.build()))).unwrap();
        reader.source.reads.clear();
        reader
    }

    /// For each of `blocks`, how many of the reads `reader` logged hold the whole block.
    fn reads_holding(
        reader: &BigWigReader<ReadLog<Cursor<Vec<u8>>>>,
        blocks: &[LeafEntry],
    ) -> Vec<usize> {
        let holds = |read: &Range<u64>, block: &LeafEntry| {
            read.start <= block.data_offset && block.data_offset + block.data_size <= read.end
        };
        blocks
            .iter()
            .map(|block| reader.source.reads.iter().filter(|read| holds(read, block)).count())
            .collect()
    }

    /// Whether `result` is an [`Error::InvalidRegion`].
    fn is_invalid_region<T>(result: &Result<T>) -> bool {
        matches!(result, Err(Error::InvalidRegion { .. }))
    }

    /// Whether `result` is an [`Error::Corrupt`].
    fn is_corrupt<T>(result: &Result<T>) -> bool {
        matches!(result, Err(Error::Corrupt { .. }))
    }

    /// Whether `result` is an [`Error::NotBigWig`].
    fn is_not_bigwig<T>(result: &Result<T>) -> bool {
        matches!(result, Err(Error::NotBigWig { .. }))
    }

    /// `per_chrom` bedGraph sections of three 10-base items each, 100 bases apart, for each of
    /// `chrom_ids` in the order given. Every item has a distinct value.
    fn tiled_sections(chrom_ids: &[u32], per_chrom: u32) -> Vec<TestSection> {
        let mut sections = Vec::new();
        for &chrom_id in chrom_ids {
            for index in 0..per_chrom {
                let base = index * 100;
                let value =
                    |k: u32| f32::from(u16::try_from(chrom_id * 1000 + index * 3 + k).unwrap());
                sections.push(bed_graph(
                    chrom_id,
                    &[
                        (base, base + 10, value(0)),
                        (base + 30, base + 40, value(1)),
                        (base + 60, base + 70, value(2)),
                    ],
                ));
            }
        }
        sections
    }

    /// The zoom record `[start, end)` on `chrom_id` summarising `summary`.
    fn record(chrom_id: u32, start: u32, end: u32, summary: Summary) -> ZoomRecord {
        ZoomRecord { chrom_id, start, end, summary }
    }

    /// The summary of `bases` bases that all have `value`.
    fn uniform(bases: u32, value: f64) -> Summary {
        let bases = f64::from(bases);
        Summary {
            covered_bases: bases,
            min: value,
            max: value,
            sum: bases * value,
            sum_squares: bases * value * value,
        }
    }

    /// Runs every kind of query over every chromosome of `bytes`, returning the first error.
    fn exercise(bytes: Vec<u8>) -> Result<()> {
        let mut reader = BigWigReader::new(Cursor::new(bytes))?;
        for item in reader.all_items()? {
            item?;
        }
        reader.compute_total_summary()?;
        let levels = reader.zoom_levels().to_vec();
        for chrom in reader.chroms().to_vec() {
            // Corruption can make a chromosome huge; a prefix is enough.
            let end = chrom.length.min(10_000);
            for item in reader.chrom_items(&chrom.name)? {
                item?;
            }
            reader.values(&chrom.name, 0, end)?;
            reader.summaries(&chrom.name, 0, end, 4, SummarySource::Raw)?;
            reader.summaries(&chrom.name, 0, end, 4, SummarySource::BestZoom)?;
            for &level in &levels {
                for record in reader.zoom_records(level, &chrom.name, 0, end)? {
                    record?;
                }
            }
        }
        Ok(())
    }

    /// Asserts that two handles hold the same chromosomes, items, values, summaries and zoom
    /// records.
    fn assert_same_contents<A: Read + Seek, B: Read + Seek>(
        a: &mut BigWigReader<A>,
        b: &mut BigWigReader<B>,
    ) {
        assert_eq!(a.chroms(), b.chroms());
        assert_eq!(collect(a.all_items().unwrap()), collect(b.all_items().unwrap()));
        let reductions = |levels: &[ZoomLevel]| -> Vec<u32> {
            levels.iter().map(ZoomLevel::reduction).collect()
        };
        assert_eq!(reductions(a.zoom_levels()), reductions(b.zoom_levels()));
        let (a_levels, b_levels) = (a.zoom_levels().to_vec(), b.zoom_levels().to_vec());
        for chrom in a.chroms().to_vec() {
            let (name, length) = (chrom.name.as_str(), chrom.length);
            assert_same_values(
                &a.values(name, 0, length).unwrap(),
                &b.values(name, 0, length).unwrap(),
            );
            for source in [SummarySource::Raw, SummarySource::BestZoom] {
                assert_eq!(
                    a.summaries(name, 0, length, 7, source).unwrap(),
                    b.summaries(name, 0, length, 7, source).unwrap()
                );
            }
            for (&a_level, &b_level) in a_levels.iter().zip(&b_levels) {
                assert_eq!(
                    collect(a.zoom_records(a_level, name, 0, length).unwrap()),
                    collect(b.zoom_records(b_level, name, 0, length).unwrap())
                );
            }
        }
    }

    /// Asserts that `reader` returns exactly `file`'s items, values and exact summaries.
    fn assert_matches_file<R: Read + Seek>(reader: &mut BigWigReader<R>, file: &TestFile) {
        let mut by_chrom_then_start = file.items();
        by_chrom_then_start.sort_by_key(|item| (item.chrom_id, item.start));
        assert_eq!(collect(reader.all_items().unwrap()), by_chrom_then_start);
        for (chrom_id, (name, length)) in (0u32..).zip(&file.chroms) {
            let expected = dense_values(file, chrom_id, 0, *length);
            assert_same_values(&reader.values(name, 0, *length).unwrap(), &expected);
            assert_eq!(
                reader.summary(name, 0, *length, SummarySource::Raw).unwrap(),
                dense_summary(&expected)
            );
        }
    }

    // ---- Opening and metadata ----------------------------------------------------------------

    #[test]
    fn info_reports_the_header_fields() {
        let reader = sample_file().open();
        let info = reader.info();
        assert_eq!(info.version, 4);
        assert_eq!(info.byte_order, ByteOrder::Little);
        assert_eq!(info.section_count, 4);
        assert!(info.uncompress_buf_size > 0);
        assert_eq!(info.zoom_levels, []);
    }

    /// bwtool #63: the header summary disagreed with the data.
    #[test]
    fn total_summary_is_reported_as_stored_even_when_it_disagrees_with_data() {
        // libBigWig's maxVal defect stores the smallest positive double as the maximum.
        let stored = Summary {
            covered_bases: 7.0,
            min: -1.0,
            max: f64::MIN_POSITIVE,
            sum: 3.0,
            sum_squares: 9.0,
        };
        let file = TestFile { total_summary: TotalSummaryStyle::Given(stored), ..sample_file() };
        assert_eq!(file.open().info().total_summary, Some(stored));
    }

    #[test]
    fn total_summary_is_none_when_its_offset_is_zero() {
        let file = TestFile { total_summary: TotalSummaryStyle::Absent, ..sample_file() };
        assert_eq!(file.open().info().total_summary, None);
    }

    /// bigtools #9, `WiggleTools` #48 and bwtool #47: tools paired chromosomes by key order rather
    /// than by `chromId`.
    #[test]
    fn chroms_are_reported_in_chrom_id_order_not_key_order() {
        // The tree's keys are byte-sorted (chr1, chr10, chr2) but the IDs run the other way.
        let file = TestFile::new(&[("chr2", 10), ("chr10", 20), ("chr1", 30)], Vec::new());
        let reader = file.open();
        assert_eq!(
            reader.chroms(),
            [ChromInfo::new("chr2", 10), ChromInfo::new("chr10", 20), ChromInfo::new("chr1", 30)]
        );
    }

    #[test]
    fn all_items_follow_chrom_id_order_not_file_order() {
        let chroms = [("chr2", 100), ("chr10", 100), ("chr1", 100)];
        let sections = vec![
            bed_graph(2, &[(0, 10, 1.0)]),
            bed_graph(1, &[(0, 10, 2.0)]),
            bed_graph(0, &[(0, 10, 3.0)]),
        ];
        let mut reader = TestFile::new(&chroms, sections).open();
        assert_eq!(
            collect(reader.all_items().unwrap()),
            [item(0, 0, 10, 3.0), item(1, 0, 10, 2.0), item(2, 0, 10, 1.0)]
        );
    }

    #[test]
    fn all_items_sort_each_chromosome_by_start() {
        let sections = vec![bed_graph(0, &[(200, 210, 2.0)]), bed_graph(0, &[(0, 10, 1.0)])];
        let mut reader = TestFile::new(&[("chr1", 1000)], sections).open();
        assert_eq!(
            collect(reader.all_items().unwrap()),
            [item(0, 0, 10, 1.0), item(0, 200, 210, 2.0)]
        );
    }

    #[test]
    fn chrom_looks_up_a_chromosome_by_name() {
        let reader = sample_file().open();
        assert_eq!(reader.chrom("chr2"), Some(&ChromInfo::new("chr2", 500)));
    }

    #[test]
    fn chrom_is_none_for_an_unknown_name() {
        let reader = sample_file().open();
        assert_eq!(reader.chrom("chrX"), None);
        assert_eq!(reader.chrom("chr"), None);
    }

    #[test]
    fn chrom_name_resolves_a_chrom_id() {
        let reader = sample_file().open();
        assert_eq!(reader.chrom_name(0), Some("chr1"));
        assert_eq!(reader.chrom_name(2), Some("chr3"));
        assert_eq!(reader.chrom_name(3), None);
    }

    #[test]
    fn duplicate_chromosome_names_are_corrupt() {
        let file = TestFile {
            chrom_tree: ChromTreeStyle::Unsorted,
            ..TestFile::new(&[("chr1", 10), ("chr1", 20)], Vec::new())
        };
        assert!(is_corrupt(&BigWigReader::new(Cursor::new(file.build()))));
    }

    #[test]
    fn short_text_file_is_not_a_bigwig() {
        let text = b"chr1\t0\t10\t1.5\n".to_vec();
        assert!(matches!(BigWigReader::new(Cursor::new(text)), Err(Error::NotBigWig { .. })));
    }

    #[test]
    fn empty_input_is_not_a_bigwig() {
        assert!(is_not_bigwig(&BigWigReader::new(Cursor::new(Vec::new()))));
    }

    #[test]
    fn one_byte_input_is_not_a_bigwig() {
        assert!(is_not_bigwig(&BigWigReader::new(Cursor::new(vec![0x26]))));
    }

    #[test]
    fn first_three_bytes_of_the_magic_are_not_a_bigwig() {
        let bytes = BIGWIG_MAGIC.to_le_bytes()[..3].to_vec();
        assert!(is_not_bigwig(&BigWigReader::new(Cursor::new(bytes))));
    }

    #[test]
    fn bigbed_file_is_rejected() {
        let mut bytes = sample_file().build();
        bytes[..4].copy_from_slice(&0x8789_F2EB_u32.to_le_bytes());
        assert!(matches!(BigWigReader::new(Cursor::new(bytes)), Err(Error::IsBigBed)));
    }

    #[test]
    fn version_five_is_unsupported() {
        let mut bytes = sample_file().build();
        bytes[4..6].copy_from_slice(&5u16.to_le_bytes());
        assert!(matches!(BigWigReader::new(Cursor::new(bytes)), Err(Error::UnsupportedVersion(5))));
    }

    /// The little-endian file `bytes` with the u64 at `at` set to `value`.
    fn with_u64_at(mut bytes: Vec<u8>, at: usize, value: u64) -> Vec<u8> {
        bytes[at..at + 8].copy_from_slice(&value.to_le_bytes());
        bytes
    }

    #[test]
    fn chrom_tree_offset_past_the_end_of_the_file_is_corrupt() {
        let bytes = sample_file().build();
        let past_end = bytes.len() as u64 + 1;
        assert!(is_corrupt(&BigWigReader::new(Cursor::new(with_u64_at(bytes, 8, past_end)))));
    }

    #[test]
    fn data_offset_past_the_end_of_the_file_is_corrupt() {
        let bytes = sample_file().build();
        let past_end = bytes.len() as u64;
        assert!(is_corrupt(&BigWigReader::new(Cursor::new(with_u64_at(bytes, 16, past_end)))));
    }

    #[test]
    fn data_index_offset_past_the_end_of_the_file_is_corrupt() {
        let bytes = sample_file().build();
        let past_end = bytes.len() as u64 + 100;
        assert!(is_corrupt(&BigWigReader::new(Cursor::new(with_u64_at(bytes, 24, past_end)))));
    }

    #[test]
    fn total_summary_offset_past_the_end_of_the_file_is_corrupt() {
        let bytes = sample_file().build();
        let past_end = bytes.len() as u64 - 10;
        assert!(is_corrupt(&BigWigReader::new(Cursor::new(with_u64_at(bytes, 44, past_end)))));
    }

    #[test]
    fn zoom_data_offset_past_the_end_of_the_file_is_corrupt() {
        let bytes = TestFile { zoom_reductions: vec![10], ..sample_file() }.build();
        let past_end = bytes.len() as u64 + 1;
        let zoom_data_offset_field = FIXED_HEADER_SIZE + 8;
        let bytes = with_u64_at(bytes, zoom_data_offset_field, past_end);
        assert!(is_corrupt(&BigWigReader::new(Cursor::new(bytes))));
    }

    #[test]
    fn zoom_index_offset_past_the_end_of_the_file_is_corrupt() {
        let bytes = TestFile { zoom_reductions: vec![10], ..sample_file() }.build();
        let past_end = bytes.len() as u64 + 1;
        let zoom_index_offset_field = FIXED_HEADER_SIZE + 16;
        let bytes = with_u64_at(bytes, zoom_index_offset_field, past_end);
        assert!(is_corrupt(&BigWigReader::new(Cursor::new(bytes))));
    }

    // ---- Handles -----------------------------------------------------------------------------

    #[test]
    fn handles_are_send() {
        fn assert_send<T: Send>() {}
        assert_send::<BigWigReader<BufReader<File>>>();
        assert_send::<BigWigReader<Cursor<Vec<u8>>>>();
    }

    #[test]
    fn with_source_shares_metadata_and_queries_independently() {
        let file = sample_file();
        let mut reader = file.open();
        let mut other = reader.with_source(Cursor::new(file.build()));
        assert_eq!(other.chroms(), reader.chroms());
        assert_eq!(other.info(), reader.info());
        let first = collect(reader.items("chr1", 0, 100).unwrap());
        assert_eq!(collect(other.items("chr2", 0, 500).unwrap()), expected_items(&file, 1, 0, 500));
        assert_eq!(first, expected_items(&file, 0, 0, 100));
    }

    #[test]
    fn into_inner_returns_the_source() {
        let bytes = sample_file().build();
        let reader = BigWigReader::new(Cursor::new(bytes.clone())).unwrap();
        assert_eq!(reader.into_inner().into_inner(), bytes);
    }

    #[test]
    fn open_reads_a_file_by_path() {
        let file = sample_file();
        let temp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), file.build()).unwrap();
        let mut reader = BigWigReader::open(temp.path()).unwrap();
        assert_matches_file(&mut reader, &file);
    }

    #[test]
    fn open_of_a_missing_file_is_an_io_error() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("missing.bw");
        assert!(matches!(BigWigReader::open(missing), Err(Error::Io(_))));
    }

    #[test]
    fn open_error_names_the_path() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("missing.bw");
        let error = BigWigReader::open(&missing).unwrap_err();
        assert!(error.to_string().contains(&missing.display().to_string()), "{error}");
    }

    #[test]
    fn open_error_keeps_the_io_error_kind() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("missing.bw");
        let Err(Error::Io(error)) = BigWigReader::open(missing) else {
            panic!("expected an I/O error");
        };
        assert_eq!(error.kind(), std::io::ErrorKind::NotFound);
    }

    #[test]
    fn reopen_error_names_the_path_once_the_file_is_gone() {
        let temp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), sample_file().build()).unwrap();
        let path = temp.path().to_path_buf();
        let reader = BigWigReader::open(&path).unwrap();
        drop(temp);
        let error = reader.reopen().unwrap_err();
        assert!(error.to_string().contains(&path.display().to_string()), "{error}");
    }

    /// bigtools #88 and `JetBrains` big #39: concurrent readers of one file interfered.
    #[test]
    fn independent_reopened_handles_can_query_concurrently() {
        let file = TestFile {
            rtree_block_size: 2,
            ..TestFile::new(&[("chr1", 5000), ("chr2", 5000)], tiled_sections(&[0, 1], 40))
        };
        let temp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), file.build()).unwrap();
        let reader = BigWigReader::open(temp.path()).unwrap();
        std::thread::scope(|scope| {
            for thread in 0u32..4 {
                let mut handle = reader.reopen().unwrap();
                let file = &file;
                scope.spawn(move || {
                    for round in 0u32..20 {
                        let chrom_id = (thread + round) % 2;
                        let name = &file.chroms[usize::try_from(chrom_id).unwrap()].0;
                        let start = (round * 197) % 4000;
                        let values = handle.values(name, start, start + 1000).unwrap();
                        assert_same_values(
                            &values,
                            &dense_values(file, chrom_id, start, start + 1000),
                        );
                    }
                });
            }
        });
    }

    #[test]
    fn reopen_needs_a_handle_opened_by_path() {
        let temp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), sample_file().build()).unwrap();
        let source = BufReader::new(File::open(temp.path()).unwrap());
        let reader = BigWigReader::new(source).unwrap();
        assert!(matches!(reader.reopen(), Err(Error::Io(_))));
    }

    #[test]
    fn reopen_fails_when_the_file_has_changed_length() {
        use std::io::Write;
        let temp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(temp.path(), sample_file().build()).unwrap();
        let reader = BigWigReader::open(temp.path()).unwrap();
        let mut appended = std::fs::OpenOptions::new().append(true).open(temp.path()).unwrap();
        appended.write_all(b"more").unwrap();
        assert!(matches!(reader.reopen(), Err(Error::Io(_))));
    }

    // ---- Items -------------------------------------------------------------------------------

    #[test]
    fn items_returns_overlapping_items_unclipped() {
        let mut reader = sample_file().open();
        assert_eq!(
            collect(reader.items("chr1", 15, 45).unwrap()),
            [item(0, 10, 20, 1.0), item(0, 20, 25, 2.0), item(0, 40, 60, -3.5)]
        );
    }

    /// megadepth #9: an item larger than the query on both sides was missed.
    #[test]
    fn item_strictly_containing_the_query_region_is_returned() {
        let mut reader = sample_file().open();
        assert_eq!(collect(reader.items("chr1", 45, 50).unwrap()), [item(0, 40, 60, -3.5)]);
    }

    #[test]
    fn items_exclude_items_that_only_touch_the_region() {
        let mut reader = sample_file().open();
        assert_eq!(collect(reader.items("chr1", 25, 40).unwrap()), []);
    }

    #[test]
    fn variable_step_item_extents_come_from_the_span() {
        let mut reader = sample_file().open();
        assert_eq!(
            collect(reader.items("chr1", 100, 250).unwrap()),
            [item(0, 100, 105, 4.0), item(0, 110, 115, 5.0), item(0, 200, 205, 6.5)]
        );
    }

    #[test]
    fn fixed_step_item_extents_come_from_the_step_and_span() {
        let mut reader = sample_file().open();
        assert_eq!(
            collect(reader.items("chr1", 300, 400).unwrap()),
            [
                item(0, 300, 305, 7.0),
                item(0, 310, 315, 8.0),
                item(0, 320, 325, 9.0),
                item(0, 330, 335, 10.0)
            ]
        );
    }

    /// pyBigWig #124 and bwtool #57: a zero-length query failed or crashed.
    #[test]
    fn zero_length_query_returns_empty_not_an_error() {
        let mut reader = sample_file().open();
        assert_eq!(collect(reader.items("chr1", 15, 15).unwrap()), []);
        assert_eq!(reader.values("chr1", 15, 15).unwrap(), []);
        assert_eq!(reader.summary("chr1", 15, 15, SummarySource::Raw).unwrap(), None);
    }

    #[test]
    fn items_on_a_chromosome_without_data_are_empty() {
        let mut reader = sample_file().open();
        assert_eq!(collect(reader.items("chr3", 0, 100).unwrap()), []);
        assert_eq!(collect(reader.chrom_items("chr3").unwrap()), []);
    }

    #[test]
    fn items_come_in_file_order() {
        let sections = vec![bed_graph(0, &[(200, 210, 2.0)]), bed_graph(0, &[(0, 10, 1.0)])];
        let mut reader = TestFile::new(&[("chr1", 1000)], sections).open();
        assert_eq!(
            collect(reader.items("chr1", 0, 1000).unwrap()),
            [item(0, 200, 210, 2.0), item(0, 0, 10, 1.0)]
        );
    }

    #[test]
    fn chrom_items_returns_every_item_of_one_chromosome() {
        let file = sample_file();
        let mut reader = file.open();
        assert_eq!(collect(reader.chrom_items("chr1").unwrap()), expected_items(&file, 0, 0, 1000));
        assert_eq!(collect(reader.chrom_items("chr2").unwrap()), expected_items(&file, 1, 0, 500));
    }

    #[test]
    fn chrom_items_of_an_unknown_chromosome_is_unknown_chrom() {
        let mut reader = sample_file().open();
        assert!(
            matches!(reader.chrom_items("chrX"), Err(Error::UnknownChrom(name)) if name == "chrX")
        );
    }

    #[test]
    fn all_items_returns_every_item() {
        let file = sample_file();
        assert_eq!(collect(file.open().all_items().unwrap()), file.items());
    }

    /// bigtools #18: iteration never ended when some chromosomes had no data.
    #[test]
    fn iteration_terminates_when_some_chromosomes_have_no_data() {
        let chroms = [("a", 100), ("b", 100), ("c", 100), ("d", 100), ("e", 100)];
        let sections = vec![bed_graph(1, &[(0, 10, 1.0)]), bed_graph(3, &[(5, 15, 2.0)])];
        let mut reader = TestFile::new(&chroms, sections).open();
        assert_eq!(
            collect(reader.all_items().unwrap()),
            [item(1, 0, 10, 1.0), item(3, 5, 15, 2.0)]
        );
        for name in ["a", "c", "e"] {
            assert_eq!(collect(reader.chrom_items(name).unwrap()), []);
        }
    }

    #[test]
    fn item_running_past_the_chromosome_end_is_clipped_to_it() {
        let file = TestFile::new(&[("chr1", 100)], vec![bed_graph(0, &[(90, 120, 1.0)])]);
        let mut reader = file.open();
        assert_eq!(collect(reader.items("chr1", 0, 100).unwrap()), [item(0, 90, 100, 1.0)]);
        assert_eq!(collect(reader.all_items().unwrap()), [item(0, 90, 100, 1.0)]);
        assert_eq!(reader.values("chr1", 95, 100).unwrap(), [1.0; 5]);
    }

    /// bigtools #78: writers without clipping store items past the chromosome's end.
    #[test]
    fn item_starting_at_or_past_the_chromosome_end_is_dropped() {
        let sections = vec![bed_graph(0, &[(10, 20, 1.0), (100, 110, 2.0), (150, 160, 3.0)])];
        let mut reader = TestFile::new(&[("chr1", 100)], sections).open();
        assert_eq!(collect(reader.all_items().unwrap()), [item(0, 10, 20, 1.0)]);
        assert_eq!(collect(reader.items("chr1", 0, 100).unwrap()), [item(0, 10, 20, 1.0)]);
        let mut values = [f32::NAN; 100];
        values[10..20].fill(1.0);
        assert_same_values(&reader.values("chr1", 0, 100).unwrap(), &values);
        let only_kept_item = Some(uniform(10, 1.0));
        assert_eq!(reader.summary("chr1", 0, 100, SummarySource::Raw).unwrap(), only_kept_item);
        assert_eq!(reader.compute_total_summary().unwrap(), uniform(10, 1.0));
    }

    /// A bedGraph section on `chr1` with the zero-length item `[30, 30)` between two others, as
    /// kent's `bedGraphToBigWig` stores it (a modENCODE file has 3 among 17.1 million items).
    fn file_with_a_zero_length_item() -> TestFile {
        let sections = vec![bed_graph(0, &[(10, 20, 1.0), (30, 30, 9.0), (40, 50, 2.0)])];
        TestFile::new(&[("chr1", 100)], sections)
    }

    /// The items of [`file_with_a_zero_length_item`] once the zero-length one is skipped.
    fn items_around_the_zero_length_item() -> [DataItem; 2] {
        [item(0, 10, 20, 1.0), item(0, 40, 50, 2.0)]
    }

    #[test]
    fn all_items_skip_a_zero_length_item_and_keep_its_neighbours() {
        let mut reader = file_with_a_zero_length_item().open();
        assert_eq!(collect(reader.all_items().unwrap()), items_around_the_zero_length_item());
    }

    #[test]
    fn items_skip_a_zero_length_item_and_keep_its_neighbours() {
        let mut reader = file_with_a_zero_length_item().open();
        let items = collect(reader.items("chr1", 0, 100).unwrap());
        assert_eq!(items, items_around_the_zero_length_item());
        assert_eq!(collect(reader.items("chr1", 25, 35).unwrap()), []);
    }

    #[test]
    fn values_ignore_a_zero_length_item() {
        let mut reader = file_with_a_zero_length_item().open();
        let mut expected = [f32::NAN; 100];
        expected[10..20].fill(1.0);
        expected[40..50].fill(2.0);
        assert_same_values(&reader.values("chr1", 0, 100).unwrap(), &expected);
    }

    #[test]
    fn summary_ignores_a_zero_length_item() {
        let mut reader = file_with_a_zero_length_item().open();
        let summary = reader.summary("chr1", 0, 100, SummarySource::Raw).unwrap().unwrap();
        let expected =
            Summary { covered_bases: 20.0, min: 1.0, max: 2.0, sum: 30.0, sum_squares: 50.0 };
        assert_eq!(summary, expected);
    }

    #[test]
    fn compute_total_summary_ignores_a_zero_length_item() {
        let mut reader = file_with_a_zero_length_item().open();
        let total = reader.compute_total_summary().unwrap();
        assert_eq!((total.covered_bases, total.sum, total.max), (20.0, 30.0, 2.0));
    }

    #[test]
    fn section_with_zero_span_reads_as_empty() {
        let sections =
            vec![variable_step(0, 0, &[(10, 1.0), (20, 2.0)]), bed_graph(0, &[(50, 60, 3.0)])];
        let mut reader = TestFile::new(&[("chr1", 100)], sections).open();
        assert_eq!(collect(reader.all_items().unwrap()), [item(0, 50, 60, 3.0)]);
    }

    #[test]
    fn section_whose_chrom_id_is_not_in_the_tree_is_corrupt() {
        let sections = vec![bed_graph(0, &[(0, 10, 1.0)]), bed_graph(5, &[(0, 10, 2.0)])];
        let mut reader = TestFile::new(&[("chr1", 100)], sections).open();
        let results: Vec<Result<DataItem>> = reader.all_items().unwrap().collect();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].as_ref().unwrap(), &item(0, 0, 10, 1.0));
        assert!(is_corrupt(&results[1]));
    }

    /// igv.js #324: a query spanning several R-tree leaves returned the first leaf's items twice.
    #[test]
    fn query_spanning_multiple_rtree_leaves_returns_each_item_once() {
        let file = TestFile {
            rtree_block_size: 2,
            ..TestFile::new(&[("chr1", 10_000)], tiled_sections(&[0], 20))
        };
        let mut reader = file.open();
        assert_eq!(collect(reader.items("chr1", 0, 10_000).unwrap()), file.items());
        assert_eq!(
            collect(reader.items("chr1", 150, 1250).unwrap()),
            expected_items(&file, 0, 150, 1250)
        );
    }

    /// The sample with every data and zoom leaf entry repeated 50 times, as a hostile file might.
    fn sample_with_repeated_leaf_entries() -> TestFile {
        TestFile { zoom_reductions: vec![10], leaf_entry_copies: 50, ..sample_file() }
    }

    /// The distinct data blocks of `file`, which lie where they would with one entry each.
    fn distinct_data_blocks(file: &TestFile) -> Vec<LeafEntry> {
        data_blocks(&TestFile { leaf_entry_copies: 1, ..file.clone() }.build())
    }

    #[test]
    fn repeated_leaf_entries_yield_each_item_once_from_items() {
        let file = sample_with_repeated_leaf_entries();
        let mut reader = file.open();
        assert_eq!(
            collect(reader.items("chr1", 0, 1000).unwrap()),
            expected_items(&file, 0, 0, 1000)
        );
    }

    #[test]
    fn repeated_leaf_entries_yield_each_item_once_from_all_items() {
        let file = sample_with_repeated_leaf_entries();
        assert_eq!(collect(file.open().all_items().unwrap()), file.items());
    }

    #[test]
    fn repeated_leaf_entries_leave_raw_summaries_exact() {
        let file = sample_with_repeated_leaf_entries();
        let mut reader = file.open();
        assert_eq!(
            reader.summary("chr1", 0, 1000, SummarySource::Raw).unwrap(),
            dense_summary(&dense_values(&file, 0, 0, 1000))
        );
    }

    #[test]
    fn repeated_zoom_leaf_entries_yield_each_record_once() {
        let file = sample_with_repeated_leaf_entries();
        let mut reader = file.open();
        let level = reader.zoom_levels()[0];
        let records = collect(reader.zoom_records(level, "chr1", 0, 1000).unwrap());
        let mut once = TestFile { leaf_entry_copies: 1, ..file }.open();
        assert_eq!(records, collect(once.zoom_records(level, "chr1", 0, 1000).unwrap()));
    }

    #[test]
    fn repeated_leaf_entries_are_read_once_by_items() {
        let file = sample_with_repeated_leaf_entries();
        let chr1_blocks: Vec<LeafEntry> = distinct_data_blocks(&file)
            .into_iter()
            .filter(|block| block.start_chrom_ix == 0)
            .collect();
        let mut reader = open_logged(&file);
        assert_eq!(collect(reader.items("chr1", 0, 1000).unwrap()).len(), 10);
        assert_eq!(reads_holding(&reader, &chr1_blocks), [1, 1, 1]);
    }

    #[test]
    fn repeated_leaf_entries_are_read_once_by_values() {
        let file = sample_with_repeated_leaf_entries();
        let chr2_block = distinct_data_blocks(&file)[3];
        let mut reader = open_logged(&file);
        reader.values("chr2", 0, 500).unwrap();
        assert_eq!(reads_holding(&reader, &[chr2_block]), [1]);
    }

    #[test]
    fn repeated_leaf_entries_are_read_once_by_all_items() {
        let file = sample_with_repeated_leaf_entries();
        let blocks = distinct_data_blocks(&file);
        let mut reader = open_logged(&file);
        assert_eq!(collect(reader.all_items().unwrap()).len(), 12);
        assert_eq!(reads_holding(&reader, &blocks), [1, 1, 1, 1]);
    }

    /// A block at `data_offset` of `data_size` bytes, with bounds that play no part.
    fn block_at(data_offset: u64, data_size: u64) -> LeafEntry {
        LeafEntry {
            start_chrom_ix: 0,
            start_base: 0,
            end_chrom_ix: 0,
            end_base: 1,
            data_offset,
            data_size,
        }
    }

    #[test]
    fn dedupe_sorts_blocks_by_offset() {
        let mut blocks = vec![block_at(300, 10), block_at(100, 10), block_at(200, 10)];
        sort_and_dedupe_blocks(&mut blocks);
        assert_eq!(blocks, [block_at(100, 10), block_at(200, 10), block_at(300, 10)]);
    }

    #[test]
    fn dedupe_keeps_one_copy_of_a_repeated_block() {
        let mut blocks = vec![block_at(100, 10), block_at(200, 10), block_at(100, 10)];
        sort_and_dedupe_blocks(&mut blocks);
        assert_eq!(blocks, [block_at(100, 10), block_at(200, 10)]);
    }

    #[test]
    fn dedupe_keeps_the_largest_of_blocks_sharing_an_offset() {
        let mut blocks = vec![block_at(100, 5), block_at(100, 10), block_at(100, 7)];
        sort_and_dedupe_blocks(&mut blocks);
        assert_eq!(blocks, [block_at(100, 10)]);
    }

    #[test]
    fn dedupe_drops_a_block_inside_an_earlier_one() {
        let mut blocks = vec![block_at(100, 50), block_at(120, 10), block_at(150, 10)];
        sort_and_dedupe_blocks(&mut blocks);
        assert_eq!(blocks, [block_at(100, 50), block_at(150, 10)]);
    }

    #[test]
    fn dedupe_keeps_blocks_that_only_partly_overlap() {
        let mut blocks = vec![block_at(100, 50), block_at(140, 20)];
        sort_and_dedupe_blocks(&mut blocks);
        assert_eq!(blocks, [block_at(100, 50), block_at(140, 20)]);
    }

    #[test]
    fn blocks_separated_by_gaps_are_read_correctly() {
        // One gap small enough to read through and one large enough to need a separate read.
        for section_gap in [100, 10_000] {
            let file = TestFile {
                section_gap,
                ..TestFile::new(&[("chr1", 2000)], tiled_sections(&[0], 6))
            };
            assert_matches_file(&mut file.open(), &file);
        }
    }

    #[test]
    fn data_larger_than_one_merged_read_is_read_in_several() {
        let sections = (0u32..6)
            .map(|index| {
                let values: Vec<f32> =
                    (0u16..60_000).map(|value| f32::from(value % 1000)).collect();
                fixed_step(0, index * 60_000, 1, 1, &values)
            })
            .collect();
        let file = TestFile { compressed: false, ..TestFile::new(&[("chr1", 400_000)], sections) };
        assert_matches_file(&mut file.open(), &file);
    }

    // ---- Values ------------------------------------------------------------------------------

    /// bigtools #71: data was lost when NaN marked missing values.
    #[test]
    fn values_with_nan_as_missing_value_still_returns_data() {
        let file = sample_file();
        let mut reader = file.open();
        assert_same_values(
            &reader.values("chr1", 0, 400).unwrap(),
            &dense_values(&file, 0, 0, 400),
        );
    }

    #[test]
    fn values_into_fills_missing_bases_and_counts_covered_ones() {
        let mut reader = sample_file().open();
        let mut out = [0.0; 50];
        let covered = reader.values_into("chr1", 0, 50, &mut out, -1.0).unwrap();
        assert_eq!(covered, 10 + 5 + 10);
        assert_eq!(out[..10], [-1.0; 10]);
        assert_eq!(out[10..20], [1.0; 10]);
        assert_eq!(out[20..25], [2.0; 5]);
        assert_eq!(out[25..40], [-1.0; 15]);
        assert_eq!(out[40..], [-3.5; 10]);
    }

    #[test]
    fn values_into_counts_a_stored_value_equal_to_missing_as_covered() {
        let mut reader = sample_file().open();
        let mut out = [9.0; 10];
        assert_eq!(reader.values_into("chr1", 10, 20, &mut out, 1.0).unwrap(), 10);
        assert_eq!(out, [1.0; 10]);
    }

    #[test]
    fn values_into_counts_covered_bases_across_a_long_region() {
        let file = TestFile::new(&[("chr1", 1000)], vec![bed_graph(0, &[(3, 197, 2.5)])]);
        let mut out = vec![0.0; 200];
        assert_eq!(file.open().values_into("chr1", 0, 200, &mut out, 0.0).unwrap(), 194);
    }

    #[test]
    fn values_into_reuses_its_buffer_without_stale_values() {
        let mut reader = sample_file().open();
        let mut out = [0.0; 10];
        reader.values_into("chr1", 10, 20, &mut out, f32::NAN).unwrap();
        let covered = reader.values_into("chr1", 25, 35, &mut out, f32::NAN).unwrap();
        assert_eq!(covered, 0);
        assert!(out.iter().all(|value| value.is_nan()));
    }

    #[test]
    fn later_item_in_file_order_wins_where_items_overlap() {
        // Lax writers (libBigWig) accept overlapping items.
        let sections = vec![bed_graph(0, &[(0, 10, 1.0)]), bed_graph(0, &[(5, 15, 2.0)])];
        let mut reader = TestFile::new(&[("chr1", 100)], sections).open();
        let mut out = [0.0; 15];
        assert_eq!(reader.values_into("chr1", 0, 15, &mut out, f32::NAN).unwrap(), 15);
        assert_eq!(out[..5], [1.0; 5]);
        assert_eq!(out[5..], [2.0; 10]);
    }

    #[test]
    fn values_into_rejects_a_buffer_of_the_wrong_length() {
        let mut reader = sample_file().open();
        let mut out = [0.0; 9];
        assert!(is_invalid_region(&reader.values_into("chr1", 10, 20, &mut out, 0.0)));
    }

    /// pyBigWig #158: values from a file without data were garbage.
    #[test]
    fn values_on_a_file_with_no_data_are_all_missing() {
        let mut reader = TestFile::new(&[("chr1", 100)], Vec::new()).open();
        assert!(reader.values("chr1", 0, 100).unwrap().iter().all(|value| value.is_nan()));
    }

    /// pyBigWig #152, bbi-js #59 and `JetBrains` big #35: with name-sorted `chromId`s, data in
    /// numeric order and a multi-level R-tree, readers that stop early missed blocks.
    #[test]
    fn values_correct_when_rtree_leaf_order_is_not_monotonic_in_chrom_id() {
        // kent numbers chr1, chr10, chr2 as 0, 1, 2; the data arrive as chr1, chr2, chr10.
        let file = TestFile {
            rtree_block_size: 2,
            ..TestFile::new(
                &[("chr1", 1000), ("chr10", 1000), ("chr2", 1000)],
                tiled_sections(&[0, 2, 1], 5),
            )
        };
        let mut reader = file.open();
        for (chrom_id, name) in [(0, "chr1"), (1, "chr10"), (2, "chr2")] {
            let values = reader.values(name, 0, 1000).unwrap();
            assert_same_values(&values, &dense_values(&file, chrom_id, 0, 1000));
            assert_eq!(
                collect(reader.items(name, 0, 1000).unwrap()),
                expected_items(&file, chrom_id, 0, 1000)
            );
        }
    }

    // ---- Summaries ---------------------------------------------------------------------------

    #[test]
    fn raw_summary_is_exact() {
        let file = sample_file();
        let mut reader = file.open();
        assert_eq!(
            reader.summary("chr1", 0, 400, SummarySource::Raw).unwrap(),
            dense_summary(&dense_values(&file, 0, 0, 400))
        );
    }

    /// pyBigWig #104 and jbrowse-components #685: a region without data gave 0 or NaN rather
    /// than "no data".
    #[test]
    fn summary_of_region_without_data_is_none() {
        let mut reader = sample_file().open();
        assert_eq!(reader.summary("chr1", 25, 40, SummarySource::Raw).unwrap(), None);
        assert_eq!(reader.summary("chr3", 0, 100, SummarySource::BestZoom).unwrap(), None);
    }

    /// pyBigWig #161: the mean divided by the region's length rather than the covered bases.
    #[test]
    fn summary_mean_divides_by_covered_bases_not_region_length() {
        let mut reader = sample_file().open();
        let summary = reader.summary("chr1", 0, 100, SummarySource::Raw).unwrap().unwrap();
        assert_eq!(summary.covered_bases, 35.0);
        assert_eq!(summary.mean(), -50.0 / 35.0);
        assert_eq!(summary.mean0(100), -0.5);
    }

    #[test]
    fn raw_summaries_match_brute_force_in_every_bin() {
        let file = sample_file();
        let mut reader = file.open();
        let values = dense_values(&file, 0, 0, 400);
        let summaries = reader.summaries("chr1", 0, 400, 7, SummarySource::Raw).unwrap();
        for (index, summary) in (0u32..).zip(summaries) {
            let (start, end) = bin_bounds(0, 400, 7, index);
            let bin = &values[usize::try_from(start).unwrap()..usize::try_from(end).unwrap()];
            assert_eq!(summary, dense_summary(bin), "bin {index}");
        }
    }

    /// pybbi #59 and `JetBrains` big #16: bins left gaps or overlapped.
    #[test]
    fn bins_partition_the_region_with_no_gap_or_overlap() {
        let mut reader = sample_file().open();
        let whole = reader.summary("chr1", 5, 333, SummarySource::Raw).unwrap().unwrap();
        for n_bins in 1..=13 {
            let bins = reader.summaries("chr1", 5, 333, n_bins, SummarySource::Raw).unwrap();
            let covered: f64 = bins.iter().flatten().map(|bin| bin.covered_bases).sum();
            let sum: f64 = bins.iter().flatten().map(|bin| bin.sum).sum();
            assert_eq!((covered, sum), (whole.covered_bases, whole.sum), "{n_bins} bins");
        }
    }

    /// pybbi #26 and `JetBrains` big #15: more bins than bases failed.
    #[test]
    fn summaries_with_more_bins_than_bases_leave_empty_bins_as_none() {
        let mut reader = sample_file().open();
        let bins = reader.summaries("chr1", 10, 13, 6, SummarySource::Raw).unwrap();
        let one_base = Some(uniform(1, 1.0));
        assert_eq!(bins, [None, one_base, None, one_base, None, one_base]);
    }

    #[test]
    fn zero_bins_is_an_invalid_region() {
        let mut reader = sample_file().open();
        assert!(is_invalid_region(&reader.summaries("chr1", 0, 100, 0, SummarySource::Raw)));
    }

    #[test]
    fn u32_max_bins_on_a_short_region_is_an_invalid_region() {
        let mut reader = sample_file().open();
        let result = reader.summaries("chr1", 0, 100, u32::MAX, SummarySource::Raw);
        assert!(is_invalid_region(&result));
    }

    #[test]
    fn one_bin_past_the_minimum_bin_limit_on_a_short_region_is_an_invalid_region() {
        let mut reader = sample_file().open();
        let result = reader.summaries("chr1", 0, 100, MIN_BIN_LIMIT + 1, SummarySource::BestZoom);
        assert!(is_invalid_region(&result));
    }

    #[test]
    fn minimum_bin_limit_on_a_short_region_is_allowed_with_extra_bins_empty() {
        let mut reader = sample_file().open();
        let bins = reader.summaries("chr1", 10, 20, MIN_BIN_LIMIT, SummarySource::Raw).unwrap();
        assert_eq!(bins.len(), 1 << 20);
        assert_eq!(bins.iter().flatten().count(), 10);
    }

    /// The bin limit keeps queries far below any count an allocator refuses, so this checks the
    /// refusal path directly.
    #[test]
    fn bin_storage_the_allocator_refuses_is_an_invalid_region_on_the_query_chromosome() {
        let reader = sample_file().open();
        let region = Region { chrom_id: 1, start: 0, end: 10 };
        let result = reader.bin_storage::<u64>(region, usize::MAX);
        assert!(matches!(
            result,
            Err(Error::InvalidRegion { chrom, start: 0, end: 10, reason: "too many bins to allocate" })
                if chrom == "chr2"
        ));
    }

    #[test]
    fn best_zoom_uses_raw_data_when_the_file_has_no_zoom_levels() {
        let mut reader = sample_file().open();
        assert_eq!(
            reader.summaries("chr1", 0, 1000, 4, SummarySource::BestZoom).unwrap(),
            reader.summaries("chr1", 0, 1000, 4, SummarySource::Raw).unwrap()
        );
    }

    #[test]
    fn best_zoom_reads_zoom_records_where_a_level_is_fine_enough() {
        // The zoom records claim every base is 5 although the data say 1, so the source shows.
        let records = (0..10).map(|bin| record(0, bin * 100, (bin + 1) * 100, uniform(100, 5.0)));
        let file = TestFile {
            zoom_levels: vec![(100, records.collect())],
            ..TestFile::new(&[("chr1", 1000)], vec![bed_graph(0, &[(0, 1000, 1.0)])])
        };
        let mut reader = file.open();
        let zoom = reader.summary("chr1", 0, 1000, SummarySource::BestZoom).unwrap().unwrap();
        let raw = reader.summary("chr1", 0, 1000, SummarySource::Raw).unwrap().unwrap();
        assert_eq!((zoom.mean(), raw.mean()), (5.0, 1.0));
    }

    #[test]
    fn best_zoom_uses_raw_data_where_no_level_is_fine_enough() {
        let records = (0..10).map(|bin| record(0, bin * 100, (bin + 1) * 100, uniform(100, 5.0)));
        let file = TestFile {
            zoom_levels: vec![(100, records.collect())],
            ..TestFile::new(&[("chr1", 1000)], vec![bed_graph(0, &[(0, 1000, 1.0)])])
        };
        let mut reader = file.open();
        // Bins 200 wide want a reduction below 100, and the only level is 100.
        let bins = reader.summaries("chr1", 0, 1000, 5, SummarySource::BestZoom).unwrap();
        assert!(bins.iter().all(|bin| bin.unwrap().mean() == 1.0));
    }

    #[test]
    fn best_zoom_summaries_of_bins_aligned_to_records_match_raw() {
        let sections = (0u32..10)
            .map(|index| {
                let values: Vec<f32> = (0u16..1000).map(|value| f32::from(value % 7)).collect();
                fixed_step(0, index * 1000, 1, 1, &values)
            })
            .collect();
        let file = TestFile {
            zoom_reductions: vec![10, 100, 1000],
            ..TestFile::new(&[("chr1", 10_000)], sections)
        };
        let mut reader = file.open();
        assert_eq!(reader.best_zoom_for(1000).map(|level| level.reduction()), Some(100));
        assert_eq!(
            reader.summaries("chr1", 0, 10_000, 10, SummarySource::BestZoom).unwrap(),
            reader.summaries("chr1", 0, 10_000, 10, SummarySource::Raw).unwrap()
        );
    }

    /// pyBigWig #12: zoom-based statistics differed from exact ones, which callers took for a
    /// bug; both are available and only `Raw` is exact.
    #[test]
    fn raw_summary_is_exact_where_zoom_summary_is_approximate() {
        let squares: Vec<f32> = (0u16..100).map(|base| f32::from(base * base)).collect();
        let file = TestFile {
            zoom_reductions: vec![10],
            ..TestFile::new(&[("chr1", 100)], vec![fixed_step(0, 0, 1, 1, &squares)])
        };
        let mut reader = file.open();
        // 22 bases want a reduction below 11, so level 10, whose records [0, 10) and [20, 30)
        // only partly overlap [5, 27).
        let raw = reader.summary("chr1", 5, 27, SummarySource::Raw).unwrap().unwrap();
        let zoom = reader.summary("chr1", 5, 27, SummarySource::BestZoom).unwrap().unwrap();
        assert_eq!(raw.sum, (5..27).map(|base: u32| f64::from(base * base)).sum());
        assert_eq!(zoom.sum, 0.5 * 285.0 + 2185.0 + 0.7 * 6085.0);
        assert_eq!(zoom.covered_bases, 22.0);
    }

    /// pyBigWig #146: libBigWig truncated the overlap fraction when summing zoom records.
    #[test]
    fn zoom_summary_sum_prorates_by_overlap_fraction() {
        let records =
            (0..5).map(|bin| record(0, bin * 10, (bin + 1) * 10, uniform(10, f64::from(bin))));
        let file = TestFile {
            zoom_levels: vec![(10, records.collect())],
            ..TestFile::new(&[("chr1", 100)], Vec::new())
        };
        let mut reader = file.open();
        let summary = reader.summary("chr1", 5, 45, SummarySource::BestZoom).unwrap().unwrap();
        assert_eq!(summary.covered_bases, 5.0 + 10.0 + 10.0 + 10.0 + 5.0);
        assert_eq!(summary.sum, 0.0 + 10.0 + 20.0 + 30.0 + 0.5 * 40.0);
        assert_eq!((summary.min, summary.max), (0.0, 4.0));
    }

    /// pyBigWig #153: a zoom maximum above the exact one is expected, not a bug.
    #[test]
    fn zoom_summary_max_may_exceed_raw_max_for_partially_overlapping_records() {
        let sections = vec![bed_graph(0, &[(0, 27, 1.0), (27, 28, 9.0), (28, 30, 1.0)])];
        let file =
            TestFile { zoom_reductions: vec![10], ..TestFile::new(&[("chr1", 100)], sections) };
        let mut reader = file.open();
        let raw = reader.summary("chr1", 0, 25, SummarySource::Raw).unwrap().unwrap();
        let zoom = reader.summary("chr1", 0, 25, SummarySource::BestZoom).unwrap().unwrap();
        assert_eq!((raw.max, zoom.max), (1.0, 9.0));
    }

    #[test]
    fn best_zoom_skips_records_with_no_covered_bases() {
        let empty =
            Summary { covered_bases: 0.0, min: -100.0, max: 100.0, sum: 0.0, sum_squares: 0.0 };
        let records = vec![record(0, 0, 10, empty), record(0, 10, 20, uniform(10, 2.0))];
        let file = TestFile {
            zoom_levels: vec![(10, records)],
            ..TestFile::new(&[("chr1", 100)], Vec::new())
        };
        let mut reader = file.open();
        let summary = reader.summary("chr1", 0, 30, SummarySource::BestZoom).unwrap();
        assert_eq!(summary, Some(uniform(10, 2.0)));
    }

    #[test]
    fn compute_total_summary_recomputes_from_raw_data() {
        let file =
            TestFile { total_summary: TotalSummaryStyle::Given(ZERO_SUMMARY), ..sample_file() };
        assert_eq!(file.open().compute_total_summary().unwrap(), file.summary_of_data());
    }

    #[test]
    fn compute_total_summary_of_a_file_without_data_is_all_zeros() {
        let mut reader = TestFile::new(&[("chr1", 100)], Vec::new()).open();
        assert_eq!(reader.compute_total_summary().unwrap(), ZERO_SUMMARY);
    }

    // ---- Zoom levels -------------------------------------------------------------------------

    #[test]
    fn zoom_levels_are_reported_in_header_order() {
        let file = TestFile {
            zoom_levels: vec![(100, Vec::new()), (10, Vec::new()), (1000, Vec::new())],
            ..sample_file()
        };
        let reader = file.open();
        let reductions: Vec<u32> = reader.zoom_levels().iter().map(ZoomLevel::reduction).collect();
        assert_eq!(reductions, [100, 10, 1000]);
        assert_eq!(reader.info().zoom_levels, reader.zoom_levels());
    }

    #[test]
    fn best_zoom_for_takes_the_largest_reduction_strictly_below_half_the_bin_width() {
        let file = TestFile {
            zoom_levels: vec![(100, Vec::new()), (10, Vec::new()), (1000, Vec::new())],
            ..sample_file()
        };
        let reader = file.open();
        let best = |bin_width| reader.best_zoom_for(bin_width).map(|level| level.reduction());
        assert_eq!(best(2000), Some(100));
        assert_eq!(best(2002), Some(1000));
        assert_eq!(best(22), Some(10));
        assert_eq!(best(21), None);
        assert_eq!(best(3), None);
    }

    #[test]
    fn zoom_records_returns_overlapping_records_unclipped() {
        let records: Vec<ZoomRecord> =
            (0..3).map(|bin| record(0, bin * 10, (bin + 1) * 10, uniform(10, 1.0))).collect();
        let file = TestFile { zoom_levels: vec![(10, records.clone())], ..sample_file() };
        let mut reader = file.open();
        let level = reader.zoom_levels()[0];
        assert_eq!(collect(reader.zoom_records(level, "chr1", 15, 25).unwrap()), records[1..]);
        assert_eq!(collect(reader.zoom_records(level, "chr1", 30, 40).unwrap()), []);
        assert_eq!(collect(reader.zoom_records(level, "chr1", 5, 5).unwrap()), []);
    }

    #[test]
    fn zoom_records_of_other_chromosomes_in_the_same_block_are_left_out() {
        // kent packs records of several chromosomes into one block.
        let records = vec![
            record(0, 990, 1000, uniform(10, 1.0)),
            record(1, 0, 10, uniform(10, 2.0)),
            record(1, 10, 20, uniform(10, 3.0)),
        ];
        let file = TestFile { zoom_levels: vec![(10, records.clone())], ..sample_file() };
        let mut reader = file.open();
        let level = reader.zoom_levels()[0];
        assert_eq!(collect(reader.zoom_records(level, "chr2", 0, 500).unwrap()), records[1..]);
        assert_eq!(collect(reader.zoom_records(level, "chr1", 0, 1000).unwrap()), records[..1]);
    }

    #[test]
    fn zoom_level_from_another_file_is_foreign() {
        let mut other = TestFile { zoom_reductions: vec![10], ..sample_file() }.open();
        let foreign = other.zoom_levels()[0];
        let mut reader = TestFile {
            zoom_reductions: vec![10],
            ..TestFile::new(&SAMPLE_CHROMS, tiled_sections(&[0], 3))
        }
        .open();
        assert!(matches!(
            reader.zoom_records(foreign, "chr1", 0, 100),
            Err(Error::ForeignZoomLevel { reduction: 10 })
        ));
        assert!(other.zoom_records(foreign, "chr1", 0, 100).is_ok());
    }

    #[test]
    fn zoom_block_with_a_trailing_partial_record_yields_its_whole_records() {
        let records: Vec<ZoomRecord> =
            (0..4).map(|bin| record(0, bin * 10, (bin + 1) * 10, uniform(10, 1.0))).collect();
        let file = TestFile {
            compressed: false,
            zoom_levels: vec![(10, records.clone())],
            ..sample_file()
        };
        let mut bytes = file.build();
        let zoom_index_offset = u64::from_le_bytes(bytes[80..88].try_into().unwrap());
        let first_leaf_size_field =
            usize::try_from(zoom_index_offset).unwrap() + RTREE_HEADER_SIZE + 4 + 24;
        let size = u64::from_le_bytes(bytes[first_leaf_size_field..][..8].try_into().unwrap());
        bytes[first_leaf_size_field..][..8].copy_from_slice(&(size - 1).to_le_bytes());
        let mut reader = BigWigReader::new(Cursor::new(bytes)).unwrap();
        let level = reader.zoom_levels()[0];
        // The block now ends one byte into its fourth record.
        assert_eq!(collect(reader.zoom_records(level, "chr1", 0, 100).unwrap()), records[..3]);
    }

    #[test]
    fn zoom_record_whose_chrom_id_is_not_in_the_tree_is_corrupt() {
        let records = vec![record(0, 0, 10, uniform(10, 1.0)), record(7, 0, 10, uniform(10, 1.0))];
        let file = TestFile { zoom_levels: vec![(10, records)], ..sample_file() };
        let mut reader = file.open();
        let level = reader.zoom_levels()[0];
        let results: Vec<_> = reader.zoom_records(level, "chr1", 0, 100).unwrap().collect();
        assert!(is_corrupt(&results[results.len() - 1]));
    }

    // ---- Regions -----------------------------------------------------------------------------

    #[test]
    fn start_after_end_is_an_invalid_region() {
        let mut reader = TestFile { zoom_reductions: vec![10], ..sample_file() }.open();
        let level = reader.zoom_levels()[0];
        assert!(is_invalid_region(&reader.items("chr1", 20, 10)));
        assert!(is_invalid_region(&reader.values("chr1", 20, 10)));
        assert!(is_invalid_region(&reader.summary("chr1", 20, 10, SummarySource::Raw)));
        assert!(is_invalid_region(&reader.zoom_records(level, "chr1", 20, 10)));
    }

    /// pyBigWig #49 and bigtools #62: a region past the chromosome's end gave a misleading error
    /// or a panic.
    #[test]
    fn region_past_the_chromosome_end_is_an_invalid_region() {
        let mut reader = sample_file().open();
        assert!(is_invalid_region(&reader.items("chr1", 990, 1001)));
        assert!(is_invalid_region(&reader.values("chr1", 990, 1010)));
        assert!(is_invalid_region(&reader.summaries("chr1", 0, 2000, 3, SummarySource::BestZoom)));
        let mut out = [0.0; 10];
        assert!(is_invalid_region(&reader.values_into("chr2", 500, 510, &mut out, 0.0)));
    }

    #[test]
    fn region_ending_exactly_at_the_chromosome_end_is_valid() {
        let file = sample_file();
        let mut reader = file.open();
        assert_eq!(
            collect(reader.items("chr2", 450, 500).unwrap()),
            expected_items(&file, 1, 450, 500)
        );
        assert_eq!(reader.values("chr2", 499, 500).unwrap(), [1.5]);
    }

    #[test]
    fn unknown_chromosome_is_reported_by_name() {
        let mut reader = sample_file().open();
        assert!(
            matches!(reader.items("chrX", 0, 10), Err(Error::UnknownChrom(name)) if name == "chrX")
        );
        assert!(matches!(reader.values("chrX", 0, 10), Err(Error::UnknownChrom(_))));
        assert!(matches!(
            reader.summary("chrX", 0, 10, SummarySource::Raw),
            Err(Error::UnknownChrom(_))
        ));
    }

    // ---- Writers' variations -----------------------------------------------------------------

    /// libBigWig #61, igv.js #786 and megadepth #14: big-endian files were misread.
    #[test]
    fn big_endian_file_reads_identically_to_little_endian() {
        let little = TestFile { zoom_reductions: vec![10, 100], ..sample_file() };
        let big = TestFile { byte_order: ByteOrder::Big, ..little.clone() };
        let mut big_reader = big.open();
        assert_eq!(big_reader.info().byte_order, ByteOrder::Big);
        assert_matches_file(&mut big_reader, &big);
        assert_same_contents(&mut big_reader, &mut little.open());
    }

    #[test]
    fn undersized_uncompress_buf_size_is_retried_with_a_larger_buffer() {
        // libBigWig declares a constant 32,768, which can be too small.
        let file = TestFile { uncompress_buf_size: Some(16), ..sample_file() };
        assert_matches_file(&mut file.open(), &file);
    }

    /// bbi-js #22 and #82, and igv.js: uncompressed files were misread when a query touched
    /// several blocks.
    #[test]
    fn uncompressed_file_with_multi_block_query_reads_correctly() {
        let file = TestFile { compressed: false, zoom_reductions: vec![10], ..sample_file() };
        let mut reader = file.open();
        assert_eq!(reader.info().uncompress_buf_size, 0);
        assert_matches_file(&mut reader, &file);
        assert_same_contents(&mut reader, &mut TestFile { compressed: true, ..file }.open());
    }

    #[test]
    fn version_3_file_without_trailing_magic_opens() {
        let file = TestFile { version: 3, trailing_magic: false, ..sample_file() };
        let mut reader = file.open();
        assert_eq!(reader.info().version, 3);
        assert_matches_file(&mut reader, &file);
    }

    #[test]
    fn file_with_zero_zoom_levels_summarises_from_raw_data() {
        let mut reader = sample_file().open();
        assert_eq!(reader.zoom_levels(), []);
        assert_eq!(reader.best_zoom_for(1_000_000), None);
        assert_eq!(
            reader.summary("chr1", 0, 1000, SummarySource::BestZoom).unwrap(),
            reader.summary("chr1", 0, 1000, SummarySource::Raw).unwrap()
        );
    }

    #[test]
    fn unused_zoom_header_slots_are_ignored() {
        // bedGraphToBigWig reserves ten zoom headers whatever it writes.
        let file =
            TestFile { zoom_reductions: vec![10, 100], zoom_header_slots: 10, ..sample_file() };
        let mut reader = file.open();
        assert_eq!(reader.zoom_levels().len(), 2);
        assert_same_contents(&mut reader, &mut TestFile { zoom_header_slots: 0, ..file }.open());
    }

    #[test]
    fn bigtools_layout_with_chrom_tree_after_data_and_no_zoom_count_reads() {
        let kent = TestFile { zoom_reductions: vec![10, 100], ..sample_file() };
        let bigtools = TestFile { layout: Layout::Bigtools, ..kent.clone() };
        let mut reader = bigtools.open();
        assert_matches_file(&mut reader, &bigtools);
        assert_same_contents(&mut reader, &mut kent.open());
    }

    #[test]
    fn libbigwig_unsorted_single_leaf_chrom_tree_reads() {
        let chroms = [("chr2", 300), ("chr10", 200), ("chr1", 100)];
        let sections = vec![
            bed_graph(0, &[(0, 10, 1.0)]),
            bed_graph(1, &[(0, 10, 2.0)]),
            bed_graph(2, &[(0, 10, 3.0)]),
        ];
        let file =
            TestFile { chrom_tree: ChromTreeStyle::Unsorted, ..TestFile::new(&chroms, sections) };
        let mut reader = file.open();
        assert_eq!(reader.chrom("chr1"), Some(&ChromInfo::new("chr1", 100)));
        assert_eq!(reader.chrom_name(0), Some("chr2"));
        assert_eq!(collect(reader.items("chr10", 0, 200).unwrap()), [item(1, 0, 10, 2.0)]);
        assert_matches_file(&mut reader, &file);
    }

    /// pyBigWig #166: libBigWig wrote fixedStep section ends six steps too far.
    #[test]
    fn section_end_six_steps_too_large_does_not_change_items() {
        let sections = vec![fixed_step(0, 100, 10, 5, &[1.0, 2.0, 3.0, 4.0])];
        let file =
            TestFile { section_end_excess: 6 * 10, ..TestFile::new(&[("chr1", 1000)], sections) };
        let mut reader = file.open();
        assert_eq!(
            collect(reader.items("chr1", 0, 1000).unwrap()),
            [
                item(0, 100, 105, 1.0),
                item(0, 110, 115, 2.0),
                item(0, 120, 125, 3.0),
                item(0, 130, 135, 4.0)
            ]
        );
        assert_matches_file(&mut reader, &file);
    }

    // ---- Issue-derived robustness ------------------------------------------------------------

    /// igv.js #317 and `JetBrains` big #28: offsets past 2³¹ or 2³² were read into 32-bit integers.
    #[test]
    fn reader_handles_file_offsets_above_4_gib() {
        let file =
            TestFile { zoom_reductions: vec![10, 100], data_start: 5 << 30, ..sample_file() };
        let (head, tail_offset, tail) = file.build_parts();
        let mut reader =
            BigWigReader::new(SparseSource { head, tail_offset, tail, position: 0 }).unwrap();
        assert!(reader.zoom_levels().iter().all(|level| level.index_offset() > 1 << 32));
        assert_matches_file(&mut reader, &file);
        assert_same_contents(&mut reader, &mut TestFile { data_start: 0, ..file }.open());
    }

    /// pyBigWig #111: a file whose header summary was all zeros was read as empty.
    #[test]
    fn data_is_readable_when_header_summary_is_all_zero() {
        let file =
            TestFile { total_summary: TotalSummaryStyle::Given(ZERO_SUMMARY), ..sample_file() };
        let mut reader = file.open();
        assert_eq!(reader.info().total_summary, Some(ZERO_SUMMARY));
        assert_matches_file(&mut reader, &file);
    }

    /// kent #25, pybbi #14 and jbrowse-components #5502: ENCODE files carry an R-tree leaf size
    /// of 0xffffffff8855ee2f, which kent tried to allocate.
    #[test]
    fn oversized_block_size_in_rtree_leaf_is_a_corrupt_error_not_an_allocation() {
        let mut bytes = sample_file().build();
        let full_index_offset =
            usize::try_from(u64::from_le_bytes(bytes[24..32].try_into().unwrap())).unwrap();
        let second_leaf_size_field = full_index_offset + RTREE_HEADER_SIZE + 4 + 32 + 24;
        bytes[second_leaf_size_field..][..8]
            .copy_from_slice(&0xFFFF_FFFF_8855_EE2F_u64.to_le_bytes());
        let mut reader = BigWigReader::new(Cursor::new(bytes)).unwrap();
        assert!(is_corrupt(&reader.items("chr1", 0, 1000)));
        assert!(is_corrupt(&reader.values("chr1", 0, 1000)));
        // Blocks the query does not touch are not read.
        assert_eq!(collect(reader.items("chr2", 0, 500).unwrap()).len(), 2);
    }

    /// libBigWig #69: a corrupt block mid-file crashed the reader. A bad deflate stream is
    /// reported as [`Error::Decompress`], and the handle keeps working.
    #[test]
    fn corrupt_deflate_stream_in_a_later_block_returns_decompress_error() {
        let file = sample_file();
        let mut bytes = file.build();
        let second = data_blocks(&bytes)[1];
        let block = usize::try_from(second.data_offset).unwrap();
        bytes[block..block + usize::try_from(second.data_size).unwrap()].fill(0xFF);
        let mut reader = BigWigReader::new(Cursor::new(bytes)).unwrap();

        let mut items = reader.items("chr1", 0, 1000).unwrap();
        for expected in &file.sections[0].items() {
            assert_eq!(&items.next().unwrap().unwrap(), expected);
        }
        let error = items.next().unwrap();
        assert!(
            matches!(error, Err(Error::Decompress { offset, .. }) if offset == second.data_offset)
        );
        assert!(items.next().is_none());

        assert_eq!(collect(reader.items("chr2", 0, 500).unwrap()), file.sections[3].items());
        assert_eq!(collect(reader.items("chr1", 0, 50).unwrap()), file.sections[0].items());
        assert!(matches!(reader.values("chr1", 0, 1000), Err(Error::Decompress { .. })));
    }

    /// pyBigWig #105 and libBigWig #69: truncated files crashed readers.
    #[test]
    fn truncated_file_returns_an_error_at_every_cut_point() {
        // Every node is full, so every byte but the trailing magic is read by some query.
        let sections = vec![
            bed_graph(0, &[(0, 5, 1.0), (20, 25, 2.0)]),
            bed_graph(0, &[(40, 45, 3.0), (60, 65, 4.0)]),
        ];
        let file = TestFile {
            zoom_reductions: vec![10],
            zoom_block_records: 1,
            rtree_block_size: 2,
            ..TestFile::new(&[("chr1", 100)], sections)
        };
        let bytes = file.build();
        exercise(bytes.clone()).unwrap();
        let magic_start = bytes.len() - 4;
        for cut in 0..4 {
            let result = exercise(bytes[..cut].to_vec());
            assert!(is_not_bigwig(&result), "cut at {cut}: {result:?}");
        }
        for cut in 4..magic_start {
            let result = exercise(bytes[..cut].to_vec());
            assert!(is_corrupt(&result), "cut at {cut}: {result:?}");
        }
        for cut in magic_start..bytes.len() {
            exercise(bytes[..cut].to_vec()).unwrap();
        }
    }

    /// libBigWig #35: a stale buffer made results depend on earlier queries.
    #[test]
    fn region_queries_in_any_order_match_a_fresh_handle() {
        let file = TestFile {
            rtree_block_size: 2,
            zoom_reductions: vec![10],
            ..TestFile::new(&[("chr1", 5000), ("chr2", 3000)], tiled_sections(&[0, 1], 30))
        };
        let bytes = file.build();
        let queries = [
            ("chr1", 0, 100),
            ("chr2", 1000, 2000),
            ("chr1", 2900, 3000),
            ("chr1", 50, 60),
            ("chr2", 0, 3000),
            ("chr1", 2500, 2600),
        ];
        let answer = |reader: &mut BigWigReader<Cursor<Vec<u8>>>,
                      (chrom, start, end): (&str, u32, u32)| {
            let items = collect(reader.items(chrom, start, end).unwrap());
            let values: Vec<u32> = reader
                .values(chrom, start, end)
                .unwrap()
                .iter()
                .map(|value| value.to_bits())
                .collect();
            let summaries =
                reader.summaries(chrom, start, end, 3, SummarySource::BestZoom).unwrap();
            (items, values, summaries)
        };
        let fresh: Vec<_> = queries
            .iter()
            .map(|&query| {
                answer(&mut BigWigReader::new(Cursor::new(bytes.clone())).unwrap(), query)
            })
            .collect();
        for order in [[0, 1, 2, 3, 4, 5], [5, 4, 3, 2, 1, 0], [3, 0, 5, 1, 4, 2]] {
            let mut reader = BigWigReader::new(Cursor::new(bytes.clone())).unwrap();
            for index in order {
                assert_eq!(answer(&mut reader, queries[index]), fresh[index], "query {index}");
            }
        }
    }

    /// bigtools #95: libBigWig writes empty files with `fullIndexOffset` 0 and no zoom levels,
    /// which other readers rejected.
    #[test]
    fn file_without_rtree_opens_and_yields_empty_results() {
        // libBigWig's empty-file total summary, with its DBL_MAX minimum and DBL_MIN maximum.
        let stored = Summary {
            covered_bases: 0.0,
            min: f64::MAX,
            max: f64::MIN_POSITIVE,
            sum: 0.0,
            sum_squares: 0.0,
        };
        let file = TestFile {
            data_index: DataIndexStyle::OffsetZero,
            total_summary: TotalSummaryStyle::Given(stored),
            ..TestFile::new(&[("chr1", 100), ("chr2", 50)], Vec::new())
        };
        let mut reader = file.open();
        assert_eq!(reader.info().total_summary, Some(stored));
        assert_eq!(collect(reader.all_items().unwrap()), []);
        assert_eq!(collect(reader.items("chr1", 0, 100).unwrap()), []);
        assert!(reader.values("chr2", 0, 50).unwrap().iter().all(|value| value.is_nan()));
        assert_eq!(reader.summary("chr1", 0, 100, SummarySource::BestZoom).unwrap(), None);
    }

    /// deepTools #598: kent's tools write a header-only file for empty input: a chromosome tree
    /// and an index with no root node, and no total summary.
    #[test]
    fn file_written_by_kent_for_empty_input_opens_as_empty() {
        let file = TestFile {
            data_index: DataIndexStyle::HeaderOnly,
            total_summary: TotalSummaryStyle::Absent,
            ..TestFile::new(&[], Vec::new())
        };
        let mut reader = file.open();
        assert_eq!(reader.chroms(), []);
        assert_eq!(reader.info().total_summary, None);
        assert_eq!(reader.info().section_count, 0);
        assert_eq!(collect(reader.all_items().unwrap()), []);
        assert_eq!(reader.compute_total_summary().unwrap(), ZERO_SUMMARY);
    }

    /// bbi-js #74: a chromosome tree with no items was rejected.
    #[test]
    fn reader_opens_a_chrom_tree_with_zero_items() {
        let mut reader = TestFile::new(&[], Vec::new()).open();
        assert_eq!(reader.chroms(), []);
        assert_eq!(collect(reader.all_items().unwrap()), []);
        assert!(matches!(reader.items("chr1", 0, 10), Err(Error::UnknownChrom(_))));
    }

    /// pyBigWig #85 and #158: a data index with no items was misread.
    #[test]
    fn reader_opens_a_data_index_with_zero_items() {
        let mut reader = TestFile::new(&[("chr1", 100)], Vec::new()).open();
        assert_eq!(collect(reader.items("chr1", 0, 100).unwrap()), []);
        assert!(reader.values("chr1", 0, 100).unwrap().iter().all(|value| value.is_nan()));
        assert_eq!(reader.summary("chr1", 0, 100, SummarySource::Raw).unwrap(), None);
    }

    proptest::proptest! {
        /// Arbitrary corruption gives errors, never a panic or unbounded work.
        #[test]
        fn reading_mutated_files_never_panics(
            mutations in proptest::collection::vec(
                (proptest::prelude::any::<proptest::sample::Index>(), proptest::prelude::any::<u8>()),
                1..8,
            ),
            keep in proptest::prelude::any::<proptest::sample::Index>(),
        ) {
            let file = TestFile {
                zoom_reductions: vec![10, 100],
                rtree_block_size: 2,
                ..sample_file()
            };
            let mut bytes = file.build();
            for (index, byte) in mutations {
                let at = index.index(bytes.len());
                bytes[at] = byte;
            }
            bytes.truncate(keep.index(bytes.len()) + 1);
            let _ = exercise(bytes);
        }
    }
}
