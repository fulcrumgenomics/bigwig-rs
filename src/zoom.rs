//! Zoom levels and summaries: the [`Summary`] aggregate, the [`ZoomRecord`] and its on-disk form,
//! the [`ZoomPolicy`] that chooses a file's reductions, grid-aligned zoom-record accumulation when
//! writing, and exact and zoom-based summaries when reading.

use crate::bbi::{ByteOrder, ByteReader};
use crate::error::{Error, Result};

/// Largest reduction any level may have: IGV reads `reductionLevel` as a signed 32-bit int, and
/// levels above it broke IGV (bigtools #85).
pub(crate) const MAX_REDUCTION: u32 = i32::MAX.unsigned_abs();
/// Number of leading items the `Auto` policy measures the data's pitch from.
pub(crate) const PITCH_SAMPLE_ITEMS: usize = 65_536;
/// Most levels an `Auto` ladder has, as kent and bigtools write at most.
pub(crate) const AUTO_MAX_LEVELS: usize = 10;
/// Ratio between successive `Auto` levels (kent's `bbiResIncrement`).
const AUTO_LEVEL_FACTOR: u64 = 4;
/// The `Auto` first level is at least this many pitches (kent's "10 × average resolution").
const AUTO_PITCH_MULTIPLIER: u64 = 10;
/// `min_first_reduction` of the default policy. On per-base data kent's `wigToBigWig` starts near
/// 35 bp, a level over half the size of the data itself; a ladder from 256 bp adds about 11%.
const DEFAULT_MIN_FIRST_REDUCTION: u32 = 256;

/// How the writer chooses zoom levels, each given by its reduction (bases per record).
///
/// Records sit on a fixed grid: bin `k` of a level with reduction `r` is `[k·r, (k+1)·r)`, clipped
/// to the chromosome, and a record is written for every bin with at least one covered base, over
/// the tight extent of those bases. Whatever the policy, the writer keeps levels in increasing
/// order, drops levels longer than the longest chromosome (unless that would leave none, as the
/// variants below describe), and stops at the first level whose record count is not smaller than
/// the previous level's. Except with [`ZoomPolicy::None`], at least one level is always written,
/// even for a file with no data. No reduction exceeds 2³¹ − 1, because IGV reads reductions as
/// signed 32-bit integers.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ZoomPolicy {
    /// A ×4 ladder derived from the data. The default is `min_first_reduction: 256`.
    ///
    /// The writer measures the mean pitch `p` of the first 65,536 items: their covered extent
    /// (first start to last end on each chromosome, summed) divided by the item count. Per-base
    /// data has `p = 1`, 1 kb bins `p = 1000`. The first reduction is
    /// `max(min_first_reduction, ⌈10·p⌉)`, where `⌈10·p⌉` is capped at the longest chromosome's
    /// length and at 2³¹ − 1. Each further level is 4× the previous one, up to 10 levels, stopping
    /// before any level longer than the longest chromosome or than 2³¹ − 1.
    ///
    /// A file with no data gets the single level `min_first_reduction`. So does a file whose
    /// chromosomes are all shorter than `min_first_reduction`: that level is kept although it is
    /// longer than every chromosome, because there must be at least one.
    ///
    /// Lower `min_first_reduction` (down to 1) for a finer first level, e.g. 32 for kent-like
    /// levels on per-base data.
    Auto {
        /// Smallest first reduction, in bases: 1 to 2³¹ − 1.
        min_first_reduction: u32,
    },
    /// Exactly these reductions, in bases: non-empty, strictly increasing, each from 1 to 2³¹ − 1.
    ///
    /// Levels longer than the longest chromosome are dropped. If every level is, the smallest is
    /// kept as the only level. A level whose reduction is a multiple of a finer level's is
    /// computed from that level's records; any other (e.g. 65 in kent's `30, 65, 130, …`) is
    /// accumulated from the raw data, at about one extra update per value.
    Explicit(Vec<u32>),
    /// No zoom levels. Smaller files, but browsers must read raw data at every scale, and
    /// IGV 2.17.0–2.19.4 throw an `ArrayIndexOutOfBoundsException` on a file with no zoom levels.
    /// Not recommended for anything a browser will display.
    None,
}

impl ZoomPolicy {
    /// Checks the policy's parameters, returning [`Error::InvalidOption`] naming the problem.
    pub(crate) fn validate(&self) -> Result<()> {
        match self {
            Self::Auto { min_first_reduction } => {
                if *min_first_reduction == 0 {
                    return Err(invalid_zoom("Auto min_first_reduction must be at least 1"));
                }
                if *min_first_reduction > MAX_REDUCTION {
                    return Err(invalid_zoom(&format!(
                        "Auto min_first_reduction must be at most {MAX_REDUCTION}, got \
                         {min_first_reduction}"
                    )));
                }
            }
            Self::Explicit(reductions) => {
                if reductions.is_empty() {
                    return Err(invalid_zoom("Explicit reductions must not be empty"));
                }
                if reductions.contains(&0) {
                    return Err(invalid_zoom("Explicit reductions must each be at least 1"));
                }
                if let Some(too_large) = reductions.iter().find(|&&r| r > MAX_REDUCTION) {
                    return Err(invalid_zoom(&format!(
                        "Explicit reductions must each be at most {MAX_REDUCTION}, got {too_large}"
                    )));
                }
                if let Some(&[previous, next]) =
                    reductions.windows(2).find(|pair| pair[1] <= pair[0])
                {
                    return Err(invalid_zoom(&format!(
                        "Explicit reductions must be strictly increasing, but {next} follows \
                         {previous}"
                    )));
                }
            }
            Self::None => {}
        }
        Ok(())
    }

    /// Whether the reductions depend on the data, so the writer must offer the leading items to
    /// a [`PitchSampler`] (and keep them for replay) before zoom accumulation can start.
    pub(crate) fn samples_pitch(&self) -> bool {
        matches!(self, Self::Auto { .. })
    }

    /// The reductions to accumulate, strictly increasing: empty for `None`, at least one level
    /// otherwise.
    ///
    /// `sample` is read only by `Auto`, and must hold every item offered so far; an empty sample
    /// means the file has no data. `longest_chrom_length` is the length of the longest chromosome
    /// in the file's chromosome list. The policy must have passed [`ZoomPolicy::validate`].
    pub(crate) fn reductions(&self, sample: &PitchSampler, longest_chrom_length: u32) -> Vec<u32> {
        match self {
            Self::Auto { min_first_reduction } => {
                auto_ladder(*min_first_reduction, sample, longest_chrom_length)
            }
            Self::Explicit(reductions) => {
                let limit = longest_chrom_length.min(MAX_REDUCTION);
                let fitting: Vec<u32> = reductions
                    .iter()
                    .copied()
                    .take_while(|&reduction| reduction <= limit)
                    .collect();
                if fitting.is_empty() {
                    reductions.iter().copied().take(1).collect()
                } else {
                    fitting
                }
            }
            Self::None => Vec::new(),
        }
    }
}

impl Default for ZoomPolicy {
    fn default() -> Self {
        Self::Auto { min_first_reduction: DEFAULT_MIN_FIRST_REDUCTION }
    }
}

/// Measures the mean pitch of the first [`PITCH_SAMPLE_ITEMS`] items, for the `Auto` policy.
///
/// The pitch is the covered extent (first start to last end on each chromosome, summed) divided
/// by the number of items. The sampler only measures: the writer keeps the items it offers and
/// replays them through a [`ZoomBuilder`] once [`ZoomPolicy::reductions`] has fixed the levels.
/// Items must be offered in the order they are written: sorted within each chromosome, with each
/// chromosome's items contiguous.
#[derive(Clone, Debug, Default)]
pub(crate) struct PitchSampler {
    /// Items sampled so far, at most [`PITCH_SAMPLE_ITEMS`].
    items: usize,
    /// Summed extents of the chromosomes sampled before the current one.
    finished_extent: u64,
    /// Extent of the chromosome currently being sampled.
    current: Option<ChromExtent>,
}

impl PitchSampler {
    /// Creates an empty sampler.
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Offers one item covering `[start, end)`. Returns whether it was sampled, which is false
    /// once the sample is full.
    pub(crate) fn offer_item(&mut self, chrom_id: u32, start: u32, end: u32) -> bool {
        if self.is_full() {
            return false;
        }
        self.extend(chrom_id, u64::from(start), u64::from(end));
        self.items += 1;
        true
    }

    /// Offers a fixedStep run of `count` items, item `i` covering
    /// `[start + i·step, start + i·step + span)`. Returns how many leading items were sampled.
    pub(crate) fn offer_fixed_step(
        &mut self,
        chrom_id: u32,
        start: u32,
        step: u32,
        span: u32,
        count: usize,
    ) -> usize {
        let taken = count.min(PITCH_SAMPLE_ITEMS - self.items);
        if taken > 0 {
            let last_start = u64::from(start) + (taken as u64 - 1) * u64::from(step);
            self.extend(chrom_id, u64::from(start), last_start + u64::from(span));
            self.items += taken;
        }
        taken
    }

    /// Whether [`PITCH_SAMPLE_ITEMS`] items have been sampled, so the pitch is final.
    pub(crate) fn is_full(&self) -> bool {
        self.items >= PITCH_SAMPLE_ITEMS
    }

    /// `⌈10·p⌉` for the sampled pitch `p`, in exact integer arithmetic; `None` with no items.
    fn ten_pitches(&self) -> Option<u64> {
        let items = self.items as u64;
        let extent = self.finished_extent + self.current.map_or(0, ChromExtent::extent);
        (items > 0).then(|| (AUTO_PITCH_MULTIPLIER * extent).div_ceil(items))
    }

    /// Grows the current chromosome's extent to `end`, or starts a new chromosome at `start`.
    fn extend(&mut self, chrom_id: u32, start: u64, end: u64) {
        match &mut self.current {
            Some(current) if current.chrom_id == chrom_id => current.end = end,
            current => {
                if let Some(previous) = current.replace(ChromExtent { chrom_id, start, end }) {
                    self.finished_extent += previous.extent();
                }
            }
        }
    }
}

/// The span of one chromosome's sampled items, from the first start to the last end.
#[derive(Clone, Copy, Debug)]
struct ChromExtent {
    /// Chromosome the items are on.
    chrom_id: u32,
    /// Start of the first item.
    start: u64,
    /// End of the last item.
    end: u64,
}

impl ChromExtent {
    /// Length of the extent in bases.
    fn extent(self) -> u64 {
        self.end - self.start
    }
}

/// Accumulates the grid-aligned zoom records of every level of a file as its data streams in.
///
/// Items must arrive sorted and non-overlapping within each chromosome, with each chromosome's
/// items contiguous, ending within the chromosome; the writer validates all of this before
/// pushing. Because records have tight extents, a record never extends past the items it
/// summarises, so the builder needs no chromosome lengths.
///
/// Completed records leave through the `emit` callback each push takes, as
/// `emit(level_index, record)`, the moment their bin can receive no more data. Each level's
/// records arrive in (chromosome, start) order, so the writer can pack them straight into blocks,
/// starting a new block when `record.chrom_id` changes. Memory is one open record per level,
/// whatever the volume of data. An error from `emit` is returned from the push; it leaves the
/// builder partly updated, which is acceptable only because such errors make the writer unusable.
///
/// The finest level is accumulated from the raw items. A coarser level whose reduction is a
/// multiple of a finer level's merges that level's records, which is exact because the grid bins
/// nest; any other level is accumulated from the raw items too.
#[derive(Debug)]
pub(crate) struct ZoomBuilder {
    /// One accumulator per level, in increasing reduction order.
    levels: Vec<LevelBuilder>,
    /// Chromosome of the records currently open; `None` before the first item and after `close`.
    chrom_id: Option<u32>,
}

impl ZoomBuilder {
    /// Creates a builder for `reductions`, which must be strictly increasing and each at least 1,
    /// as [`ZoomPolicy::reductions`] returns them. The emitted level indices index this slice.
    pub(crate) fn new(reductions: &[u32]) -> Self {
        debug_assert!(reductions.first().is_none_or(|&first| first >= 1));
        debug_assert!(reductions.windows(2).all(|pair| pair[0] < pair[1]));
        let levels = reductions
            .iter()
            .enumerate()
            .map(|(index, &reduction)| {
                // The coarsest dividing level has the fewest records to merge.
                let source = (0..index).rev().find(|&finer| reduction % reductions[finer] == 0);
                LevelBuilder::new(reduction, source)
            })
            .collect();
        Self { levels, chrom_id: None }
    }

    /// Adds one item: `value` on every base of `[start, end)` of chromosome `chrom_id`.
    ///
    /// An item on a different chromosome from the previous one first closes every open record.
    pub(crate) fn push_item<E>(
        &mut self,
        chrom_id: u32,
        start: u32,
        end: u32,
        value: f32,
        emit: &mut E,
    ) -> Result<()>
    where
        E: FnMut(usize, ZoomRecord) -> Result<()>,
    {
        debug_assert!(start < end);
        self.enter_chrom(chrom_id, emit)?;
        for level_index in 0..self.levels.len() {
            if self.levels[level_index].source.is_none() {
                self.add_item_to_level(level_index, chrom_id, start, end, value, emit)?;
            }
        }
        Ok(())
    }

    /// Adds a fixedStep run on chromosome `chrom_id`: `values[i]` on every base of
    /// `[start + i·step, start + i·step + span)`, with `1 ≤ span ≤ step`.
    ///
    /// Values falling wholly inside one bin are added in a single pass over the slice, so dense
    /// per-base data costs a few floating-point operations per value. An item on a different
    /// chromosome from the previous one first closes every open record.
    pub(crate) fn push_fixed_step<E>(
        &mut self,
        chrom_id: u32,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
        emit: &mut E,
    ) -> Result<()>
    where
        E: FnMut(usize, ZoomRecord) -> Result<()>,
    {
        debug_assert!(span >= 1 && span <= step);
        if values.is_empty() {
            return Ok(());
        }
        self.enter_chrom(chrom_id, emit)?;
        for level_index in 0..self.levels.len() {
            if self.levels[level_index].source.is_none() {
                self.add_fixed_step_to_level(
                    level_index,
                    chrom_id,
                    start,
                    step,
                    span,
                    values,
                    emit,
                )?;
            }
        }
        Ok(())
    }

    /// Completes and emits every open record, finest level first so that each record also
    /// reaches the levels computed from it. Call at the end of the data; the writer may also call
    /// it on a chromosome change, which pushes do implicitly.
    pub(crate) fn close<E>(&mut self, emit: &mut E) -> Result<()>
    where
        E: FnMut(usize, ZoomRecord) -> Result<()>,
    {
        for level_index in 0..self.levels.len() {
            if let Some(record) = self.levels[level_index].take_record() {
                self.complete(level_index, record, emit)?;
            }
        }
        self.chrom_id = None;
        Ok(())
    }

    /// Closes every open record if `chrom_id` differs from the chromosome being accumulated.
    fn enter_chrom<E>(&mut self, chrom_id: u32, emit: &mut E) -> Result<()>
    where
        E: FnMut(usize, ZoomRecord) -> Result<()>,
    {
        if self.chrom_id != Some(chrom_id) {
            self.close(emit)?;
            self.chrom_id = Some(chrom_id);
        }
        Ok(())
    }

    /// Adds a fixedStep run to the raw-sourced level `level_index`.
    #[expect(clippy::too_many_arguments, reason = "a fixedStep run plus the level and callback")]
    fn add_fixed_step_to_level<E>(
        &mut self,
        level_index: usize,
        chrom_id: u32,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
        emit: &mut E,
    ) -> Result<()>
    where
        E: FnMut(usize, ZoomRecord) -> Result<()>,
    {
        let (step, span64) = (u64::from(step), u64::from(span));
        let mut index = 0;
        while index < values.len() {
            let item_start = position_u32(u64::from(start) + index as u64 * step);
            if let Some(record) = self.levels[level_index].enter_bin(chrom_id, item_start) {
                self.complete(level_index, record, emit)?;
            }
            let level = &mut self.levels[level_index];
            let item_end = u64::from(item_start) + span64;
            if item_end <= level.bin_end {
                // Items index.. whose ends stay inside the bin: the j-th ends at item_end + j·step.
                let in_bin = usize::try_from((level.bin_end - item_end) / step + 1)
                    .map_or(values.len() - index, |in_bin| in_bin.min(values.len() - index));
                let last_end = item_end + (in_bin as u64 - 1) * step;
                level.add_values(&values[index..index + in_bin], span, position_u32(last_end));
                index += in_bin;
            } else {
                let item_end = position_u32(item_end);
                self.add_item_to_level(
                    level_index,
                    chrom_id,
                    item_start,
                    item_end,
                    values[index],
                    emit,
                )?;
                index += 1;
            }
        }
        Ok(())
    }

    /// Adds one item to the raw-sourced level `level_index`, splitting it at bin edges.
    fn add_item_to_level<E>(
        &mut self,
        level_index: usize,
        chrom_id: u32,
        start: u32,
        end: u32,
        value: f32,
        emit: &mut E,
    ) -> Result<()>
    where
        E: FnMut(usize, ZoomRecord) -> Result<()>,
    {
        let mut position = start;
        while position < end {
            if let Some(record) = self.levels[level_index].enter_bin(chrom_id, position) {
                self.complete(level_index, record, emit)?;
            }
            let level = &mut self.levels[level_index];
            // A bin end that doesn't fit in u32 lies beyond every item end.
            let piece_end = u32::try_from(level.bin_end).map_or(end, |bin_end| bin_end.min(end));
            level.add_run(piece_end - position, value, piece_end);
            position = piece_end;
        }
        Ok(())
    }

    /// Emits a completed record of level `level_index`, then merges it into each level computed
    /// from that level, completing their records in turn.
    fn complete<E>(&mut self, level_index: usize, record: ZoomRecord, emit: &mut E) -> Result<()>
    where
        E: FnMut(usize, ZoomRecord) -> Result<()>,
    {
        emit(level_index, record)?;
        for coarser_index in level_index + 1..self.levels.len() {
            if self.levels[coarser_index].source != Some(level_index) {
                continue;
            }
            if let Some(completed) = self.levels[coarser_index].merge(record) {
                self.complete(coarser_index, completed, emit)?;
            }
        }
        Ok(())
    }
}

/// One level of a [`ZoomBuilder`]: its reduction, where its records come from, and its open
/// record.
#[derive(Debug)]
struct LevelBuilder {
    /// Bases per bin.
    reduction: u32,
    /// The finer level this level merges records from, or `None` for raw items.
    source: Option<usize>,
    /// End of the open record's bin (unclipped, so it may pass `u32::MAX`); 0 when none is open.
    bin_end: u64,
    /// Chromosome of the open record.
    chrom_id: u32,
    /// First covered base of the open record.
    start: u32,
    /// One past the last covered base of the open record.
    end: u32,
    /// Aggregate of the open record's covered bases.
    running: RunningSummary,
}

impl LevelBuilder {
    /// Creates a level with no open record.
    fn new(reduction: u32, source: Option<usize>) -> Self {
        Self {
            reduction,
            source,
            bin_end: 0,
            chrom_id: 0,
            start: 0,
            end: 0,
            running: RunningSummary::EMPTY,
        }
    }

    /// Merges a completed record of the source level, whose bin nests inside one of this
    /// level's bins. Returns this level's previous record if the new one starts a new bin.
    fn merge(&mut self, record: ZoomRecord) -> Option<ZoomRecord> {
        let completed = self.enter_bin(record.chrom_id, record.start);
        debug_assert!(u64::from(record.end) <= self.bin_end);
        self.running.add_weighted(&record.summary, 1.0);
        self.end = record.end;
        completed
    }

    /// Makes the bin containing `position` the open one, starting its record at `position`.
    /// Returns the previous record if it was in another bin; positions must not decrease.
    fn enter_bin(&mut self, chrom_id: u32, position: u32) -> Option<ZoomRecord> {
        if u64::from(position) < self.bin_end {
            return None;
        }
        let completed = self.take_record();
        let reduction = u64::from(self.reduction);
        self.bin_end = (u64::from(position) / reduction + 1) * reduction;
        self.chrom_id = chrom_id;
        self.start = position;
        self.end = position;
        completed
    }

    /// Adds `bases` bases of `value` ending at `end` to the open record.
    fn add_run(&mut self, bases: u32, value: f32, end: u32) {
        self.running.add_run(bases, value);
        self.end = end;
    }

    /// Adds fixedStep values of `span` bases each, the last ending at `end`, to the open record.
    fn add_values(&mut self, values: &[f32], span: u32, end: u32) {
        self.running.add_values(values, span);
        self.end = end;
    }

    /// Closes the open record, returning it unless it covered nothing.
    fn take_record(&mut self) -> Option<ZoomRecord> {
        let summary = self.running.to_summary();
        self.running = RunningSummary::EMPTY;
        self.bin_end = 0;
        summary.map(|summary| ZoomRecord {
            chrom_id: self.chrom_id,
            start: self.start,
            end: self.end,
            summary,
        })
    }
}

/// One zoom record: a summary of the covered bases in `[start, end)` on chromosome `chrom_id`.
///
/// Records written by this crate each cover one bin of their level's grid, and `[start, end)` is
/// the tight extent of the bin's covered bases. Other writers' records may extend over uncovered
/// bases: kent's run to `start + reduction`. On disk the summary's values are f32, and
/// `covered_bases` is the integer `validCount`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZoomRecord {
    /// Index into the file's chromosome list (the file's chromId).
    pub chrom_id: u32,
    /// Start of the record's extent (0-based).
    pub start: u32,
    /// One past the end of the record's extent.
    pub end: u32,
    /// Aggregate over the covered bases in `[start, end)`.
    pub summary: Summary,
}

impl ZoomRecord {
    /// Size of an encoded record in bytes.
    pub(crate) const ENCODED_LEN: usize = 32;

    /// Appends the 32-byte on-disk form: chromId, start, end, validCount (u32 each), then
    /// minVal, maxVal, sumData and sumSquares (f32 each), in `byte_order`.
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "the format stores f32 values and an integer validCount; records this crate \
                  writes have integral covered_bases of at most the reduction"
    )]
    pub(crate) fn encode(&self, byte_order: ByteOrder, buf: &mut Vec<u8>) {
        let summary = &self.summary;
        byte_order.put_u32(buf, self.chrom_id);
        byte_order.put_u32(buf, self.start);
        byte_order.put_u32(buf, self.end);
        byte_order.put_u32(buf, summary.covered_bases as u32);
        byte_order.put_f32(buf, summary.min as f32);
        byte_order.put_f32(buf, summary.max as f32);
        byte_order.put_f32(buf, summary.sum as f32);
        byte_order.put_f32(buf, summary.sum_squares as f32);
    }

    /// Decodes a record from its 32-byte on-disk form; a truncated record is
    /// [`Error::Corrupt`].
    pub(crate) fn decode(reader: &mut ByteReader<'_>) -> Result<Self> {
        let chrom_id = reader.u32()?;
        let start = reader.u32()?;
        let end = reader.u32()?;
        let valid_count = reader.u32()?;
        let min = reader.f32()?;
        let max = reader.f32()?;
        let sum = reader.f32()?;
        let sum_squares = reader.f32()?;
        let summary = Summary {
            covered_bases: f64::from(valid_count),
            min: f64::from(min),
            max: f64::from(max),
            sum: f64::from(sum),
            sum_squares: f64::from(sum_squares),
        };
        Ok(Self { chrom_id, start, end, summary })
    }
}

/// Summary of one region from zoom records, pro-rated as kent's `bbiSummarySlice` does but
/// without its final `normalizeCount` step.
///
/// A record overlapping the region contributes the fraction
/// `f = |record ∩ region| / (record.end − record.start)` of its `covered_bases`, `sum` and
/// `sum_squares`, and its `min` and `max` unweighted. Records with no covered bases are skipped.
/// The result is an approximation unless records lie wholly inside the region.
///
/// kent then rounds the pro-rated coverage `c` up to `⌈c⌉` and scales `sum` and `sum_squares` by
/// `⌈c⌉ / c`. Here `c` is kept as it is, so `covered_bases`, `sum` and `sum_squares` can be smaller
/// than kent's by that factor; the mean, `min` and `max` agree.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ZoomSummaryAccumulator {
    /// Start of the region.
    start: u32,
    /// End of the region (exclusive).
    end: u32,
    /// Aggregate of the contributions so far.
    running: RunningSummary,
}

impl ZoomSummaryAccumulator {
    /// Creates an accumulator for the region `[start, end)`.
    pub(crate) fn new(start: u32, end: u32) -> Self {
        Self { start, end, running: RunningSummary::EMPTY }
    }

    /// Adds the part of `record` that overlaps the region. The caller filters by chromosome.
    pub(crate) fn add_record(&mut self, record: &ZoomRecord) {
        let overlap_start = record.start.max(self.start);
        let overlap_end = record.end.min(self.end);
        // A positive overlap implies record.end > record.start, so the division is safe.
        if overlap_end <= overlap_start || record.summary.covered_bases <= 0.0 {
            return;
        }
        let fraction =
            f64::from(overlap_end - overlap_start) / f64::from(record.end - record.start);
        self.running.add_weighted(&record.summary, fraction);
    }

    /// The pro-rated summary, or `None` if no record contributed.
    pub(crate) fn summary(&self) -> Option<Summary> {
        self.running.to_summary()
    }
}

/// Exact summary of the raw items within one region, accumulated in f64.
///
/// Serves both the reader's raw summaries, with a region, and the writer's total summary,
/// [`unbounded`](Self::unbounded). Items are clipped to the region; the caller filters by
/// chromosome.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RawSummaryAccumulator {
    /// Start of the region.
    start: u64,
    /// End of the region (exclusive).
    end: u64,
    /// Aggregate of the covered bases so far.
    running: RunningSummary,
}

impl RawSummaryAccumulator {
    /// Creates an accumulator for the region `[start, end)`.
    pub(crate) fn new(start: u32, end: u32) -> Self {
        Self { start: u64::from(start), end: u64::from(end), running: RunningSummary::EMPTY }
    }

    /// Creates an accumulator that counts every base of every item, for a total summary.
    pub(crate) fn unbounded() -> Self {
        Self { start: 0, end: u64::MAX, running: RunningSummary::EMPTY }
    }

    /// Adds the bases of the item `[start, end)` with `value` that lie in the region.
    pub(crate) fn add_item(&mut self, start: u32, end: u32, value: f32) {
        self.add_clipped(u64::from(start), u64::from(end), value);
    }

    /// Adds the bases in the region of a fixedStep run: `values[i]` on every base of
    /// `[start + i·step, start + i·step + span)`, with `1 ≤ span ≤ step`.
    pub(crate) fn add_fixed_step(&mut self, start: u32, step: u32, span: u32, values: &[f32]) {
        let (start, step, span64) = (u64::from(start), u64::from(step), u64::from(span));
        let Some(last_index) = (values.len() as u64).checked_sub(1) else { return };
        if start >= self.start && start + last_index * step + span64 <= self.end {
            self.running.add_values(values, span);
        } else {
            for (index, &value) in values.iter().enumerate() {
                let item_start = start + index as u64 * step;
                self.add_clipped(item_start, item_start + span64, value);
            }
        }
    }

    /// The exact summary, or `None` if no base in the region is covered.
    pub(crate) fn summary(&self) -> Option<Summary> {
        self.running.to_summary()
    }

    /// Adds the bases of `[start, end)` with `value` that lie in the region.
    #[expect(
        clippy::cast_possible_truncation,
        reason = "an item clipped to a region is no longer than the item, whose ends are u32"
    )]
    fn add_clipped(&mut self, start: u64, end: u64, value: f32) {
        let (start, end) = (start.max(self.start), end.min(self.end));
        if start < end {
            self.running.add_run((end - start) as u32, value);
        }
    }
}

/// Aggregate over the covered bases of a region, a zoom record, or the whole file.
///
/// For raw data `covered_bases` is exact (integral). For summaries computed from zoom records it
/// is pro-rated by overlap and may be fractional.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Summary {
    /// Bases with data. Exact (integral) for raw data; pro-rated (possibly fractional) for zoom.
    pub covered_bases: f64,
    /// Smallest value over the covered bases.
    pub min: f64,
    /// Largest value over the covered bases.
    pub max: f64,
    /// Sum of the values over the covered bases.
    pub sum: f64,
    /// Sum of the squared values over the covered bases.
    pub sum_squares: f64,
}

impl Summary {
    /// Mean over the covered bases: `sum / covered_bases`. NaN when no base is covered.
    pub fn mean(&self) -> f64 {
        self.sum / self.covered_bases
    }

    /// Sample standard deviation over the covered bases:
    /// `sqrt((sum_squares - sum²/n) / (n - 1))`, as kent's `calcStdFromSums`, libBigWig and
    /// pyBigWig compute it. It is 0 when `n <= 1`, and a variance made negative by rounding is
    /// treated as 0.
    pub fn std(&self) -> f64 {
        let n = self.covered_bases;
        if n <= 1.0 {
            return 0.0;
        }
        let variance = (self.sum_squares - self.sum * self.sum / n) / (n - 1.0);
        variance.max(0.0).sqrt()
    }

    /// Fraction of a region of `region_length` bases that is covered, in `[0, 1]`.
    pub fn coverage(&self, region_length: u64) -> f64 {
        self.covered_bases / length_as_f64(region_length)
    }

    /// Mean over a region of `region_length` bases with uncovered bases counted as 0 (kent's
    /// "mean0"): `sum / region_length`.
    pub fn mean0(&self, region_length: u64) -> f64 {
        self.sum / length_as_f64(region_length)
    }
}

/// A [`Summary`] being accumulated in f64, with `min` and `max` updated independently (unlike
/// libBigWig's `else if`, which loses the first value's maximum).
#[derive(Clone, Copy, Debug)]
struct RunningSummary {
    /// Bases covered so far (possibly fractional when pro-rating).
    covered_bases: f64,
    /// Smallest value so far; +∞ before the first.
    min: f64,
    /// Largest value so far; −∞ before the first.
    max: f64,
    /// Sum of the values over the covered bases.
    sum: f64,
    /// Sum of the squared values over the covered bases.
    sum_squares: f64,
}

impl RunningSummary {
    /// Nothing covered yet.
    const EMPTY: Self = Self {
        covered_bases: 0.0,
        min: f64::INFINITY,
        max: f64::NEG_INFINITY,
        sum: 0.0,
        sum_squares: 0.0,
    };

    /// Adds `bases` bases that all have `value`.
    fn add_run(&mut self, bases: u32, value: f32) {
        let (bases, value) = (f64::from(bases), f64::from(value));
        self.covered_bases += bases;
        self.sum += bases * value;
        self.sum_squares += bases * value * value;
        self.min = self.min.min(value);
        self.max = self.max.max(value);
    }

    /// Adds `values`, each covering `span` bases. This is the per-value path for dense data, so
    /// it sums into locals and scales by `span` once.
    fn add_values(&mut self, values: &[f32], span: u32) {
        let (mut sum, mut sum_squares) = (0.0, 0.0);
        let (mut min, mut max) = (self.min, self.max);
        for &value in values {
            let value = f64::from(value);
            sum += value;
            sum_squares += value * value;
            min = min.min(value);
            max = max.max(value);
        }
        let span = f64::from(span);
        self.covered_bases += span * length_as_f64(values.len() as u64);
        self.sum += span * sum;
        self.sum_squares += span * sum_squares;
        self.min = min;
        self.max = max;
    }

    /// Adds `weight` × the counts and sums of `summary`, and its `min` and `max` unweighted.
    /// A weight of 1 merges a whole summary exactly.
    fn add_weighted(&mut self, summary: &Summary, weight: f64) {
        self.covered_bases += weight * summary.covered_bases;
        self.sum += weight * summary.sum;
        self.sum_squares += weight * summary.sum_squares;
        self.min = self.min.min(summary.min);
        self.max = self.max.max(summary.max);
    }

    /// The summary so far, or `None` if nothing is covered.
    fn to_summary(self) -> Option<Summary> {
        (self.covered_bases > 0.0).then_some(Summary {
            covered_bases: self.covered_bases,
            min: self.min,
            max: self.max,
            sum: self.sum,
            sum_squares: self.sum_squares,
        })
    }
}

/// The `Auto` ladder: with `limit = min(longest chromosome, MAX_REDUCTION)`,
/// `r₀ = max(min_first, min(⌈10·p⌉, limit))`, then ×4 up to [`AUTO_MAX_LEVELS`] levels while a
/// level is at most `limit`. With no data, the single level `min_first`.
///
/// Capping the pitch term at `limit` means `r₀` exceeds it only when `min_first` does, and
/// validation keeps `min_first` within [`MAX_REDUCTION`]. The arithmetic is in u64, so the ×4
/// steps cannot wrap (bigtools #79).
fn auto_ladder(min_first: u32, sample: &PitchSampler, longest_chrom_length: u32) -> Vec<u32> {
    let Some(ten_pitches) = sample.ten_pitches() else { return vec![min_first] };
    let limit = u64::from(longest_chrom_length.min(MAX_REDUCTION));
    let first = u64::from(min_first).max(ten_pitches.min(limit));
    let mut ladder = vec![first];
    while ladder.len() < AUTO_MAX_LEVELS {
        let next = ladder[ladder.len() - 1] * AUTO_LEVEL_FACTOR;
        if next > limit {
            break;
        }
        ladder.push(next);
    }
    // Every level is at most max(min_first, limit), both u32.
    ladder.into_iter().map(|reduction| u32::try_from(reduction).unwrap_or(u32::MAX)).collect()
}

/// Number of leading levels to keep at `finish()`, given each level's record count in increasing
/// reduction order.
///
/// Levels are kept while each has strictly fewer records than the one before; the first level
/// that does not shrink is dropped with every level after it, as kent and bigtools stop there.
/// The first level is always kept. For nested ladders (each reduction a multiple of the previous)
/// counts can only fall, so this drops exactly the non-shrinking tail. A non-nested `Explicit`
/// ladder can shrink again after a non-shrinking level; those later levels are dropped too.
pub(crate) fn kept_level_count(record_counts: &[u64]) -> usize {
    if record_counts.is_empty() {
        return 0;
    }
    1 + record_counts.windows(2).take_while(|pair| pair[1] < pair[0]).count()
}

/// Index of the level a zoom-based summary of bins `bin_width` bases wide should use, given the
/// file's reductions in header order (not assumed sorted), or `None` for raw data.
///
/// This is kent's rule: with `desired = ⌊bin_width / 2⌋`, the largest reduction strictly less
/// than `desired`, and `None` if there is none or `desired ≤ 1`. libBigWig uses `≤`, so the two
/// differ when a reduction is exactly half the bin width. Among equal reductions the first wins.
pub(crate) fn best_zoom(
    reductions: impl IntoIterator<Item = u32>,
    bin_width: u32,
) -> Option<usize> {
    let desired = bin_width / 2;
    if desired <= 1 {
        return None;
    }
    let mut best: Option<(usize, u32)> = None;
    for (index, reduction) in reductions.into_iter().enumerate() {
        if reduction < desired && best.is_none_or(|(_, best_reduction)| reduction > best_reduction)
        {
            best = Some((index, reduction));
        }
    }
    best.map(|(index, _)| index)
}

/// Bounds of bin `index` of `count` equal bins over `[start, end)`:
/// `[start + ⌊L·index/count⌋, start + ⌊L·(index+1)/count⌋)` with `L = end − start`, as kent and
/// libBigWig compute them. When `count > L` some bins are empty (`start == end`).
///
/// Requires `start ≤ end` and `index < count`. The products are at most `(2³² − 1)²`, so they fit
/// in u64 and the result is exact.
#[expect(
    clippy::cast_possible_truncation,
    reason = "each offset is at most L, so start + offset is at most end"
)]
pub(crate) fn bin_bounds(start: u32, end: u32, count: u32, index: u32) -> (u32, u32) {
    debug_assert!(start <= end && index < count);
    let length = u64::from(end - start);
    let count = u64::from(count);
    let offset = |index: u64| (length * index / count) as u32;
    (start + offset(u64::from(index)), start + offset(u64::from(index) + 1))
}

/// Builds the [`Error::InvalidOption`] for an invalid zoom policy.
fn invalid_zoom(reason: &str) -> Error {
    Error::InvalidOption(format!("zoom: {reason}"))
}

/// Narrows a position the writer has checked lies within a chromosome, so within u32.
fn position_u32(position: u64) -> u32 {
    u32::try_from(position).expect("the writer checks that items end within the chromosome")
}

/// Converts a region length to f64 for division.
#[allow(clippy::cast_precision_loss, reason = "region lengths are far below 2^53, so exact")]
fn length_as_f64(length: u64) -> f64 {
    length as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- Test data -----------------------------------------------------------------------

    /// Deterministic pseudo-random numbers for generated test data.
    struct Lcg(u64);

    impl Lcg {
        /// A number in `0..bound`.
        fn below(&mut self, bound: u32) -> u32 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            u32::try_from((self.0 >> 33) % u64::from(bound)).unwrap()
        }

        /// A small integer value in `-4..=4`, so every f64 sum in the tests is exact.
        fn value(&mut self) -> f32 {
            f32::from(u8::try_from(self.below(9)).unwrap()) - 4.0
        }
    }

    /// Sorted, non-overlapping `(start, end, value)` items with gaps of 0 to `max_gap` bases and
    /// lengths of 1 to `max_length`, filling at most `chrom_length` bases.
    fn random_items(
        rng: &mut Lcg,
        chrom_length: u32,
        max_gap: u32,
        max_length: u32,
    ) -> Vec<(u32, u32, f32)> {
        let mut items = Vec::new();
        let mut position = rng.below(max_gap + 1);
        while position < chrom_length {
            let end = (position + 1 + rng.below(max_length)).min(chrom_length);
            items.push((position, end, rng.value()));
            position = end + rng.below(max_gap + 1);
        }
        items
    }

    /// The value of each base of a chromosome, `None` where no item covers it.
    fn dense_model(chrom_length: u32, items: &[(u32, u32, f32)]) -> Vec<Option<f32>> {
        let mut bases = vec![None; usize::try_from(chrom_length).unwrap()];
        for &(start, end, value) in items {
            for base in start..end {
                bases[usize::try_from(base).unwrap()] = Some(value);
            }
        }
        bases
    }

    /// The records a level with `reduction` must hold for a chromosome, computed base by base
    /// from its dense model, independently of the builder.
    fn expected_records(chrom_id: u32, bases: &[Option<f32>], reduction: u32) -> Vec<ZoomRecord> {
        let mut records = Vec::new();
        let bin_length = usize::try_from(reduction).unwrap();
        for (bin_index, bin) in bases.chunks(bin_length).enumerate() {
            let bin_start = u32::try_from(bin_index).unwrap() * reduction;
            let covered: Vec<(u32, f64)> = bin
                .iter()
                .enumerate()
                .filter_map(|(offset, value)| {
                    value
                        .map(|value| (bin_start + u32::try_from(offset).unwrap(), f64::from(value)))
                })
                .collect();
            let (Some(&(first, _)), Some(&(last, _))) = (covered.first(), covered.last()) else {
                continue;
            };
            let mut summary = Summary {
                covered_bases: 0.0,
                min: f64::INFINITY,
                max: f64::NEG_INFINITY,
                sum: 0.0,
                sum_squares: 0.0,
            };
            for &(_, value) in &covered {
                summary.covered_bases += 1.0;
                summary.min = summary.min.min(value);
                summary.max = summary.max.max(value);
                summary.sum += value;
                summary.sum_squares += value * value;
            }
            records.push(ZoomRecord { chrom_id, start: first, end: last + 1, summary });
        }
        records
    }

    /// Pushes `(chrom_id, start, end, value)` items through a builder for `reductions`, closes
    /// it, and returns each level's records.
    fn records_from_items(
        reductions: &[u32],
        items: &[(u32, u32, u32, f32)],
    ) -> Vec<Vec<ZoomRecord>> {
        let mut builder = ZoomBuilder::new(reductions);
        let mut levels = vec![Vec::new(); reductions.len()];
        let mut emit = |level: usize, record: ZoomRecord| -> Result<()> {
            levels[level].push(record);
            Ok(())
        };
        for &(chrom_id, start, end, value) in items {
            builder.push_item(chrom_id, start, end, value, &mut emit).unwrap();
        }
        builder.close(&mut emit).unwrap();
        levels
    }

    /// Pushes one fixedStep run on chromosome 0 through a builder for `reductions`, closes it,
    /// and returns each level's records.
    fn records_from_fixed_step(
        reductions: &[u32],
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Vec<Vec<ZoomRecord>> {
        let mut builder = ZoomBuilder::new(reductions);
        let mut levels = vec![Vec::new(); reductions.len()];
        let mut emit = |level: usize, record: ZoomRecord| -> Result<()> {
            levels[level].push(record);
            Ok(())
        };
        builder.push_fixed_step(0, start, step, span, values, &mut emit).unwrap();
        builder.close(&mut emit).unwrap();
        levels
    }

    /// Items `(start, end, value)` all on chromosome 0.
    fn on_chrom_0(items: &[(u32, u32, f32)]) -> Vec<(u32, u32, u32, f32)> {
        items.iter().map(|&(start, end, value)| (0, start, end, value)).collect()
    }

    /// A record on chromosome 0 summarising `covered` bases all of value `value`.
    fn uniform_record(start: u32, end: u32, covered: u32, value: f64) -> ZoomRecord {
        let covered_bases = f64::from(covered);
        ZoomRecord {
            chrom_id: 0,
            start,
            end,
            summary: Summary {
                covered_bases,
                min: value,
                max: value,
                sum: covered_bases * value,
                sum_squares: covered_bases * value * value,
            },
        }
    }

    /// The `(start, end)` extents of a level's records.
    fn extents(records: &[ZoomRecord]) -> Vec<(u32, u32)> {
        records.iter().map(|record| (record.start, record.end)).collect()
    }

    /// A sampler that has seen `count` per-base values on one chromosome.
    fn per_base_sample(count: usize) -> PitchSampler {
        let mut sampler = PitchSampler::new();
        sampler.offer_fixed_step(0, 0, 1, 1, count);
        sampler
    }

    // ---- Summary ---------------------------------------------------------------------------

    /// Summary of the values 1, 2, 3, 4 (one base each).
    fn one_to_four() -> Summary {
        Summary { covered_bases: 4.0, min: 1.0, max: 4.0, sum: 10.0, sum_squares: 30.0 }
    }

    #[test]
    fn mean_is_sum_over_covered_bases() {
        assert_eq!(one_to_four().mean(), 2.5);
    }

    #[test]
    fn std_is_the_sample_standard_deviation() {
        // Sample variance of 1..=4 is (30 - 100/4) / 3 = 5/3.
        assert!((one_to_four().std() - (5.0f64 / 3.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn std_is_zero_for_a_single_covered_base() {
        let single =
            Summary { covered_bases: 1.0, min: 7.0, max: 7.0, sum: 7.0, sum_squares: 49.0 };
        assert_eq!(single.std(), 0.0);
    }

    #[test]
    fn std_is_zero_for_fractional_coverage_below_one() {
        let partial =
            Summary { covered_bases: 0.5, min: 2.0, max: 2.0, sum: 1.0, sum_squares: 2.0 };
        assert_eq!(partial.std(), 0.0);
    }

    #[test]
    fn std_clamps_a_rounding_negative_variance_to_zero() {
        // Three bases of value 1 whose sum of squares came out slightly below the exact 3.
        let constant =
            Summary { covered_bases: 3.0, min: 1.0, max: 1.0, sum: 3.0, sum_squares: 3.0 - 1e-12 };
        assert_eq!(constant.std(), 0.0);
    }

    #[test]
    fn coverage_is_covered_bases_over_region_length() {
        assert_eq!(one_to_four().coverage(8), 0.5);
    }

    #[test]
    fn mean0_counts_uncovered_bases_as_zero() {
        assert_eq!(one_to_four().mean0(20), 0.5);
    }

    // ---- ZoomPolicy validation -------------------------------------------------------------

    #[test]
    fn default_policy_is_auto_with_a_256_bp_minimum_first_reduction() {
        assert_eq!(ZoomPolicy::default(), ZoomPolicy::Auto { min_first_reduction: 256 });
    }

    #[test]
    fn auto_policy_rejects_a_zero_min_first_reduction() {
        let policy = ZoomPolicy::Auto { min_first_reduction: 0 };
        assert!(matches!(policy.validate(), Err(Error::InvalidOption(_))));
    }

    #[test]
    fn auto_policy_accepts_a_min_first_reduction_of_one() {
        assert!(ZoomPolicy::Auto { min_first_reduction: 1 }.validate().is_ok());
    }

    #[test]
    fn auto_policy_rejects_a_min_first_reduction_above_i32_max() {
        let policy = ZoomPolicy::Auto { min_first_reduction: MAX_REDUCTION + 1 };
        assert!(matches!(policy.validate(), Err(Error::InvalidOption(_))));
        assert!(ZoomPolicy::Auto { min_first_reduction: MAX_REDUCTION }.validate().is_ok());
    }

    #[test]
    fn explicit_policy_rejects_a_reduction_above_i32_max() {
        let policy = ZoomPolicy::Explicit(vec![256, MAX_REDUCTION + 1]);
        assert!(matches!(policy.validate(), Err(Error::InvalidOption(_))));
        assert!(ZoomPolicy::Explicit(vec![256, MAX_REDUCTION]).validate().is_ok());
    }

    #[test]
    fn explicit_policy_rejects_an_empty_list() {
        assert!(matches!(ZoomPolicy::Explicit(vec![]).validate(), Err(Error::InvalidOption(_))));
    }

    #[test]
    fn explicit_policy_rejects_a_zero_reduction() {
        let policy = ZoomPolicy::Explicit(vec![0, 10]);
        assert!(matches!(policy.validate(), Err(Error::InvalidOption(_))));
    }

    #[test]
    fn explicit_policy_rejects_a_repeated_reduction() {
        let policy = ZoomPolicy::Explicit(vec![10, 40, 40]);
        assert!(matches!(policy.validate(), Err(Error::InvalidOption(_))));
    }

    #[test]
    fn explicit_policy_rejects_decreasing_reductions() {
        let policy = ZoomPolicy::Explicit(vec![40, 10]);
        assert!(matches!(policy.validate(), Err(Error::InvalidOption(_))));
    }

    #[test]
    fn explicit_policy_accepts_increasing_reductions_that_are_not_multiples() {
        assert!(ZoomPolicy::Explicit(vec![30, 65, 130]).validate().is_ok());
    }

    #[test]
    fn none_policy_is_valid() {
        assert!(ZoomPolicy::None.validate().is_ok());
    }

    #[test]
    fn validation_errors_name_the_zoom_option() {
        let error = ZoomPolicy::Explicit(vec![40, 10]).validate().unwrap_err();
        assert!(error.to_string().contains("zoom"), "{error}");
        assert!(error.to_string().contains("10 follows 40"), "{error}");
    }

    #[test]
    fn only_the_auto_policy_samples_the_pitch() {
        assert!(ZoomPolicy::default().samples_pitch());
        assert!(!ZoomPolicy::Explicit(vec![10]).samples_pitch());
        assert!(!ZoomPolicy::None.samples_pitch());
    }

    // ---- PitchSampler ----------------------------------------------------------------------

    #[test]
    fn sampler_counts_each_fixed_step_value_as_an_item() {
        let mut sampler = PitchSampler::new();
        assert_eq!(sampler.offer_fixed_step(0, 100, 1, 1, 1000), 1000);
        assert_eq!(sampler.ten_pitches(), Some(10));
    }

    #[test]
    fn sampler_takes_only_the_leading_items_that_fit_the_sample() {
        let mut sampler = PitchSampler::new();
        assert_eq!(sampler.offer_fixed_step(0, 0, 1, 1, 60_000), 60_000);
        assert_eq!(sampler.offer_fixed_step(0, 60_000, 1, 1, 10_000), 5_536);
        assert!(sampler.is_full());
        assert_eq!(sampler.offer_fixed_step(0, 70_000, 1, 1, 10), 0);
        assert!(!sampler.offer_item(0, 80_000, 80_001));
    }

    #[test]
    fn sampler_is_not_full_below_the_sample_size() {
        assert!(!per_base_sample(PITCH_SAMPLE_ITEMS - 1).is_full());
        assert!(per_base_sample(PITCH_SAMPLE_ITEMS).is_full());
    }

    #[test]
    fn sampler_extent_includes_gaps_within_a_chromosome() {
        let mut sampler = PitchSampler::new();
        assert!(sampler.offer_item(0, 0, 10));
        assert!(sampler.offer_item(0, 90, 100));
        // Extent 100 over 2 items: p = 50.
        assert_eq!(sampler.ten_pitches(), Some(500));
    }

    #[test]
    fn sampler_extent_excludes_the_space_between_chromosomes() {
        let mut sampler = PitchSampler::new();
        sampler.offer_item(0, 1000, 1010);
        sampler.offer_item(1, 0, 10);
        sampler.offer_item(1, 10, 20);
        // Extents 10 + 20 over 3 items: p = 10.
        assert_eq!(sampler.ten_pitches(), Some(100));
    }

    #[test]
    fn sampler_pitch_rounds_up() {
        let mut sampler = PitchSampler::new();
        sampler.offer_fixed_step(0, 0, 3, 1, 3);
        // Extent 7 over 3 items: 10·p = 23.3, rounded up.
        assert_eq!(sampler.ten_pitches(), Some(24));
    }

    #[test]
    fn sampler_extent_of_a_fixed_step_run_ends_at_its_last_span() {
        let mut sampler = PitchSampler::new();
        sampler.offer_fixed_step(0, 0, 10, 1, 2);
        // Items [0, 1) and [10, 11): extent 11 over 2 items.
        assert_eq!(sampler.ten_pitches(), Some(55));
    }

    // ---- Auto ladder -----------------------------------------------------------------------

    #[test]
    fn auto_ladder_for_per_base_data_starts_at_256_and_grows_by_four() {
        let reductions = ZoomPolicy::default().reductions(&per_base_sample(1000), 248_956_422);
        let expected: Vec<u32> = (0..10).map(|k| 256 * 4u32.pow(k)).collect();
        assert_eq!(reductions, expected);
    }

    #[test]
    fn auto_ladder_has_at_most_ten_levels() {
        let policy = ZoomPolicy::Auto { min_first_reduction: 1 };
        let reductions = policy.reductions(&per_base_sample(1000), u32::MAX);
        assert_eq!(reductions.len(), 10);
        assert_eq!(reductions[0], 10);
        assert_eq!(reductions[9], 10 * 4u32.pow(9));
    }

    #[test]
    fn auto_ladder_honours_a_lower_min_first_reduction() {
        let policy = ZoomPolicy::Auto { min_first_reduction: 32 };
        // Per-base data: 10·p = 10 < 32.
        let reductions = policy.reductions(&per_base_sample(1000), 100_000_000);
        assert_eq!(&reductions[..3], &[32, 128, 512]);
    }

    #[test]
    fn auto_ladder_for_1_kb_bins_starts_at_10_kb() {
        let mut sampler = PitchSampler::new();
        sampler.offer_fixed_step(0, 0, 1000, 1000, 5000);
        let reductions = ZoomPolicy::default().reductions(&sampler, 248_956_422);
        assert_eq!(&reductions[..3], &[10_000, 40_000, 160_000]);
    }

    #[test]
    fn auto_ladder_for_sparse_sites_starts_at_ten_times_their_spacing() {
        // Single-base sites every 100 bp: extent 99,901 over 1000 items, 10·p = 999.01.
        let mut sampler = PitchSampler::new();
        for site in 0..1000 {
            sampler.offer_item(0, site * 100, site * 100 + 1);
        }
        let reductions = ZoomPolicy::default().reductions(&sampler, 248_956_422);
        assert_eq!(&reductions[..2], &[1000, 4000]);
    }

    #[test]
    fn auto_ladder_stops_before_levels_longer_than_the_longest_chromosome() {
        let reductions = ZoomPolicy::default().reductions(&per_base_sample(1000), 100_000);
        assert_eq!(reductions, [256, 1024, 4096, 16_384, 65_536]);
    }

    #[test]
    fn auto_ladder_keeps_a_level_equal_to_the_longest_chromosome() {
        let reductions = ZoomPolicy::default().reductions(&per_base_sample(1000), 4096);
        assert_eq!(reductions, [256, 1024, 4096]);
    }

    #[test]
    fn auto_ladder_for_a_file_with_no_data_is_one_level_at_min_first_reduction() {
        let policy = ZoomPolicy::Auto { min_first_reduction: 100 };
        assert_eq!(policy.reductions(&PitchSampler::new(), 248_956_422), [100]);
    }

    /// libBigWig #26: a chromosome holding one whole-chromosome interval crashed the writer, and
    /// the interim fix of writing no zoom levels broke IGV.
    #[test]
    fn auto_ladder_caps_the_pitch_term_at_the_longest_chromosome() {
        // One item covering a whole 5,000-bp chromosome: 10·p = 50,000.
        let mut sampler = PitchSampler::new();
        sampler.offer_item(0, 0, 5000);
        assert_eq!(ZoomPolicy::default().reductions(&sampler, 5000), [5000]);
    }

    #[test]
    fn auto_ladder_keeps_a_min_first_reduction_longer_than_every_chromosome() {
        assert_eq!(ZoomPolicy::default().reductions(&per_base_sample(50), 100), [256]);
    }

    #[test]
    fn auto_ladder_caps_the_pitch_term_at_i32_max_on_longer_chromosomes() {
        // Two items covering a 4 Gb chromosome: 10·p is far above 2³¹ − 1.
        let mut sampler = PitchSampler::new();
        sampler.offer_item(0, 0, 2_000_000_000);
        sampler.offer_item(0, 2_000_000_000, 4_000_000_000);
        assert_eq!(ZoomPolicy::default().reductions(&sampler, 4_000_000_000), [MAX_REDUCTION]);
    }

    #[test]
    fn auto_ladder_stops_before_a_level_above_i32_max_on_longer_chromosomes() {
        // The next level, 2³¹, would fit the chromosome but not IGV's signed int.
        let policy = ZoomPolicy::Auto { min_first_reduction: 1 << 29 };
        assert_eq!(policy.reductions(&per_base_sample(10), u32::MAX), [1 << 29]);
    }

    #[test]
    fn explicit_reductions_are_capped_at_i32_max_on_longer_chromosomes() {
        let policy = ZoomPolicy::Explicit(vec![1 << 20, 1 << 30, 1 << 31]);
        assert_eq!(policy.reductions(&PitchSampler::new(), u32::MAX), [1 << 20, 1 << 30]);
    }

    #[test]
    fn auto_ladder_does_not_wrap_near_the_u32_limit() {
        let policy = ZoomPolicy::Auto { min_first_reduction: 1 << 30 };
        assert_eq!(policy.reductions(&per_base_sample(10), u32::MAX), [1 << 30]);
    }

    // ---- Explicit and None reductions ------------------------------------------------------

    #[test]
    fn explicit_reductions_are_used_as_given() {
        let policy = ZoomPolicy::Explicit(vec![30, 65, 130]);
        assert_eq!(policy.reductions(&PitchSampler::new(), 1_000_000), [30, 65, 130]);
    }

    #[test]
    fn explicit_reductions_drop_levels_longer_than_the_longest_chromosome() {
        let policy = ZoomPolicy::Explicit(vec![10, 100, 1000, 10_000]);
        assert_eq!(policy.reductions(&PitchSampler::new(), 1000), [10, 100, 1000]);
    }

    #[test]
    fn explicit_reductions_keep_the_smallest_when_every_level_is_too_long() {
        let policy = ZoomPolicy::Explicit(vec![500, 2000]);
        assert_eq!(policy.reductions(&PitchSampler::new(), 100), [500]);
    }

    #[test]
    fn none_policy_has_no_levels() {
        assert!(ZoomPolicy::None.reductions(&per_base_sample(10), 1000).is_empty());
    }

    // ---- Shapes that broke other writers ---------------------------------------------------

    /// Runs `(chrom_id, start, end, value)` items through the `Auto` path the writer follows:
    /// sample them, choose the default ladder, then accumulate. Returns the reductions and each
    /// level's records.
    fn auto_zoom(
        items: &[(u32, u32, u32, f32)],
        longest_chrom_length: u32,
    ) -> (Vec<u32>, Vec<Vec<ZoomRecord>>) {
        let mut sampler = PitchSampler::new();
        for &(chrom_id, start, end, _) in items {
            sampler.offer_item(chrom_id, start, end);
        }
        let reductions = ZoomPolicy::default().reductions(&sampler, longest_chrom_length);
        let levels = records_from_items(&reductions, items);
        (reductions, levels)
    }

    /// Asserts what readers need of a ladder: at least one level, strictly increasing, and each
    /// at most the longest chromosome and 2³¹ − 1.
    fn assert_readable_ladder(reductions: &[u32], longest_chrom_length: u32) {
        assert!(!reductions.is_empty());
        assert!(reductions.windows(2).all(|pair| pair[0] < pair[1]), "{reductions:?}");
        let limit = longest_chrom_length.min(MAX_REDUCTION);
        assert!(reductions.iter().all(|&reduction| reduction <= limit), "{reductions:?}");
    }

    /// Asserts that every level's records lie within their chromosomes and together cover
    /// exactly `covered_bases`.
    fn assert_levels_cover(levels: &[Vec<ZoomRecord>], chrom_lengths: &[u32], covered_bases: f64) {
        for level in levels {
            for record in level {
                let chrom_length = chrom_lengths[usize::try_from(record.chrom_id).unwrap()];
                assert!(record.start < record.end && record.end <= chrom_length, "{record:?}");
            }
            let total: f64 = level.iter().map(|record| record.summary.covered_bases).sum();
            assert_eq!(total, covered_bases);
        }
    }

    /// libBigWig #18: when every interval spanned tens of megabases, the zoom bin size grew past
    /// the longest chromosome and overflowed a u32, and `bwFinalize` looped forever.
    #[test]
    fn zoom_ladder_never_overflows_when_every_interval_spans_most_of_a_chromosome() {
        let chrom_lengths = [248_956_422, 4_000_000_000];
        let items = [
            (0, 1_000_000, 200_000_000, 1.0),
            (0, 200_000_000, 248_000_000, 2.0),
            (1, 100, 3_900_000_000, 3.0),
        ];
        let (reductions, levels) = auto_zoom(&items, chrom_lengths[1]);
        assert_readable_ladder(&reductions, chrom_lengths[1]);
        assert_levels_cover(
            &levels,
            &chrom_lengths,
            199_000_000.0 + 48_000_000.0 + 3_899_999_900.0,
        );
    }

    /// bigtools #79: a bedGraph starting with a ~123 Mb interval followed by tiny ones made the
    /// zoom sizes wrap past `u32::MAX`, and `bedgraphtobigwig` hung.
    #[test]
    fn zoom_ladder_does_not_overflow_when_a_huge_interval_precedes_dense_data() {
        let chrom_length = 248_956_422;
        let mut rng = Lcg(9);
        let mut items = vec![(0, 0, 123_554_979, 0.0)];
        let mut position = 123_554_979;
        for _ in 0..20_000 {
            let end = position + 2 + rng.below(79);
            items.push((0, position, end, 1e-8));
            position = end;
        }
        let (reductions, levels) = auto_zoom(&items, chrom_length);
        assert_readable_ladder(&reductions, chrom_length);
        assert_levels_cover(&levels, &[chrom_length], f64::from(position));
    }

    /// bigtools #85: sparse single-base items got zoom reductions above 2³¹ − 1, which IGV read
    /// as negative, so it showed no data (IGV #1701).
    #[test]
    fn sparse_single_base_items_get_zoom_reductions_below_i32_max() {
        let chrom_length = 248_956_422;
        let items: Vec<(u32, u32, u32, f32)> =
            (0..10).map(|site| (0, site * 24_000_000, site * 24_000_000 + 1, 1.0)).collect();
        let (reductions, levels) = auto_zoom(&items, chrom_length);
        assert_readable_ladder(&reductions, chrom_length);
        assert_levels_cover(&levels, &[chrom_length], 10.0);
    }

    // ---- Dropping non-shrinking levels -----------------------------------------------------

    #[test]
    fn kept_levels_include_every_strictly_shrinking_level() {
        assert_eq!(kept_level_count(&[1000, 250, 60, 15]), 4);
    }

    #[test]
    fn kept_levels_stop_at_the_first_level_that_does_not_shrink() {
        assert_eq!(kept_level_count(&[1000, 250, 250, 60]), 2);
    }

    #[test]
    fn kept_levels_drop_a_level_with_more_records_than_the_previous() {
        assert_eq!(kept_level_count(&[10, 12]), 1);
    }

    #[test]
    fn kept_levels_always_include_the_first_level() {
        assert_eq!(kept_level_count(&[0, 0, 0]), 1);
        assert_eq!(kept_level_count(&[7]), 1);
    }

    #[test]
    fn kept_levels_of_no_levels_is_zero() {
        assert_eq!(kept_level_count(&[]), 0);
    }

    // ---- ZoomBuilder: record semantics -----------------------------------------------------

    #[test]
    fn contiguous_data_gives_one_full_record_per_bin() {
        let levels = records_from_items(&[10], &on_chrom_0(&[(0, 30, 2.0)]));
        assert_eq!(
            levels[0],
            [
                uniform_record(0, 10, 10, 2.0),
                uniform_record(10, 20, 10, 2.0),
                uniform_record(20, 30, 10, 2.0)
            ]
        );
    }

    #[test]
    fn record_fields_aggregate_the_covered_bases_exactly() {
        let levels =
            records_from_items(&[100], &on_chrom_0(&[(0, 2, 1.0), (2, 3, 4.0), (3, 5, -2.0)]));
        let summary = Summary {
            covered_bases: 5.0,
            min: -2.0,
            max: 4.0,
            sum: 2.0 + 4.0 - 4.0,
            sum_squares: 2.0 + 16.0 + 8.0,
        };
        assert_eq!(levels[0], [ZoomRecord { chrom_id: 0, start: 0, end: 5, summary }]);
    }

    #[test]
    fn record_min_and_max_are_updated_independently() {
        // libBigWig's `else if` never lets the first value reach the maximum.
        let decreasing =
            records_from_items(&[10], &on_chrom_0(&[(0, 1, 5.0), (1, 2, 3.0), (2, 3, 1.0)]));
        assert_eq!((decreasing[0][0].summary.min, decreasing[0][0].summary.max), (1.0, 5.0));
        let negative = records_from_items(&[10], &on_chrom_0(&[(0, 1, -1.0), (1, 2, -2.0)]));
        assert_eq!((negative[0][0].summary.min, negative[0][0].summary.max), (-2.0, -1.0));
    }

    #[test]
    fn record_accumulation_is_in_f64() {
        // 2^24 + 1 is not representable in f32, so an f32 sum of these would lose the 1.
        let levels = records_from_items(&[10], &on_chrom_0(&[(0, 1, 16_777_216.0), (1, 2, 1.0)]));
        assert_eq!(levels[0][0].summary.sum, 16_777_217.0);
    }

    #[test]
    fn gap_shorter_than_a_bin_stays_in_one_record_with_a_tight_extent() {
        let levels = records_from_items(&[10], &on_chrom_0(&[(1, 3, 1.0), (6, 8, 1.0)]));
        assert_eq!(levels[0], [uniform_record(1, 8, 4, 1.0)]);
    }

    #[test]
    fn gap_equal_to_an_aligned_bin_leaves_that_bin_without_a_record() {
        let levels = records_from_items(&[10], &on_chrom_0(&[(0, 10, 1.0), (20, 30, 1.0)]));
        assert_eq!(extents(&levels[0]), [(0, 10), (20, 30)]);
    }

    #[test]
    fn gap_equal_to_an_unaligned_bin_gives_records_only_where_there_is_data() {
        let levels = records_from_items(&[10], &on_chrom_0(&[(0, 5, 1.0), (15, 25, 1.0)]));
        assert_eq!(
            levels[0],
            [
                uniform_record(0, 5, 5, 1.0),
                uniform_record(15, 20, 5, 1.0),
                uniform_record(20, 25, 5, 1.0)
            ]
        );
    }

    #[test]
    fn gap_longer_than_a_bin_produces_no_phantom_records() {
        // bigtools PR #103: a cursor moved back into the gap emitted a record per empty bin.
        let levels = records_from_items(&[10], &on_chrom_0(&[(0, 5, 1.0), (100, 105, 3.0)]));
        assert_eq!(levels[0], [uniform_record(0, 5, 5, 1.0), uniform_record(100, 105, 5, 3.0)]);
    }

    #[test]
    fn item_starting_exactly_on_a_bin_edge_opens_a_new_record() {
        let levels = records_from_items(&[10], &on_chrom_0(&[(4, 10, 1.0), (10, 12, 2.0)]));
        assert_eq!(levels[0], [uniform_record(4, 10, 6, 1.0), uniform_record(10, 12, 2, 2.0)]);
    }

    #[test]
    fn item_ending_exactly_on_a_bin_edge_touches_only_its_own_bin() {
        let levels = records_from_items(&[256], &on_chrom_0(&[(0, 256, 1.0)]));
        assert_eq!(levels[0], [uniform_record(0, 256, 256, 1.0)]);
    }

    #[test]
    fn item_spanning_several_bins_contributes_its_overlap_to_each() {
        let levels = records_from_items(&[10], &on_chrom_0(&[(5, 37, 2.0)]));
        assert_eq!(
            levels[0],
            [
                uniform_record(5, 10, 5, 2.0),
                uniform_record(10, 20, 10, 2.0),
                uniform_record(20, 30, 10, 2.0),
                uniform_record(30, 37, 7, 2.0),
            ]
        );
    }

    #[test]
    fn last_bin_is_clipped_at_the_chromosome_end() {
        // A 25-bp chromosome fully covered: the last bin [20, 30) holds only 5 bases.
        let levels = records_from_items(&[10, 20], &on_chrom_0(&[(0, 25, 1.0)]));
        assert_eq!(levels[0].last(), Some(&uniform_record(20, 25, 5, 1.0)));
        assert_eq!(levels[1].last(), Some(&uniform_record(20, 25, 5, 1.0)));
    }

    #[test]
    fn records_never_have_zero_covered_bases() {
        let mut rng = Lcg(7);
        let items = on_chrom_0(&random_items(&mut rng, 5000, 50, 20));
        for level in records_from_items(&[3, 12, 50, 200], &items) {
            assert!(level.iter().all(|record| record.summary.covered_bases >= 1.0));
        }
    }

    #[test]
    fn records_on_different_chromosomes_are_never_merged() {
        let levels = records_from_items(&[10, 40], &[(0, 0, 5, 1.0), (1, 3, 8, 2.0)]);
        for level in &levels {
            assert_eq!(level.len(), 2);
            assert_eq!((level[0].chrom_id, level[0].start, level[0].end), (0, 0, 5));
            assert_eq!((level[1].chrom_id, level[1].start, level[1].end), (1, 3, 8));
        }
    }

    // ---- ZoomBuilder: emission -------------------------------------------------------------

    #[test]
    fn a_record_is_emitted_once_data_arrives_beyond_its_bin() {
        let mut builder = ZoomBuilder::new(&[10]);
        let mut emitted = Vec::new();
        let mut emit = |level: usize, record: ZoomRecord| -> Result<()> {
            emitted.push((level, record));
            Ok(())
        };
        builder.push_item(0, 0, 5, 1.0, &mut emit).unwrap();
        builder.push_item(0, 12, 13, 1.0, &mut emit).unwrap();
        assert_eq!(emitted, [(0, uniform_record(0, 5, 5, 1.0))]);
    }

    #[test]
    fn a_chromosome_change_emits_every_open_record_of_the_previous_chromosome() {
        let mut builder = ZoomBuilder::new(&[10, 40, 160]);
        let mut emitted = Vec::new();
        let mut emit = |level: usize, record: ZoomRecord| -> Result<()> {
            emitted.push((level, record.chrom_id));
            Ok(())
        };
        builder.push_item(0, 0, 5, 1.0, &mut emit).unwrap();
        builder.push_item(1, 0, 5, 1.0, &mut emit).unwrap();
        assert_eq!(emitted, [(0, 0), (1, 0), (2, 0)]);
    }

    #[test]
    fn close_emits_every_open_record() {
        let mut builder = ZoomBuilder::new(&[10, 40]);
        let mut emitted = Vec::new();
        let mut collect = |level: usize, record: ZoomRecord| -> Result<()> {
            emitted.push((level, record));
            Ok(())
        };
        builder.push_item(0, 0, 5, 1.0, &mut collect).unwrap();
        let emitted_by_push = emitted.len();
        let mut collect = |level: usize, record: ZoomRecord| -> Result<()> {
            emitted.push((level, record));
            Ok(())
        };
        builder.close(&mut collect).unwrap();
        assert_eq!(emitted_by_push, 0);
        assert_eq!(emitted, [(0, uniform_record(0, 5, 5, 1.0)), (1, uniform_record(0, 5, 5, 1.0))]);
    }

    #[test]
    fn close_without_data_emits_nothing() {
        let mut builder = ZoomBuilder::new(&[10]);
        let mut emitted = 0;
        builder
            .close(&mut |_, _| {
                emitted += 1;
                Ok(())
            })
            .unwrap();
        assert_eq!(emitted, 0);
    }

    #[test]
    fn an_emit_error_is_returned_from_the_push() {
        let mut builder = ZoomBuilder::new(&[10]);
        let mut fail = |_, _| -> Result<()> { Err(Error::WriterFailed) };
        builder.push_item(0, 0, 5, 1.0, &mut fail).unwrap();
        assert!(matches!(builder.push_item(0, 20, 25, 1.0, &mut fail), Err(Error::WriterFailed)));
    }

    // ---- ZoomBuilder: fixedStep runs -------------------------------------------------------

    #[test]
    fn fixed_step_per_base_run_gives_full_bins() {
        let levels = records_from_fixed_step(&[4], 0, 1, 1, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(extents(&levels[0]), [(0, 4), (4, 6)]);
        assert_eq!(levels[0][0].summary.sum, 10.0);
        assert_eq!(levels[0][1].summary.sum_squares, 61.0);
    }

    #[test]
    fn fixed_step_values_with_span_below_step_count_only_their_span() {
        // Items [0, 3), [10, 13), [20, 23) at r = 100.
        let levels = records_from_fixed_step(&[100], 0, 10, 3, &[1.0, 2.0, 3.0]);
        let summary = Summary {
            covered_bases: 9.0,
            min: 1.0,
            max: 3.0,
            sum: 3.0 * 6.0,
            sum_squares: 3.0 * 14.0,
        };
        assert_eq!(levels[0], [ZoomRecord { chrom_id: 0, start: 0, end: 23, summary }]);
    }

    #[test]
    fn fixed_step_value_straddling_a_bin_edge_is_split() {
        // Items [8, 13) and [13, 18) at r = 10: the first is split at 10.
        let levels = records_from_fixed_step(&[10], 8, 5, 5, &[1.0, 2.0]);
        assert_eq!(levels[0][0], uniform_record(8, 10, 2, 1.0));
        assert_eq!(levels[0][1].start, 10);
        assert_eq!(levels[0][1].end, 18);
        assert_eq!(levels[0][1].summary.sum, 3.0 + 10.0);
    }

    #[test]
    fn fixed_step_run_with_step_longer_than_a_bin_skips_empty_bins() {
        let levels = records_from_fixed_step(&[10], 5, 100, 2, &[1.0, 1.0, 1.0]);
        assert_eq!(extents(&levels[0]), [(5, 7), (105, 107), (205, 207)]);
    }

    #[test]
    fn fixed_step_run_matches_the_same_items_pushed_one_by_one() {
        let mut rng = Lcg(11);
        let values: Vec<f32> = (0..3000).map(|_| rng.value()).collect();
        let (start, step, span) = (17, 7, 4);
        let items: Vec<(u32, u32, u32, f32)> = values
            .iter()
            .zip(0u32..)
            .map(|(&value, index)| (0, start + index * step, start + index * step + span, value))
            .collect();
        let reductions = [5, 20, 33, 66];
        assert_eq!(
            records_from_fixed_step(&reductions, start, step, span, &values),
            records_from_items(&reductions, &items)
        );
    }

    #[test]
    fn fixed_step_runs_continue_the_open_record() {
        let mut builder = ZoomBuilder::new(&[10]);
        let mut records = Vec::new();
        let mut emit = |_, record: ZoomRecord| -> Result<()> {
            records.push(record);
            Ok(())
        };
        builder.push_fixed_step(0, 0, 1, 1, &[1.0; 4], &mut emit).unwrap();
        builder.push_fixed_step(0, 4, 1, 1, &[1.0; 8], &mut emit).unwrap();
        builder.close(&mut emit).unwrap();
        assert_eq!(records, [uniform_record(0, 10, 10, 1.0), uniform_record(10, 12, 2, 1.0)]);
    }

    // ---- ZoomBuilder: against the dense model, and the cascade -----------------------------

    #[test]
    fn nested_ladder_matches_the_dense_model_on_gappy_data() {
        let mut rng = Lcg(1);
        let items = random_items(&mut rng, 20_000, 40, 30);
        let bases = dense_model(20_000, &items);
        let reductions = [4, 16, 64, 256];
        let levels = records_from_items(&reductions, &on_chrom_0(&items));
        for (level, &reduction) in levels.iter().zip(&reductions) {
            assert_eq!(*level, expected_records(0, &bases, reduction), "reduction {reduction}");
        }
    }

    #[test]
    fn non_nested_ladder_matches_the_dense_model_on_gappy_data() {
        let mut rng = Lcg(2);
        let items = random_items(&mut rng, 20_000, 100, 80);
        let bases = dense_model(20_000, &items);
        let reductions = [30, 65, 130, 390];
        let levels = records_from_items(&reductions, &on_chrom_0(&items));
        for (level, &reduction) in levels.iter().zip(&reductions) {
            assert_eq!(*level, expected_records(0, &bases, reduction), "reduction {reduction}");
        }
    }

    #[test]
    fn ladder_matches_the_dense_model_across_chromosomes_of_unaligned_lengths() {
        let mut rng = Lcg(3);
        let chrom_lengths = [1_237, 5_001, 999];
        let reductions = [7, 21, 84, 100];
        let mut items = Vec::new();
        let mut expected = vec![Vec::new(); reductions.len()];
        for (chrom_id, &length) in (0u32..).zip(&chrom_lengths) {
            let chrom_items = random_items(&mut rng, length, 25, 60);
            let bases = dense_model(length, &chrom_items);
            for (level, &reduction) in expected.iter_mut().zip(&reductions) {
                level.extend(expected_records(chrom_id, &bases, reduction));
            }
            items.extend(
                chrom_items.iter().map(|&(start, end, value)| (chrom_id, start, end, value)),
            );
        }
        assert_eq!(records_from_items(&reductions, &items), expected);
    }

    #[test]
    fn fixed_step_per_base_data_matches_the_dense_model() {
        let mut rng = Lcg(4);
        let values: Vec<f32> = (0..10_000).map(|_| rng.value()).collect();
        let items: Vec<(u32, u32, f32)> =
            values.iter().zip(3u32..).map(|(&value, base)| (base, base + 1, value)).collect();
        let bases = dense_model(10_003, &items);
        let reductions = [1, 4, 100, 1000];
        let levels = records_from_fixed_step(&reductions, 3, 1, 1, &values);
        for (level, &reduction) in levels.iter().zip(&reductions) {
            assert_eq!(*level, expected_records(0, &bases, reduction), "reduction {reduction}");
        }
    }

    #[test]
    fn cascaded_nested_level_equals_direct_accumulation_field_by_field() {
        let mut rng = Lcg(5);
        let items = on_chrom_0(&random_items(&mut rng, 30_000, 60, 45));
        let cascaded = records_from_items(&[10, 40, 160], &items);
        let direct = records_from_items(&[160], &items);
        assert_eq!(cascaded[2], direct[0]);
    }

    #[test]
    fn cascaded_level_of_a_non_nested_ladder_equals_direct_accumulation_field_by_field() {
        // 65 is accumulated from raw items alongside 30; 130 merges 65's records.
        let mut rng = Lcg(6);
        let items = on_chrom_0(&random_items(&mut rng, 30_000, 60, 45));
        let ladder = records_from_items(&[30, 65, 130], &items);
        assert_eq!(ladder[1], records_from_items(&[65], &items)[0]);
        assert_eq!(ladder[2], records_from_items(&[130], &items)[0]);
    }

    // ---- ZoomRecord encoding ---------------------------------------------------------------

    /// A record with distinct values in every field.
    fn sample_record() -> ZoomRecord {
        ZoomRecord {
            chrom_id: 3,
            start: 1024,
            end: 1280,
            summary: Summary {
                covered_bases: 200.0,
                min: -1.5,
                max: 8.25,
                sum: 412.5,
                sum_squares: 9001.0,
            },
        }
    }

    #[test]
    fn zoom_record_encodes_to_32_bytes_in_field_order() {
        let mut buf = Vec::new();
        sample_record().encode(ByteOrder::Little, &mut buf);
        assert_eq!(buf.len(), ZoomRecord::ENCODED_LEN);
        let mut expected = Vec::new();
        for field in [3u32, 1024, 1280, 200] {
            expected.extend_from_slice(&field.to_le_bytes());
        }
        for field in [-1.5f32, 8.25, 412.5, 9001.0] {
            expected.extend_from_slice(&field.to_le_bytes());
        }
        assert_eq!(buf, expected);
    }

    #[test]
    fn zoom_record_round_trips_in_both_byte_orders() {
        for byte_order in [ByteOrder::Little, ByteOrder::Big] {
            let mut buf = Vec::new();
            sample_record().encode(byte_order, &mut buf);
            let mut reader = ByteReader::new(&buf, byte_order, 0);
            assert_eq!(ZoomRecord::decode(&mut reader).unwrap(), sample_record());
            assert_eq!(reader.remaining(), 0);
        }
    }

    #[test]
    fn zoom_record_byte_orders_differ_on_disk() {
        let (mut little, mut big) = (Vec::new(), Vec::new());
        sample_record().encode(ByteOrder::Little, &mut little);
        sample_record().encode(ByteOrder::Big, &mut big);
        assert_eq!(&big[..4], &3u32.to_be_bytes());
        assert_ne!(little, big);
    }

    #[test]
    fn zoom_record_stores_f32_rounded_values() {
        let mut record = sample_record();
        record.summary.sum = 0.1;
        let mut buf = Vec::new();
        record.encode(ByteOrder::Little, &mut buf);
        let decoded = ZoomRecord::decode(&mut ByteReader::new(&buf, ByteOrder::Little, 0)).unwrap();
        assert_eq!(decoded.summary.sum, f64::from(0.1f32));
    }

    #[test]
    fn truncated_zoom_record_is_corrupt() {
        let mut buf = Vec::new();
        sample_record().encode(ByteOrder::Little, &mut buf);
        buf.truncate(31);
        let result = ZoomRecord::decode(&mut ByteReader::new(&buf, ByteOrder::Little, 500));
        assert!(matches!(result, Err(Error::Corrupt { offset: 528, .. })));
    }

    // ---- Best zoom -------------------------------------------------------------------------

    #[test]
    fn best_zoom_picks_the_largest_reduction_below_half_the_bin() {
        assert_eq!(best_zoom([256, 1024, 4096], 5000), Some(1));
    }

    #[test]
    fn best_zoom_excludes_a_reduction_equal_to_half_the_bin() {
        // libBigWig's `≤` would pick 1024 here; kent's `<` picks 256.
        assert_eq!(best_zoom([256, 1024, 4096], 2048), Some(0));
        assert_eq!(best_zoom([256, 1024, 4096], 2050), Some(1));
    }

    #[test]
    fn best_zoom_is_none_when_every_reduction_is_too_coarse() {
        assert_eq!(best_zoom([256, 1024], 512), None);
    }

    #[test]
    fn best_zoom_is_none_for_bins_narrower_than_four_bases() {
        assert_eq!(best_zoom([1], 3), None);
        assert_eq!(best_zoom([1], 4), Some(0));
    }

    #[test]
    fn best_zoom_does_not_assume_sorted_reductions() {
        assert_eq!(best_zoom([4096, 256, 1024], 5000), Some(2));
    }

    #[test]
    fn best_zoom_takes_the_first_of_equal_reductions() {
        assert_eq!(best_zoom([256, 1024, 1024], 5000), Some(1));
    }

    #[test]
    fn best_zoom_without_levels_is_none() {
        assert_eq!(best_zoom([], 1_000_000), None);
    }

    // ---- Zoom-based summaries --------------------------------------------------------------

    #[test]
    fn zoom_summary_of_a_contained_record_is_the_record() {
        let record = uniform_record(10, 20, 10, 3.0);
        let mut accumulator = ZoomSummaryAccumulator::new(0, 100);
        accumulator.add_record(&record);
        assert_eq!(accumulator.summary(), Some(record.summary));
    }

    #[test]
    fn zoom_summary_pro_rates_a_partial_overlap_by_record_extent() {
        let record = ZoomRecord {
            chrom_id: 0,
            start: 0,
            end: 100,
            summary: Summary {
                covered_bases: 50.0,
                min: 1.0,
                max: 9.0,
                sum: 200.0,
                sum_squares: 1000.0,
            },
        };
        let mut accumulator = ZoomSummaryAccumulator::new(75, 200);
        accumulator.add_record(&record);
        let expected =
            Summary { covered_bases: 12.5, min: 1.0, max: 9.0, sum: 50.0, sum_squares: 250.0 };
        assert_eq!(accumulator.summary(), Some(expected));
    }

    #[test]
    fn zoom_summary_takes_min_and_max_over_overlapping_records_unweighted() {
        let mut accumulator = ZoomSummaryAccumulator::new(5, 15);
        accumulator.add_record(&uniform_record(0, 10, 10, -3.0));
        accumulator.add_record(&uniform_record(10, 20, 10, 7.0));
        let summary = accumulator.summary().unwrap();
        assert_eq!((summary.min, summary.max), (-3.0, 7.0));
        assert_eq!(summary.covered_bases, 10.0);
        assert_eq!(summary.sum, -15.0 + 35.0);
    }

    #[test]
    fn zoom_summary_skips_records_with_no_covered_bases() {
        let mut empty = uniform_record(0, 10, 0, 0.0);
        empty.summary.min = -1000.0;
        let mut accumulator = ZoomSummaryAccumulator::new(0, 10);
        accumulator.add_record(&empty);
        assert_eq!(accumulator.summary(), None);
        accumulator.add_record(&uniform_record(0, 10, 10, 2.0));
        assert_eq!(accumulator.summary().unwrap().min, 2.0);
    }

    #[test]
    fn zoom_summary_ignores_records_that_only_touch_the_region() {
        let mut accumulator = ZoomSummaryAccumulator::new(10, 20);
        accumulator.add_record(&uniform_record(0, 10, 10, 1.0));
        accumulator.add_record(&uniform_record(20, 30, 10, 1.0));
        assert_eq!(accumulator.summary(), None);
    }

    #[test]
    fn zoom_summary_skips_a_record_with_an_empty_extent() {
        let mut accumulator = ZoomSummaryAccumulator::new(0, 100);
        accumulator.add_record(&uniform_record(50, 50, 4, 1.0));
        accumulator.add_record(&uniform_record(60, 40, 4, 1.0));
        assert_eq!(accumulator.summary(), None);
    }

    #[test]
    fn zoom_summary_over_whole_bins_equals_the_exact_summary() {
        let mut rng = Lcg(8);
        let items = random_items(&mut rng, 10_000, 30, 30);
        let records = records_from_items(&[100], &on_chrom_0(&items)).remove(0);
        let (start, end) = (2000, 7000);
        let mut zoom = ZoomSummaryAccumulator::new(start, end);
        for record in &records {
            zoom.add_record(record);
        }
        let mut raw = RawSummaryAccumulator::new(start, end);
        for &(item_start, item_end, value) in &items {
            raw.add_item(item_start, item_end, value);
        }
        assert_eq!(zoom.summary(), raw.summary());
    }

    // ---- Exact raw summaries ---------------------------------------------------------------

    #[test]
    fn raw_summary_clips_items_to_the_region() {
        let mut accumulator = RawSummaryAccumulator::new(10, 20);
        accumulator.add_item(5, 12, 2.0);
        accumulator.add_item(18, 30, 3.0);
        let expected = Summary {
            covered_bases: 4.0,
            min: 2.0,
            max: 3.0,
            sum: 4.0 + 6.0,
            sum_squares: 8.0 + 18.0,
        };
        assert_eq!(accumulator.summary(), Some(expected));
    }

    #[test]
    fn raw_summary_is_none_when_nothing_in_the_region_is_covered() {
        let mut accumulator = RawSummaryAccumulator::new(10, 20);
        assert_eq!(accumulator.summary(), None);
        accumulator.add_item(0, 10, 1.0);
        accumulator.add_item(20, 25, 1.0);
        assert_eq!(accumulator.summary(), None);
    }

    #[test]
    fn raw_summary_min_and_max_are_updated_independently() {
        let mut accumulator = RawSummaryAccumulator::new(0, 10);
        accumulator.add_item(0, 1, 5.0);
        accumulator.add_item(1, 2, 3.0);
        let summary = accumulator.summary().unwrap();
        assert_eq!((summary.min, summary.max), (3.0, 5.0));
    }

    #[test]
    fn raw_summary_of_a_fixed_step_run_equals_its_items() {
        let values = [1.0, -2.0, 3.0, 4.0, 0.5];
        for (region_start, region_end) in [(0, 1000), (12, 33), (40, 41)] {
            let mut run = RawSummaryAccumulator::new(region_start, region_end);
            run.add_fixed_step(10, 7, 3, &values);
            let mut items = RawSummaryAccumulator::new(region_start, region_end);
            for (&value, index) in values.iter().zip(0u32..) {
                items.add_item(10 + index * 7, 13 + index * 7, value);
            }
            assert_eq!(run.summary(), items.summary(), "region {region_start}-{region_end}");
        }
    }

    #[test]
    fn raw_summary_of_an_empty_fixed_step_run_is_none() {
        let mut accumulator = RawSummaryAccumulator::new(0, 10);
        accumulator.add_fixed_step(0, 1, 1, &[]);
        assert_eq!(accumulator.summary(), None);
    }

    #[test]
    fn unbounded_raw_summary_counts_every_base() {
        let mut accumulator = RawSummaryAccumulator::unbounded();
        accumulator.add_item(0, 10, 1.0);
        accumulator.add_item(u32::MAX - 10, u32::MAX, 2.0);
        accumulator.add_fixed_step(100, 1, 1, &[3.0; 5]);
        let summary = accumulator.summary().unwrap();
        assert_eq!(summary.covered_bases, 25.0);
        assert_eq!(summary.sum, 10.0 + 20.0 + 15.0);
    }

    // ---- Bin boundaries --------------------------------------------------------------------

    #[test]
    fn bin_bounds_follow_the_floor_formula() {
        let bins: Vec<(u32, u32)> = (0..3).map(|index| bin_bounds(100, 110, 3, index)).collect();
        assert_eq!(bins, [(100, 103), (103, 106), (106, 110)]);
    }

    #[test]
    fn bins_partition_the_region_exactly() {
        let (start, end, count) = (1_000, 1_997, 64);
        let bins: Vec<(u32, u32)> =
            (0..count).map(|index| bin_bounds(start, end, count, index)).collect();
        assert_eq!(bins[0].0, start);
        assert_eq!(bins[63].1, end);
        assert!(bins.windows(2).all(|pair| pair[0].1 == pair[1].0));
    }

    #[test]
    fn bins_beyond_the_region_length_are_empty() {
        let bins: Vec<(u32, u32)> = (0..4).map(|index| bin_bounds(0, 2, 4, index)).collect();
        assert_eq!(bins, [(0, 0), (0, 1), (1, 1), (1, 2)]);
    }

    #[test]
    fn bin_bounds_do_not_overflow_at_the_u32_limit() {
        assert_eq!(bin_bounds(0, u32::MAX, u32::MAX, u32::MAX - 1), (u32::MAX - 1, u32::MAX));
        assert_eq!(bin_bounds(0, u32::MAX, 2, 1), (u32::MAX / 2, u32::MAX));
    }
}
