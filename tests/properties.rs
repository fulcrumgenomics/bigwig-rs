//! Round-trip property tests: random chromosome lists, data and writer options are written with
//! `BigWigWriter`, read back with `BigWigReader`, and every answer is compared with a brute-force
//! model built from the pushed items.
//!
//! A case is a chromosome list in random order (with prefix collisions such as `chr1`, `chr10`
//! and `chr1_alt`), sorted non-overlapping bedGraph, variableStep and fixedStep runs with gaps on
//! each chromosome, random options, a random chromosome push order and a random push granularity.
//! Data may touch the chromosome end, sit just below `u32::MAX`, or (with `clip`) run past the
//! end. Each property runs 64 cases by default; set `PROPTEST_CASES` for a longer soak. Failures
//! are not persisted to disk: a failing case is printed shrunk, and each bug it exposes gets a
//! deterministic regression test at the end of this file.

#![allow(
    clippy::float_cmp,
    reason = "stored values, counts, min and max round-trip exactly, so they are compared exactly"
)]
#![allow(clippy::cast_precision_loss, reason = "base counts here are far below 2^52")]

use std::io::{Cursor, Seek, Write};

use bigwig::{
    BigWigReader, BigWigWriter, ByteOrder, ChromIdOrder, ChromInfo, Compression, DataItem, Summary,
    SummarySource, WriterOptions, ZoomLevel, ZoomPolicy, ZoomRecord,
};
use proptest::collection::{btree_set, vec};
use proptest::prelude::*;
use proptest::test_runner::TestCaseError;

/// Cases per property when `PROPTEST_CASES` is not set.
const DEFAULT_CASES: u32 = 64;
/// Largest reduction a zoom level may have: IGV reads reductions as signed 32-bit integers.
const MAX_REDUCTION: u32 = i32::MAX.unsigned_abs();
/// Leading items the `Auto` zoom policy measures the data's pitch from.
const PITCH_SAMPLE_ITEMS: usize = 65_536;
/// Most levels an `Auto` ladder has.
const AUTO_MAX_LEVELS: usize = 10;
/// Ratio between successive `Auto` levels.
const AUTO_LEVEL_FACTOR: u64 = 4;
/// Exact (raw) summaries differ from brute force only in the order of f64 additions.
const RAW_TOLERANCE: Tolerance = Tolerance { relative: 1e-9, f32_storage: false };
/// Summaries made from whole zoom records carry each record's sums rounded to f32.
const ZOOM_TOLERANCE: Tolerance = Tolerance { relative: 1e-6, f32_storage: true };
/// Chromosomes with data whose summaries each case checks; 64 bin counts each cost a full read.
const SUMMARY_CHROMS: usize = 3;
/// Longest region whose per-base values are checked against a dense model.
const MAX_DENSE_REGION: u32 = 1 << 20;
/// The chromosome a [`GappedRun`] is written to.
const GAPPED_CHROM: &str = "chr1";

/// A reader over a file written to memory.
type Reader = BigWigReader<Cursor<Vec<u8>>>;
/// Outcome of checking one property on one case.
type Check = Result<(), TestCaseError>;

/// One generated input: a chromosome list with its data, the writer options, the chromosome push
/// order, and seeds for the query regions.
#[derive(Clone, Debug)]
struct Case {
    /// The chromosome list passed to the writer, in list order.
    chroms: Vec<ChromCase>,
    /// List positions in the order chromosomes receive data under [`ChromIdOrder::Arrival`].
    /// Under [`ChromIdOrder::Listed`] data is pushed in list order instead.
    write_order: Vec<usize>,
    options: Options,
    /// Seeds from which each chromosome's query regions are derived.
    query_seeds: Vec<(u32, u32)>,
}

impl Case {
    /// Places each chromosome's planned runs one after another from position 0, stopping a
    /// chromosome's data at `profile`'s extent and the whole case's at its item and base budgets,
    /// ends the chromosome after the data or, with clipping, possibly within it, then shifts the
    /// data by the plan's offset as far as the chromosome and the data can still end by
    /// `u32::MAX`.
    fn lay_out(
        profile: Profile,
        names: Vec<String>,
        plans: Vec<ChromPlan>,
        write_order: Vec<usize>,
        options: Options,
        query_seeds: Vec<(u32, u32)>,
    ) -> Self {
        let mut budget = Budget { items: profile.item_budget(), bases: profile.base_budget() };
        let chroms = names
            .into_iter()
            .zip(plans)
            .map(|(name, plan)| {
                let mut runs = Vec::new();
                let mut cursor = 0;
                for unit in plan.units {
                    let start = cursor + u64::from(unit.gap);
                    let Some(data) = unit.shape.place(start, profile.max_extent(), &mut budget)
                    else {
                        break;
                    };
                    cursor = u64::from(data.end());
                    runs.push(Run { data, pushes: unit.pushes });
                }
                if runs.is_empty() {
                    return ChromCase { name, length: plan.empty_length, runs };
                }
                let unshifted_length = if options.clip && plan.overhang > 0 {
                    to_u32(cursor.saturating_sub(u64::from(plan.overhang)).max(1))
                } else {
                    to_u32(cursor + u64::from(plan.tail))
                };
                let offset =
                    plan.offset.min(u32::MAX - to_u32(cursor.max(u64::from(unshifted_length))));
                for run in &mut runs {
                    run.data.shift(offset);
                }
                ChromCase { name, length: unshifted_length + offset, runs }
            })
            .collect();
        Self { chroms, write_order, options, query_seeds }
    }

    /// List positions of the chromosomes with runs, in the order the runs are pushed.
    fn push_order(&self) -> Vec<usize> {
        let order: Vec<usize> = if self.options.chrom_id_order == ChromIdOrder::Listed {
            (0..self.chroms.len()).collect()
        } else {
            self.write_order.clone()
        };
        order.into_iter().filter(|&position| !self.chroms[position].runs.is_empty()).collect()
    }

    /// Number of data sections the writer must cut: a section ends when it holds
    /// `items_per_section` items, and before an item on another chromosome, of another section
    /// type, step or span, or, for fixedStep, one that doesn't continue the section's run.
    fn expected_section_count(&self) -> u64 {
        // The open section: its chromosome, type (with step and span), item count (0 when none
        // is open), and where a fixedStep run's next item must start to continue it.
        let (mut open_chrom, mut open_kind, mut open_items, mut next_start) = (0, (0, 0, 0), 0, 0);
        let mut sections = 0;
        for position in self.push_order() {
            for run in &self.chroms[position].runs {
                let mut items = Vec::new();
                run.data.expand(0, &mut items);
                for item in clip_items(items, self.chroms[position].length) {
                    // A clipped item's span differs from its run's, so it has a kind of its own.
                    let (kind, fixed_step) = match run.data {
                        RunData::Intervals(_) => ((1, 0, 0), None),
                        RunData::VariableStep { .. } => ((2, 0, item.end - item.start), None),
                        RunData::FixedStep { step, .. } => {
                            ((3, step, item.end - item.start), Some(step))
                        }
                    };
                    let continues = open_items > 0
                        && open_chrom == position
                        && open_kind == kind
                        && (fixed_step.is_none() || next_start == item.start);
                    if continues {
                        open_items += 1;
                    } else {
                        (open_chrom, open_kind, open_items) = (position, kind, 1);
                        sections += 1;
                    }
                    next_start = item.start.wrapping_add(fixed_step.unwrap_or(0));
                    if open_items == self.options.items_per_section {
                        open_items = 0;
                    }
                }
            }
        }
        sections
    }

    /// The writer's chromosome list.
    fn chrom_infos(&self) -> Vec<ChromInfo> {
        self.chroms.iter().map(|chrom| ChromInfo::new(chrom.name.clone(), chrom.length)).collect()
    }

    /// Length of the longest chromosome.
    fn longest_length(&self) -> u32 {
        self.chroms.iter().map(|chrom| chrom.length).max().unwrap_or(0)
    }

    /// Query regions on a chromosome of `length` with `items`: the whole chromosome, empty
    /// regions at both ends, the span of the data, a random region per seed, and regions at the
    /// edges of a random item.
    fn regions(&self, length: u32, items: &[DataItem]) -> Vec<(u32, u32)> {
        let mut regions = vec![(0, length), (0, 0), (length, length)];
        if let (Some(first), Some(last)) = (items.first(), items.last()) {
            regions.push((first.start, last.end));
        }
        for &seed in &self.query_seeds {
            regions.push(seeded_region(length, seed));
            if !items.is_empty() {
                let item = items[seed.0 as usize % items.len()];
                regions.push((item.start, item.end));
                regions.push((item.start.saturating_sub(1), item.start));
                regions.push((item.end - 1, item.end.saturating_add(1).min(length)));
                regions.push((item.end, item.end.saturating_add(seed.1 % 3).min(length)));
            }
        }
        regions
    }
}

/// One chromosome of a case: its name, its length and its data, in position order.
#[derive(Clone, Debug)]
struct ChromCase {
    name: String,
    length: u32,
    runs: Vec<Run>,
}

impl ChromCase {
    /// The chromosome's items as a reader returns them, with fixedStep and variableStep runs
    /// expanded and clipped to the chromosome, tagged with `chrom_id`.
    fn items(&self, chrom_id: u32) -> Vec<DataItem> {
        let mut items = Vec::new();
        for run in &self.runs {
            run.data.expand(chrom_id, &mut items);
        }
        clip_items(items, self.length)
    }
}

/// One run of data and how it is split into pushes.
#[derive(Clone, Debug)]
struct Run {
    data: RunData,
    pushes: Pushes,
}

/// A run of one section type, in absolute coordinates.
#[derive(Clone, Debug)]
enum RunData {
    /// bedGraph intervals `(start, end, value)`.
    Intervals(Vec<(u32, u32, f32)>),
    /// variableStep items `(start, value)`, each covering `span` bases.
    VariableStep { span: u32, items: Vec<(u32, f32)> },
    /// `count` fixedStep values cycling through `pattern`: value `i` covers
    /// `[start + i·step, start + i·step + span)`.
    FixedStep { start: u32, step: u32, span: u32, pattern: Vec<f32>, count: u32 },
}

impl RunData {
    /// Appends the run's items, tagged with `chrom_id`, to `out`.
    fn expand(&self, chrom_id: u32, out: &mut Vec<DataItem>) {
        let item = |start: u32, end: u32, value: f32| DataItem { chrom_id, start, end, value };
        match self {
            Self::Intervals(intervals) => {
                out.extend(intervals.iter().map(|&(start, end, value)| item(start, end, value)));
            }
            Self::VariableStep { span, items } => {
                out.extend(items.iter().map(|&(start, value)| item(start, start + span, value)));
            }
            Self::FixedStep { start, step, span, .. } => {
                for (index, value) in (0u32..).zip(self.fixed_values()) {
                    let item_start = start + index * step;
                    out.push(item(item_start, item_start + span, value));
                }
            }
        }
    }

    /// The values of a fixedStep run; empty for other runs.
    fn fixed_values(&self) -> Vec<f32> {
        match self {
            Self::FixedStep { pattern, count, .. } => {
                pattern.iter().copied().cycle().take(*count as usize).collect()
            }
            Self::Intervals(_) | Self::VariableStep { .. } => Vec::new(),
        }
    }

    /// Moves the run `by` bases towards the chromosome end.
    fn shift(&mut self, by: u32) {
        match self {
            Self::Intervals(intervals) => {
                for (start, end, _) in intervals {
                    *start += by;
                    *end += by;
                }
            }
            Self::VariableStep { items, .. } => {
                for (start, _) in items {
                    *start += by;
                }
            }
            Self::FixedStep { start, .. } => *start += by,
        }
    }

    /// End of the run's last item.
    fn end(&self) -> u32 {
        match self {
            Self::Intervals(intervals) => intervals.last().map_or(0, |&(_, end, _)| end),
            Self::VariableStep { span, items } => {
                items.last().map_or(0, |&(start, _)| start + span)
            }
            Self::FixedStep { start, step, span, count, .. } => start + (count - 1) * step + span,
        }
    }

    /// Pushes the run to `chrom`, `pushes` deciding how it is sliced.
    fn push<W: Write + Seek>(
        &self,
        writer: &mut BigWigWriter<W>,
        chrom: &str,
        pushes: Pushes,
    ) -> bigwig::Result<()> {
        let chunk_len = match pushes {
            Pushes::Whole => usize::MAX,
            Pushes::Chunks(len) => len,
            Pushes::Single => 1,
        };
        match self {
            Self::Intervals(intervals) => {
                for chunk in intervals.chunks(chunk_len) {
                    if let [(start, end, value)] = *chunk {
                        writer.push_interval(chrom, start, end, value)?;
                    } else {
                        writer.push_intervals(chrom, chunk)?;
                    }
                }
            }
            Self::VariableStep { span, items } => {
                for chunk in items.chunks(chunk_len) {
                    writer.push_variable_step(chrom, *span, chunk)?;
                }
            }
            Self::FixedStep { start, step, span, .. } => {
                let values = self.fixed_values();
                // In u64: the start after the last chunk may lie past u32::MAX.
                let mut chunk_start = u64::from(*start);
                for chunk in values.chunks(chunk_len) {
                    writer.push_fixed_step(chrom, to_u32(chunk_start), *step, *span, chunk)?;
                    chunk_start += chunk.len() as u64 * u64::from(*step);
                }
            }
        }
        Ok(())
    }
}

/// How a run is sliced into push calls; the file must not depend on it.
#[derive(Clone, Copy, Debug)]
enum Pushes {
    /// One call with the whole run.
    Whole,
    /// Calls of this many items each (the last may have fewer).
    Chunks(usize),
    /// One call per item (`push_interval` for bedGraph).
    Single,
}

/// The writer options of a case.
#[derive(Clone, Debug)]
struct Options {
    items_per_section: u16,
    zoom: ZoomPolicy,
    compression: Compression,
    byte_order: ByteOrder,
    chrom_id_order: ChromIdOrder,
    zoom_items_per_block: u32,
    index_block_size: u32,
    /// Bytes of compressed zoom blocks held in memory before they spill to a temporary file;
    /// `None` for the default, which these cases never reach.
    zoom_memory_limit: Option<usize>,
    clip: bool,
}

impl Options {
    /// The [`WriterOptions`] these describe.
    fn writer_options(&self) -> WriterOptions {
        let options = WriterOptions::new()
            .items_per_section(self.items_per_section)
            .zoom(self.zoom.clone())
            .compression(self.compression)
            .byte_order(self.byte_order)
            .chrom_id_order(self.chrom_id_order)
            .zoom_items_per_block(self.zoom_items_per_block)
            .index_block_size(self.index_block_size)
            .clip(self.clip);
        match self.zoom_memory_limit {
            Some(limit) => options.zoom_memory_limit(limit),
            None => options,
        }
    }
}

/// Size limits on generated data, so the dense models and the many queries stay cheap.
#[derive(Clone, Copy, Debug)]
enum Profile {
    /// Modest data, for the properties that query it many ways.
    Typical,
    /// Room for per-base runs longer than the `Auto` policy's pitch sample, so the sample fills
    /// partway through the data.
    LongRuns,
}

impl Profile {
    /// Most items in one case.
    fn item_budget(self) -> u64 {
        match self {
            Self::Typical => 8_000,
            Self::LongRuns => 150_000,
        }
    }

    /// Most covered bases in one case, which bounds the records of a reduction-1 zoom level.
    fn base_budget(self) -> u64 {
        match self {
            Self::Typical => 100_000,
            Self::LongRuns => 300_000,
        }
    }

    /// Positions past which no item on a chromosome may end.
    fn max_extent(self) -> u64 {
        match self {
            Self::Typical => 200_000,
            Self::LongRuns => 400_000,
        }
    }
}

/// Items and covered bases a case may still add.
#[derive(Clone, Copy, Debug)]
struct Budget {
    items: u64,
    bases: u64,
}

impl Budget {
    /// Takes one item of `bases` bases if the budget allows.
    fn take(&mut self, bases: u64) -> bool {
        let fits = self.items >= 1 && self.bases >= bases;
        if fits {
            self.items -= 1;
            self.bases -= bases;
        }
        fits
    }
}

/// A chromosome's planned data before it is laid out.
#[derive(Clone, Debug)]
struct ChromPlan {
    units: Vec<UnitPlan>,
    /// Bases after the last item; 0 makes the data touch the chromosome end.
    tail: u32,
    /// With clipping, how far the data runs past the chromosome end (0 for not at all).
    overhang: u32,
    /// How far to shift the data, at most to where the chromosome ends at `u32::MAX`; large
    /// shifts test coordinates near the top of the u32 range.
    offset: u32,
    /// Length of the chromosome if it ends up with no data.
    empty_length: u32,
}

/// One planned run: the gap before it, its shape and how it is pushed.
#[derive(Clone, Debug)]
struct UnitPlan {
    gap: u32,
    shape: UnitShape,
    pushes: Pushes,
}

/// A run's shape in relative coordinates.
#[derive(Clone, Debug)]
enum UnitShape {
    /// bedGraph intervals, each `(gap before it, length, value)`.
    Intervals(Vec<(u32, u32, f32)>),
    /// variableStep items of `span` bases, each `(gap before it, value)`.
    VariableStep { span: u32, items: Vec<(u32, f32)> },
    /// `count` fixedStep values cycling through `pattern`.
    FixedStep { step: u32, span: u32, pattern: Vec<f32>, count: u32 },
}

impl UnitShape {
    /// Places the run at `start`, keeping the leading items that end by `max_extent` and fit the
    /// budget. `None` if no item does.
    fn place(&self, start: u64, max_extent: u64, budget: &mut Budget) -> Option<RunData> {
        let data = match self {
            Self::Intervals(shape) => {
                let mut intervals = Vec::new();
                let mut position = start;
                for &(gap, length, value) in shape {
                    let item_start = position + u64::from(gap);
                    let item_end = item_start + u64::from(length);
                    if item_end > max_extent || !budget.take(u64::from(length)) {
                        break;
                    }
                    intervals.push((to_u32(item_start), to_u32(item_end), value));
                    position = item_end;
                }
                (!intervals.is_empty()).then_some(RunData::Intervals(intervals))?
            }
            Self::VariableStep { span, items: shape } => {
                let mut items = Vec::new();
                let mut position = start;
                for &(gap, value) in shape {
                    let item_start = position + u64::from(gap);
                    let item_end = item_start + u64::from(*span);
                    if item_end > max_extent || !budget.take(u64::from(*span)) {
                        break;
                    }
                    items.push((to_u32(item_start), value));
                    position = item_end;
                }
                (!items.is_empty()).then_some(RunData::VariableStep { span: *span, items })?
            }
            Self::FixedStep { step, span, pattern, count } => {
                let (step64, span64) = (u64::from(*step), u64::from(*span));
                let fitting = if start + span64 > max_extent {
                    0
                } else {
                    (max_extent - start - span64) / step64 + 1
                };
                let mut kept = 0;
                while kept < u64::from(*count).min(fitting) && budget.take(span64) {
                    kept += 1;
                }
                (kept > 0).then_some(RunData::FixedStep {
                    start: to_u32(start),
                    step: *step,
                    span: *span,
                    pattern: pattern.clone(),
                    count: to_u32(kept),
                })?
            }
        };
        Some(data)
    }
}

/// What a file written from a case must contain, derived from the case alone.
#[derive(Debug)]
struct Model {
    /// List positions of the chromosomes, in `chromId` order.
    by_chrom_id: Vec<usize>,
    /// Each chromosome's items, by list position, tagged with its expected `chromId`.
    items: Vec<Vec<DataItem>>,
    /// List positions of the chromosomes that receive data (some items survive clipping), in
    /// push order.
    data_order: Vec<usize>,
}

impl Model {
    /// The model of `case`: `chromId`s by list position under [`ChromIdOrder::Listed`], and
    /// otherwise by first data, then the chromosomes without data in list order.
    fn new(case: &Case) -> Self {
        let count = case.chroms.len();
        let has_data: Vec<bool> =
            case.chroms.iter().map(|chrom| !chrom.items(0).is_empty()).collect();
        let data_order: Vec<usize> =
            case.push_order().into_iter().filter(|&position| has_data[position]).collect();
        let by_chrom_id: Vec<usize> = if case.options.chrom_id_order == ChromIdOrder::Listed {
            (0..count).collect()
        } else {
            let without_data = (0..count).filter(|&position| !has_data[position]);
            data_order.iter().copied().chain(without_data).collect()
        };
        let mut chrom_ids = vec![0; count];
        for (chrom_id, &position) in (0u32..).zip(&by_chrom_id) {
            chrom_ids[position] = chrom_id;
        }
        let items = case
            .chroms
            .iter()
            .zip(&chrom_ids)
            .map(|(chrom, &chrom_id)| chrom.items(chrom_id))
            .collect();
        Self { by_chrom_id, items, data_order }
    }

    /// Every item, chromosome by chromosome in `chromId` order.
    fn all_items(&self) -> Vec<DataItem> {
        self.by_chrom_id.iter().flat_map(|&position| self.items[position].iter().copied()).collect()
    }

    /// List positions of the chromosomes whose summaries are checked: the first few with data,
    /// and the first without.
    fn summary_chroms(&self) -> Vec<usize> {
        let with_data = (0..self.items.len()).filter(|&p| !self.items[p].is_empty());
        let without_data = (0..self.items.len()).filter(|&p| self.items[p].is_empty());
        with_data.take(SUMMARY_CHROMS).chain(without_data.take(1)).collect()
    }

    /// The reductions of the zoom levels the writer must keep: the policy's ladder, cut at the
    /// first level that doesn't have fewer records than the level before it.
    fn expected_reductions(&self, case: &Case) -> Vec<u32> {
        let limit = case.longest_length().min(MAX_REDUCTION);
        let ladder = match &case.options.zoom {
            ZoomPolicy::Auto { min_first_reduction } => {
                self.auto_ladder(*min_first_reduction, limit)
            }
            ZoomPolicy::Explicit(reductions) => {
                let fitting: Vec<u32> = reductions
                    .iter()
                    .copied()
                    .take_while(|&reduction| reduction <= limit)
                    .collect();
                if fitting.is_empty() { reductions[..1].to_vec() } else { fitting }
            }
            _ => Vec::new(),
        };
        let record_counts: Vec<u64> = ladder
            .iter()
            .map(|&reduction| {
                self.items.iter().map(|items| grid_record_count(items, reduction)).sum()
            })
            .collect();
        let kept = if ladder.is_empty() {
            0
        } else {
            1 + record_counts.windows(2).take_while(|pair| pair[1] < pair[0]).count()
        };
        ladder[..kept].to_vec()
    }

    /// The `Auto` ladder before levels are dropped: the first level is
    /// `max(min_first, min(⌈10·p⌉, limit))` for the mean pitch `p` of the first
    /// [`PITCH_SAMPLE_ITEMS`] items pushed, then ×4 while a level is at most `limit`.
    fn auto_ladder(&self, min_first: u32, limit: u32) -> Vec<u32> {
        let (mut sampled, mut extent) = (0, 0);
        for &position in &self.data_order {
            let items = &self.items[position];
            let taken = items.len().min(PITCH_SAMPLE_ITEMS - sampled);
            if taken == 0 {
                break;
            }
            extent += u64::from(items[taken - 1].end - items[0].start);
            sampled += taken;
        }
        if sampled == 0 {
            return vec![min_first];
        }
        let ten_pitches = (10 * extent).div_ceil(sampled as u64);
        let mut ladder = vec![u64::from(min_first).max(ten_pitches.min(u64::from(limit)))];
        while ladder.len() < AUTO_MAX_LEVELS {
            let next = ladder[ladder.len() - 1] * AUTO_LEVEL_FACTOR;
            if next > u64::from(limit) {
                break;
            }
            ladder.push(next);
        }
        ladder.into_iter().map(to_u32).collect()
    }
}

/// A brute-force summary of covered bases, accumulated in f64.
#[derive(Clone, Copy, Debug)]
struct Exact {
    covered: u64,
    min: f64,
    max: f64,
    sum: f64,
    sum_squares: f64,
    /// Σ|value| over the covered bases: the scale of `sum`'s rounding error, since `sum` itself
    /// may cancel to near 0.
    abs_sum: f64,
}

impl Exact {
    /// Nothing covered yet.
    const EMPTY: Self = Self {
        covered: 0,
        min: f64::INFINITY,
        max: f64::NEG_INFINITY,
        sum: 0.0,
        sum_squares: 0.0,
        abs_sum: 0.0,
    };

    /// The summary of the bases of `items` (sorted, non-overlapping) in `[start, end)`, or `None`
    /// when none is covered.
    fn of(items: &[DataItem], start: u32, end: u32) -> Option<Self> {
        let mut exact = Self::EMPTY;
        for item in overlapping(items, start, end) {
            let bases = item.end.min(end) - item.start.max(start);
            exact.add(u64::from(bases), item.value);
        }
        exact.nonempty()
    }

    /// Adds `bases` bases of `value`.
    fn add(&mut self, bases: u64, value: f32) {
        let (count, value) = (bases as f64, f64::from(value));
        self.covered += bases;
        self.min = self.min.min(value);
        self.max = self.max.max(value);
        self.sum += count * value;
        self.sum_squares += count * value * value;
        self.abs_sum += count * value.abs();
    }

    /// The summary, or `None` if nothing is covered.
    fn nonempty(self) -> Option<Self> {
        (self.covered > 0).then_some(self)
    }
}

/// How far a summary's sums may stray from brute force.
#[derive(Clone, Copy, Debug)]
struct Tolerance {
    /// Allowed error relative to the larger of the sum and its scale.
    relative: f64,
    /// Whether the sums went through f32 storage, which also loses values below f32's normal
    /// range (up to [`f32::MIN_POSITIVE`] per record, and a record covers at least one base).
    f32_storage: bool,
}

impl Tolerance {
    /// Whether `actual` is close enough to `expected`, a sum over `covered` bases whose terms'
    /// magnitudes sum to `scale`.
    fn accepts(self, actual: f64, expected: f64, scale: f64, covered: u64) -> bool {
        let underflow =
            if self.f32_storage { f64::from(f32::MIN_POSITIVE) * covered as f64 } else { 0.0 };
        (actual - expected).abs() <= self.relative * expected.abs().max(scale) + underflow
    }
}

/// Dense fixedStep values with NaN marking missing data, for `push_fixed_step_skip`, alone on a
/// chromosome that ends after them.
#[derive(Clone, Debug)]
struct GappedRun {
    start: u32,
    step: u32,
    span: u32,
    /// Value `i` covers `[start + i·step, start + i·step + span)`, unless it is NaN.
    values: Vec<f32>,
    /// Ascending indices strictly inside `values` at which it may be split into pushes.
    splits: Vec<usize>,
    /// Bases of the chromosome after the last value's step.
    tail: u32,
    options: Options,
}

impl GappedRun {
    /// Length of the chromosome.
    fn chrom_length(&self) -> u32 {
        self.start + to_u32(self.values.len() as u64) * self.step + self.tail
    }

    /// Writes the values to memory with one `push_fixed_step_skip` per slice between `splits`.
    fn write(&self, splits: &[usize]) -> bigwig::Result<Vec<u8>> {
        let chroms = [ChromInfo::new(GAPPED_CHROM, self.chrom_length())];
        let mut writer =
            self.options.writer_options().create_with(Cursor::new(Vec::new()), &chroms)?;
        let mut from = 0;
        for &to in splits.iter().chain([&self.values.len()]) {
            let chunk_start = self.start + to_u32(from as u64) * self.step;
            let chunk = &self.values[from..to];
            writer.push_fixed_step_skip(GAPPED_CHROM, chunk_start, self.step, self.span, chunk)?;
            from = to;
        }
        Ok(writer.finish()?.into_inner())
    }

    /// The value of every base of the chromosome, NaN where no value covers it.
    fn dense_model(&self) -> Vec<f32> {
        let mut bases = vec![f32::NAN; self.chrom_length() as usize];
        for (index, &value) in (0u32..).zip(&self.values) {
            let item_start = (self.start + index * self.step) as usize;
            bases[item_start..item_start + self.span as usize].fill(value);
        }
        bases
    }
}

/// One expected zoom record: the tight extent of a grid bin's covered bases and their summary.
#[derive(Clone, Copy, Debug)]
struct GridRecord {
    start: u32,
    end: u32,
    summary: Exact,
}

/// A proptest that writes each case, opens the file, checks the invariants every file must meet,
/// then runs `$check(case, reader, model)`.
macro_rules! read_property {
    ($(#[$meta:meta])* $name:ident, $profile:expr, $check:path) => {
        proptest! {
            #![proptest_config(config($profile))]
            $(#[$meta])*
            #[test]
            fn $name(case in case_strategy($profile)) {
                let (mut reader, model) = round_trip(&case)?;
                $check(&case, &mut reader, &model)?;
            }
        }
    };
}

/// A proptest that runs `$check(case)`, which writes the case itself.
macro_rules! write_property {
    ($(#[$meta:meta])* $name:ident, $profile:expr, $check:path) => {
        proptest! {
            #![proptest_config(config($profile))]
            $(#[$meta])*
            #[test]
            fn $name(case in case_strategy($profile)) {
                $check(&case)?;
            }
        }
    };
}

read_property!(
    /// Property 1: every chromosome's items, and the whole file's, read back exactly as pushed.
    items_round_trip,
    Profile::Typical,
    check_items_round_trip
);
read_property!(
    /// Property 2: a region query returns exactly the items overlapping the region.
    region_items_match_brute_force,
    Profile::Typical,
    check_region_items
);
read_property!(
    /// Property 3: per-base values equal a dense model, and `values_into` counts covered bases.
    values_match_dense_model,
    Profile::Typical,
    check_values
);
read_property!(
    /// Property 4: exact summaries over 1 to 64 bins equal brute force.
    raw_summaries_match_brute_force,
    Profile::Typical,
    check_raw_summaries
);
read_property!(
    /// Property 5: the kept zoom levels are the policy's ladder, cut where records stop shrinking.
    zoom_levels_follow_the_policy,
    Profile::Typical,
    check_zoom_levels
);
read_property!(
    /// Property 5: each level's records are the grid-aligned records of the data.
    zoom_records_match_grid_model,
    Profile::Typical,
    check_zoom_records
);
read_property!(
    /// Property 5: zoom-based summaries over whole zoom bins equal exact summaries.
    best_zoom_summaries_over_whole_bins_are_exact,
    Profile::Typical,
    check_best_zoom_whole_bins
);
read_property!(
    /// Property 5: zoom-based summaries over arbitrary regions stay within the exact bounds.
    best_zoom_summaries_over_any_region_are_bounded,
    Profile::Typical,
    check_best_zoom_any_region
);
read_property!(
    /// Property 6: the header's total summary equals brute force over all the data.
    total_summary_matches_brute_force,
    Profile::Typical,
    check_total_summary
);
write_property!(
    /// Property 7: writing the same input twice gives identical bytes.
    writing_twice_gives_identical_files,
    Profile::Typical,
    check_deterministic
);
write_property!(
    /// The writer's contract: slicing the same data into different pushes gives identical bytes.
    push_granularity_does_not_change_the_file,
    Profile::Typical,
    check_push_granularity
);
write_property!(
    /// The writer's contract: spilling zoom blocks to a temporary file gives identical bytes.
    zoom_spill_does_not_change_the_file,
    Profile::Typical,
    check_zoom_spill
);
read_property!(
    /// Properties 1, 5 and 6 with enough items that the `Auto` pitch sample fills partway through.
    long_runs_round_trip_with_a_full_pitch_sample,
    Profile::LongRuns,
    check_long_runs
);

proptest! {
    #![proptest_config(config(Profile::Typical))]
    /// `push_fixed_step_skip`: values split into pushes anywhere give the bytes of one push, and
    /// every base reads back as the value covering it, or as missing in a gap.
    #[test]
    fn gapped_fixed_step_pushes_round_trip_however_split(run in gapped_run_strategy()) {
        check_gapped_run(&run)?;
    }
}

/// The proptest configuration: `PROPTEST_CASES` cases if set, otherwise [`DEFAULT_CASES`] (a
/// quarter as many with only the zlib-rs backend, which deflates far more slowly in debug builds;
/// the properties don't depend on the backend), a quarter as many again for the costly
/// [`Profile::LongRuns`], and no regression files.
fn config(profile: Profile) -> ProptestConfig {
    let mut config = ProptestConfig::default();
    if std::env::var_os("PROPTEST_CASES").is_none() {
        config.cases = if cfg!(feature = "libdeflate") { DEFAULT_CASES } else { DEFAULT_CASES / 4 };
    }
    if matches!(profile, Profile::LongRuns) {
        config.cases = (config.cases / 4).max(1);
    }
    config.failure_persistence = None;
    config
}

/// Random cases of `profile`.
fn case_strategy(profile: Profile) -> impl Strategy<Value = Case> {
    chrom_count(profile)
        .prop_flat_map(move |count| {
            (
                vec(chrom_name(), count).prop_map(unique_names),
                vec(chrom_plan(profile), count),
                Just((0..count).collect::<Vec<_>>()).prop_shuffle(),
                options_strategy(),
                vec((any::<u32>(), any::<u32>()), 3),
            )
        })
        .prop_map(move |(names, plans, write_order, options, query_seeds)| {
            Case::lay_out(profile, names, plans, write_order, options, query_seeds)
        })
}

/// Chromosomes per case: 1 to 50, mostly few.
fn chrom_count(profile: Profile) -> BoxedStrategy<usize> {
    match profile {
        Profile::Typical => {
            prop_oneof![6 => 1..=4usize, 3 => 5..=12usize, 1 => 13..=50usize].boxed()
        }
        Profile::LongRuns => (1..=4usize).boxed(),
    }
}

/// Chromosome names, often prefixes of one another (`chr1`, `chr10`, `chr1_alt`).
fn chrom_name() -> impl Strategy<Value = String> {
    let common = ["chr1", "chr10", "chr1_alt", "chr1_random", "chr11", "chr2", "chrX", "1", "10"];
    prop_oneof![
        3 => prop::sample::select(common.to_vec()).prop_map(str::to_owned),
        2 => "chr[0-9]{1,2}(_alt|_random)?",
        1 => "[!-~]{1,16}",
    ]
}

/// `names` with each repeat suffixed by its position, which keeps them unique and adds prefix
/// collisions (`chr1`, `chr1_3`).
fn unique_names(names: Vec<String>) -> Vec<String> {
    let mut seen = std::collections::HashSet::new();
    names
        .into_iter()
        .enumerate()
        .map(|(position, name)| {
            let mut unique = name;
            while !seen.insert(unique.clone()) {
                unique = format!("{unique}_{position}");
            }
            unique
        })
        .collect()
}

/// One chromosome's planned runs, tail, offset and length if empty. Offsets and lengths reach
/// `u32::MAX` only in [`Profile::Typical`].
fn chrom_plan(profile: Profile) -> impl Strategy<Value = ChromPlan> {
    let (offset, empty_length) = match profile {
        Profile::Typical => (
            prop_oneof![6 => Just(0u32), 1 => Just(u32::MAX), 1 => any::<u32>()].boxed(),
            prop_oneof![3 => 1..=10u32, 3 => 1..=200_000u32, 1 => 1..=u32::MAX].boxed(),
        ),
        Profile::LongRuns => (Just(0u32).boxed(), (1..=200_000u32).boxed()),
    };
    (
        vec(unit_plan(profile), 0..=5),
        prop_oneof![2 => Just(0u32), 1 => 0..=10u32, 1 => 0..=50_000u32],
        prop_oneof![2 => Just(0u32), 1 => 1..=10u32, 1 => 1..=5_000u32],
        offset,
        empty_length,
    )
        .prop_map(|(units, tail, overhang, offset, empty_length)| ChromPlan {
            units,
            tail,
            overhang,
            offset,
            empty_length,
        })
}

/// One planned run.
fn unit_plan(profile: Profile) -> impl Strategy<Value = UnitPlan> {
    let gap = prop_oneof![3 => 0..=3u32, 2 => 0..=100u32, 1 => 0..=20_000u32];
    let pushes = prop_oneof![
        2 => Just(Pushes::Whole),
        1 => (1..=64usize).prop_map(Pushes::Chunks),
        1 => Just(Pushes::Single),
    ];
    (gap, unit_shape(profile), pushes).prop_map(|(gap, shape, pushes)| UnitPlan {
        gap,
        shape,
        pushes,
    })
}

/// A run's shape: bedGraph intervals, a variableStep run or a fixedStep run.
fn unit_shape(profile: Profile) -> BoxedStrategy<UnitShape> {
    let small_gap = || prop_oneof![3 => Just(0u32), 2 => 0..=5u32, 1 => 0..=1_000u32];
    let interval_length = prop_oneof![3 => 1..=10u32, 2 => 1..=1_000u32, 1 => 1..=20_000u32];
    let intervals =
        vec((small_gap(), interval_length, value()), 1..=40).prop_map(UnitShape::Intervals);
    let span = prop_oneof![2 => Just(1u32), 2 => 1..=20u32, 1 => 1..=2_000u32];
    let variable_step = (span, vec((small_gap(), value()), 1..=40))
        .prop_map(|(span, items)| UnitShape::VariableStep { span, items });
    let step = prop_oneof![3 => Just(1u32), 2 => 1..=10u32, 1 => 1..=2_000u32];
    let count = prop_oneof![3 => 1..=32u32, 2 => 1..=500u32, 1 => 1..=3_000u32];
    let fixed_step = (
        step.prop_flat_map(|step| (Just(step), prop_oneof![Just(step), 1..=step])),
        vec(value(), 1..=16),
        count,
    )
        .prop_map(|((step, span), pattern, count)| UnitShape::FixedStep {
            step,
            span,
            pattern,
            count,
        });
    match profile {
        Profile::Typical => prop_oneof![intervals, variable_step, fixed_step].boxed(),
        Profile::LongRuns => {
            let per_base =
                (vec(value(), 1..=16), 40_000..=80_000u32).prop_map(|(pattern, count)| {
                    UnitShape::FixedStep { step: 1, span: 1, pattern, count }
                });
            prop_oneof![1 => intervals, 1 => variable_step, 1 => fixed_step, 2 => per_base].boxed()
        }
    }
}

/// Finite values: small integers, zeros of both signs, and a spread of magnitudes small enough
/// that zoom sums of squares stay finite in f32.
fn value() -> impl Strategy<Value = f32> {
    prop_oneof![
        3 => (-8i8..=8).prop_map(f32::from),
        3 => -1_000.0f32..1_000.0,
        1 => Just(0.0f32),
        1 => Just(-0.0f32),
        1 => -1e9f32..1e9,
        1 => prop::sample::select(vec![1e-30f32, -1e-30, f32::MIN_POSITIVE, 1e12, -1e12]),
    ]
}

/// Random [`GappedRun`]s whose values are NaN at a rate of 0, 1 or anything between, with up to
/// seven split points and random writer options.
fn gapped_run_strategy() -> impl Strategy<Value = GappedRun> {
    let nan_rate = prop_oneof![Just(0.0), Just(1.0), 0.0..=1.0f64];
    (1..=10u32, 1..=2_000usize, nan_rate)
        .prop_flat_map(|(step, count, nan_rate)| {
            let values = vec((value(), 0.0..1.0f64), count).prop_map(move |draws| {
                let nan_or =
                    |(value, draw): (f32, f64)| if draw < nan_rate { f32::NAN } else { value };
                draws.into_iter().map(nan_or).collect::<Vec<f32>>()
            });
            (
                Just(step),
                1..=step,
                0..1_000u32,
                values,
                btree_set(1..count.max(2), 0..8),
                0..100u32,
                options_strategy(),
            )
        })
        .prop_map(|(step, span, start, values, splits, tail, options)| {
            let splits = splits.into_iter().filter(|&split| split < values.len()).collect();
            GappedRun { start, step, span, values, splits, tail, options }
        })
}

/// Random writer options.
fn options_strategy() -> impl Strategy<Value = Options> {
    let items_per_section = prop_oneof![3 => 1..=8u16, 2 => 1..=256u16, 1 => 1..=u16::MAX];
    let compression = prop_oneof![
        Just(Compression::None),
        Just(Compression::Deflate { level: 1 }),
        Just(Compression::Deflate { level: 6 }),
    ];
    let byte_order = prop_oneof![Just(ByteOrder::Little), Just(ByteOrder::Big)];
    let chrom_id_order = prop_oneof![Just(ChromIdOrder::Arrival), Just(ChromIdOrder::Listed)];
    let zoom_items_per_block = prop_oneof![2 => 1..=3u32, 1 => 1..=2_048u32];
    let index_block_size = prop_oneof![2 => 2..=4u32, 1 => 2..=512u32];
    (
        items_per_section,
        zoom_policy(),
        compression,
        byte_order,
        chrom_id_order,
        zoom_items_per_block,
        index_block_size,
        prop_oneof![8 => Just(None), 1 => Just(Some(0)), 1 => (1..=4_096usize).prop_map(Some)],
        prop::bool::weighted(0.25),
    )
        .prop_map(
            |(
                items_per_section,
                zoom,
                compression,
                byte_order,
                chrom_id_order,
                zoom_items_per_block,
                index_block_size,
                zoom_memory_limit,
                clip,
            )| Options {
                items_per_section,
                zoom,
                compression,
                byte_order,
                chrom_id_order,
                zoom_items_per_block,
                index_block_size,
                zoom_memory_limit,
                clip,
            },
        )
}

/// `Auto` with a random first reduction (often 1), `Explicit` with random increasing reductions
/// that need not divide one another or with a geometric ladder of up to 14 levels (more than the
/// 10 header slots reserved by default), or occasionally `None`.
fn zoom_policy() -> impl Strategy<Value = ZoomPolicy> {
    let min_first_reduction =
        prop_oneof![Just(1u32), 1..=64u32, 1..=4_096u32, 1..=300_000u32, Just(MAX_REDUCTION)];
    let reduction = prop_oneof![1..=16u32, 1..=1_000u32, 1..=300_000u32, Just(MAX_REDUCTION)];
    prop_oneof![
        4 => min_first_reduction
            .prop_map(|min_first_reduction| ZoomPolicy::Auto { min_first_reduction }),
        4 => btree_set(reduction, 1..=6)
            .prop_map(|reductions| ZoomPolicy::Explicit(reductions.into_iter().collect())),
        1 => (1..=8u32, 2..=3u32, 1..=14u32).prop_map(|(first, factor, levels)| {
            ZoomPolicy::Explicit((0..levels).map(|level| first * factor.pow(level)).collect())
        }),
        1 => Just(ZoomPolicy::None),
    ]
}

/// Writes `case` to memory, slicing every run as `pushes` says, or as the case says if `None`.
fn write(case: &Case, pushes: Option<Pushes>) -> bigwig::Result<Vec<u8>> {
    let mut writer =
        case.options.writer_options().create_with(Cursor::new(Vec::new()), &case.chrom_infos())?;
    for position in case.push_order() {
        let chrom = &case.chroms[position];
        for run in &chrom.runs {
            run.data.push(&mut writer, &chrom.name, pushes.unwrap_or(run.pushes))?;
        }
    }
    Ok(writer.finish()?.into_inner())
}

/// Writes `case`, opens the file and checks the invariants every file must meet.
fn round_trip(case: &Case) -> Result<(Reader, Model), TestCaseError> {
    let reader = BigWigReader::new(Cursor::new(write(case, None)?))?;
    let model = Model::new(case);
    check_invariants(&reader, case, &model)?;
    Ok((reader, model))
}

/// Every chromosome is listed under its dense `chromId`; zoom reductions strictly increase and
/// stay below 2³¹; the header matches the options.
fn check_invariants(reader: &Reader, case: &Case, model: &Model) -> Check {
    let chroms = reader.chroms();
    prop_assert_eq!(chroms.len(), case.chroms.len());
    for (chrom_id, &position) in (0u32..).zip(&model.by_chrom_id) {
        let expected = &case.chroms[position];
        let listed = &chroms[chrom_id as usize];
        prop_assert_eq!(&listed.name, &expected.name, "name of chromId {}", chrom_id);
        prop_assert_eq!(listed.length, expected.length, "length of {}", expected.name);
        prop_assert_eq!(reader.chrom_name(chrom_id), Some(expected.name.as_str()));
        prop_assert_eq!(reader.chrom(&expected.name), Some(listed));
    }
    prop_assert_eq!(reader.chrom_name(to_u32(chroms.len() as u64)), None);

    let reductions: Vec<u32> = reader.zoom_levels().iter().map(ZoomLevel::reduction).collect();
    prop_assert!(
        reductions.windows(2).all(|pair| pair[0] < pair[1]),
        "zoom reductions {:?} are not strictly increasing",
        reductions
    );
    prop_assert!(
        reductions.iter().all(|&reduction| (1..=MAX_REDUCTION).contains(&reduction)),
        "zoom reductions {:?} include one outside 1 to 2^31 - 1",
        reductions
    );
    prop_assert_eq!(reader.zoom_levels().is_empty(), case.options.zoom == ZoomPolicy::None);

    let info = reader.info();
    prop_assert_eq!(info.version, 4);
    prop_assert_eq!(info.byte_order, case.options.byte_order);
    if case.options.compression == Compression::None {
        prop_assert_eq!(info.uncompress_buf_size, 0);
    }
    Ok(())
}

/// Property 1.
fn check_items_round_trip(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    for (chrom, expected) in case.chroms.iter().zip(&model.items) {
        let actual = reader.chrom_items(&chrom.name)?.collect::<bigwig::Result<Vec<_>>>()?;
        check_items(&actual, expected, &format!("chrom_items({})", chrom.name))?;
    }
    let actual = reader.all_items()?.collect::<bigwig::Result<Vec<_>>>()?;
    check_items(&actual, &model.all_items(), "all_items()")?;
    prop_assert_eq!(reader.info().section_count, case.expected_section_count());
    Ok(())
}

/// Property 2.
fn check_region_items(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    for (chrom, items) in case.chroms.iter().zip(&model.items) {
        for (start, end) in case.regions(chrom.length, items) {
            let actual =
                reader.items(&chrom.name, start, end)?.collect::<bigwig::Result<Vec<_>>>()?;
            let expected = overlapping(items, start, end);
            check_items(&actual, expected, &format!("items({}, {start}, {end})", chrom.name))?;
        }
    }
    Ok(())
}

/// Property 3.
fn check_values(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    for (chrom, items) in case.chroms.iter().zip(&model.items) {
        for (start, end) in case.regions(chrom.length, items) {
            if end - start > MAX_DENSE_REGION {
                continue;
            }
            let mut expected = vec![f32::NAN; (end - start) as usize];
            for item in overlapping(items, start, end) {
                let from = item.start.max(start) - start;
                let to = item.end.min(end) - start;
                expected[from as usize..to as usize].fill(item.value);
            }
            let actual = reader.values(&chrom.name, start, end)?;
            check_values_equal(&actual, &expected, start, &chrom.name)?;

            let covered = expected.iter().filter(|value| !value.is_nan()).count() as u64;
            for missing in [0.0, -7.5] {
                let mut out = vec![1.25; expected.len()];
                let count = reader.values_into(&chrom.name, start, end, &mut out, missing)?;
                prop_assert_eq!(
                    count,
                    covered,
                    "values_into({}, {}, {}) covered count with missing {}",
                    chrom.name,
                    start,
                    end,
                    missing
                );
                let with_missing: Vec<f32> = expected
                    .iter()
                    .map(|&value| if value.is_nan() { missing } else { value })
                    .collect();
                check_values_equal(&out, &with_missing, start, &chrom.name)?;
            }
        }
    }
    Ok(())
}

/// Property 4.
fn check_raw_summaries(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    for position in model.summary_chroms() {
        let (chrom, items) = (&case.chroms[position], &model.items[position]);
        let mut regions = vec![(0, chrom.length), seeded_region(chrom.length, case.query_seeds[0])];
        if !items.is_empty() {
            // A short region, so some bin counts exceed its length.
            let item = items[case.query_seeds[0].1 as usize % items.len()];
            regions
                .push((item.start.saturating_sub(1), item.end.saturating_add(1).min(chrom.length)));
        }
        for (start, end) in regions {
            let single = reader.summary(&chrom.name, start, end, SummarySource::Raw)?;
            check_summary(single, Exact::of(items, start, end), RAW_TOLERANCE, &|| {
                format!("summary({}, {start}, {end})", chrom.name)
            })?;
            for n_bins in 1..=64 {
                let bins = reader.summaries(&chrom.name, start, end, n_bins, SummarySource::Raw)?;
                prop_assert_eq!(bins.len(), n_bins as usize);
                for (index, &actual) in (0u32..).zip(&bins) {
                    let (bin_start, bin_end) = bin_bounds(start, end, n_bins, index);
                    let expected = Exact::of(items, bin_start, bin_end);
                    check_summary(actual, expected, RAW_TOLERANCE, &|| {
                        format!(
                            "bin {index} [{bin_start}, {bin_end}) of {n_bins} over {}:{start}-{end}",
                            chrom.name
                        )
                    })?;
                }
            }
        }
    }
    Ok(())
}

/// Property 5: the kept levels.
fn check_zoom_levels(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    let reductions: Vec<u32> = reader.zoom_levels().iter().map(ZoomLevel::reduction).collect();
    prop_assert_eq!(&reductions, &model.expected_reductions(case));

    // `best_zoom_for` follows kent's rule: the largest reduction strictly below
    // ⌊width / 2⌋, and raw data when ⌊width / 2⌋ ≤ 1. Widths around 2r and 4r probe the edge.
    let seeded = case.query_seeds.iter().map(|&(a, _)| u64::from(a));
    let near_levels = reductions
        .iter()
        .map(|&r| u64::from(r))
        .flat_map(|r| [2 * r, 2 * r + 1, 2 * r + 2, 2 * r + 3, 4 * r, 4 * r + 1]);
    let widths = (0..=6).chain(near_levels).chain(seeded).filter_map(|w| u32::try_from(w).ok());
    for width in widths {
        let desired = width / 2;
        let expected = reductions
            .iter()
            .copied()
            .filter(|&reduction| desired > 1 && reduction < desired)
            .max();
        let actual = reader.best_zoom_for(width).map(|level| level.reduction());
        prop_assert_eq!(actual, expected, "best_zoom_for({}) with levels {:?}", width, reductions);
    }
    Ok(())
}

/// Property 5: every level's records, whole-chromosome and in sub-regions.
fn check_zoom_records(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    for level in reader.zoom_levels().to_vec() {
        for (position, chrom) in case.chroms.iter().enumerate() {
            let items = &model.items[position];
            let chrom_id = model.by_chrom_id.iter().position(|&p| p == position).unwrap_or(0);
            let expected = grid_records(items, level.reduction());
            for (start, end) in case.regions(chrom.length, items) {
                let actual = reader
                    .zoom_records(level, &chrom.name, start, end)?
                    .collect::<bigwig::Result<Vec<_>>>()?;
                let wanted: Vec<GridRecord> = expected
                    .iter()
                    .copied()
                    .filter(|record| start < end && record.start < end && record.end > start)
                    .collect();
                let what = || {
                    format!(
                        "zoom_records(level {}, {}, {start}, {end})",
                        level.reduction(),
                        chrom.name
                    )
                };
                prop_assert_eq!(actual.len(), wanted.len(), "record count of {}", what());
                for (record, expected) in actual.iter().zip(&wanted) {
                    check_record(record, expected, to_u32(chrom_id as u64), &what)?;
                }
            }
        }
    }
    Ok(())
}

/// Property 5: `BestZoom` over regions made of whole zoom bins, wherever a bin width selects the
/// level, equals the exact summary up to f32 rounding.
fn check_best_zoom_whole_bins(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    let (seed_a, seed_b) = case.query_seeds[1];
    for level in reader.zoom_levels().to_vec() {
        let reduction = u64::from(level.reduction());
        for position in model.summary_chroms() {
            let (chrom, items) = (&case.chroms[position], &model.items[position]);
            let length = u64::from(chrom.length);
            let bins_in_chrom = length / reduction;
            // Bin widths of whole zoom bins for which kent's rule picks this level.
            let widths: Vec<u64> = (2..=64)
                .map(|multiple| multiple * reduction)
                .filter(|&width| width <= length)
                .filter(|&width| reader.best_zoom_for(to_u32(width)) == Some(level))
                .collect();
            let mut regions = Vec::new();
            if let Some(&width) = widths.get(seed_a as usize % widths.len().max(1)) {
                let n_bins = u64::from(seed_b % 4 + 1).min(length / width);
                let span = n_bins * width;
                let first_bin = u64::from(seed_a) % (bins_in_chrom + 1);
                let start = (first_bin * reduction).min((length - span) / reduction * reduction);
                regions.push((start, start + span, n_bins));
            }
            // From a bin start to the chromosome end, where the last zoom bin is clipped.
            if bins_in_chrom > 0 {
                let start = u64::from(seed_b) % bins_in_chrom * reduction;
                if reader.best_zoom_for(to_u32(length - start)) == Some(level) {
                    regions.push((start, length, 1));
                }
            }
            for (start, end, n_bins) in regions {
                let (start, end, n_bins) = (to_u32(start), to_u32(end), to_u32(n_bins));
                let bins =
                    reader.summaries(&chrom.name, start, end, n_bins, SummarySource::BestZoom)?;
                for (index, &actual) in (0u32..).zip(&bins) {
                    let (bin_start, bin_end) = bin_bounds(start, end, n_bins, index);
                    let expected = Exact::of(items, bin_start, bin_end);
                    check_summary(actual, expected, ZOOM_TOLERANCE, &|| {
                        format!(
                            "BestZoom bin [{bin_start}, {bin_end}) of {}:{start}-{end} at level {}",
                            chrom.name, reduction
                        )
                    })?;
                }
            }
        }
    }
    Ok(())
}

/// Property 5: `BestZoom` over arbitrary regions. With a zoom level, `min` and `max` lie between
/// the exact values over the bin and over the bin widened to whole zoom bins, and
/// `covered_bases` is within one zoom bin per partial edge; without one, the result is exact.
fn check_best_zoom_any_region(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    for position in model.summary_chroms() {
        let (chrom, items) = (&case.chroms[position], &model.items[position]);
        for (start, end) in case.regions(chrom.length, items) {
            for n_bins in [1, 2, 5] {
                let level = reader.best_zoom_for((end - start) / n_bins);
                let bins =
                    reader.summaries(&chrom.name, start, end, n_bins, SummarySource::BestZoom)?;
                for (index, &actual) in (0u32..).zip(&bins) {
                    let (bin_start, bin_end) = bin_bounds(start, end, n_bins, index);
                    let exact = Exact::of(items, bin_start, bin_end);
                    let what = || {
                        format!(
                            "BestZoom bin [{bin_start}, {bin_end}) of {n_bins} over {}:{start}-{end}",
                            chrom.name
                        )
                    };
                    match level {
                        None => check_summary(actual, exact, RAW_TOLERANCE, &what)?,
                        Some(level) => check_zoom_bounds(
                            actual,
                            items,
                            (bin_start, bin_end),
                            chrom.length,
                            level.reduction(),
                            &what,
                        )?,
                    }
                }
            }
        }
    }
    Ok(())
}

/// Property 6.
fn check_total_summary(_case: &Case, reader: &mut Reader, model: &Model) -> Check {
    let mut expected = Exact::EMPTY;
    for item in model.all_items() {
        expected.add(u64::from(item.end - item.start), item.value);
    }
    let stored = reader.info().total_summary;
    let computed = reader.compute_total_summary()?;
    for (summary, what) in [(stored, "header total summary"), (Some(computed), "computed")] {
        let Some(summary) = summary else {
            return Err(TestCaseError::fail("the header has no total summary"));
        };
        if expected.covered > 0 {
            check_summary(Some(summary), Some(expected), RAW_TOLERANCE, &|| what.to_owned())?;
        } else {
            let zero =
                Summary { covered_bases: 0.0, min: 0.0, max: 0.0, sum: 0.0, sum_squares: 0.0 };
            prop_assert_eq!(summary, zero, "{} of a file without data", what);
        }
    }
    Ok(())
}

/// Property 7, without a parallel executor.
fn check_deterministic(case: &Case) -> Check {
    let first = write(case, None)?;
    let second = write(case, None)?;
    prop_assert!(first == second, "two writes of the same case differ");
    Ok(())
}

/// Slicing every run whole, item by item, or as generated gives the same bytes.
fn check_push_granularity(case: &Case) -> Check {
    let generated = write(case, None)?;
    let whole = write(case, Some(Pushes::Whole))?;
    let single = write(case, Some(Pushes::Single))?;
    prop_assert!(generated == whole, "pushing runs as generated differs from pushing them whole");
    prop_assert!(single == whole, "pushing items one by one differs from pushing runs whole");
    Ok(())
}

/// Spilling no zoom block, every block, or those beyond a small limit gives the same bytes.
fn check_zoom_spill(case: &Case) -> Check {
    let with_limit = |limit| {
        let mut case = case.clone();
        case.options.zoom_memory_limit = limit;
        write(&case, None)
    };
    let in_memory = with_limit(None)?;
    prop_assert!(with_limit(Some(0))? == in_memory, "spilling every zoom block changes the file");
    prop_assert!(
        with_limit(Some(1_000))? == in_memory,
        "spilling the zoom blocks beyond 1,000 bytes changes the file"
    );
    Ok(())
}

/// One push of a [`GappedRun`] and pushes split at its split points give the same bytes, and
/// `values_into` with NaN for missing data reads back its dense model, counting the covered
/// bases.
fn check_gapped_run(run: &GappedRun) -> Check {
    let whole = run.write(&[])?;
    prop_assert!(
        run.write(&run.splits)? == whole,
        "splitting the pushes at {:?} changes the file",
        run.splits
    );
    let mut reader = BigWigReader::new(Cursor::new(whole))?;
    let expected = run.dense_model();
    let mut actual = vec![0.0; expected.len()];
    let covered = reader.values_into(GAPPED_CHROM, 0, run.chrom_length(), &mut actual, f32::NAN)?;
    check_values_equal(&actual, &expected, 0, GAPPED_CHROM)?;
    prop_assert_eq!(covered, expected.iter().filter(|value| !value.is_nan()).count() as u64);
    Ok(())
}

/// Properties 1, 5 and 6 on [`Profile::LongRuns`] cases.
fn check_long_runs(case: &Case, reader: &mut Reader, model: &Model) -> Check {
    check_items_round_trip(case, reader, model)?;
    check_zoom_levels(case, reader, model)?;
    check_zoom_records(case, reader, model)?;
    check_total_summary(case, reader, model)
}

/// Checks that `actual` items equal `expected`, value bits included, reporting the first
/// difference.
fn check_items(actual: &[DataItem], expected: &[DataItem], what: &str) -> Check {
    let key = |item: &DataItem| (item.chrom_id, item.start, item.end, item.value.to_bits());
    if let Some(index) = (0..actual.len().min(expected.len()))
        .find(|&index| key(&actual[index]) != key(&expected[index]))
    {
        return Err(TestCaseError::fail(format!(
            "{what}: item {index} is {:?}, expected {:?}",
            actual[index], expected[index]
        )));
    }
    prop_assert_eq!(
        actual.len(),
        expected.len(),
        "{}: item counts differ; first extra item {:?}",
        what,
        actual.get(expected.len()).or(expected.get(actual.len()))
    );
    Ok(())
}

/// Checks that per-base values starting at base `start` equal `expected`, value bits included,
/// with any NaN matching any NaN.
fn check_values_equal(actual: &[f32], expected: &[f32], start: u32, chrom: &str) -> Check {
    prop_assert_eq!(actual.len(), expected.len());
    let same = |a: f32, b: f32| (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits();
    if let Some(offset) = (0..actual.len()).find(|&i| !same(actual[i], expected[i])) {
        return Err(TestCaseError::fail(format!(
            "value of {chrom}:{} is {}, expected {}",
            start as usize + offset,
            actual[offset],
            expected[offset]
        )));
    }
    Ok(())
}

/// Checks a summary against brute force: both `None`, or equal counts, `min` and `max`, and sums
/// within `tolerance`.
fn check_summary(
    actual: Option<Summary>,
    expected: Option<Exact>,
    tolerance: Tolerance,
    what: &dyn Fn() -> String,
) -> Check {
    match (actual, expected) {
        (None, None) => Ok(()),
        (Some(actual), Some(expected)) => {
            prop_assert!(
                actual.covered_bases == expected.covered as f64
                    && actual.min == expected.min
                    && actual.max == expected.max
                    && tolerance.accepts(
                        actual.sum,
                        expected.sum,
                        expected.abs_sum,
                        expected.covered
                    )
                    && tolerance.accepts(
                        actual.sum_squares,
                        expected.sum_squares,
                        0.0,
                        expected.covered
                    ),
                "{}: {:?}, expected {:?}",
                what(),
                actual,
                expected
            );
            Ok(())
        }
        (actual, expected) => {
            Err(TestCaseError::fail(format!("{}: {actual:?}, expected {expected:?}", what())))
        }
    }
}

/// Checks a zoom record against its expected grid record.
fn check_record(
    record: &ZoomRecord,
    expected: &GridRecord,
    chrom_id: u32,
    what: &dyn Fn() -> String,
) -> Check {
    prop_assert!(
        record.chrom_id == chrom_id && record.start == expected.start && record.end == expected.end,
        "{}: record {:?}, expected chromId {} over [{}, {})",
        what(),
        record,
        chrom_id,
        expected.start,
        expected.end
    );
    check_summary(Some(record.summary), Some(expected.summary), ZOOM_TOLERANCE, &|| {
        format!("{} record [{}, {})", what(), expected.start, expected.end)
    })
}

/// Checks a zoom-based summary of `bin` against the bounds pro-rating implies (see
/// [`check_best_zoom_any_region`]).
fn check_zoom_bounds(
    actual: Option<Summary>,
    items: &[DataItem],
    (bin_start, bin_end): (u32, u32),
    chrom_length: u32,
    reduction: u32,
    what: &dyn Fn() -> String,
) -> Check {
    let exact = Exact::of(items, bin_start, bin_end);
    let reduction64 = u64::from(reduction);
    let widened_start = to_u32(u64::from(bin_start) / reduction64 * reduction64);
    let widened_end = to_u32(
        (u64::from(bin_end).div_ceil(reduction64) * reduction64).min(u64::from(chrom_length)),
    );
    let widened = Exact::of(items, widened_start, widened_end);
    let Some(actual) = actual else {
        prop_assert!(exact.is_none(), "{}: None, but the bin has data {:?}", what(), exact);
        return Ok(());
    };
    let Some(widened) = widened else {
        return Err(TestCaseError::fail(format!(
            "{}: {actual:?} where no zoom record can be",
            what()
        )));
    };
    let (exact_min, exact_max) =
        exact.map_or((f64::INFINITY, f64::NEG_INFINITY), |e| (e.min, e.max));
    prop_assert!(
        widened.min <= actual.min && actual.min <= exact_min,
        "{}: min {} outside [{}, {}]",
        what(),
        actual.min,
        widened.min,
        exact_min
    );
    prop_assert!(
        exact_max <= actual.max && actual.max <= widened.max,
        "{}: max {} outside [{}, {}]",
        what(),
        actual.max,
        exact_max,
        widened.max
    );
    let partial_first = bin_start % reduction != 0;
    let partial_last = bin_end % reduction != 0;
    let partial_edges = if bin_start / reduction == (bin_end - 1) / reduction {
        u32::from(partial_first || partial_last)
    } else {
        u32::from(partial_first) + u32::from(partial_last)
    };
    let allowance = f64::from(partial_edges) * f64::from(reduction) + 1e-6;
    let exact_covered = exact.map_or(0.0, |e| e.covered as f64);
    prop_assert!(
        (actual.covered_bases - exact_covered).abs() <= allowance,
        "{}: covered {} is more than {} from the exact {}",
        what(),
        actual.covered_bases,
        allowance,
        exact_covered
    );
    Ok(())
}

/// The items of `items` (sorted, non-overlapping) that share a base with `[start, end)`; none for
/// an empty region.
fn overlapping(items: &[DataItem], start: u32, end: u32) -> &[DataItem] {
    if start >= end {
        return &[];
    }
    let first = items.partition_point(|item| item.end <= start);
    let last = items.partition_point(|item| item.start < end);
    &items[first..last.max(first)]
}

/// The grid-aligned zoom records of one chromosome's `items` at `reduction`: one per bin
/// `[k·r, (k+1)·r)` with a covered base, over the tight extent of its covered bases.
fn grid_records(items: &[DataItem], reduction: u32) -> Vec<GridRecord> {
    let reduction = u64::from(reduction);
    let mut records: Vec<GridRecord> = Vec::new();
    let mut open: Option<(u64, GridRecord)> = None;
    for item in items {
        let (mut position, end) = (u64::from(item.start), u64::from(item.end));
        while position < end {
            let bin = position / reduction;
            let piece_end = end.min((bin + 1) * reduction);
            match &mut open {
                Some((open_bin, record)) if *open_bin == bin => record.end = to_u32(piece_end),
                _ => {
                    records.extend(open.take().map(|(_, record)| record));
                    let record = GridRecord {
                        start: to_u32(position),
                        end: to_u32(piece_end),
                        summary: Exact::EMPTY,
                    };
                    open = Some((bin, record));
                }
            }
            if let Some((_, record)) = &mut open {
                record.summary.add(piece_end - position, item.value);
            }
            position = piece_end;
        }
    }
    records.extend(open.map(|(_, record)| record));
    records
}

/// Number of grid bins at `reduction` holding at least one base of `items` (sorted,
/// non-overlapping): the number of records [`grid_records`] gives, without building them.
fn grid_record_count(items: &[DataItem], reduction: u32) -> u64 {
    let reduction = u64::from(reduction);
    let mut count = 0;
    let mut previous_last_bin = None;
    for item in items {
        let first_bin = u64::from(item.start) / reduction;
        let last_bin = (u64::from(item.end) - 1) / reduction;
        count += last_bin - first_bin + 1;
        if previous_last_bin == Some(first_bin) {
            count -= 1;
        }
        previous_last_bin = Some(last_bin);
    }
    count
}

/// Bin `index` of `count` equal bins over `[start, end)`, as kent computes them.
fn bin_bounds(start: u32, end: u32, count: u32, index: u32) -> (u32, u32) {
    let length = u64::from(end - start);
    let offset = |index: u32| to_u32(length * u64::from(index) / u64::from(count));
    (start + offset(index), start + offset(index + 1))
}

/// A region of a chromosome of `length` derived from `seed`: any start, then any end after it.
fn seeded_region(length: u32, (a, b): (u32, u32)) -> (u32, u32) {
    let length = u64::from(length);
    let start = u64::from(a) % (length + 1);
    let end = start + u64::from(b) % (length - start + 1);
    (to_u32(start), to_u32(end))
}

/// `items` clipped to a chromosome of `length`: those starting at or past the end are dropped and
/// those running past it are truncated, as the writer's `clip(true)` does.
fn clip_items(items: Vec<DataItem>, length: u32) -> Vec<DataItem> {
    items
        .into_iter()
        .filter(|item| item.start < length)
        .map(|item| DataItem { end: item.end.min(length), ..item })
        .collect()
}

/// Narrows a position or count known to fit in u32.
fn to_u32(value: u64) -> u32 {
    u32::try_from(value).expect("generated positions and counts fit in u32")
}
