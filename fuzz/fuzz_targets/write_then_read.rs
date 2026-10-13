#![no_main]

//! Structured fuzzing of the writer: a random chromosome list, options and push sequence (valid
//! and invalid) is replayed against a simple model of the documented validation rules. Every
//! push must be accepted or rejected exactly as the model says, and after `finish` the file must
//! read back as exactly the accepted items.

use std::io::Cursor;

use arbitrary::Arbitrary;
use bigwig::{
    BigWigReader, BigWigWriter, ByteOrder, ChromIdOrder, ChromInfo, Compression, SummarySource,
    WriterOptions, ZoomPolicy,
};
use libfuzzer_sys::fuzz_target;

/// Longest chromosome generated; keeps zoom work per accepted item bounded.
const MAX_CHROM_LENGTH: u32 = 1 << 14;
/// Chromosomes at most this long also have `values()` checked against the model.
const MAX_VALUES_CHECK_LENGTH: u32 = 1 << 14;
const MAX_CHROMS: usize = 6;
const MAX_OPS: usize = 40;
const MAX_ITEMS_PER_OP: usize = 64;

#[derive(Arbitrary, Debug)]
struct Plan {
    chrom_lengths: Vec<u32>,
    options: OptionsSpec,
    ops: Vec<Op>,
}

#[derive(Arbitrary, Debug)]
struct OptionsSpec {
    items_per_section: u8,
    zoom: ZoomSpec,
    zoom_items_per_block: u8,
    index_block_size: u8,
    compression: CompressionSpec,
    spill_zoom_blocks: bool,
    listed_order: bool,
    clip: bool,
    big_endian: bool,
}

#[derive(Arbitrary, Debug)]
enum ZoomSpec {
    Default,
    Auto(u16),
    Explicit(Vec<u16>),
    None,
}

#[derive(Arbitrary, Debug)]
enum CompressionSpec {
    None,
    Deflate(u8),
}

#[derive(Arbitrary, Debug)]
struct Op {
    chrom: u8,
    kind: Kind,
}

#[derive(Arbitrary, Debug)]
enum Kind {
    Intervals(Vec<(Coord, u16, Val)>),
    Variable { span: u32, items: Vec<(Coord, Val)> },
    Fixed { start: Coord, step: u32, span: u32, values: Vec<Val> },
    FixedSkip { start: Coord, step: u32, span: u32, values: Vec<Option<Val>> },
}

/// A coordinate, usually near where the chromosome's data currently ends so most pushes are valid.
#[derive(Arbitrary, Debug, Clone, Copy)]
enum Coord {
    After(u16),
    Before(u16),
    Absolute(u32),
}

#[derive(Arbitrary, Debug, Clone, Copy)]
enum Val {
    Quarters(i8),
    Raw(f32),
}

impl Val {
    fn get(self) -> f32 {
        match self {
            Val::Quarters(quarters) => f32::from(quarters) / 4.0,
            Val::Raw(value) => value,
        }
    }
}

impl Coord {
    fn resolve(self, base: u64) -> u32 {
        let clamp = |wide: u64| u32::try_from(wide).unwrap_or(u32::MAX);
        match self {
            Coord::After(gap) => clamp(base + u64::from(gap)),
            Coord::Before(gap) => clamp(base.saturating_sub(u64::from(gap))),
            Coord::Absolute(value) => value,
        }
    }
}

/// Mostly small, sometimes raw, so both ordinary and overflow-prone values are exercised.
fn small_or_raw(value: u32) -> u32 {
    if value % 4 == 0 { value } else { value % 64 }
}

/// An item as stored: `[start, end)` and the bit pattern of its value.
type Item = (u32, u32, u32);

/// The writer's documented validation rules, as a model of which pushes it accepts and what the
/// file then holds.
struct Model {
    lengths: Vec<u32>,
    listed: bool,
    clip: bool,
    has_data: Vec<bool>,
    /// Each chromosome's end of its last accepted item before clipping, which orders later
    /// items even when clipping dropped it.
    last_end: Vec<u64>,
    /// The chromosome receiving data.
    current: Option<usize>,
    /// Chromosomes in the order they first received data.
    arrival: Vec<usize>,
    items: Vec<Vec<Item>>,
}

impl Model {
    fn new(lengths: Vec<u32>, listed: bool, clip: bool) -> Self {
        let n = lengths.len();
        Self {
            lengths,
            listed,
            clip,
            has_data: vec![false; n],
            last_end: vec![0; n],
            current: None,
            arrival: Vec::new(),
            items: vec![Vec::new(); n],
        }
    }

    /// The chromosome a push may go to and the end of its last accepted item, or `None` if the
    /// push is rejected on account of the chromosome.
    fn target(&self, index: Option<usize>) -> Option<(usize, u64)> {
        let index = index?;
        if let Some(current) = self.current {
            if current == index {
                return Some((index, self.last_end[index]));
            }
            if self.listed && index < current {
                return None;
            }
        }
        if self.has_data[index] { None } else { Some((index, self.last_end[index])) }
    }

    /// Records accepted items `[start, end, value)` (not yet clipped): every one orders later
    /// pushes, and those starting inside the chromosome are written.
    fn commit(&mut self, index: usize, accepted: &[(u64, u64, f32)]) {
        let Some(&(_, last_end, _)) = accepted.last() else { return };
        self.last_end[index] = last_end;
        let length = u64::from(self.lengths[index]);
        let kept: Vec<_> = accepted.iter().filter(|item| item.0 < length).collect();
        if kept.is_empty() {
            return;
        }
        if !self.has_data[index] {
            self.has_data[index] = true;
            self.arrival.push(index);
        }
        for &&(start, end, value) in &kept {
            self.items[index].push((start as u32, end.min(length) as u32, value.to_bits()));
        }
        self.current = Some(index);
    }

    /// Whether the items `[start, end, value)` in order pass the per-item rules: sorted, within
    /// the chromosome (unless clipping) and finite.
    fn items_valid(&self, index: usize, mut previous_end: u64, items: &[(u64, u64, f32)]) -> bool {
        let length = u64::from(self.lengths[index]);
        for &(start, end, value) in items {
            if start < previous_end || (end > length && !self.clip) || !value.is_finite() {
                return false;
            }
            previous_end = end;
        }
        true
    }
}

fn chrom_name(index: usize) -> String {
    format!("c{index}")
}

fn build_options(spec: &OptionsSpec) -> WriterOptions {
    let mut options = WriterOptions::new()
        .clip(spec.clip)
        .byte_order(if spec.big_endian { ByteOrder::Big } else { ByteOrder::Little })
        .chrom_id_order(if spec.listed_order {
            ChromIdOrder::Listed
        } else {
            ChromIdOrder::Arrival
        })
        .compression(match spec.compression {
            CompressionSpec::None => Compression::None,
            CompressionSpec::Deflate(level) => Compression::Deflate { level: level % 10 },
        });
    if spec.items_per_section < 240 {
        options = options.items_per_section(u16::from(spec.items_per_section % 16) + 1);
    }
    if spec.zoom_items_per_block < 200 {
        options = options.zoom_items_per_block(u32::from(spec.zoom_items_per_block % 8) + 1);
    }
    if spec.index_block_size < 200 {
        options = options.index_block_size(u32::from(spec.index_block_size % 15) + 2);
    }
    if spec.spill_zoom_blocks {
        options = options.zoom_memory_limit(0);
    }
    match &spec.zoom {
        ZoomSpec::Default => {}
        ZoomSpec::Auto(min) => {
            options =
                options.zoom(ZoomPolicy::Auto { min_first_reduction: u32::from(*min).max(1) });
        }
        ZoomSpec::Explicit(reductions) => {
            let mut reductions: Vec<u32> =
                reductions.iter().map(|&r| u32::from(r).max(1)).collect();
            reductions.sort_unstable();
            reductions.dedup();
            if reductions.is_empty() {
                reductions.push(1);
            }
            options = options.zoom(ZoomPolicy::Explicit(reductions));
        }
        ZoomSpec::None => options = options.zoom(ZoomPolicy::None),
    }
    options
}

/// Applies one op to the writer and the model, asserting they agree on acceptance.
fn apply(writer: &mut BigWigWriter<Cursor<Vec<u8>>>, model: &mut Model, op: &Op) {
    let n = model.lengths.len();
    let pick = usize::from(op.chrom) % (n + 2);
    let (name, index) = match pick {
        p if p < n => (chrom_name(p), Some(p)),
        p if p == n => ("unknown".to_owned(), None),
        _ => {
            let index = model.current.unwrap_or(0);
            (chrom_name(index), Some(index))
        }
    };
    let target = model.target(index);
    let base = target.map_or(0, |(_, end)| end);

    let (result, expected) = match &op.kind {
        Kind::Intervals(raw) => {
            let mut previous_end = base;
            let mut items = Vec::new();
            let mut intervals = Vec::new();
            for &(coord, len, value) in raw.iter().take(MAX_ITEMS_PER_OP) {
                let start = coord.resolve(previous_end);
                let end = u32::try_from(u64::from(start) + u64::from(len)).unwrap_or(u32::MAX);
                previous_end = u64::from(end);
                intervals.push((start, end, value.get()));
                items.push((u64::from(start), u64::from(end), value.get()));
            }
            let valid = !items.is_empty()
                && items.iter().all(|item| item.0 < item.1)
                && target.is_some_and(|(i, prev)| model.items_valid(i, prev, &items));
            (writer.push_intervals(&name, &intervals), valid.then_some(items))
        }
        Kind::Variable { span, items: raw } => {
            let span = small_or_raw(*span);
            let mut previous_end = base;
            let mut items = Vec::new();
            let mut pairs = Vec::new();
            for &(coord, value) in raw.iter().take(MAX_ITEMS_PER_OP) {
                let start = coord.resolve(previous_end);
                previous_end = u64::from(start) + u64::from(span);
                pairs.push((start, value.get()));
                items.push((u64::from(start), previous_end, value.get()));
            }
            let valid = !items.is_empty()
                && span >= 1
                && target.is_some_and(|(i, prev)| model.items_valid(i, prev, &items));
            (writer.push_variable_step(&name, span, &pairs), valid.then_some(items))
        }
        Kind::Fixed { start, step, span, values: raw } => {
            let start = start.resolve(base);
            let step = small_or_raw(*step);
            let span = if *span % 3 == 0 { step } else { small_or_raw(*span) };
            let values: Vec<f32> = raw.iter().take(MAX_ITEMS_PER_OP).map(|v| v.get()).collect();
            let items: Vec<(u64, u64, f32)> = values
                .iter()
                .enumerate()
                .map(|(i, &value)| {
                    let item_start = u64::from(start) + i as u64 * u64::from(step);
                    (item_start, item_start + u64::from(span), value)
                })
                .collect();
            // Non-finite values are rejected anywhere in the push, even among items that
            // clipping would drop. Overlap cannot arise because `span <= step`.
            let valid = !items.is_empty()
                && step >= 1
                && span >= 1
                && span <= step
                && values.iter().all(|v| v.is_finite())
                && target.is_some_and(|(i, prev)| {
                    u64::from(start) >= prev && model.items_valid_bounds_only(i, &items)
                });
            (writer.push_fixed_step(&name, start, step, span, &values), valid.then_some(items))
        }
        Kind::FixedSkip { start, step, span, values: raw } => {
            let start = start.resolve(base);
            let step = small_or_raw(*step);
            let span = if *span % 3 == 0 { step } else { small_or_raw(*span) };
            let values: Vec<f32> = raw
                .iter()
                .take(MAX_ITEMS_PER_OP)
                .map(|value| value.map_or(f32::NAN, Val::get))
                .collect();
            // `None` values are NaN, which `push_fixed_step_skip` takes as gaps: only the other
            // values are items.
            let items: Vec<(u64, u64, f32)> = values
                .iter()
                .enumerate()
                .filter(|(_, value)| !value.is_nan())
                .map(|(i, &value)| {
                    let item_start = u64::from(start) + i as u64 * u64::from(step);
                    (item_start, item_start + u64::from(span), value)
                })
                .collect();
            let arguments_valid = !values.is_empty() && step >= 1 && span >= 1 && span <= step;
            // Gaps alone are checked only for the chromosome name and the arguments, and are no
            // data. Otherwise the first item is ordered, the items are bounded, and ±infinity is
            // rejected anywhere.
            let valid = arguments_valid
                && match items.first() {
                    None => index.is_some(),
                    Some(&(first_start, _, _)) => {
                        items.iter().all(|item| item.2.is_finite())
                            && target.is_some_and(|(i, prev)| {
                                first_start >= prev && model.items_valid_bounds_only(i, &items)
                            })
                    }
                };
            (writer.push_fixed_step_skip(&name, start, step, span, &values), valid.then_some(items))
        }
    };

    assert_eq!(
        result.is_ok(),
        expected.is_some(),
        "writer and model disagree on {op:?} -> {result:?}; model state: current {:?}",
        model.current
    );
    if let (Some(items), Some((index, _))) = (expected, target) {
        model.commit(index, &items);
    }
}

impl Model {
    /// The bounds rule alone: without clipping every item must end within the chromosome.
    fn items_valid_bounds_only(&self, index: usize, items: &[(u64, u64, f32)]) -> bool {
        let length = u64::from(self.lengths[index]);
        self.clip || items.iter().all(|item| item.1 <= length)
    }
}

fuzz_target!(|plan: Plan| {
    let mut lengths: Vec<u32> = plan
        .chrom_lengths
        .iter()
        .take(MAX_CHROMS)
        .map(
            |&raw| if raw & 0x8000_0000 != 0 { 1 + raw % MAX_CHROM_LENGTH } else { 1 + raw % 4096 },
        )
        .collect();
    if lengths.is_empty() {
        lengths.push(100);
    }
    let chroms: Vec<ChromInfo> = lengths
        .iter()
        .enumerate()
        .map(|(i, &length)| ChromInfo::new(chrom_name(i), length))
        .collect();

    let options = build_options(&plan.options);
    let mut writer = options
        .create_with(Cursor::new(Vec::new()), &chroms)
        .expect("generated options and chromosomes are valid");
    let mut model = Model::new(lengths.clone(), plan.options.listed_order, plan.options.clip);
    for op in plan.ops.iter().take(MAX_OPS) {
        apply(&mut writer, &mut model, op);
    }
    let bytes = writer.finish().expect("finish after only validation errors").into_inner();

    let mut reader = BigWigReader::new(Cursor::new(bytes)).expect("a finished file opens");

    let mut read_chroms: Vec<(String, u32)> =
        reader.chroms().iter().map(|c| (c.name.clone(), c.length)).collect();
    read_chroms.sort();
    let mut written_chroms: Vec<(String, u32)> =
        chroms.iter().map(|c| (c.name.clone(), c.length)).collect();
    written_chroms.sort();
    assert_eq!(read_chroms, written_chroms, "chromosome list");

    // all_items yields chromosomes in chromId order, which is first-data order (Arrival) or list
    // order (Listed); both are the order chromosomes first received data.
    let expected: Vec<(String, u32, u32, u32)> = model
        .arrival
        .iter()
        .flat_map(|&i| model.items[i].iter().map(move |&(s, e, v)| (chrom_name(i), s, e, v)))
        .collect();
    let names_by_id: Vec<String> = (0..lengths.len() as u32)
        .map(|id| reader.chrom_name(id).expect("dense chromIds").to_owned())
        .collect();
    let actual: Vec<(String, u32, u32, u32)> = reader
        .all_items()
        .expect("all_items")
        .map(|item| {
            let item = item.expect("item reads back");
            let name = names_by_id[item.chrom_id as usize].clone();
            (name, item.start, item.end, item.value.to_bits())
        })
        .collect();
    assert_eq!(actual, expected, "items read back differ from the accepted pushes");

    for (i, &length) in lengths.iter().enumerate() {
        let name = chrom_name(i);
        let items = &model.items[i];

        let by_chrom: Vec<_> = reader
            .chrom_items(&name)
            .expect("chrom_items")
            .map(|item| item.map(|d| (d.start, d.end, d.value.to_bits())))
            .collect::<Result<_, _>>()
            .expect("chrom items read");
        assert_eq!(&by_chrom, items, "chrom_items of {name}");

        if length <= MAX_VALUES_CHECK_LENGTH {
            let mut expected_values = vec![f32::NAN; length as usize];
            for &(start, end, bits) in items {
                expected_values[start as usize..end as usize].fill(f32::from_bits(bits));
            }
            let values = reader.values(&name, 0, length).expect("values");
            let same = values.len() == expected_values.len()
                && values
                    .iter()
                    .zip(&expected_values)
                    .all(|(a, b)| a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()));
            assert!(same, "values of {name} differ from the model");
        }

        let summary = reader.summary(&name, 0, length, SummarySource::Raw).expect("raw summary");
        if items.is_empty() {
            assert!(summary.is_none(), "summary of {name} with no data");
        } else {
            let summary = summary.expect("summary of {name}");
            let covered: u64 = items.iter().map(|&(s, e, _)| u64::from(e - s)).sum();
            let min = items.iter().map(|i| f64::from(f32::from_bits(i.2))).fold(f64::MAX, f64::min);
            let max = items.iter().map(|i| f64::from(f32::from_bits(i.2))).fold(f64::MIN, f64::max);
            assert_eq!(summary.covered_bases, covered as f64, "covered bases of {name}");
            assert_eq!(summary.min, min, "min of {name}");
            assert_eq!(summary.max, max, "max of {name}");
        }
        reader.summary(&name, 0, length, SummarySource::BestZoom).expect("zoom summary");
        for level in reader.zoom_levels().to_vec() {
            for record in reader.zoom_records(level, &name, 0, length).expect("zoom_records") {
                record.expect("zoom record reads");
            }
        }
    }
});
