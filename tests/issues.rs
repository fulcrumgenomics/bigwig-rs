//! Regression tests derived from bugs reported against other wig and bigWig implementations.
//! Each test is named after the behaviour it asserts and cites the issues it guards against. The
//! unit-test modules under `src/` hold the cases that need hand-built files or crate internals;
//! the tests here round-trip data through the public writer and reader.
//!
//! Files are written to memory. Where a test checks the on-disk layout, the helpers below decode
//! the few structures it needs (header offsets, the chromosome B+ tree, R-tree nodes, uncompressed
//! sections and zoom blocks) straight from the bytes, which the writer emits little-endian.

// Values read back must equal the stored f32 values bit for bit.
#![allow(clippy::float_cmp)]

use std::io::Cursor;

use bigwig::{
    BigWigReader, BigWigWriter, ChromIdOrder, ChromInfo, Compression, DataItem, Summary,
    SummarySource, WriterOptions, ZoomLevel, ZoomPolicy, ZoomRecord,
};

type MemoryWriter = BigWigWriter<Cursor<Vec<u8>>>;
type MemoryReader = BigWigReader<Cursor<Vec<u8>>>;

/// Size of the fixed header at the start of every file; the zoom headers follow it.
const FIXED_HEADER_SIZE: u64 = 64;
/// Size of one zoom header.
const ZOOM_HEADER_SIZE: u64 = 24;
/// Size of the chromosome B+ tree header; the root node follows it.
const CHROM_TREE_HEADER_SIZE: u64 = 32;
/// Size of an R-tree header; the root node follows it.
const RTREE_HEADER_SIZE: u64 = 48;
/// Size of a data section header; the items follow it.
const SECTION_HEADER_SIZE: u64 = 24;
/// Size of one bedGraph item in a data section.
const BED_GRAPH_ITEM_SIZE: u64 = 12;
/// Size of one zoom record in a zoom block.
const ZOOM_RECORD_SIZE: u64 = 32;
/// Section `type` of a bedGraph section.
const BED_GRAPH_SECTION: u8 = 1;
/// Section `type` of a fixedStep section.
const FIXED_STEP_SECTION: u8 = 3;
/// Most items a chromosome B+ tree node may hold for readers that read the count as a signed
/// 16-bit integer (IGV, per pyBigWig #36).
const MAX_SIGNED_NODE_COUNT: usize = 32_767;
/// The total summary of a file without data.
const ZERO_SUMMARY: Summary =
    Summary { covered_bases: 0.0, min: 0.0, max: 0.0, sum: 0.0, sum_squares: 0.0 };

// ---- Writing and reading ----------------------------------------------------------------------

/// Writes a file over `chroms` with `options`, adding data with `pushes`, and returns its bytes.
fn write_file(
    options: &WriterOptions,
    chroms: &[ChromInfo],
    pushes: impl FnOnce(&mut MemoryWriter),
) -> Vec<u8> {
    let mut writer = options.create_with(Cursor::new(Vec::new()), chroms).unwrap();
    pushes(&mut writer);
    writer.finish().unwrap().into_inner()
}

/// Opens `bytes` as a bigWig.
fn open(bytes: Vec<u8>) -> MemoryReader {
    BigWigReader::new(Cursor::new(bytes)).unwrap()
}

/// Collects an item or record iterator, failing the test on any error.
fn collect<T>(iter: impl Iterator<Item = bigwig::Result<T>>) -> Vec<T> {
    iter.map(Result::unwrap).collect()
}

/// The `(start, end, value)` of each item.
fn extents(items: &[DataItem]) -> Vec<(u32, u32, f32)> {
    items.iter().map(|item| (item.start, item.end, item.value)).collect()
}

/// The reductions of the file's zoom levels, in header order.
fn reductions(reader: &MemoryReader) -> Vec<u32> {
    reader.zoom_levels().iter().map(ZoomLevel::reduction).collect()
}

/// Asserts that `actual` holds exactly `expected`, where a NaN in `expected` (no data) must be the
/// NaN the reader fills gaps with.
fn assert_values(actual: &[f32], expected: &[f32]) {
    let bits = |values: &[f32]| values.iter().map(|value| value.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(actual), bits(expected), "{actual:?} != {expected:?}");
}

// ---- Decoding the on-disk layout --------------------------------------------------------------

/// The `N` bytes at `offset`.
fn field<const N: usize>(bytes: &[u8], offset: u64) -> [u8; N] {
    let start = usize::try_from(offset).unwrap();
    bytes[start..start + N].try_into().unwrap()
}

/// The little-endian u16 at `offset`.
fn u16_at(bytes: &[u8], offset: u64) -> u16 {
    u16::from_le_bytes(field(bytes, offset))
}

/// The little-endian u32 at `offset`.
fn u32_at(bytes: &[u8], offset: u64) -> u32 {
    u32::from_le_bytes(field(bytes, offset))
}

/// The little-endian u64 at `offset`.
fn u64_at(bytes: &[u8], offset: u64) -> u64 {
    u64::from_le_bytes(field(bytes, offset))
}

/// Overwrites the little-endian u32 at `offset` with `value`.
fn put_u32_at(bytes: &mut [u8], offset: u64, value: u32) {
    let start = usize::try_from(offset).unwrap();
    bytes[start..start + 4].copy_from_slice(&value.to_le_bytes());
}

/// `fullIndexOffset`, from the fixed header.
fn data_index_offset(bytes: &[u8]) -> u64 {
    u64_at(bytes, 24)
}

/// `indexOffset` of zoom level `level`, from its zoom header.
fn zoom_index_offset(bytes: &[u8], level: u64) -> u64 {
    u64_at(bytes, FIXED_HEADER_SIZE + level * ZOOM_HEADER_SIZE + 16)
}

/// The keys of the chromosome B+ tree, one list per leaf, in tree order.
fn chrom_tree_leaf_keys(bytes: &[u8]) -> Vec<Vec<String>> {
    let tree = u64_at(bytes, 8);
    let key_size = u64::from(u32_at(bytes, tree + 8));
    let mut leaves = Vec::new();
    collect_chrom_tree_leaves(bytes, tree + CHROM_TREE_HEADER_SIZE, key_size, &mut leaves);
    leaves
}

/// Appends the keys of every leaf under the B+ tree node at `node` to `leaves`. Leaf slots hold a
/// key then `chromId` and `chromSize`; non-leaf slots a key then a child offset: 8 bytes either
/// way.
fn collect_chrom_tree_leaves(
    bytes: &[u8],
    node: u64,
    key_size: u64,
    leaves: &mut Vec<Vec<String>>,
) {
    let is_leaf = bytes[usize::try_from(node).unwrap()] == 1;
    let slots =
        (0..u64::from(u16_at(bytes, node + 2))).map(|index| node + 4 + index * (key_size + 8));
    if is_leaf {
        let key = |slot: u64| {
            let start = usize::try_from(slot).unwrap();
            let raw = &bytes[start..start + usize::try_from(key_size).unwrap()];
            String::from_utf8(raw.iter().copied().take_while(|&byte| byte != 0).collect()).unwrap()
        };
        leaves.push(slots.map(key).collect());
    } else {
        for slot in slots {
            collect_chrom_tree_leaves(bytes, u64_at(bytes, slot + key_size), key_size, leaves);
        }
    }
}

/// An R-tree slot's bounds, from `(startChromIx, startBase)` to `(endChromIx, endBase)`, compared
/// as `(chromId, base)` pairs as kent compares them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Bounds {
    start: (u32, u32),
    end: (u32, u32),
}

impl Bounds {
    /// The bounds of an item or zoom record covering `[start, end)` on chromosome `chrom_id`.
    fn of_extent(chrom_id: u32, start: u32, end: u32) -> Self {
        Self { start: (chrom_id, start), end: (chrom_id, end) }
    }

    /// The four bounds fields at the start of the slot at `slot`.
    fn decode(bytes: &[u8], slot: u64) -> Self {
        Self {
            start: (u32_at(bytes, slot), u32_at(bytes, slot + 4)),
            end: (u32_at(bytes, slot + 8), u32_at(bytes, slot + 12)),
        }
    }

    /// Whether `other` lies within these bounds.
    fn contains(self, other: Self) -> bool {
        self.start <= other.start && other.end <= self.end
    }

    /// The smallest bounds containing both.
    fn union(self, other: Self) -> Self {
        Self { start: self.start.min(other.start), end: self.end.max(other.end) }
    }
}

/// One R-tree leaf entry: the bounds of one block and where the block lies.
#[derive(Clone, Copy, Debug)]
struct LeafEntry {
    bounds: Bounds,
    data_offset: u64,
    data_size: u64,
}

impl LeafEntry {
    /// The block this entry points at.
    fn block(self, bytes: &[u8]) -> &[u8] {
        let start = usize::try_from(self.data_offset).unwrap();
        &bytes[start..start + usize::try_from(self.data_size).unwrap()]
    }
}

/// The leaf entries of the R-tree at `index_offset`, in tree order, asserting along the way that
/// every non-leaf slot's bounds enclose all the slots below it.
fn rtree_leaves(bytes: &[u8], index_offset: u64) -> Vec<LeafEntry> {
    let mut leaves = Vec::new();
    walk_rtree_node(bytes, index_offset + RTREE_HEADER_SIZE, &mut leaves);
    leaves
}

/// Appends the leaf entries under the R-tree node at `node` to `leaves`, and returns the union of
/// the node's slot bounds, or `None` for an empty node. Leaf slots are 32 bytes (bounds, offset,
/// size); non-leaf slots 24 (bounds, child offset).
fn walk_rtree_node(bytes: &[u8], node: u64, leaves: &mut Vec<LeafEntry>) -> Option<Bounds> {
    let is_leaf = bytes[usize::try_from(node).unwrap()] == 1;
    let mut node_bounds: Option<Bounds> = None;
    for index in 0..u64::from(u16_at(bytes, node + 2)) {
        let slot = node + 4 + index * if is_leaf { 32 } else { 24 };
        let bounds = Bounds::decode(bytes, slot);
        if is_leaf {
            let (data_offset, data_size) = (u64_at(bytes, slot + 16), u64_at(bytes, slot + 24));
            leaves.push(LeafEntry { bounds, data_offset, data_size });
        } else {
            let below = walk_rtree_node(bytes, u64_at(bytes, slot + 16), leaves).unwrap();
            assert!(bounds.contains(below), "node bounds {bounds:?} do not enclose {below:?}");
        }
        node_bounds = Some(node_bounds.map_or(bounds, |union| union.union(bounds)));
    }
    node_bounds
}

/// The header of an uncompressed data section.
#[derive(Clone, Copy, Debug)]
struct SectionHeader {
    chrom_id: u32,
    start: u32,
    end: u32,
    step: u32,
    span: u32,
    kind: u8,
    item_count: u16,
}

impl SectionHeader {
    /// Decodes the header at the start of the uncompressed section `block`.
    fn decode(block: &[u8]) -> Self {
        Self {
            chrom_id: u32_at(block, 0),
            start: u32_at(block, 4),
            end: u32_at(block, 8),
            step: u32_at(block, 12),
            span: u32_at(block, 16),
            kind: block[20],
            item_count: u16_at(block, 22),
        }
    }
}

/// The bounds of every item of the uncompressed bedGraph section `block`.
fn bed_graph_item_bounds(block: &[u8]) -> Vec<Bounds> {
    let header = SectionHeader::decode(block);
    assert_eq!(header.kind, BED_GRAPH_SECTION);
    (0..u64::from(header.item_count))
        .map(|index| SECTION_HEADER_SIZE + index * BED_GRAPH_ITEM_SIZE)
        .map(|item| {
            Bounds::of_extent(header.chrom_id, u32_at(block, item), u32_at(block, item + 4))
        })
        .collect()
}

/// The bounds of every record of the uncompressed zoom block `block`.
fn zoom_record_bounds(block: &[u8]) -> Vec<Bounds> {
    (0..u64::try_from(block.len()).unwrap() / ZOOM_RECORD_SIZE)
        .map(|index| index * ZOOM_RECORD_SIZE)
        .map(|record| {
            Bounds::of_extent(
                u32_at(block, record),
                u32_at(block, record + 4),
                u32_at(block, record + 8),
            )
        })
        .collect()
}

// ---- Empty files and zoom levels --------------------------------------------------------------

/// pyBigWig #38 (<https://github.com/deeptools/pyBigWig/issues/38>), libBigWig PR #53
/// (<https://github.com/dpryan79/libBigWig/pull/53>) and bigtools #67
/// (<https://github.com/jackh726/bigtools/issues/67>): writers could not produce a file without
/// data, or hung trying, and readers rejected the empty files other writers produced. Those
/// other writers' variants are hand-built in `reader.rs`.
#[test]
fn empty_file_written_by_us_is_valid_and_readable() {
    let chroms = [ChromInfo::new("chr1", 1_000), ChromInfo::new("chr2", 500)];
    let mut reader = open(write_file(&WriterOptions::new(), &chroms, |_| {}));
    assert_eq!(reader.chroms(), chroms);
    assert_eq!(reader.info().section_count, 0);
    assert_eq!(reader.info().total_summary, Some(ZERO_SUMMARY));
    assert_eq!(collect(reader.all_items().unwrap()), []);
    assert_eq!(collect(reader.items("chr1", 0, 1_000).unwrap()), []);
    assert!(reader.values("chr2", 0, 500).unwrap().iter().all(|value| value.is_nan()));
    assert_eq!(reader.summary("chr1", 0, 1_000, SummarySource::Raw).unwrap(), None);
    assert_eq!(reader.summary("chr1", 0, 1_000, SummarySource::BestZoom).unwrap(), None);
    assert_eq!(reader.compute_total_summary().unwrap(), ZERO_SUMMARY);
    let level = *reader.zoom_levels().first().expect("an empty file still has a zoom level");
    assert_eq!(collect(reader.zoom_records(level, "chr1", 0, 1_000).unwrap()), []);
}

/// libBigWig #26 (<https://github.com/dpryan79/libBigWig/issues/26>): a contig holding a single
/// whole-contig interval crashed libBigWig's writer, and the interim fix wrote no zoom levels,
/// which IGV cannot open. Reading files without zoom levels is covered by
/// `file_with_zero_zoom_levels_summarises_from_raw_data` in `reader.rs`.
#[test]
fn single_whole_chromosome_interval_still_gets_one_zoom_level() {
    let length = 248_956_422;
    let chroms = [ChromInfo::new("chr1", length)];
    let mut reader = open(write_file(&WriterOptions::new(), &chroms, |writer| {
        writer.push_interval("chr1", 0, length, 0.0).unwrap();
    }));
    // The default ladder's first reduction is ten times the pitch, capped at the chromosome.
    assert_eq!(reductions(&reader), [length]);
    let level = reader.zoom_levels()[0];
    let records = collect(reader.zoom_records(level, "chr1", 0, length).unwrap());
    let summary = Summary {
        covered_bases: f64::from(length),
        min: 0.0,
        max: 0.0,
        sum: 0.0,
        sum_squares: 0.0,
    };
    assert_eq!(records, [ZoomRecord { chrom_id: 0, start: 0, end: length, summary }]);
}

/// pyBigWig #89 (<https://github.com/deeptools/pyBigWig/issues/89>): IGV could not read
/// pyBigWig files that had no zoom levels.
#[test]
fn file_has_at_least_one_zoom_level_by_default_for_tiny_data() {
    // One base of data on a chromosome shorter than the smallest default reduction.
    let chroms = [ChromInfo::new("chrM", 10)];
    let reader = open(write_file(&WriterOptions::new(), &chroms, |writer| {
        writer.push_interval("chrM", 3, 4, 1.0).unwrap();
    }));
    assert_eq!(reader.zoom_levels().len(), 1);
}

/// pyBigWig #89 (<https://github.com/deeptools/pyBigWig/issues/89>): IGV could not read
/// pyBigWig files that had no zoom levels.
#[test]
fn file_has_at_least_one_zoom_level_by_default_for_large_data() {
    let length = 1_000_000;
    let depths: Vec<f32> =
        (0..length).map(|base| f32::from(u8::try_from(base % 50).unwrap())).collect();
    let reader =
        open(write_file(&WriterOptions::new(), &[ChromInfo::new("chr1", length)], |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &depths).unwrap();
        }));
    let reductions = reductions(&reader);
    assert!(reductions.len() > 1, "{reductions:?}");
    assert!(reductions.windows(2).all(|pair| pair[0] < pair[1]), "{reductions:?}");
}

/// bigtools #3 (<https://github.com/jackh726/bigtools/pull/3>): a reader decoded zoom blocks
/// with the wrong record size. Writes `count` one-base items ten bases apart, so the single zoom
/// level (reduction 10) holds one record per item, packed 1,024 records to a block, and asserts
/// that the reader returns exactly those records.
fn assert_zoom_level_reads_back_every_record(count: u16) {
    let length = u32::from(count) * 10 + 10;
    let items: Vec<(u32, u32, f32)> = (0..count)
        .map(|index| (u32::from(index) * 10, u32::from(index) * 10 + 1, f32::from(index)))
        .collect();
    let options =
        WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![10])).zoom_items_per_block(1_024);
    let mut reader = open(write_file(&options, &[ChromInfo::new("chr1", length)], |writer| {
        writer.push_intervals("chr1", &items).unwrap();
    }));
    let level = reader.zoom_levels()[0];
    assert_eq!(level.reduction(), 10);
    let expected: Vec<ZoomRecord> = items
        .iter()
        .map(|&(start, end, value)| {
            let value = f64::from(value);
            let summary = Summary {
                covered_bases: 1.0,
                min: value,
                max: value,
                sum: value,
                sum_squares: value * value,
            };
            ZoomRecord { chrom_id: 0, start, end, summary }
        })
        .collect();
    assert_eq!(collect(reader.zoom_records(level, "chr1", 0, length).unwrap()), expected);
}

/// bigtools #3 (<https://github.com/jackh726/bigtools/pull/3>); see
/// `assert_zoom_level_reads_back_every_record`.
#[test]
fn zoom_level_of_one_record_reads_back_exactly_one() {
    assert_zoom_level_reads_back_every_record(1);
}

/// bigtools #3 (<https://github.com/jackh726/bigtools/pull/3>); see
/// `assert_zoom_level_reads_back_every_record`.
#[test]
fn zoom_block_one_record_short_of_full_reads_back_every_record() {
    assert_zoom_level_reads_back_every_record(1_023);
}

/// bigtools #3 (<https://github.com/jackh726/bigtools/pull/3>); see
/// `assert_zoom_level_reads_back_every_record`.
#[test]
fn full_zoom_block_reads_back_exactly_its_records() {
    assert_zoom_level_reads_back_every_record(1_024);
}

/// bigtools #3 (<https://github.com/jackh726/bigtools/pull/3>); see
/// `assert_zoom_level_reads_back_every_record`.
#[test]
fn zoom_level_spilling_one_record_into_a_second_block_reads_back_every_record() {
    assert_zoom_level_reads_back_every_record(1_025);
}

/// bigtools PR #103 (<https://github.com/jackh726/bigtools/pull/103>): moving the zoom cursor
/// back into a gap wrote a phantom record per empty bin, which browsers drew as signal.
/// `gap_longer_than_a_bin_produces_no_phantom_records` in `zoom.rs` checks the records as built;
/// this checks what a reader sees.
#[test]
fn zoom_records_never_span_gaps_longer_than_a_bin() {
    let options = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![1_000, 100_000]));
    let mut reader = open(write_file(&options, &[ChromInfo::new("chr1", 1_000_000)], |writer| {
        writer.push_fixed_step("chr1", 10_000, 1, 1, &[5.0; 3_000]).unwrap();
        writer.push_fixed_step("chr1", 900_000, 1, 1, &[7.0; 3_000]).unwrap();
    }));
    assert_eq!(reductions(&reader), [1_000, 100_000]);
    for level in reader.zoom_levels().to_vec() {
        let in_gap = collect(reader.zoom_records(level, "chr1", 13_000, 900_000).unwrap());
        assert_eq!(in_gap, [], "level {}", level.reduction());
    }
    // Bins of 100 kb read the 1 kb level, and one 600 kb bin reads the 100 kb level.
    assert_eq!(reader.best_zoom_for(100_000).map(|level| level.reduction()), Some(1_000));
    assert_eq!(reader.best_zoom_for(600_000).map(|level| level.reduction()), Some(100_000));
    for n_bins in [6, 1] {
        let bins = reader.summaries("chr1", 200_000, 800_000, n_bins, SummarySource::BestZoom);
        assert!(bins.unwrap().iter().all(Option::is_none), "{n_bins} bins");
    }
}

// ---- Chromosomes ------------------------------------------------------------------------------

/// pyBigWig #36 (<https://github.com/deeptools/pyBigWig/issues/36>) and IGV #294
/// (<https://github.com/igvteam/igv/issues/294>): libBigWig could not write tens of thousands of
/// contigs (it never built a non-leaf node, and its fix had to cap leaves at 32,767 for IGV), and
/// IGV lost scaffolds of a 138,607-scaffold genome.
#[test]
fn chrom_tree_above_leaf_capacity_round_trips_with_multi_level_tree() {
    let names: Vec<String> = (0..40_000_u16).map(|index| format!("scaffold_{index}")).collect();
    let chroms: Vec<ChromInfo> =
        names.iter().map(|name| ChromInfo::new(name.as_str(), 100)).collect();
    // Raw blocks: compressing 80,000 tiny blocks is slow in unoptimised builds, and irrelevant.
    let options =
        WriterOptions::new().compression(Compression::None).zoom(ZoomPolicy::Explicit(vec![10]));
    let bytes = write_file(&options, &chroms, |writer| {
        for (index, name) in (0_u16..).zip(&names) {
            writer.push_interval(name, 0, 1, f32::from(index)).unwrap();
        }
    });

    let leaves = chrom_tree_leaf_keys(&bytes);
    assert!(leaves.len() > 1);
    assert!(leaves.iter().all(|keys| keys.len() <= MAX_SIGNED_NODE_COUNT));
    assert_eq!(leaves.iter().map(Vec::len).sum::<usize>(), names.len());

    let mut reader = open(bytes);
    assert_eq!(reader.chroms().len(), names.len());
    for name in &names {
        assert_eq!(reader.chrom(name).map(|chrom| chrom.length), Some(100), "{name}");
    }
    let items = collect(reader.all_items().unwrap());
    assert_eq!(items.len(), names.len());
    for (index, item) in (0_u16..).zip(&items) {
        assert_eq!(reader.chrom_name(item.chrom_id), Some(names[usize::from(index)].as_str()));
        assert_eq!((item.start, item.end, item.value), (0, 1, f32::from(index)));
    }
}

/// pyBigWig #14 (<https://github.com/deeptools/pyBigWig/issues/14>): libBigWig crashed opening
/// Drosophila r6 bigWigs, whose sequence names range from two characters to dozens.
#[test]
fn reader_opens_files_with_short_and_long_scaffold_names_mixed() {
    let long_name = |suffix: &str| format!("Unmapped_Scaffold_{suffix:0>22}");
    let names = [
        "2L".to_string(),
        "X".to_string(),
        "4".to_string(),
        "chr1".to_string(),
        "chr10".to_string(),
        "chr1_alt".to_string(),
        "211000022278279".to_string(),
        "mitochondrion_genome".to_string(),
        long_name("1"),
        long_name("2"),
    ];
    assert_eq!(names[8].len(), 40);
    let chroms: Vec<ChromInfo> = (0_u16..)
        .zip(&names)
        .map(|(index, name)| ChromInfo::new(name.as_str(), 1_000 + u32::from(index)))
        .collect();
    let bytes = write_file(&WriterOptions::new(), &chroms, |writer| {
        for (index, name) in (0_u16..).zip(&names) {
            let start = u32::from(index);
            writer.push_interval(name, start, start + 10, f32::from(index)).unwrap();
        }
    });

    let mut reader = open(bytes);
    for (index, chrom) in (0_u16..).zip(&chroms) {
        let start = u32::from(index);
        assert_eq!(reader.chrom(&chrom.name), Some(chrom));
        assert_values(
            &reader.values(&chrom.name, start, start + 10).unwrap(),
            &[f32::from(index); 10],
        );
    }
    assert_eq!(collect(reader.all_items().unwrap()).len(), names.len());
}

/// bbi-js #87 (<https://github.com/GMOD/bbi-js/pull/87>): name lookups compared keys with
/// `localeCompare`, which orders mixed case and underscores differently from the bytewise order
/// kent sorts keys in, so some names were never found. Our writer's tree is a single leaf below
/// 32,768 names; the multi-level trees above that are searched kent's way in `bptree.rs`.
#[test]
fn chrom_lookup_uses_byte_order_for_mixed_case_and_underscore_names() {
    let names = ["chr1", "chrX", "chrx", "chr_1", "chr10", "YWHAH", "ywhah", "Metazoa_SRP"];
    let chroms: Vec<ChromInfo> = (0_u16..)
        .zip(names)
        .map(|(index, name)| ChromInfo::new(name, 100 + u32::from(index)))
        .collect();
    let bytes = write_file(&WriterOptions::new(), &chroms, |writer| {
        for (index, name) in (0_u16..).zip(names) {
            writer.push_interval(name, 0, 10, f32::from(index)).unwrap();
        }
    });

    // Bytewise: upper case before `_` before lower case, and `1` before `10`.
    let byte_order = ["Metazoa_SRP", "YWHAH", "chr1", "chr10", "chrX", "chr_1", "chrx", "ywhah"];
    assert_eq!(chrom_tree_leaf_keys(&bytes), [byte_order]);

    let mut reader = open(bytes);
    for (index, chrom) in (0_u16..).zip(&chroms) {
        assert_eq!(reader.chrom(&chrom.name), Some(chrom));
        assert_values(&reader.values(&chrom.name, 0, 10).unwrap(), &[f32::from(index); 10]);
    }
}

/// pyBigWig #68 (<https://github.com/deeptools/pyBigWig/issues/68>), `deepTools` #761
/// (<https://github.com/deeptools/deepTools/issues/761>), kent #38
/// (<https://github.com/ucscGenomeBrowser/kent/issues/38>) and kent #44
/// (<https://github.com/ucscGenomeBrowser/kent/issues/44>): tools broke on chromosomes longer than
/// 2³¹ bases, which the format's unsigned 32-bit coordinates allow.
#[test]
fn chromosome_longer_than_2_pow_31_round_trips() {
    let length = u32::MAX;
    let two_pow_31 = 1_u32 << 31;
    let items = [
        (two_pow_31 - 5, two_pow_31 + 5, 1.0),
        (3_000_000_000, 3_000_000_010, 2.0),
        (length - 11, length - 1, 3.0),
        (length - 1, length, 4.0),
    ];
    let options = WriterOptions::new().zoom(ZoomPolicy::Explicit(vec![1 << 16]));
    let mut reader = open(write_file(&options, &[ChromInfo::new("huge", length)], |writer| {
        writer.push_intervals("huge", &items).unwrap();
    }));
    assert_eq!(reader.chrom("huge").map(|chrom| chrom.length), Some(length));
    assert_eq!(extents(&collect(reader.chrom_items("huge").unwrap())), items);

    let nan = f32::NAN;
    let around_two_pow_31 = reader.values("huge", two_pow_31 - 6, two_pow_31 + 6).unwrap();
    assert_values(
        &around_two_pow_31,
        &[nan, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, nan],
    );
    let at_the_end = reader.values("huge", length - 13, length).unwrap();
    assert_values(&at_the_end, &[nan, nan, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0]);

    // The 64 kb level's grid splits the first item at 2³¹.
    let level = reader.zoom_levels()[0];
    let records =
        collect(reader.zoom_records(level, "huge", two_pow_31 - 10, two_pow_31 + 10).unwrap());
    let record_extents: Vec<(u32, u32)> =
        records.iter().map(|record| (record.start, record.end)).collect();
    assert_eq!(record_extents, [(two_pow_31 - 5, two_pow_31), (two_pow_31, two_pow_31 + 5)]);

    // Over the whole chromosome, the raw summary and the zoom level agree exactly.
    let whole = Summary { covered_bases: 31.0, min: 1.0, max: 4.0, sum: 64.0, sum_squares: 156.0 };
    assert_eq!(reader.summary("huge", 0, length, SummarySource::Raw).unwrap(), Some(whole));
    assert_eq!(reader.summary("huge", 0, length, SummarySource::BestZoom).unwrap(), Some(whole));
}

/// rtracklayer #132 (<https://github.com/lawremi/rtracklayer/issues/132>): `export.bw` wrote
/// only the last chromosome on some platforms, and reading the file back crashed.
#[test]
fn multi_chromosome_write_preserves_every_chromosome() {
    let chroms: Vec<ChromInfo> =
        (1..=5).map(|number| ChromInfo::new(format!("chr{number}"), 10_000)).collect();
    // Run-length-encoded coverage, different on each chromosome.
    let runs = |offset: u16| {
        let value = f32::from(offset);
        vec![
            (0, 100, value),
            (100, 250, value + 0.5),
            (400, 1_000 + u32::from(offset), value + 1.0),
        ]
    };
    for write_order in [vec![0, 1, 2, 3, 4], vec![4, 2, 0, 3, 1]] {
        let bytes = write_file(&WriterOptions::new(), &chroms, |writer| {
            for &index in &write_order {
                writer
                    .push_intervals(&chroms[index].name, &runs(u16::try_from(index).unwrap()))
                    .unwrap();
            }
        });
        let mut reader = open(bytes);
        for (index, chrom) in chroms.iter().enumerate() {
            let items = collect(reader.chrom_items(&chrom.name).unwrap());
            assert_eq!(extents(&items), runs(u16::try_from(index).unwrap()), "{write_order:?}");
        }
        assert_eq!(collect(reader.all_items().unwrap()).len(), 15, "{write_order:?}");
    }
}

// ---- Data layout ------------------------------------------------------------------------------

/// libBigWig PR #70 (<https://github.com/dpryan79/libBigWig/pull/70>) and pyBigWig #166
/// (<https://github.com/deeptools/pyBigWig/issues/166>): libBigWig set the `end` of the last
/// fixedStep section of each `addEntries` call six steps too far. Our writer must use kent's
/// formula for every section, however the run is split across calls and sections, and our reader
/// must give the same items when `end` is too large.
#[test]
fn fixed_step_section_end_equals_start_plus_n_minus_one_steps_plus_span() {
    let (step, span) = (5, 3);
    let options = WriterOptions::new()
        .compression(Compression::None)
        .items_per_section(64)
        .zoom(ZoomPolicy::Explicit(vec![1_000]));
    let mut bytes = write_file(&options, &[ChromInfo::new("chr1", 100_000)], |writer| {
        writer.push_fixed_step("chr1", 0, step, span, &[1.0]).unwrap();
        writer.push_fixed_step("chr1", 1_000, step, span, &[2.0; 7]).unwrap();
        // One run of 1,000 values pushed in three calls, each continuing where the last ended.
        writer.push_fixed_step("chr1", 10_000, step, span, &[3.0; 300]).unwrap();
        writer.push_fixed_step("chr1", 11_500, step, span, &[4.0; 300]).unwrap();
        writer.push_fixed_step("chr1", 13_000, step, span, &[5.0; 400]).unwrap();
    });

    let sections: Vec<(LeafEntry, SectionHeader)> = rtree_leaves(&bytes, data_index_offset(&bytes))
        .into_iter()
        .map(|leaf| (leaf, SectionHeader::decode(leaf.block(&bytes))))
        .collect();
    let counts: Vec<u16> = sections.iter().map(|(_, header)| header.item_count).collect();
    let mut expected_counts = vec![1, 7];
    expected_counts.extend([64; 15]);
    expected_counts.push(40);
    assert_eq!(counts, expected_counts);
    for (_, header) in &sections {
        assert_eq!((header.kind, header.step, header.span), (FIXED_STEP_SECTION, step, span));
        assert_eq!(header.end, header.start + (u32::from(header.item_count) - 1) * step + span);
    }

    let written = collect(open(bytes.clone()).all_items().unwrap());
    for (leaf, header) in &sections {
        put_u32_at(&mut bytes, leaf.data_offset + 8, header.end + 6 * step);
    }
    assert_eq!(collect(open(bytes).all_items().unwrap()), written);
}

/// gonetics #2 (<https://github.com/pbenner/gonetics/issues/2>): its writer built a wrong R-tree
/// when data had gaps or chromosomes without data. Every leaf's bounds must contain every item
/// (or zoom record) in its block, every node's bounds its children's, and each item must be found
/// by a query on its own extent.
#[test]
fn rtree_bounds_cover_every_item_when_data_has_large_gaps() {
    let chroms = [
        ChromInfo::new("chrA", 200_000_000),
        ChromInfo::new("chrB", 1_000_000),
        ChromInfo::new("chrC", 50_000_000),
    ];
    // chrA has data at both ends with a 150 Mb gap between; chrB, between chrA and chrC in
    // chromId order, has none.
    let near =
        |base: u32| (0..20).map(move |index| (base + index * 50, base + index * 50 + 10, 1.0));
    let chr_a: Vec<(u32, u32, f32)> = near(0).chain(near(150_000_000)).collect();
    let chr_c: Vec<(u32, u32, f32)> =
        (0..10).map(|index| (1_000 + index * 1_000, 1_010 + index * 1_000, 2.0)).collect();
    let options = WriterOptions::new()
        .chrom_id_order(ChromIdOrder::Listed)
        .compression(Compression::None)
        .items_per_section(4)
        .index_block_size(2)
        .zoom(ZoomPolicy::Explicit(vec![100]))
        .zoom_items_per_block(3);
    let bytes = write_file(&options, &chroms, |writer| {
        writer.push_intervals("chrA", &chr_a).unwrap();
        writer.push_intervals("chrC", &chr_c).unwrap();
    });

    // Ten sections on chrA and three on chrC: several levels of fan-out 2.
    let data_leaves = rtree_leaves(&bytes, data_index_offset(&bytes));
    assert_eq!(data_leaves.len(), 13);
    for leaf in data_leaves {
        for item in bed_graph_item_bounds(leaf.block(&bytes)) {
            assert!(leaf.bounds.contains(item), "{leaf:?} does not contain item {item:?}");
        }
    }
    for leaf in rtree_leaves(&bytes, zoom_index_offset(&bytes, 0)) {
        for record in zoom_record_bounds(leaf.block(&bytes)) {
            assert!(leaf.bounds.contains(record), "{leaf:?} does not contain record {record:?}");
        }
    }

    let mut reader = open(bytes);
    for (chrom, items) in [("chrA", &chr_a), ("chrC", &chr_c)] {
        for &(start, end, value) in items {
            let found = collect(reader.items(chrom, start, end).unwrap());
            assert_eq!(extents(&found), [(start, end, value)], "{chrom}:{start}-{end}");
        }
    }
    assert_eq!(collect(reader.items("chrA", 1_000, 150_000_000).unwrap()), []);
    assert_eq!(collect(reader.chrom_items("chrB").unwrap()), []);
}

// ---- Values -----------------------------------------------------------------------------------

/// pyBigWig #132 (<https://github.com/deeptools/pyBigWig/issues/132>): callers holding dense
/// arrays that mark missing data with NaN wrote the NaN into their files, where readers disagree
/// on it and it poisons zoom and total summaries. `push_fixed_step` still rejects NaN, and
/// `push_fixed_step_skip` writes it as a gap.
#[test]
fn dense_values_with_nan_runs_can_be_written_as_gaps_by_a_helper() {
    let dense = [1.0, f32::NAN, f32::NAN, 3.0];
    let mut writer = WriterOptions::new()
        .create_with(Cursor::new(Vec::new()), &[ChromInfo::new("chr1", 10)])
        .unwrap();
    let error = writer.push_fixed_step("chr1", 0, 1, 1, &dense).unwrap_err();
    assert!(matches!(error, bigwig::Error::NonFiniteValue { position: 1, .. }), "{error:?}");
    writer.push_fixed_step_skip("chr1", 0, 1, 1, &dense).unwrap();

    let mut reader = open(writer.finish().unwrap().into_inner());
    let items = collect(reader.chrom_items("chr1").unwrap());
    assert_eq!(extents(&items), [(0, 1, 1.0), (3, 4, 3.0)]);
    assert_values(&reader.values("chr1", 0, 4).unwrap(), &dense);
    let summary = Summary { covered_bases: 2.0, min: 1.0, max: 3.0, sum: 4.0, sum_squares: 10.0 };
    assert_eq!(reader.info().total_summary, Some(summary));
}
