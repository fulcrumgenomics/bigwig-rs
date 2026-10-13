//! Criterion benchmarks on synthetic per-base depth from a seeded generator, so no data is checked
//! in. The first measured numbers on a quiet machine are the regression baselines.
//!
//! Write sizes default to 1 Mbp and 10 Mbp; set `BIGWIG_BENCH_BASES` to a comma-separated list
//! (e.g. `BIGWIG_BENCH_BASES=100000000`) to benchmark other sizes. Everything is written to
//! memory, so the benches measure the library and not the disk.

use std::hint::black_box;
use std::io::Cursor;

use bigwig::{BigWigReader, ChromInfo, SummarySource, WriterOptions};
use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const CHROM: &str = "chr1";
const DEFAULT_WRITE_SIZES: [u32; 2] = [1_000_000, 10_000_000];
/// Chromosome length of the file the read and summary benches query.
const QUERY_FILE_BASES: u32 = 10_000_000;
const QUERY_REGION_BASES: u32 = 10_000;
const SUMMARY_REGION_BASES: u32 = 1_000_000;
const SUMMARY_BINS: u32 = 100;
const REGION_POOL_SIZE: usize = 1_024;
const SEED: u64 = 0x5EED_B16B;

/// A `SplitMix64` generator: tiny, seedable and identical on every platform.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform value in `0..bound`.
    fn below(&mut self, bound: u32) -> u32 {
        u32::try_from(self.next_u64() % u64::from(bound)).expect("remainder is below a u32 bound")
    }
}

/// Integer-valued depth that drifts slowly and jitters per base, like read depth: the mean moves
/// to a new level every few hundred bases and each base is the mean plus up to ±3.
fn synthetic_depth(bases: u32, rng: &mut Rng) -> Vec<f32> {
    let mut depth = Vec::with_capacity(bases as usize);
    let mut mean = 30_u32;
    let mut until_shift = 0;
    for _ in 0..bases {
        if until_shift == 0 {
            mean = (mean + rng.below(9)).saturating_sub(4).max(4);
            until_shift = 100 + rng.below(400);
        }
        until_shift -= 1;
        let jittered = (mean + rng.below(7)).saturating_sub(3);
        // Depths are far below 2^24, so the conversion is exact.
        depth.push(f32::from(u16::try_from(jittered).expect("depth fits in u16")));
    }
    depth
}

/// Run-length encodes `depth` into bedGraph intervals, merging equal neighbours.
fn depth_to_intervals(depth: &[f32]) -> Vec<(u32, u32, f32)> {
    let mut intervals: Vec<(u32, u32, f32)> = Vec::new();
    for (position, &value) in (0_u32..).zip(depth) {
        match intervals.last_mut() {
            Some(last) if last.2.to_bits() == value.to_bits() => last.1 = position + 1,
            _ => intervals.push((position, position + 1, value)),
        }
    }
    intervals
}

/// Sparse items like per-CpG methylation: a position every 1 to 100 bases, a value in 0..=100.
fn sparse_items(bases: u32, rng: &mut Rng) -> Vec<(u32, f32)> {
    let mut items = Vec::new();
    let mut position = rng.below(100);
    while position < bases {
        items.push((position, f32::from(u16::try_from(rng.below(101)).expect("percentage"))));
        position += 1 + rng.below(100);
    }
    items
}

fn write_sizes() -> Vec<u32> {
    match std::env::var("BIGWIG_BENCH_BASES") {
        Ok(list) => list
            .split(',')
            .map(|size| size.trim().parse().expect("BIGWIG_BENCH_BASES: comma-separated integers"))
            .collect(),
        Err(_) => DEFAULT_WRITE_SIZES.to_vec(),
    }
}

fn chroms(bases: u32) -> [ChromInfo; 1] {
    [ChromInfo::new(CHROM, bases)]
}

fn element_count(len: usize) -> Throughput {
    Throughput::Elements(u64::try_from(len).expect("length fits in u64"))
}

/// A per-base fixedStep file of `depth` with default options, as bytes.
fn fixed_step_file(depth: &[f32]) -> Vec<u8> {
    let bases = u32::try_from(depth.len()).expect("bench sizes fit in u32");
    let mut writer = WriterOptions::new()
        .create_with(Cursor::new(Vec::new()), &chroms(bases))
        .expect("create writer");
    writer.push_fixed_step(CHROM, 0, 1, 1, depth).expect("push depth");
    writer.finish().expect("finish").into_inner()
}

fn bench_write_fixed_step(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("write_fixed_step");
    group.sample_size(10);
    for bases in write_sizes() {
        let depth = synthetic_depth(bases, &mut Rng(SEED));
        group.throughput(Throughput::Elements(u64::from(bases)));
        group.bench_function(format!("{bases}_bp"), |bencher| {
            bencher.iter(|| black_box(fixed_step_file(&depth)));
        });
    }
    group.finish();
}

fn bench_write_bedgraph(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("write_bedgraph");
    group.sample_size(10);
    let bases = DEFAULT_WRITE_SIZES[0];
    let intervals = depth_to_intervals(&synthetic_depth(bases, &mut Rng(SEED)));
    group.throughput(element_count(intervals.len()));
    group.bench_function(format!("{bases}_bp_runs"), |bencher| {
        bencher.iter(|| {
            let mut writer = WriterOptions::new()
                .create_with(Cursor::new(Vec::new()), &chroms(bases))
                .expect("create writer");
            writer.push_intervals(CHROM, &intervals).expect("push intervals");
            black_box(writer.finish().expect("finish").into_inner())
        });
    });
    group.finish();
}

fn bench_write_sparse_variable_step(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("write_variable_step");
    group.sample_size(10);
    let bases = DEFAULT_WRITE_SIZES[1];
    let items = sparse_items(bases, &mut Rng(SEED));
    group.throughput(element_count(items.len()));
    group.bench_function(format!("sparse_{bases}_bp"), |bencher| {
        bencher.iter(|| {
            let mut writer = WriterOptions::new()
                .create_with(Cursor::new(Vec::new()), &chroms(bases))
                .expect("create writer");
            writer.push_variable_step(CHROM, 1, &items).expect("push items");
            black_box(writer.finish().expect("finish").into_inner())
        });
    });
    group.finish();
}

fn bench_read_values(criterion: &mut Criterion) {
    let mut rng = Rng(SEED);
    let bytes = fixed_step_file(&synthetic_depth(QUERY_FILE_BASES, &mut rng));
    let regions: Vec<u32> =
        (0..REGION_POOL_SIZE).map(|_| rng.below(QUERY_FILE_BASES - QUERY_REGION_BASES)).collect();
    let mut reader = BigWigReader::new(Cursor::new(bytes)).expect("open file");
    let mut out = vec![0.0_f32; QUERY_REGION_BASES as usize];

    let mut group = criterion.benchmark_group("read_values");
    group.throughput(Throughput::Elements(u64::from(QUERY_REGION_BASES)));
    let mut next = 0;
    group.bench_function(format!("random_{QUERY_REGION_BASES}_bp"), |bencher| {
        bencher.iter(|| {
            let start = regions[next % regions.len()];
            next += 1;
            let covered = reader
                .values_into(CHROM, start, start + QUERY_REGION_BASES, &mut out, f32::NAN)
                .expect("read values");
            black_box(covered)
        });
    });
    group.finish();
}

fn bench_summaries(criterion: &mut Criterion) {
    let mut rng = Rng(SEED);
    let bytes = fixed_step_file(&synthetic_depth(QUERY_FILE_BASES, &mut rng));
    let regions: Vec<u32> =
        (0..REGION_POOL_SIZE).map(|_| rng.below(QUERY_FILE_BASES - SUMMARY_REGION_BASES)).collect();
    let mut reader = BigWigReader::new(Cursor::new(bytes)).expect("open file");

    let mut group = criterion.benchmark_group("summaries");
    for (label, source) in [("best_zoom", SummarySource::BestZoom), ("raw", SummarySource::Raw)] {
        let mut next = 0;
        group.bench_function(format!("{label}_{SUMMARY_BINS}_bins"), |bencher| {
            bencher.iter(|| {
                let start = regions[next % regions.len()];
                next += 1;
                let bins = reader
                    .summaries(CHROM, start, start + SUMMARY_REGION_BASES, SUMMARY_BINS, source)
                    .expect("summaries");
                black_box(bins)
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_write_fixed_step,
    bench_write_bedgraph,
    bench_write_sparse_variable_step,
    bench_read_values,
    bench_summaries
);
criterion_main!(benches);
