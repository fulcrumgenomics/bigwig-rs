//! `generate`: every interop fixture, written deterministically from a seed.
//!
//! Our bigWigs go under `ours/`, one per fixture. Text inputs for the other tools' writers go
//! under `text/`, each with an `.items.tsv` listing exactly the items a writer must store for it,
//! taken from the generated data rather than from any parser.
//!
//! Most values are multiples of 1/8 below 512 in magnitude, which f32, kent's `%g` and every
//! text round trip represent exactly; a few arbitrary f32 values exercise full-precision
//! printing.

use std::error::Error;
use std::fs::{self, File};
use std::io::{BufWriter, Seek, Write};
use std::path::Path;

use bigwig::{
    BigWigWriter, ByteOrder, ChromIdOrder, ChromInfo, Compression, WriterOptions, ZoomPolicy,
    bedgraph, wig,
};
use fgoxide::io::DelimFile;
use serde::Serialize;

use crate::dump::ItemRow;

/// Subdirectory of our bigWig fixtures.
const OURS_DIR: &str = "ours";
/// Subdirectory of the text inputs for the other tools' writers.
const TEXT_DIR: &str = "text";
/// Values per `push_fixed_step` call for long runs, as a streaming producer would push them.
const FIXED_STEP_CHUNK: usize = 10_000;
/// Items per `push_variable_step` call for long runs.
const VARIABLE_STEP_CHUNK: usize = 2_000;
/// Primary hg38 chromosomes in BAM-header order, with their lengths.
const HG38_PRIMARY: [(&str, u32); 25] = [
    ("chr1", 248_956_422),
    ("chr2", 242_193_529),
    ("chr3", 198_295_559),
    ("chr4", 190_214_555),
    ("chr5", 181_538_259),
    ("chr6", 170_805_979),
    ("chr7", 159_345_973),
    ("chr8", 145_138_636),
    ("chr9", 138_394_717),
    ("chr10", 133_797_422),
    ("chr11", 135_086_622),
    ("chr12", 133_275_309),
    ("chr13", 114_364_328),
    ("chr14", 107_043_718),
    ("chr15", 101_991_189),
    ("chr16", 90_338_345),
    ("chr17", 83_257_441),
    ("chr18", 80_373_285),
    ("chr19", 58_617_616),
    ("chr20", 64_444_167),
    ("chr21", 46_709_983),
    ("chr22", 50_818_468),
    ("chrX", 156_040_895),
    ("chrY", 57_227_415),
    ("chrM", 16_569),
];

/// One bigWig fixture: what to write, and how.
struct Fixture {
    name: &'static str,
    description: &'static str,
    chroms: Vec<ChromInfo>,
    options: WriterOptions,
    pushes: Vec<Push>,
}

impl Fixture {
    /// Writes the fixture's bigWig to `path` by replaying its pushes.
    fn write_bigwig(&self, path: &Path) -> Result<(), Box<dyn Error>> {
        let mut writer = self.options.create(path, &self.chroms)?;
        for push in &self.pushes {
            push.write_to(&mut writer)?;
        }
        writer.finish()?;
        Ok(())
    }
}

/// One row of `manifest.tsv`.
#[derive(Debug, Serialize)]
struct ManifestRow {
    fixture: &'static str,
    file: String,
    description: &'static str,
}

/// One writer call, which is also one wig block (or run of bedGraph lines) in the text
/// fixtures.
#[derive(Clone, Debug)]
enum Push {
    /// `values[i]` covers `[start + i·step, start + i·step + span)`.
    Fixed { chrom: String, start: u32, step: u32, span: u32, values: Vec<f32> },
    /// Each `(start, value)` covers `[start, start + span)`.
    Variable { chrom: String, span: u32, items: Vec<(u32, f32)> },
    /// bedGraph intervals `(start, end, value)`.
    Intervals { chrom: String, intervals: Vec<(u32, u32, f32)> },
}

impl Push {
    /// Replays this push on `writer`, in chunks the size a streaming producer would use.
    fn write_to<W: Write + Seek>(&self, writer: &mut BigWigWriter<W>) -> Result<(), bigwig::Error> {
        match self {
            Self::Fixed { chrom, start, step, span, values } => {
                let mut chunk_start = *start;
                for chunk in values.chunks(FIXED_STEP_CHUNK) {
                    writer.push_fixed_step(chrom, chunk_start, *step, *span, chunk)?;
                    chunk_start += step * chunk.len() as u32;
                }
                Ok(())
            }
            Self::Variable { chrom, span, items } => {
                for chunk in items.chunks(VARIABLE_STEP_CHUNK) {
                    writer.push_variable_step(chrom, *span, chunk)?;
                }
                Ok(())
            }
            Self::Intervals { chrom, intervals } => writer.push_intervals(chrom, intervals),
        }
    }

    /// Writes this push as a wig block, or as bedGraph lines for intervals.
    fn write_wig<W: Write>(&self, writer: &mut wig::Writer<W>) -> Result<(), bigwig::Error> {
        match self {
            Self::Fixed { chrom, start, step, span, values } => {
                writer.fixed_step(chrom, *start, *step, *span)?;
                values.iter().try_for_each(|&value| writer.fixed_value(value))
            }
            Self::Variable { chrom, span, items } => {
                writer.variable_step(chrom, *span)?;
                items.iter().try_for_each(|&(start, value)| writer.variable_value(start, value))
            }
            Self::Intervals { chrom, intervals } => {
                intervals.iter().try_for_each(|&(start, end, value)| {
                    writer.bedgraph_record(chrom, start, end, value)
                })
            }
        }
    }

    /// The chromosome this push writes to.
    fn chrom(&self) -> &str {
        match self {
            Self::Fixed { chrom, .. }
            | Self::Variable { chrom, .. }
            | Self::Intervals { chrom, .. } => chrom,
        }
    }

    /// The items this push stores, as `(start, end, value)`, in order.
    fn items(&self) -> Vec<(u32, u32, f32)> {
        match self {
            Self::Fixed { start, step, span, values, .. } => (0u32..)
                .zip(values)
                .map(|(index, &value)| {
                    let item_start = start + index * step;
                    (item_start, item_start + span, value)
                })
                .collect(),
            Self::Variable { span, items, .. } => {
                items.iter().map(|&(start, value)| (start, start + span, value)).collect()
            }
            Self::Intervals { intervals, .. } => intervals.clone(),
        }
    }
}

/// splitmix64 (Steele, Lea and Flood 2014): a tiny, well-known generator, ample for
/// reproducible test data and needing no dependency.
pub(crate) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// A generator for the stream called `name` under `seed`, so each fixture's data is
    /// independent of which other fixtures exist.
    pub(crate) fn new(seed: u64, name: &str) -> Self {
        // FNV-1a of the name, mixed into the seed.
        let name_hash = name.bytes().fold(0xCBF2_9CE4_8422_2325_u64, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01B3)
        });
        Self { state: seed ^ name_hash }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A value in `[low, high)`; the modulo bias is irrelevant here.
    fn range(&mut self, low: u32, high: u32) -> u32 {
        assert!(low < high, "empty range {low}..{high}");
        low + (self.next_u64() % u64::from(high - low)) as u32
    }

    /// True with probability `1 / n`.
    fn one_in(&mut self, n: u32) -> bool {
        self.range(0, n) == 0
    }

    /// A signed signal value: usually a multiple of 1/8 in `[-256, 256)`, sometimes an arbitrary
    /// f32 in `[-1000, 1000)` that needs all nine significant digits.
    fn signal(&mut self) -> f32 {
        if self.one_in(8) {
            let unit = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
            (unit * 2_000.0 - 1_000.0) as f32
        } else {
            (f64::from(self.range(0, 4_096)) / 8.0 - 256.0) as f32
        }
    }

    /// A non-negative multiple of 1/8 below `max`.
    fn eighths(&mut self, max: u32) -> f32 {
        self.range(0, max * 8) as f32 / 8.0
    }
}

/// Writes every fixture under `out_dir`, plus `manifest.tsv` describing each bigWig.
pub fn generate(seed: u64, out_dir: &Path) -> Result<(), Box<dyn Error>> {
    let ours = out_dir.join(OURS_DIR);
    fs::create_dir_all(&ours)?;
    let mut manifest = Vec::new();
    for fixture in bigwig_fixtures(seed) {
        let file = format!("{OURS_DIR}/{}.bw", fixture.name);
        fixture.write_bigwig(&out_dir.join(&file))?;
        manifest.push(ManifestRow {
            fixture: fixture.name,
            file,
            description: fixture.description,
        });
    }
    DelimFile::default().write_tsv(out_dir.join("manifest.tsv"), manifest)?;

    let text = out_dir.join(TEXT_DIR);
    fs::create_dir_all(&text)?;
    write_text_fixtures(seed, &text)
}

/// Every bigWig fixture.
fn bigwig_fixtures(seed: u64) -> Vec<Fixture> {
    let sections_chroms = vec![
        ChromInfo::new("chr1", 400_000),
        ChromInfo::new("chr2", 250_000),
        ChromInfo::new("chr3", 180_000),
        ChromInfo::new("chrM", 16_569),
    ];
    let sections = mixed_sections(&mut SplitMix64::new(seed, "sections"), &sections_chroms);
    let riker_chroms = vec![
        ChromInfo::new("chr1", 300_000),
        ChromInfo::new("chr2", 220_000),
        ChromInfo::new("chr10", 150_000),
        ChromInfo::new("chrX", 120_000),
        ChromInfo::new("chrM", 16_569),
    ];
    let listed = WriterOptions::new().chrom_id_order(ChromIdOrder::Listed);
    let cpg_chroms = vec![ChromInfo::new("chr1", 600_000), ChromInfo::new("chr2", 300_000)];

    vec![
        Fixture {
            name: "sections",
            description: "fixedStep, variableStep and bedGraph sections mixed on each chromosome; \
                          default options (Auto zoom)",
            chroms: sections_chroms.clone(),
            options: WriterOptions::new(),
            pushes: sections.clone(),
        },
        Fixture {
            name: "zoom_explicit",
            description: "the sections data with Explicit zooms 50, 130, 1000, 4000, 20000 and \
                          1,000 items per section",
            chroms: sections_chroms.clone(),
            options: WriterOptions::new()
                .zoom(ZoomPolicy::Explicit(vec![50, 130, 1_000, 4_000, 20_000]))
                .items_per_section(1_000),
            pushes: sections.clone(),
        },
        Fixture {
            name: "uncompressed",
            description: "the sections data stored uncompressed",
            chroms: sections_chroms.clone(),
            options: WriterOptions::new().compression(Compression::None),
            pushes: sections.clone(),
        },
        Fixture {
            name: "big_endian",
            description: "the sections data written big-endian",
            chroms: sections_chroms,
            options: WriterOptions::new().byte_order(ByteOrder::Big),
            pushes: sections,
        },
        Fixture {
            name: "riker_per_base",
            description: "riker per-base depth: step-1 fixedStep runs in chunks, BAM-order \
                          chromosomes with Listed IDs and a zero-data chromosome",
            chroms: riker_chroms.clone(),
            options: listed.clone(),
            pushes: per_base_depth(&mut SplitMix64::new(seed, "riker_per_base"), &riker_chroms),
        },
        Fixture {
            name: "riker_binned",
            description: "riker binned depth: 100-base fixedStep bins over whole chromosomes",
            chroms: riker_chroms.clone(),
            options: listed.clone(),
            pushes: binned_depth(&mut SplitMix64::new(seed, "riker_binned"), &riker_chroms, 100),
        },
        Fixture {
            name: "riker_quantised",
            description: "riker quantised depth: abutting bedGraph runs of quantised values",
            chroms: riker_chroms.clone(),
            options: listed.clone(),
            pushes: quantised_runs(&mut SplitMix64::new(seed, "riker_quantised"), &riker_chroms),
        },
        contigs_3000(seed),
        contigs_40000(seed),
        locale_names(seed),
        Fixture {
            name: "cpg_sparse",
            description: "sparse per-CpG methylation fractions as one-base variableStep items",
            pushes: cpg_sites(&mut SplitMix64::new(seed, "cpg_sparse"), &cpg_chroms),
            chroms: cpg_chroms,
            options: WriterOptions::new(),
        },
        Fixture {
            name: "empty",
            description: "two chromosomes and no data",
            chroms: vec![ChromInfo::new("chr1", 1_000_000), ChromInfo::new("chr2", 500_000)],
            options: WriterOptions::new(),
            pushes: Vec::new(),
        },
    ]
}

/// 3,000 contigs in BAM-header order (primary chromosomes, then alts and unplaced contigs), with
/// Listed IDs and data on six of them, including the last.
fn contigs_3000(seed: u64) -> Fixture {
    let mut rng = SplitMix64::new(seed, "contigs_3000");
    let mut chroms: Vec<ChromInfo> =
        HG38_PRIMARY.iter().map(|&(name, length)| ChromInfo::new(name, length)).collect();
    let mut index = 0u32;
    while chroms.len() < 3_000 {
        let name = match index % 3 {
            0 => format!("chr{}_KI{:06}v1_random", 1 + index % 22, 270_000 + index),
            1 => format!("chrUn_GL{:06}v1", 1_000 + index),
            _ => format!("HLA-A*{:02}:{:02}:{:02}", index / 900, (index / 30) % 30, index % 30),
        };
        chroms.push(ChromInfo::new(name, rng.range(1_000, 200_000)));
        index += 1;
    }
    let with_data = [0, 6, 22, 24, 1_500, 2_999];
    let pushes = with_data
        .iter()
        .flat_map(|&position| {
            let chrom = &chroms[position];
            let window = ChromInfo::new(chrom.name.clone(), chrom.length.min(250_000));
            quantised_runs(&mut rng, &[window])
        })
        .collect();
    Fixture {
        name: "contigs_3000",
        description: "3,000 contigs in BAM-header (non-lexicographic) order with Listed IDs, \
                      data on six",
        chroms,
        options: WriterOptions::new().chrom_id_order(ChromIdOrder::Listed),
        pushes,
    }
}

/// 40,000 contigs named `scaffold_1` .. `scaffold_40000` (so not in byte order), which needs a
/// multi-level chromosome B+ tree, with data on twelve pushed in an order unrelated to the list.
fn contigs_40000(seed: u64) -> Fixture {
    let mut rng = SplitMix64::new(seed, "contigs_40000");
    let chroms: Vec<ChromInfo> = (1..=40_000)
        .map(|number| ChromInfo::new(format!("scaffold_{number}"), rng.range(500, 50_000)))
        .collect();
    let with_data = [39_999, 0, 32_767, 9, 19_999, 32_766, 99, 32_768, 1, 9_999, 999, 39_998];
    let pushes = with_data
        .iter()
        .flat_map(|&position| cpg_sites(&mut rng, &chroms[position..=position]))
        .collect();
    Fixture {
        name: "contigs_40000",
        description: "40,000 contigs (multi-level B+ tree), data on twelve in arrival order",
        chroms,
        options: WriterOptions::new(),
        pushes,
    }
}

/// Names whose byte order differs from JavaScript `localeCompare` order, which igv.js uses to
/// search the chromosome tree: case (`chrZ`/`chra`, `Chr3`/`chr1`), punctuation
/// (`chr_alt`/`chrEBV`) and mixed-case suffixes.
fn locale_names(seed: u64) -> Fixture {
    let names = [
        "chr1",
        "chr10",
        "chr1_KI270706v1_random",
        "chr1_gl000191_random",
        "chr2",
        "Chr3",
        "chrUn_KI270302v1",
        "chrUn_gl000220",
        "chrEBV",
        "chr_alt",
        "chrZ",
        "chra",
        "chrM",
        "HLA-A*01:01:01:01",
        "hla_lower",
        "chrX",
    ];
    let mut rng = SplitMix64::new(seed, "locale_names");
    let chroms: Vec<ChromInfo> =
        names.iter().map(|&name| ChromInfo::new(name, rng.range(20_000, 60_000))).collect();
    let pushes = quantised_runs(&mut rng, &chroms);
    Fixture {
        name: "locale_names",
        description: "contig names that localeCompare orders differently from bytes",
        chroms,
        options: WriterOptions::new(),
        pushes,
    }
}

/// Runs of all three section types along each chromosome, with gaps, abutting runs, abutting
/// equal-valued intervals, and fixedStep runs whose span is less than their step.
fn mixed_sections(rng: &mut SplitMix64, chroms: &[ChromInfo]) -> Vec<Push> {
    let mut pushes = Vec::new();
    for chrom in chroms {
        let name = chrom.name.clone();
        let mut position = rng.range(0, 1_000);
        // Stop well short of the end so every run fits.
        while position + 20_000 < chrom.length {
            let push = match rng.range(0, 3) {
                0 => {
                    let step = [1, 1, 5, 20, 50][rng.range(0, 5) as usize];
                    let span = rng.range(1, step + 1);
                    let count = rng.range(50, 3_000).min((chrom.length - position - span) / step);
                    let values = (0..count).map(|_| rng.signal()).collect();
                    Push::Fixed { chrom: name.clone(), start: position, step, span, values }
                }
                1 => {
                    let span = [1, 3, 10][rng.range(0, 3) as usize];
                    let mut start = position;
                    let mut items = Vec::new();
                    for _ in 0..rng.range(20, 1_500) {
                        if start + span > chrom.length {
                            break;
                        }
                        items.push((start, rng.signal()));
                        start += span + rng.range(0, 40);
                    }
                    Push::Variable { chrom: name.clone(), span, items }
                }
                _ => {
                    let mut start = position;
                    let mut intervals: Vec<(u32, u32, f32)> = Vec::new();
                    for _ in 0..rng.range(20, 1_500) {
                        if start >= chrom.length {
                            break;
                        }
                        let end = (start + rng.range(1, 500)).min(chrom.length);
                        let value = match intervals.last() {
                            Some(&(_, previous_end, previous))
                                if previous_end == start && rng.one_in(3) =>
                            {
                                previous
                            }
                            _ => rng.signal(),
                        };
                        intervals.push((start, end, value));
                        start = if rng.one_in(3) { end } else { end + rng.range(1, 200) };
                    }
                    Push::Intervals { chrom: name.clone(), intervals }
                }
            };
            let end = push.items().last().map_or(position, |&(_, end, _)| end);
            pushes.push(push);
            position = if rng.one_in(4) { end } else { end + rng.range(1, 2_000) };
        }
    }
    pushes
}

/// Per-base depth as riker writes it: covered stretches of step-1 fixedStep values from a
/// random walk (so long runs of equal values), separated by uncovered gaps. The fourth
/// chromosome gets no data.
fn per_base_depth(rng: &mut SplitMix64, chroms: &[ChromInfo]) -> Vec<Push> {
    let mut pushes = Vec::new();
    for (index, chrom) in chroms.iter().enumerate() {
        if index == 3 {
            continue;
        }
        let mut position = rng.range(0, 5_000);
        while position < chrom.length {
            let count = rng.range(2_000, 60_000).min(chrom.length - position);
            let mut depth = rng.range(0, 60);
            let values = (0..count)
                .map(|_| {
                    if rng.one_in(8) {
                        depth = (depth + rng.range(0, 7)).saturating_sub(3);
                    }
                    depth as f32
                })
                .collect();
            pushes.push(Push::Fixed {
                chrom: chrom.name.clone(),
                start: position,
                step: 1,
                span: 1,
                values,
            });
            position += count + rng.range(100, 20_000);
        }
    }
    pushes
}

/// Binned depth: one fixedStep run of `bin`-wide bins over each whole chromosome, with mean
/// depths in eighths and stretches of zero.
fn binned_depth(rng: &mut SplitMix64, chroms: &[ChromInfo], bin: u32) -> Vec<Push> {
    chroms
        .iter()
        .map(|chrom| {
            let mut zero_run = 0;
            let values = (0..chrom.length / bin)
                .map(|_| {
                    if zero_run == 0 && rng.one_in(200) {
                        zero_run = rng.range(1, 50);
                    }
                    if zero_run > 0 {
                        zero_run -= 1;
                        0.0
                    } else {
                        rng.eighths(80)
                    }
                })
                .collect();
            Push::Fixed { chrom: chrom.name.clone(), start: 0, step: bin, span: bin, values }
        })
        .collect()
}

/// Quantised depth: abutting bedGraph runs whose value changes at every run boundary, with
/// occasional gaps.
fn quantised_runs(rng: &mut SplitMix64, chroms: &[ChromInfo]) -> Vec<Push> {
    let mut pushes = Vec::new();
    for chrom in chroms {
        let mut intervals: Vec<(u32, u32, f32)> = Vec::new();
        let mut start = rng.range(0, 2_000).min(chrom.length - 1);
        while start < chrom.length {
            let end = (start + rng.range(1, 3_000)).min(chrom.length);
            let mut value = rng.range(0, 129) as f32 / 2.0;
            if let Some(&(_, previous_end, previous)) = intervals.last()
                && previous_end == start
                && previous == value
            {
                value += 0.5;
            }
            intervals.push((start, end, value));
            start = if rng.one_in(10) { end + rng.range(1, 5_000) } else { end };
        }
        pushes.push(Push::Intervals { chrom: chrom.name.clone(), intervals });
    }
    pushes
}

/// Sparse per-CpG data: one-base variableStep items a few to a few hundred bases apart, valued
/// at methylation fractions `k / n`.
fn cpg_sites(rng: &mut SplitMix64, chroms: &[ChromInfo]) -> Vec<Push> {
    let mut pushes = Vec::new();
    for chrom in chroms {
        let mut items = Vec::new();
        let mut position = rng.range(0, 400);
        while position < chrom.length {
            let reads = rng.range(1, 31);
            items.push((position, rng.range(0, reads + 1) as f32 / reads as f32));
            position += rng.range(2, 400);
        }
        pushes.push(Push::Variable { chrom: chrom.name.clone(), span: 1, items });
    }
    pushes
}

/// Writes the text inputs for the other tools' writers, each with its `.items.tsv`, and the
/// shared `chrom.sizes`. Chromosomes are in BAM order (`chr10` after `chr2`), which is not byte
/// order.
fn write_text_fixtures(seed: u64, dir: &Path) -> Result<(), Box<dyn Error>> {
    let chroms = [
        ChromInfo::new("chr1", 200_000),
        ChromInfo::new("chr2", 150_000),
        ChromInfo::new("chr10", 100_000),
        ChromInfo::new("chrX", 80_000),
    ];
    let mut sizes = BufWriter::new(File::create(dir.join("chrom.sizes"))?);
    for chrom in &chroms {
        writeln!(sizes, "{}\t{}", chrom.name, chrom.length)?;
    }
    sizes.flush()?;

    let mut rng = SplitMix64::new(seed, "text");
    let bedgraph_pushes = text_intervals(&mut rng, &chroms);
    let mut bedgraph_writer =
        bedgraph::Writer::new(BufWriter::new(File::create(dir.join("input.bedGraph"))?));
    for push in &bedgraph_pushes {
        for (start, end, value) in push.items() {
            bedgraph_writer.write_record(push.chrom(), start, end, value)?;
        }
    }
    bedgraph_writer.into_inner()?.flush()?;
    write_expected_items(&dir.join("input.items.tsv"), &bedgraph_pushes)?;

    let fixed = text_blocks(&mut rng, &chroms, &[BlockKind::Fixed]);
    let variable = text_blocks(&mut rng, &chroms, &[BlockKind::Variable]);
    let mixed = text_blocks(
        &mut rng,
        &chroms,
        &[BlockKind::Fixed, BlockKind::Variable, BlockKind::BedGraph],
    );
    for (name, pushes) in [("fixed", fixed), ("variable", variable), ("mixed", mixed)] {
        let path = dir.join(format!("{name}.wig"));
        let mut writer = wig::Writer::new(BufWriter::new(File::create(&path)?));
        for push in &pushes {
            push.write_wig(&mut writer)?;
        }
        writer.into_inner()?.flush()?;
        write_expected_items(&dir.join(format!("{name}.items.tsv")), &pushes)?;
    }
    Ok(())
}

/// The kind of block `text_blocks` lays down.
#[derive(Clone, Copy, Debug)]
enum BlockKind {
    Fixed,
    Variable,
    BedGraph,
}

/// bedGraph intervals for the text input: quantised runs with some abutting, a few
/// non-dyadic values (`0.1`, `1/3`) and some full-precision ones.
fn text_intervals(rng: &mut SplitMix64, chroms: &[ChromInfo]) -> Vec<Push> {
    let mut pushes = Vec::new();
    for chrom in chroms {
        let mut intervals = Vec::new();
        let mut start = rng.range(0, 1_000);
        while start + 2_000 < chrom.length {
            let end = start + rng.range(1, 2_000);
            let value = match rng.range(0, 20) {
                0 => 0.1,
                1 => 1.0 / 3.0,
                2 => rng.signal(),
                _ => rng.eighths(100),
            };
            intervals.push((start, end, value));
            start = if rng.one_in(3) { end } else { end + rng.range(1, 3_000) };
        }
        pushes.push(Push::Intervals { chrom: chrom.name.clone(), intervals });
    }
    pushes
}

/// Wig blocks along each chromosome, cycling through `kinds` at random: fixedStep with
/// `(step, span)` of (1, 1), (10, 10), (25, 5) or (100, 50); variableStep with span 1, 5 or 20;
/// or a run of bedGraph lines.
fn text_blocks(rng: &mut SplitMix64, chroms: &[ChromInfo], kinds: &[BlockKind]) -> Vec<Push> {
    let mut pushes = Vec::new();
    for chrom in chroms {
        let name = chrom.name.clone();
        let mut position = rng.range(0, 1_000);
        while position + 30_000 < chrom.length {
            let push = match kinds[rng.range(0, kinds.len() as u32) as usize] {
                BlockKind::Fixed => {
                    let (step, span) =
                        [(1, 1), (10, 10), (25, 5), (100, 50)][rng.range(0, 4) as usize];
                    let values = (0..rng.range(10, 250)).map(|_| rng.eighths(64)).collect();
                    Push::Fixed { chrom: name.clone(), start: position, step, span, values }
                }
                BlockKind::Variable => {
                    let span = [1, 5, 20][rng.range(0, 3) as usize];
                    let mut start = position;
                    let items = (0..rng.range(10, 250))
                        .map(|_| {
                            let item = (start, rng.eighths(64));
                            start += span + rng.range(0, 30);
                            item
                        })
                        .collect();
                    Push::Variable { chrom: name.clone(), span, items }
                }
                BlockKind::BedGraph => {
                    let mut start = position;
                    let intervals = (0..rng.range(5, 100))
                        .map(|_| {
                            let end = start + rng.range(1, 100);
                            let interval = (start, end, rng.eighths(64));
                            start = end + rng.range(0, 50);
                            interval
                        })
                        .collect();
                    Push::Intervals { chrom: name.clone(), intervals }
                }
            };
            // pyBigWig takes a fixedStep block to end at `start + count·step`, past its last item
            // when span < step, and rejects a later block that starts sooner (pyBigWig.c
            // addEntries), so the next block starts after that.
            let end = match &push {
                Push::Fixed { start, step, values, .. } => start + step * values.len() as u32,
                _ => push.items().last().map_or(position, |&(_, end, _)| end),
            };
            pushes.push(push);
            position = end + rng.range(1, 2_000);
        }
    }
    pushes
}

/// Writes the items `pushes` store, in order, to `path`.
fn write_expected_items(path: &Path, pushes: &[Push]) -> Result<(), Box<dyn Error>> {
    let mut rows = Vec::new();
    for push in pushes {
        for (start, end, value) in push.items() {
            rows.push(ItemRow { chrom: push.chrom(), start, end, value: f64::from(value) });
        }
    }
    DelimFile::default().write_tsv(path, rows)?;
    Ok(())
}
