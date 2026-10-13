//! Writes a dozen varied, valid bigWig files into a directory (default `corpus/read_bigwig`) as
//! the seed corpus for the `read_bigwig` target. The files are generated, never checked in.

use std::fs;
use std::io::Cursor;
use std::path::PathBuf;

use bigwig::{
    BigWigWriter, ByteOrder, ChromIdOrder, ChromInfo, Compression, WriterOptions, ZoomPolicy,
};

type Writer = BigWigWriter<Cursor<Vec<u8>>>;

/// A small deterministic generator, so seeds are the same on every run.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u32 {
        self.0 =
            self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) as u32
    }

    fn value(&mut self) -> f32 {
        (self.next() % 400) as f32 / 8.0 - 20.0
    }
}

fn chroms(count: usize, length: u32) -> Vec<ChromInfo> {
    (0..count).map(|i| ChromInfo::new(format!("chr{i}"), length)).collect()
}

fn bedgraph(writer: &mut Writer, chrom: &str, length: u32, rng: &mut Lcg) {
    let mut position = 0;
    while position + 20 < length {
        let start = position + rng.next() % 10;
        let end = start + 1 + rng.next() % 9;
        writer.push_interval(chrom, start, end, rng.value()).unwrap();
        position = end;
    }
}

fn variable(writer: &mut Writer, chrom: &str, length: u32, span: u32, rng: &mut Lcg) {
    let mut items = Vec::new();
    let mut position = 0;
    while position + span + 10 < length {
        let start = position + rng.next() % 10;
        items.push((start, rng.value()));
        position = start + span;
    }
    writer.push_variable_step(chrom, span, &items).unwrap();
}

fn fixed(writer: &mut Writer, chrom: &str, length: u32, step: u32, rng: &mut Lcg) {
    let count = (length / step).min(500) as usize;
    let values: Vec<f32> = (0..count).map(|_| rng.value()).collect();
    writer.push_fixed_step(chrom, 0, step, step, &values).unwrap();
}

fn finish(writer: Writer) -> Vec<u8> {
    writer.finish().unwrap().into_inner()
}

fn main() {
    let dir = PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "corpus/read_bigwig".into()));
    fs::create_dir_all(&dir).unwrap();
    let mut rng = Lcg(1);
    let mut seeds: Vec<(&str, Vec<u8>)> = Vec::new();
    let new = |options: WriterOptions, chroms: &[ChromInfo]| {
        options.create_with(Cursor::new(Vec::new()), chroms).unwrap()
    };

    let list = chroms(2, 2_000);
    let mut w = new(WriterOptions::new(), &list);
    bedgraph(&mut w, "chr0", 2_000, &mut rng);
    bedgraph(&mut w, "chr1", 2_000, &mut rng);
    seeds.push(("bedgraph-default", finish(w)));

    let mut w = new(WriterOptions::new().items_per_section(7), &list);
    variable(&mut w, "chr0", 2_000, 3, &mut rng);
    variable(&mut w, "chr1", 2_000, 1, &mut rng);
    seeds.push(("variable-step-small-sections", finish(w)));

    let mut w = new(WriterOptions::new().items_per_section(50), &list);
    fixed(&mut w, "chr0", 2_000, 5, &mut rng);
    fixed(&mut w, "chr1", 2_000, 1, &mut rng);
    seeds.push(("fixed-step", finish(w)));

    let mut w = new(WriterOptions::new().items_per_section(9), &list);
    bedgraph(&mut w, "chr0", 1_000, &mut rng);
    seeds.push(("mixed-bedgraph-then-variable", {
        variable(&mut w, "chr1", 2_000, 2, &mut rng);
        finish(w)
    }));

    let options = WriterOptions::new().compression(Compression::None).zoom(ZoomPolicy::None);
    let mut w = new(options, &list);
    bedgraph(&mut w, "chr0", 2_000, &mut rng);
    seeds.push(("uncompressed-no-zoom", finish(w)));

    let options = WriterOptions::new()
        .byte_order(ByteOrder::Big)
        .zoom(ZoomPolicy::Auto { min_first_reduction: 4 });
    let mut w = new(options, &list);
    bedgraph(&mut w, "chr0", 2_000, &mut rng);
    fixed(&mut w, "chr1", 2_000, 4, &mut rng);
    seeds.push(("big-endian-fine-zoom", finish(w)));

    let options = WriterOptions::new()
        .zoom(ZoomPolicy::Explicit(vec![4, 16, 65, 200]))
        .zoom_items_per_block(5)
        .compression(Compression::Deflate { level: 9 });
    let mut w = new(options, &list);
    bedgraph(&mut w, "chr0", 2_000, &mut rng);
    bedgraph(&mut w, "chr1", 2_000, &mut rng);
    seeds.push(("explicit-zoom-small-blocks", finish(w)));

    seeds.push(("empty", finish(new(WriterOptions::new(), &list))));
    seeds.push((
        "empty-no-zoom-big-endian",
        finish(new(
            WriterOptions::new().zoom(ZoomPolicy::None).byte_order(ByteOrder::Big),
            &chroms(1, 10),
        )),
    ));

    let many = chroms(300, 100);
    let options = WriterOptions::new().index_block_size(2).items_per_section(2);
    let mut w = new(options, &many);
    for i in (0..300).step_by(7) {
        bedgraph(&mut w, &format!("chr{i}"), 100, &mut rng);
    }
    seeds.push(("many-chroms-binary-rtree", finish(w)));

    let options = WriterOptions::new().chrom_id_order(ChromIdOrder::Listed).clip(true);
    let mut w = new(options, &chroms(3, 500));
    bedgraph(&mut w, "chr0", 500, &mut rng);
    w.push_interval("chr2", 490, 900, 3.5).unwrap();
    seeds.push(("listed-order-clipped", finish(w)));

    let options =
        WriterOptions::new().zoom_memory_limit(0).zoom(ZoomPolicy::Auto { min_first_reduction: 1 });
    let mut w = new(options, &chroms(1, 5_000));
    fixed(&mut w, "chr0", 5_000, 1, &mut rng);
    seeds.push(("spilled-zoom-blocks", finish(w)));

    for (name, bytes) in &seeds {
        fs::write(dir.join(format!("{name}.bw")), bytes).unwrap();
    }
    println!("wrote {} seeds to {}", seeds.len(), dir.display());
}
