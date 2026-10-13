//! `dump`: this crate's reading of a bigWig, on stdout.
//!
//! `bedgraph` and `wig` print what `bigwig::convert` writes. Every other subcommand writes TSV:
//! one header line and then one row per record, through fgoxide's delimited writer. f32 values
//! are widened to f64 before printing, so the printed shortest round-trip decimal parses back
//! (as a double) to exactly the stored value.

use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;

use bigwig::convert::{self, Selection};
use bigwig::{BigWigReader, ByteOrder, DataItem, Summary, SummarySource};
use fgoxide::io::DelimFileWriter;
use serde::Serialize;

/// Header fields and the stored total summary of a file. The summary columns are empty when the
/// file stores none.
#[derive(Debug, Serialize)]
struct InfoRow {
    version: u16,
    byte_order: &'static str,
    /// Zoom reductions in header order, comma-separated.
    zoom_reductions: String,
    uncompress_buf_size: u32,
    section_count: u64,
    chrom_count: usize,
    bases_covered: Option<f64>,
    min: Option<f64>,
    max: Option<f64>,
    sum: Option<f64>,
    sum_squares: Option<f64>,
    mean: Option<f64>,
    std: Option<f64>,
}

/// One chromosome.
#[derive(Debug, Serialize)]
struct ChromRow<'a> {
    name: &'a str,
    chrom_id: u32,
    length: u32,
}

/// One stored item, or one merged run of items. Also the format of the expected items that
/// `generate` writes next to each text input.
#[derive(Debug, Serialize)]
pub struct ItemRow<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub value: f64,
}

/// One base of a `values` query.
#[derive(Debug, Serialize)]
struct ValueRow {
    position: u32,
    value: f64,
}

/// One bin of a `summaries` query. `source` is `raw` or `zoom:<reduction>`; the statistics are
/// empty when no base of the bin has data.
#[derive(Debug, Serialize)]
struct SummaryRow {
    start: u32,
    end: u32,
    source: String,
    covered_bases: Option<f64>,
    min: Option<f64>,
    max: Option<f64>,
    sum: Option<f64>,
    sum_squares: Option<f64>,
    mean: Option<f64>,
    std: Option<f64>,
}

/// One zoom record, with its fields as stored.
#[derive(Debug, Serialize)]
struct ZoomRow<'a> {
    chrom: &'a str,
    start: u32,
    end: u32,
    covered_bases: f64,
    min: f64,
    max: f64,
    sum: f64,
    sum_squares: f64,
}

/// Prints the header fields and stored total summary.
pub fn info(path: &Path) -> Result<(), Box<dyn Error>> {
    let reader = BigWigReader::open(path)?;
    let info = reader.info();
    let summary = info.total_summary;
    let zoom_reductions: Vec<String> =
        info.zoom_levels.iter().map(|level| level.reduction().to_string()).collect();
    let row = InfoRow {
        version: info.version,
        byte_order: match info.byte_order {
            ByteOrder::Little => "little",
            ByteOrder::Big => "big",
        },
        zoom_reductions: zoom_reductions.join(","),
        uncompress_buf_size: info.uncompress_buf_size,
        section_count: info.section_count,
        chrom_count: reader.chroms().len(),
        bases_covered: summary.map(|s| s.covered_bases),
        min: summary.map(|s| s.min),
        max: summary.map(|s| s.max),
        sum: summary.map(|s| s.sum),
        sum_squares: summary.map(|s| s.sum_squares),
        mean: summary.map(|s| s.mean()),
        std: summary.map(|s| s.std()),
    };
    write_rows([row])
}

/// Prints the chromosomes in `chromId` order.
pub fn chroms(path: &Path) -> Result<(), Box<dyn Error>> {
    let reader = BigWigReader::open(path)?;
    let mut rows = Vec::with_capacity(reader.chroms().len());
    for (chrom_id, chrom) in (0u32..).zip(reader.chroms()) {
        rows.push(ChromRow { name: &chrom.name, chrom_id, length: chrom.length });
    }
    write_rows(rows)
}

/// Prints the stored items of the selection, unclipped: those overlapping it.
pub fn items(path: &Path, selection: &Selection) -> Result<(), Box<dyn Error>> {
    let mut reader = BigWigReader::open(path)?;
    let names = chrom_names(&reader);
    let items = match selection {
        Selection::All => reader.all_items()?,
        Selection::Chrom(chrom) => reader.chrom_items(chrom)?,
        Selection::Region { chrom, start, end } => reader.items(chrom, *start, *end)?,
    };
    let mut writer = stdout_writer();
    for item in items {
        let item = item?;
        writer.write(&item_row(&names, &item))?;
    }
    writer.flush()?;
    Ok(())
}

/// Prints the selection as bedGraph text with abutting equal values merged, as kent's
/// `bigWigToBedGraph` does: `convert::bigwig_to_bedgraph`'s output.
pub fn bedgraph(path: &Path, selection: &Selection) -> Result<(), Box<dyn Error>> {
    let mut reader = BigWigReader::open(path)?;
    convert::bigwig_to_bedgraph(&mut reader, selection, true, io::stdout().lock())?;
    Ok(())
}

/// Prints the selection as wig text, `convert::bigwig_to_wig`'s output.
pub fn wig(path: &Path, selection: &Selection) -> Result<(), Box<dyn Error>> {
    let mut reader = BigWigReader::open(path)?;
    convert::bigwig_to_wig(&mut reader, selection, io::stdout().lock())?;
    Ok(())
}

/// Prints the value of every base of `[start, end)` on `chrom`.
pub fn values(path: &Path, chrom: &str, start: u32, end: u32) -> Result<(), Box<dyn Error>> {
    let mut reader = BigWigReader::open(path)?;
    let values = reader.values(chrom, start, end)?;
    let rows = (start..)
        .zip(values)
        .map(|(position, value)| ValueRow { position, value: f64::from(value) });
    write_rows(rows)
}

/// Prints the summaries of `bins` equal bins of `[start, end)` on `chrom`.
pub fn summaries(
    path: &Path,
    chrom: &str,
    start: u32,
    end: u32,
    bins: u32,
    source: SummarySource,
) -> Result<(), Box<dyn Error>> {
    let mut reader = BigWigReader::open(path)?;
    // `summaries` validates the arguments, so the bin width is only used once they are valid.
    let summaries = reader.summaries(chrom, start, end, bins, source)?;
    let used = match source {
        SummarySource::Raw => None,
        SummarySource::BestZoom => reader.best_zoom_for((end - start) / bins),
    };
    let source =
        used.map_or_else(|| "raw".to_owned(), |level| format!("zoom:{}", level.reduction()));
    let length = u64::from(end - start);
    let rows = (0..bins).zip(summaries).map(|(index, summary)| {
        // The bin bounds as `summaries` documents them.
        let bin_start = start + (length * u64::from(index) / u64::from(bins)) as u32;
        let bin_end = start + (length * u64::from(index + 1) / u64::from(bins)) as u32;
        summary_row(bin_start, bin_end, &source, summary)
    });
    write_rows(rows)
}

/// Prints the records of zoom level `level` (a header index) overlapping `[start, end)` on
/// `chrom`.
pub fn zoom_records(
    path: &Path,
    level: usize,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<(), Box<dyn Error>> {
    let mut reader = BigWigReader::open(path)?;
    let names = chrom_names(&reader);
    let Some(&zoom_level) = reader.zoom_levels().get(level) else {
        return Err(format!("{} has no zoom level {level}", path.display()).into());
    };
    let mut writer = stdout_writer();
    for record in reader.zoom_records(zoom_level, chrom, start, end)? {
        let record = record?;
        let Summary { covered_bases, min, max, sum, sum_squares } = record.summary;
        writer.write(&ZoomRow {
            chrom: &names[record.chrom_id as usize],
            start: record.start,
            end: record.end,
            covered_bases,
            min,
            max,
            sum,
            sum_squares,
        })?;
    }
    writer.flush()?;
    Ok(())
}

/// The chromosome names indexed by `chromId`, copied so items can be named while the reader is
/// borrowed by an iterator.
fn chrom_names(reader: &BigWigReader<BufReader<File>>) -> Vec<String> {
    reader.chroms().iter().map(|chrom| chrom.name.clone()).collect()
}

/// The row for one item.
fn item_row<'a>(names: &'a [String], item: &DataItem) -> ItemRow<'a> {
    ItemRow {
        chrom: &names[item.chrom_id as usize],
        start: item.start,
        end: item.end,
        value: f64::from(item.value),
    }
}

/// The row for one summary bin.
fn summary_row(start: u32, end: u32, source: &str, summary: Option<Summary>) -> SummaryRow {
    SummaryRow {
        start,
        end,
        source: source.to_owned(),
        covered_bases: summary.map(|s| s.covered_bases),
        min: summary.map(|s| s.min),
        max: summary.map(|s| s.max),
        sum: summary.map(|s| s.sum),
        sum_squares: summary.map(|s| s.sum_squares),
        mean: summary.map(|s| s.mean()),
        std: summary.map(|s| s.std()),
    }
}

/// A tab-delimited, unquoted writer on stdout.
fn stdout_writer<S: Serialize>() -> DelimFileWriter<S> {
    let stdout: Box<dyn Write + Send> = Box::new(io::stdout());
    DelimFileWriter::new(BufWriter::new(stdout), b'\t', false)
}

/// Writes `rows` to stdout.
fn write_rows<S: Serialize>(rows: impl IntoIterator<Item = S>) -> Result<(), Box<dyn Error>> {
    let mut writer = stdout_writer();
    writer.write_all(rows)?;
    writer.flush()?;
    Ok(())
}
