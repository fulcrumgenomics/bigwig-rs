//! Support for the public-file corpus: `scan`, one pass over every item of a bigWig, a recomputed
//! total summary, and a write → read round trip of the items through our writer; and `sample`,
//! random items around which to place query regions.

use std::error::Error;
use std::io;
use std::path::Path;
use std::time::Instant;

use bigwig::{BigWigReader, ByteOrder, DataItem, Summary, WriterOptions};
use fgoxide::io::DelimFileWriter;
use serde::Serialize;

use crate::dump::ItemRow;
use crate::generate::SplitMix64;

/// bedGraph intervals buffered per `push_intervals` call in the round trip.
const PUSH_BATCH: usize = 65_536;

/// The scan of one file, as one TSV row. Header summary columns are empty when the file stores
/// none; round-trip columns describe what happened when the items were rewritten.
#[derive(Debug, Serialize)]
struct ScanRow {
    version: u16,
    byte_order: &'static str,
    zoom_reductions: String,
    uncompress_buf_size: u32,
    chrom_count: usize,
    chroms_with_data: usize,
    items: u64,
    item_bases: u64,
    /// Items that start before the previous item on their chromosome ends: overlaps or
    /// unsorted data, which only lax writers produce.
    overlapping_items: u64,
    header_bases_covered: Option<f64>,
    header_min: Option<f64>,
    header_max: Option<f64>,
    header_sum: Option<f64>,
    computed_bases_covered: f64,
    computed_min: f64,
    computed_max: f64,
    computed_sum: f64,
    /// `identical`, `skipped`, `rejected: <writer error>` or `differs: <first difference>`.
    round_trip: String,
    seconds: f64,
}

/// Scans the bigWig at `path`, rewriting its items to `round_trip` when given, and prints one
/// TSV row.
pub fn scan(path: &Path, round_trip: Option<&Path>) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let mut reader = BigWigReader::open(path)?;
    let info = reader.info().clone();
    let chroms = reader.chroms().to_vec();
    let stored = info.total_summary;

    let mut items = 0u64;
    let mut item_bases = 0u64;
    let mut overlapping_items = 0u64;
    let mut chroms_with_data = 0usize;
    let mut previous: Option<DataItem> = None;
    for item in reader.all_items()? {
        let item = item?;
        items += 1;
        item_bases += u64::from(item.end - item.start);
        match previous {
            Some(last) if last.chrom_id == item.chrom_id => {
                if item.start < last.end {
                    overlapping_items += 1;
                }
            }
            _ => chroms_with_data += 1,
        }
        previous = Some(item);
    }
    let computed = reader.compute_total_summary()?;
    let round_trip = match round_trip {
        None => "skipped".to_owned(),
        Some(out) => rewrite_and_compare(&mut reader, out, &computed)?,
    };

    let zoom_reductions: Vec<String> =
        info.zoom_levels.iter().map(|level| level.reduction().to_string()).collect();
    let row = ScanRow {
        version: info.version,
        byte_order: match info.byte_order {
            ByteOrder::Little => "little",
            ByteOrder::Big => "big",
        },
        zoom_reductions: zoom_reductions.join(","),
        uncompress_buf_size: info.uncompress_buf_size,
        chrom_count: chroms.len(),
        chroms_with_data,
        items,
        item_bases,
        overlapping_items,
        header_bases_covered: stored.map(|s| s.covered_bases),
        header_min: stored.map(|s| s.min),
        header_max: stored.map(|s| s.max),
        header_sum: stored.map(|s| s.sum),
        computed_bases_covered: computed.covered_bases,
        computed_min: computed.min,
        computed_max: computed.max,
        computed_sum: computed.sum,
        round_trip,
        seconds: started.elapsed().as_secs_f64(),
    };
    let stdout: Box<dyn io::Write + Send> = Box::new(io::stdout());
    let mut writer = DelimFileWriter::new(io::BufWriter::new(stdout), b'\t', false);
    writer.write(&row)?;
    writer.flush()?;
    Ok(())
}

/// Prints `count` items chosen uniformly at random from the whole file (reservoir sampling in
/// one pass, seeded so a file always gives the same sample), as item TSV rows in file order.
pub fn sample(path: &Path, count: usize, seed: u64) -> Result<(), Box<dyn Error>> {
    let mut reader = BigWigReader::open(path)?;
    let names: Vec<String> = reader.chroms().iter().map(|chrom| chrom.name.clone()).collect();
    let mut rng = SplitMix64::new(seed, "sample");
    let mut chosen: Vec<DataItem> = Vec::with_capacity(count);
    for (seen, item) in (0u64..).zip(reader.all_items()?) {
        let item = item?;
        if chosen.len() < count {
            chosen.push(item);
        } else {
            let slot = rng.next_u64() % (seen + 1);
            if let Some(kept) = usize::try_from(slot).ok().and_then(|slot| chosen.get_mut(slot)) {
                *kept = item;
            }
        }
    }
    chosen.sort_by_key(|item| (item.chrom_id, item.start));
    let stdout: Box<dyn io::Write + Send> = Box::new(io::stdout());
    let mut writer = DelimFileWriter::new(io::BufWriter::new(stdout), b'\t', false);
    for item in &chosen {
        writer.write(&ItemRow {
            chrom: &names[item.chrom_id as usize],
            start: item.start,
            end: item.end,
            value: f64::from(item.value),
        })?;
    }
    writer.flush()?;
    Ok(())
}

/// Writes every item of `reader` to a new bigWig at `out` as bedGraph intervals, with default
/// options, then reads it back and compares it with the source item by item and with the
/// source's recomputed total summary. Returns `identical`, `rejected: ...` when the writer
/// refuses the data (overlapping items or NaN values from lax writers), or `differs: ...`.
/// The rewritten file is removed afterwards.
fn rewrite_and_compare<R: io::Read + io::Seek>(
    reader: &mut BigWigReader<R>,
    out: &Path,
    computed: &Summary,
) -> Result<String, Box<dyn Error>> {
    let chroms = reader.chroms().to_vec();
    let outcome = (|| -> Result<String, Box<dyn Error>> {
        let mut writer = match WriterOptions::new().create(out, &chroms) {
            Ok(writer) => writer,
            Err(error) => return Ok(format!("rejected: {error}")),
        };
        let mut batch: Vec<(u32, u32, f32)> = Vec::with_capacity(PUSH_BATCH);
        let mut batch_chrom: Option<u32> = None;
        for item in reader.all_items()? {
            let item = item?;
            if batch_chrom != Some(item.chrom_id) || batch.len() == PUSH_BATCH {
                if let Some(chrom_id) = batch_chrom
                    && let Err(error) =
                        writer.push_intervals(&chroms[chrom_id as usize].name, &batch)
                {
                    return Ok(format!("rejected: {error}"));
                }
                batch.clear();
                batch_chrom = Some(item.chrom_id);
            }
            batch.push((item.start, item.end, item.value));
        }
        if let Some(chrom_id) = batch_chrom
            && let Err(error) = writer.push_intervals(&chroms[chrom_id as usize].name, &batch)
        {
            return Ok(format!("rejected: {error}"));
        }
        writer.finish()?;

        let mut rewritten = BigWigReader::open(out)?;
        let written_names: Vec<String> =
            rewritten.chroms().iter().map(|chrom| chrom.name.clone()).collect();
        let mut theirs = reader.all_items()?;
        let mut ours = rewritten.all_items()?;
        let mut compared = 0u64;
        loop {
            match (theirs.next().transpose()?, ours.next().transpose()?) {
                (None, None) => break,
                (Some(source), Some(copy)) => {
                    let same = chroms[source.chrom_id as usize].name
                        == written_names[copy.chrom_id as usize]
                        && (source.start, source.end) == (copy.start, copy.end)
                        && source.value.to_bits() == copy.value.to_bits();
                    if !same {
                        return Ok(format!("differs: item {compared}: {source:?} vs {copy:?}"));
                    }
                    compared += 1;
                }
                (source, copy) => {
                    return Ok(format!("differs: after {compared} items: {source:?} vs {copy:?}"));
                }
            }
        }
        drop(ours);
        let written = rewritten.info().total_summary;
        let summary_matches = written.is_some_and(|w| {
            w.covered_bases == computed.covered_bases
                && w.min == computed.min
                && w.max == computed.max
                && (w.sum - computed.sum).abs() <= 1e-9 * computed.sum.abs().max(1.0)
        });
        if !summary_matches {
            return Ok(format!("differs: total summary {written:?} vs recomputed {computed:?}"));
        }
        Ok("identical".to_owned())
    })();
    // The rewritten file is scratch whatever happened; a failure to remove it is not a finding.
    let _ = std::fs::remove_file(out);
    outcome
}
