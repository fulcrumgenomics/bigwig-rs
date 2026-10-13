#![no_main]

//! Opens arbitrary bytes as a bigWig and, if that works, exercises every read path. Errors are
//! expected; panics, hangs and runaway allocation are bugs. Run with
//! `-rss_limit_mb=2048 -timeout=10` so the last two are reported.

use std::io::Cursor;

use arbitrary::Unstructured;
use bigwig::{BigWigReader, SummarySource};
use libfuzzer_sys::fuzz_target;

/// Items or zoom records consumed from any one iterator, so a hostile file cannot make an
/// iteration run for long.
const MAX_ITERATED: usize = 10_000;
/// Longest region whose `values` are materialised.
const MAX_VALUES_REGION: u32 = 1 << 16;
const MAX_REGIONS: usize = 4;
const MAX_ZOOM_LEVELS_QUERIED: usize = 4;
/// Trailing input bytes the regions are derived from.
const REGION_SEED_BYTES: usize = 64;

fuzz_target!(|data: &[u8]| {
    let Ok(mut reader) = BigWigReader::new(Cursor::new(data)) else { return };

    let _ = reader.info().clone();
    let chroms = reader.chroms().to_vec();
    if let Ok(items) = reader.all_items() {
        for _ in items.take(MAX_ITERATED) {}
    }
    if chroms.is_empty() {
        return;
    }

    let seed = &data[data.len().saturating_sub(REGION_SEED_BYTES)..];
    let mut unstructured = Unstructured::new(seed);
    let levels: Vec<_> =
        reader.zoom_levels().iter().copied().take(MAX_ZOOM_LEVELS_QUERIED).collect();
    for _ in 0..MAX_REGIONS {
        let Ok(chrom) = unstructured.choose(&chroms) else { return };
        let length = chrom.length;
        let start = unstructured.int_in_range(0..=length).unwrap_or(0);
        let end = unstructured.int_in_range(start..=length).unwrap_or(length);
        let n_bins = unstructured.int_in_range(1..=64u32).unwrap_or(1);
        let name = chrom.name.clone();

        if let Ok(items) = reader.items(&name, start, end) {
            for _ in items.take(MAX_ITERATED) {}
        }
        if end - start <= MAX_VALUES_REGION {
            let _ = reader.values(&name, start, end);
        }
        let _ = reader.summary(&name, start, end, SummarySource::Raw);
        let _ = reader.summary(&name, start, end, SummarySource::BestZoom);
        let _ = reader.summaries(&name, start, end, n_bins, SummarySource::Raw);
        let _ = reader.summaries(&name, start, end, n_bins, SummarySource::BestZoom);
        for &level in &levels {
            if let Ok(records) = reader.zoom_records(level, &name, start, end) {
                for _ in records.take(MAX_ITERATED) {}
            }
        }
    }
});
