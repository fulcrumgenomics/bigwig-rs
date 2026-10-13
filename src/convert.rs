//! Conversions between bigWig and the bedGraph and wig text formats.
//!
//! [`bedgraph_to_bigwig`] and [`wig_to_bigwig`] read their input once, from any [`BufRead`], with
//! the writer's bounded memory, so input must be sorted: each chromosome's data contiguous and in
//! order without overlaps. Every error the writer raises for an item is reported as an
//! [`Error::AtLine`] that wraps it with the line that holds the item.
//!
//! [`bigwig_to_bedgraph`] and [`bigwig_to_wig`] write all of a file, one chromosome or one region
//! (a [`Selection`]), chromosome by chromosome in `chromId` order.

use std::io::{BufRead, BufWriter, Read, Seek, Write};
use std::iter;

use crate::bedgraph::{self, Record, parse_record};
use crate::section::SectionKind;
use crate::text::LineReader;
use crate::wig::{self, Event};
use crate::{BigWigReader, BigWigWriter, ChromInfo, DataItem, Error, Result, WriterOptions};

/// Most fixedStep or variableStep values pushed to the writer in one call. The writer's output
/// does not depend on how a run is sliced, so this only bounds the memory a long block takes.
const BATCH_LEN: usize = 4_096;

/// What a text → bigWig conversion did.
///
/// More fields may be added; `Report::default()` gives one with every count 0.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct Report {
    /// Data items read and passed to the writer: fixedStep and variableStep values and bedGraph
    /// intervals, including those counted in `clipped` and `dropped`.
    pub items: u64,
    /// Zero-length bedGraph intervals (`start == end`), which hold no data and were skipped.
    pub skipped_empty: u64,
    /// Items extending past their chromosome's end that [`WriterOptions::clip`] truncated at the
    /// end and wrote. Without `clip` such items are errors.
    pub clipped: u64,
    /// Items starting at or past their chromosome's end, which [`WriterOptions::clip`] dropped.
    /// Without `clip` such items are errors.
    pub dropped: u64,
}

/// Which part of a bigWig to export.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Selection {
    /// Every chromosome.
    All,
    /// One whole chromosome.
    Chrom(String),
    /// The 0-based half-open range `[start, end)` of one chromosome. Items crossing either end
    /// are clipped to it.
    Region {
        /// Chromosome name.
        chrom: String,
        /// 0-based first base.
        start: u32,
        /// One past the last base; at most the chromosome's length.
        end: u32,
    },
}

impl Selection {
    /// The region to read: `None` for the whole file, otherwise a chromosome and a range in it.
    /// Fails with [`Error::UnknownChrom`] for a chromosome the file does not have.
    fn query_region<R: Read + Seek>(
        &self,
        reader: &BigWigReader<R>,
    ) -> Result<Option<(&str, u32, u32)>> {
        match self {
            Self::All => Ok(None),
            Self::Chrom(chrom) => {
                let chrom_info =
                    reader.chrom(chrom).ok_or_else(|| Error::UnknownChrom(chrom.clone()))?;
                Ok(Some((chrom, 0, chrom_info.length)))
            }
            Self::Region { chrom, start, end } => Ok(Some((chrom, *start, *end))),
        }
    }
}

/// A [`BigWigWriter`] fed from text. It counts what it passes to the writer in a [`Report`], and
/// reports each error the writer raises for an item as an [`Error::AtLine`] at the item's input
/// line.
struct ReportingWriter<W: Write + Seek> {
    writer: BigWigWriter<W>,
    /// The counts so far, except those the writer keeps.
    report: Report,
}

impl<W: Write + Seek> ReportingWriter<W> {
    /// Creates the writer, which validates `options` and `chroms`.
    fn new(options: &WriterOptions, sink: W, chroms: &[ChromInfo]) -> Result<Self> {
        let writer = options.create_with(sink, chroms)?;
        Ok(Self { writer, report: Report::default() })
    }

    /// Finishes the bigWig and returns the report, with the writer's clipping counts.
    fn finish(mut self) -> Result<Report> {
        self.report.clipped = self.writer.clipped_items();
        self.report.dropped = self.writer.dropped_items();
        self.writer.finish()?;
        Ok(self.report)
    }

    /// Pushes a bedGraph interval read from input line `line`, or skips it if it is empty and
    /// otherwise valid.
    fn push_interval(&mut self, record: &Record, line: u64) -> Result<()> {
        if record.start == record.end {
            self.writer
                .skip_empty_interval(&record.chrom, record.start, record.value)
                .map_err(|error| at_line(error, line))?;
            self.report.skipped_empty += 1;
            return Ok(());
        }
        self.writer
            .push_interval(&record.chrom, record.start, record.end, record.value)
            .map_err(|error| at_line(error, line))?;
        self.report.items += 1;
        Ok(())
    }

    /// Pushes fixedStep `values`, the first at `start`, read from input lines `lines`, one line
    /// per value. An error the writer raises names the line of the first value it rejects.
    fn push_fixed_step(
        &mut self,
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
        lines: &[u64],
    ) -> Result<()> {
        self.writer
            .push_fixed_step_indexed(chrom, start, step, span, values)
            .map_err(|rejected| at_line(rejected.error, lines[rejected.index]))?;
        self.report.items += values.len() as u64;
        Ok(())
    }

    /// Pushes variableStep `items`, each covering `span` bases, read from input lines `lines`,
    /// one line per item. An error the writer raises names the line of the first item it rejects.
    fn push_variable_step(
        &mut self,
        chrom: &str,
        span: u32,
        items: &[(u32, f32)],
        lines: &[u64],
    ) -> Result<()> {
        self.writer
            .push_variable_step_indexed(chrom, span, items)
            .map_err(|rejected| at_line(rejected.error, lines[rejected.index]))?;
        self.report.items += items.len() as u64;
        Ok(())
    }
}

/// The wig block being read, and its values not yet pushed.
#[derive(Default)]
struct WigBlock {
    header: BlockHeader,
    chrom: String,
    fixed_values: Vec<f32>,
    variable_items: Vec<(u32, f32)>,
    /// The input line of each pending value.
    lines: Vec<u64>,
}

impl WigBlock {
    /// Pushes the pending values and starts a fixedStep block.
    fn start_fixed_step<W: Write + Seek>(
        &mut self,
        output: &mut ReportingWriter<W>,
        chrom: &str,
        start: u32,
        step: u32,
        span: u32,
    ) -> Result<()> {
        self.flush(output)?;
        self.chrom.clear();
        self.chrom.push_str(chrom);
        self.header = BlockHeader::FixedStep { next_start: start, step, span };
        Ok(())
    }

    /// Pushes the pending values and starts a variableStep block.
    fn start_variable_step<W: Write + Seek>(
        &mut self,
        output: &mut ReportingWriter<W>,
        chrom: &str,
        span: u32,
    ) -> Result<()> {
        self.flush(output)?;
        self.chrom.clear();
        self.chrom.push_str(chrom);
        self.header = BlockHeader::VariableStep { span };
        Ok(())
    }

    /// Adds a value of the fixedStep block, read from input line `line`.
    fn add_fixed_value<W: Write + Seek>(
        &mut self,
        output: &mut ReportingWriter<W>,
        value: f32,
        line: u64,
    ) -> Result<()> {
        self.fixed_values.push(value);
        self.lines.push(line);
        if self.lines.len() == BATCH_LEN { self.flush(output) } else { Ok(()) }
    }

    /// Adds an item of the variableStep block, read from input line `line`.
    fn add_variable_value<W: Write + Seek>(
        &mut self,
        output: &mut ReportingWriter<W>,
        start: u32,
        value: f32,
        line: u64,
    ) -> Result<()> {
        self.variable_items.push((start, value));
        self.lines.push(line);
        if self.lines.len() == BATCH_LEN { self.flush(output) } else { Ok(()) }
    }

    /// Pushes the pending values and closes the block, as a bedGraph line or the end of input
    /// does.
    fn end<W: Write + Seek>(&mut self, output: &mut ReportingWriter<W>) -> Result<()> {
        self.flush(output)?;
        self.header = BlockHeader::None;
        Ok(())
    }

    /// Pushes the pending values, keeping the block open.
    fn flush<W: Write + Seek>(&mut self, output: &mut ReportingWriter<W>) -> Result<()> {
        // The wig reader yields values only inside a block of their own kind, so values pending
        // for the other kind never exist.
        match &mut self.header {
            BlockHeader::FixedStep { next_start, step, span } if !self.fixed_values.is_empty() => {
                output.push_fixed_step(
                    &self.chrom,
                    *next_start,
                    *step,
                    *span,
                    &self.fixed_values,
                    &self.lines,
                )?;
                *next_start = fixed_step_position(*next_start, *step, self.fixed_values.len());
            }
            BlockHeader::VariableStep { span } if !self.variable_items.is_empty() => {
                output.push_variable_step(&self.chrom, *span, &self.variable_items, &self.lines)?;
            }
            _ => {}
        }
        self.fixed_values.clear();
        self.variable_items.clear();
        self.lines.clear();
        Ok(())
    }
}

/// The declaration of the wig block being read.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum BlockHeader {
    /// No block is open: before the first declaration, or after a bedGraph line.
    #[default]
    None,
    /// A fixedStep block, whose next pending value starts at `next_start`.
    FixedStep { next_start: u32, step: u32, span: u32 },
    /// A variableStep block.
    VariableStep { span: u32 },
}

/// Writes bigWig sections as wig text, joining consecutive sections of one fixedStep or
/// variableStep run under a single header.
struct WigEmitter<W: Write> {
    writer: wig::Writer<W>,
    /// The file's chromosomes, indexed by `chromId`.
    chroms: Vec<ChromInfo>,
    /// The selected range of a single chromosome, which items are clipped to; `None` for all.
    range: Option<(u32, u32)>,
    /// The block the last value was written to.
    open: OpenBlock,
}

impl<W: Write> WigEmitter<W> {
    /// Creates an emitter writing to `out`, clipping items to `range` when it is given.
    fn new(out: W, chroms: Vec<ChromInfo>, range: Option<(u32, u32)>) -> Self {
        Self { writer: wig::Writer::new(out), chroms, range, open: OpenBlock::None }
    }

    /// Flushes the output.
    fn finish(self) -> Result<()> {
        self.writer.into_inner()?;
        Ok(())
    }

    /// Writes the items of a section of type `kind`.
    fn write_section(&mut self, kind: SectionKind, items: &[DataItem]) -> Result<()> {
        let range = self.range;
        let clipped = || items.iter().filter_map(move |&item| clip_to_range(item, range));
        match kind {
            SectionKind::BedGraph => self.write_bed_graph_section(items),
            SectionKind::VariableStep { .. } => {
                clipped().try_for_each(|item| self.write_variable(item))
            }
            SectionKind::FixedStep { step, .. } => {
                clipped().try_for_each(|item| self.write_fixed(item, step))
            }
        }
    }

    /// Writes a bedGraph section's items as bedGraph lines, after kent's comment giving their
    /// extent. Nothing is written when no item is in the selected range.
    fn write_bed_graph_section(&mut self, items: &[DataItem]) -> Result<()> {
        let range = self.range;
        let mut clipped = items.iter().filter_map(|&item| clip_to_range(item, range));
        let Some(first) = clipped.next() else { return Ok(()) };
        let last_end = items
            .iter()
            .rev()
            .find_map(|&item| clip_to_range(item, range))
            .map_or(first.end, |last| last.end);
        let chrom = chrom_name(&self.chroms, first.chrom_id);
        self.writer.comment(&format!("bedGraph section {chrom}:{}-{last_end}", first.start))?;
        for item in iter::once(first).chain(clipped) {
            self.writer.bedgraph_record(chrom, item.start, item.end, item.value)?;
        }
        // A bedGraph line ends the open block, so the next value needs a new header.
        self.open = OpenBlock::None;
        Ok(())
    }

    /// Writes a variableStep item, under a new header unless it can join the open block.
    fn write_variable(&mut self, item: DataItem) -> Result<()> {
        // A clipped item covers fewer bases than the section's span, so it gets its own header.
        let span = item.end - item.start;
        let block = OpenBlock::VariableStep { chrom_id: item.chrom_id, span };
        if self.open != block {
            self.writer.variable_step(chrom_name(&self.chroms, item.chrom_id), span)?;
            self.open = block;
        }
        self.writer.variable_value(item.start, item.value)
    }

    /// Writes a fixedStep item of a section with `step`, under a new header unless it continues
    /// the open block's run.
    fn write_fixed(&mut self, item: DataItem, step: u32) -> Result<()> {
        let span = item.end - item.start;
        let continued =
            OpenBlock::FixedStep { chrom_id: item.chrom_id, step, span, next_start: item.start };
        if self.open != continued {
            let chrom = chrom_name(&self.chroms, item.chrom_id);
            self.writer.fixed_step(chrom, item.start, step, span)?;
        }
        self.writer.fixed_value(item.value)?;
        // No item starts at u32::MAX, so a saturated start continues nothing.
        let next_start = item.start.saturating_add(step);
        self.open = OpenBlock::FixedStep { chrom_id: item.chrom_id, step, span, next_start };
        Ok(())
    }
}

/// The wig block open in the output, as what a value must match to be written under it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OpenBlock {
    /// No header yet, or bedGraph lines ended the last block.
    None,
    /// A fixedStep block whose next value would start at `next_start`.
    FixedStep { chrom_id: u32, step: u32, span: u32, next_start: u32 },
    /// A variableStep block.
    VariableStep { chrom_id: u32, span: u32 },
}

/// Streams bedGraph text into a bigWig written to `sink`, and reports what was converted.
///
/// The input is read once, so it can be a pipe. Chromosomes may come in any order, but each
/// one's intervals must be contiguous and sorted by start without overlaps, as kent
/// `bedGraphToBigWig` 2.10 and later require. Zero-length intervals hold no data: one that is
/// otherwise valid is skipped and counted in [`Report::skipped_empty`], where kent writes it, and
/// one that isn't (on an unknown or finished chromosome, starting before the previous interval
/// ends, past the chromosome end or with a non-finite value) is an error, as in kent. Later
/// intervals may not start before a skipped one. Items past a chromosome's end are errors unless
/// [`WriterOptions::clip`] is set.
///
/// `sink` must be empty (see [`WriterOptions::create_with`]); pass `&mut sink` to keep it, for
/// example a `Cursor<Vec<u8>>`. Fails with [`Error::InvalidOption`] or
/// [`Error::InvalidChromList`] before reading anything, with [`Error::Io`] for a sink that isn't
/// empty, with [`Error::Parse`] for a line the bedGraph reader cannot parse, with
/// [`Error::AtLine`] wrapping the writer's error for a line the writer rejects (unsorted or
/// overlapping, an unknown chromosome, a NaN value, ...), both carrying the 1-based line number,
/// and with [`Error::Io`]. After an error the sink holds an unfinished file that readers reject.
///
/// # Example
///
/// ```
/// use std::io::Cursor;
///
/// use bigwig::convert::bedgraph_to_bigwig;
/// use bigwig::{BigWigReader, ChromInfo, Error, WriterOptions};
///
/// let chroms = ChromInfo::read_chrom_sizes("chr1\t1000\n".as_bytes())?;
/// let text = "chr1\t0\t100\t1.5\nchr1\t100\t100\t9\nchr1\t150\t200\t2\n";
/// let mut bigwig = Cursor::new(Vec::new());
/// let report = bedgraph_to_bigwig(text.as_bytes(), &chroms, &WriterOptions::new(), &mut bigwig)?;
/// // The zero-length interval at 100 holds no data.
/// assert_eq!((report.items, report.skipped_empty), (2, 1));
/// let mut reader = BigWigReader::new(bigwig)?;
/// assert_eq!(reader.values("chr1", 99, 100)?, [1.5]);
///
/// // An error the writer raises is wrapped with the input line.
/// let unsorted = "chr1\t50\t60\t1\nchr1\t0\t10\t2\n";
/// let sink = Cursor::new(Vec::new());
/// let error =
///     bedgraph_to_bigwig(unsorted.as_bytes(), &chroms, &WriterOptions::new(), sink).unwrap_err();
/// assert_eq!(
///     error.to_string(),
///     "line 2: chr1:0 starts before the previous item ends at 60 (unsorted or overlapping)"
/// );
/// let Error::AtLine { line, error } = error else { panic!("expected AtLine") };
/// assert_eq!(line, 2);
/// assert!(matches!(*error, Error::UnsortedOrOverlapping { previous_end: 60, .. }));
/// # Ok::<(), bigwig::Error>(())
/// ```
pub fn bedgraph_to_bigwig<R: BufRead, W: Write + Seek>(
    input: R,
    chroms: &[ChromInfo],
    options: &WriterOptions,
    sink: W,
) -> Result<Report> {
    let mut output = ReportingWriter::new(options, sink, chroms)?;
    let mut lines = LineReader::new(input);
    let mut record = Record::default();
    while lines.next_line()? {
        parse_record(lines.line(), lines.line_number(), &mut record)?;
        output.push_interval(&record, lines.line_number())?;
    }
    output.finish()
}

/// Streams wig text into a bigWig written to `sink`, keeping each block's section type:
/// `fixedStep` blocks become fixedStep sections, `variableStep` blocks variableStep sections and
/// bedGraph lines bedGraph sections. Reports what was converted.
///
/// The input is read once, so it can be a pipe, and must be sorted: each chromosome's items
/// contiguous and in order without overlaps (kent `wigToBigWig` sorts in memory instead).
/// `sink` must be empty. Zero-length bedGraph lines, clipping and errors are as for
/// [`bedgraph_to_bigwig`]; an error the writer raises names the line of the value it rejected.
pub fn wig_to_bigwig<R: BufRead, W: Write + Seek>(
    input: R,
    chroms: &[ChromInfo],
    options: &WriterOptions,
    sink: W,
) -> Result<Report> {
    let mut output = ReportingWriter::new(options, sink, chroms)?;
    let mut reader = wig::Reader::new(input);
    let mut event = Event::FixedValue(0.0);
    let mut block = WigBlock::default();
    while reader.read_event(&mut event)? {
        let line = reader.line_number();
        match &event {
            Event::FixedStep { chrom, start, step, span } => {
                block.start_fixed_step(&mut output, chrom, *start, *step, *span)?;
            }
            Event::VariableStep { chrom, span } => {
                block.start_variable_step(&mut output, chrom, *span)?;
            }
            Event::FixedValue(value) => block.add_fixed_value(&mut output, *value, line)?,
            Event::VariableValue { start, value } => {
                block.add_variable_value(&mut output, *start, *value, line)?;
            }
            Event::BedGraph(record) => {
                block.end(&mut output)?;
                output.push_interval(record, line)?;
            }
        }
    }
    block.end(&mut output)?;
    output.finish()
}

/// Writes the items of `reader` in `selection` to `out` as bedGraph lines, chromosome by
/// chromosome in `chromId` order. Items crossing the ends of a [`Selection::Region`] are clipped
/// to it. With `merge_equal`, adjacent items of equal value on one chromosome are written as one
/// line, as kent `bigWigToBedGraph` always does; without it each stored item is a line.
///
/// Values are written in the shortest form that parses back to the same `f32`, where kent writes
/// six significant digits. `out` is buffered here. Fails with [`Error::UnknownChrom`] or
/// [`Error::InvalidRegion`] for a selection the file does not have, with what reading items may
/// raise ([`Error::Corrupt`], [`Error::Decompress`], [`Error::Io`]), and with
/// [`Error::InvalidItem`] for a chromosome name bedGraph cannot hold.
///
/// # Example
///
/// ```
/// use std::io::Cursor;
///
/// use bigwig::convert::{Selection, bigwig_to_bedgraph};
/// use bigwig::{BigWigReader, ChromInfo, WriterOptions};
///
/// let chroms = [ChromInfo::new("chr1", 1_000)];
/// let mut writer = WriterOptions::new().create_with(Cursor::new(Vec::new()), &chroms)?;
/// writer.push_fixed_step("chr1", 100, 10, 10, &[1.0, 1.0, 2.0])?;
/// let mut reader = BigWigReader::new(writer.finish()?)?;
///
/// let mut merged = Vec::new();
/// bigwig_to_bedgraph(&mut reader, &Selection::All, true, &mut merged)?;
/// assert_eq!(merged, b"chr1\t100\t120\t1\nchr1\t120\t130\t2\n");
///
/// let region = Selection::Region { chrom: "chr1".to_string(), start: 105, end: 125 };
/// let mut clipped = Vec::new();
/// bigwig_to_bedgraph(&mut reader, &region, false, &mut clipped)?;
/// assert_eq!(clipped, b"chr1\t105\t110\t1\nchr1\t110\t120\t1\nchr1\t120\t125\t2\n");
/// # Ok::<(), bigwig::Error>(())
/// ```
pub fn bigwig_to_bedgraph<R: Read + Seek, W: Write>(
    reader: &mut BigWigReader<R>,
    selection: &Selection,
    merge_equal: bool,
    out: W,
) -> Result<()> {
    let chroms = reader.chroms().to_vec();
    let region = selection.query_region(reader)?;
    let range = region.map(|(_, start, end)| (start, end));
    let items = match region {
        None => reader.all_items()?,
        Some((chrom, start, end)) => reader.items(chrom, start, end)?,
    };
    let mut writer = bedgraph::Writer::new(BufWriter::new(out));
    let mut write = |item: DataItem| {
        writer.write_record(chrom_name(&chroms, item.chrom_id), item.start, item.end, item.value)
    };
    // With merging, the run of equal values not yet written.
    let mut run: Option<DataItem> = None;
    for item in items {
        let Some(item) = clip_to_range(item?, range) else { continue };
        if !merge_equal {
            write(item)?;
            continue;
        }
        match &mut run {
            Some(run) if extends_run(run, &item) => run.end = item.end,
            _ => {
                if let Some(finished) = run.replace(item) {
                    write(finished)?;
                }
            }
        }
    }
    if let Some(finished) = run {
        write(finished)?;
    }
    writer.into_inner()?;
    Ok(())
}

/// Writes the items of `reader` in `selection` to `out` as wig, chromosome by chromosome in
/// `chromId` order, keeping section types.
///
/// Consecutive fixedStep or variableStep sections that continue one run are written under a
/// single header, and a header is written only when a value follows it, where kent
/// `bigWigToWig` writes one per stored section. bedGraph sections become bedGraph lines, each
/// section preceded by a `#bedGraph section chrom:start-end` comment giving the 0-based extent of
/// its lines, as kent writes. Headers carry 1-based positions and an explicit `span`.
///
/// Items crossing the ends of a [`Selection::Region`] are clipped to it. A clipped fixedStep or
/// variableStep item covers fewer bases than its neighbours, so it is written under a header of
/// its own. `out` is buffered here. Fails as [`bigwig_to_bedgraph`] does.
pub fn bigwig_to_wig<R: Read + Seek, W: Write>(
    reader: &mut BigWigReader<R>,
    selection: &Selection,
    out: W,
) -> Result<()> {
    let region = selection.query_region(reader)?;
    let range = region.map(|(_, start, end)| (start, end));
    let mut emitter = WigEmitter::new(BufWriter::new(out), reader.chroms().to_vec(), range);
    reader.for_each_section(region, |kind, items| emitter.write_section(kind, items))?;
    emitter.finish()
}

/// `error` wrapped in an [`Error::AtLine`] at input line `line` if it is the writer rejecting
/// data; otherwise unchanged. An `AtLine` is not a data-validation error, so none is wrapped twice.
fn at_line(error: Error, line: u64) -> Error {
    if error.is_data_validation() { Error::AtLine { line, error: Box::new(error) } } else { error }
}

/// Start of value `index` of a fixedStep run that starts at `start`. The wig reader rejects a
/// value whose extent passes `u32::MAX`, so saturation never reaches a value that is pushed.
fn fixed_step_position(start: u32, step: u32, index: usize) -> u32 {
    let offset = u32::try_from(index).unwrap_or(u32::MAX).saturating_mul(step);
    start.saturating_add(offset)
}

/// `item` cut to the range `[start, end)`, or `None` if none of it is inside. With no range the
/// item is returned whole.
fn clip_to_range(item: DataItem, range: Option<(u32, u32)>) -> Option<DataItem> {
    let Some((range_start, range_end)) = range else { return Some(item) };
    let start = item.start.max(range_start);
    let end = item.end.min(range_end);
    (start < end).then_some(DataItem { start, end, ..item })
}

/// Whether `item` directly follows `run` on the same chromosome with the same value.
#[allow(clippy::float_cmp, reason = "kent's bigWigToBedGraph merges values that compare equal")]
fn extends_run(run: &DataItem, item: &DataItem) -> bool {
    run.chrom_id == item.chrom_id && run.end == item.start && run.value == item.value
}

/// Name of the chromosome with `chrom_id`. Readers only return items whose `chromId` is in their
/// chromosome list; any other gets an empty name, which the text writers reject.
fn chrom_name(chroms: &[ChromInfo], chrom_id: u32) -> &str {
    usize::try_from(chrom_id)
        .ok()
        .and_then(|index| chroms.get(index))
        .map_or("", |chrom| chrom.name.as_str())
}

#[cfg(test)]
mod tests {
    use std::fmt::{Debug, Write as _};
    use std::io::{self, Cursor, SeekFrom};

    use super::*;
    use crate::{ChromIdOrder, Compression};

    /// The chromosomes every test converts against. `chr1` is long enough for several batches of
    /// per-base values, clipping tests use the short `chr2` and `chr3`, and no export test gives
    /// `chr3` data.
    fn chroms() -> Vec<ChromInfo> {
        vec![
            ChromInfo::new("chr1", 20_000),
            ChromInfo::new("chr2", 500),
            ChromInfo::new("chr3", 300),
        ]
    }

    fn bedgraph_with(text: &str, options: &WriterOptions) -> Result<(Report, Vec<u8>)> {
        let mut sink = Cursor::new(Vec::new());
        let report = bedgraph_to_bigwig(text.as_bytes(), &chroms(), options, &mut sink)?;
        Ok((report, sink.into_inner()))
    }

    fn bedgraph(text: &str) -> (Report, Vec<u8>) {
        bedgraph_with(text, &WriterOptions::new()).unwrap()
    }

    fn wig_with(text: &str, options: &WriterOptions) -> Result<(Report, Vec<u8>)> {
        let mut sink = Cursor::new(Vec::new());
        let report = wig_to_bigwig(text.as_bytes(), &chroms(), options, &mut sink)?;
        Ok((report, sink.into_inner()))
    }

    fn wig(text: &str) -> (Report, Vec<u8>) {
        wig_with(text, &WriterOptions::new()).unwrap()
    }

    /// A bigWig over [`chroms`] written by `fill`.
    fn build(
        options: &WriterOptions,
        fill: impl FnOnce(&mut BigWigWriter<Cursor<Vec<u8>>>) -> Result<()>,
    ) -> Vec<u8> {
        let mut writer = options.create_with(Cursor::new(Vec::new()), &chroms()).unwrap();
        fill(&mut writer).unwrap();
        writer.finish().unwrap().into_inner()
    }

    fn open(bytes: &[u8]) -> BigWigReader<Cursor<Vec<u8>>> {
        BigWigReader::new(Cursor::new(bytes.to_vec())).unwrap()
    }

    /// Every item of a bigWig as `(chrom, start, end, value)`.
    fn items(bytes: &[u8]) -> Vec<(String, u32, u32, f32)> {
        let mut reader = open(bytes);
        let chroms = reader.chroms().to_vec();
        reader
            .all_items()
            .unwrap()
            .map(|item| {
                let item = item.unwrap();
                let chrom = &chroms[usize::try_from(item.chrom_id).unwrap()].name;
                (chrom.clone(), item.start, item.end, item.value)
            })
            .collect()
    }

    fn item(chrom: &str, start: u32, end: u32, value: f32) -> (String, u32, u32, f32) {
        (chrom.to_string(), start, end, value)
    }

    /// The kind of every data section of a bigWig, in order.
    fn section_kinds(bytes: &[u8]) -> Vec<SectionKind> {
        let mut kinds = Vec::new();
        open(bytes)
            .for_each_section(None, |kind, _| {
                kinds.push(kind);
                Ok(())
            })
            .unwrap();
        kinds
    }

    fn try_to_bedgraph(bytes: &[u8], selection: &Selection, merge_equal: bool) -> Result<String> {
        let mut out = Vec::new();
        bigwig_to_bedgraph(&mut open(bytes), selection, merge_equal, &mut out)?;
        Ok(String::from_utf8(out).unwrap())
    }

    fn to_bedgraph(bytes: &[u8], selection: &Selection, merge_equal: bool) -> String {
        try_to_bedgraph(bytes, selection, merge_equal).unwrap()
    }

    fn try_to_wig(bytes: &[u8], selection: &Selection) -> Result<String> {
        let mut out = Vec::new();
        bigwig_to_wig(&mut open(bytes), selection, &mut out)?;
        Ok(String::from_utf8(out).unwrap())
    }

    fn to_wig(bytes: &[u8], selection: &Selection) -> String {
        try_to_wig(bytes, selection).unwrap()
    }

    fn region(chrom: &str, start: u32, end: u32) -> Selection {
        Selection::Region { chrom: chrom.to_string(), start, end }
    }

    /// The line and message of the parse error in `result`.
    fn parse_failure<T: Debug>(result: Result<T>) -> (u64, String) {
        match result {
            Err(Error::Parse { line, message }) => (line, message),
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    /// The line and the writer's error of the [`Error::AtLine`] in `result`.
    fn writer_rejection<T: Debug>(result: Result<T>) -> (u64, Error) {
        match result {
            Err(Error::AtLine { line, error }) => (line, *error),
            other => panic!("expected an error at a line, got {other:?}"),
        }
    }

    /// `count` one-base bedGraph intervals on `chr1`, starting at 0.
    fn one_base_intervals(count: u32) -> String {
        let mut text = String::new();
        for start in 0..count {
            writeln!(text, "chr1\t{start}\t{}\t1", start + 1).unwrap();
        }
        text
    }

    /// Options under which every item is a section of its own, written uncompressed, so the
    /// whole of `chr1` as one-base items overflows the writer's output buffer while pushing.
    fn options_writing_to_the_sink_while_pushing() -> WriterOptions {
        WriterOptions::new().items_per_section(1).compression(Compression::None)
    }

    /// A sink that accepts nothing.
    struct FailingSink;

    impl Write for FailingSink {
        fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("sink refuses writes"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    impl Seek for FailingSink {
        fn seek(&mut self, _pos: SeekFrom) -> io::Result<u64> {
            Ok(0)
        }
    }

    // ---- bedgraph_to_bigwig ------------------------------------------------------------------

    #[test]
    fn bedgraph_to_bigwig_writes_every_interval() {
        let (_, bytes) = bedgraph("chr1\t0\t10\t1.5\nchr1\t20\t30\t-2\nchr2\t5\t6\t3\n");
        assert_eq!(
            items(&bytes),
            [item("chr1", 0, 10, 1.5), item("chr1", 20, 30, -2.0), item("chr2", 5, 6, 3.0)]
        );
    }

    #[test]
    fn bedgraph_to_bigwig_writes_bedgraph_sections() {
        let (_, bytes) = bedgraph("chr1\t0\t10\t1\nchr2\t0\t10\t2\n");
        assert_eq!(section_kinds(&bytes), [SectionKind::BedGraph, SectionKind::BedGraph]);
    }

    #[test]
    fn bedgraph_to_bigwig_reports_the_items_written() {
        let (report, _) =
            bedgraph("track type=bedGraph\nchr1\t0\t10\t1\n# note\nchr1\t10\t20\t2\n");
        assert_eq!(report, Report { items: 2, ..Report::default() });
    }

    #[test]
    fn bedgraph_to_bigwig_accepts_chromosomes_in_any_order() {
        let (_, bytes) = bedgraph("chr2\t0\t10\t1\nchr1\t0\t10\t2\n");
        assert_eq!(items(&bytes), [item("chr2", 0, 10, 1.0), item("chr1", 0, 10, 2.0)]);
    }

    #[test]
    fn bedgraph_to_bigwig_skips_and_counts_zero_length_intervals() {
        let (report, bytes) = bedgraph("chr1\t0\t10\t1\nchr1\t10\t10\t5\nchr1\t20\t30\t2\n");
        assert_eq!(report, Report { items: 2, skipped_empty: 1, ..Report::default() });
        assert_eq!(items(&bytes), [item("chr1", 0, 10, 1.0), item("chr1", 20, 30, 2.0)]);
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_zero_length_interval_on_an_unknown_chromosome_at_its_line() {
        let (line, error) = writer_rejection(bedgraph_with(
            "chr1\t0\t10\t1\nchrX\t5\t5\t2\n",
            &WriterOptions::new(),
        ));
        assert_eq!(line, 2);
        assert!(matches!(error, Error::UnknownChrom(ref chrom) if chrom == "chrX"), "{error:?}");
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_zero_length_interval_with_a_nan_value_at_its_line() {
        let (line, error) = writer_rejection(bedgraph_with(
            "chr1\t0\t10\t1\nchr1\t20\t20\tnan\n",
            &WriterOptions::new(),
        ));
        assert_eq!(line, 2);
        assert!(
            matches!(error, Error::NonFiniteValue { value, .. } if value.is_nan()),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_zero_length_interval_inside_the_previous_one_at_its_line() {
        let (line, error) = writer_rejection(bedgraph_with(
            "chr1\t0\t10\t1\nchr1\t5\t5\t2\n",
            &WriterOptions::new(),
        ));
        assert_eq!(line, 2);
        assert!(
            matches!(error, Error::UnsortedOrOverlapping { start: 5, previous_end: 10, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_zero_length_interval_on_a_finished_chromosome_at_its_line() {
        let text = "chr1\t0\t10\t1\nchr2\t0\t10\t2\nchr1\t20\t20\t3\n";
        let (line, error) = writer_rejection(bedgraph_with(text, &WriterOptions::new()));
        assert_eq!(line, 3);
        assert!(
            matches!(error, Error::ChromAlreadyFinished(ref chrom) if chrom == "chr1"),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_zero_length_interval_past_the_chromosome_end_at_its_line() {
        let (line, error) =
            writer_rejection(bedgraph_with("chr2\t510\t510\t1\n", &WriterOptions::new()));
        assert_eq!(line, 1);
        assert!(
            matches!(error, Error::OutOfBounds { start: 510, end: 510, length: 500, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_skips_a_zero_length_interval_past_the_chromosome_end_with_clip() {
        let options = WriterOptions::new().clip(true);
        let (report, bytes) =
            bedgraph_with("chr2\t0\t10\t1\nchr2\t510\t510\t2\n", &options).unwrap();
        assert_eq!(report.skipped_empty, 1);
        assert_eq!(items(&bytes), [item("chr2", 0, 10, 1.0)]);
    }

    #[test]
    fn bedgraph_to_bigwig_rejects_an_interval_starting_before_a_skipped_zero_length_one() {
        let text = "chr1\t0\t10\t1\nchr1\t25\t25\t2\nchr1\t22\t30\t3\n";
        let (line, error) = writer_rejection(bedgraph_with(text, &WriterOptions::new()));
        assert_eq!(line, 3);
        assert!(
            matches!(error, Error::UnsortedOrOverlapping { start: 22, previous_end: 25, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_counts_an_interval_truncated_by_clip() {
        let options = WriterOptions::new().clip(true);
        let (report, bytes) = bedgraph_with("chr2\t490\t510\t1\n", &options).unwrap();
        assert_eq!(report, Report { items: 1, clipped: 1, ..Report::default() });
        assert_eq!(items(&bytes), [item("chr2", 490, 500, 1.0)]);
    }

    #[test]
    fn bedgraph_to_bigwig_counts_an_interval_dropped_by_clip() {
        let options = WriterOptions::new().clip(true);
        let (report, bytes) =
            bedgraph_with("chr2\t0\t10\t1\nchr2\t500\t510\t2\n", &options).unwrap();
        assert_eq!(report, Report { items: 2, dropped: 1, ..Report::default() });
        assert_eq!(items(&bytes), [item("chr2", 0, 10, 1.0)]);
    }

    #[test]
    fn bedgraph_to_bigwig_reports_unsorted_input_at_its_line() {
        let (line, error) = writer_rejection(bedgraph_with(
            "chr1\t50\t60\t1\nchr1\t60\t70\t2\nchr1\t0\t10\t3\n",
            &WriterOptions::new(),
        ));
        assert_eq!(line, 3);
        assert!(
            matches!(
                error,
                Error::UnsortedOrOverlapping { ref chrom, start: 0, previous_end: 70 }
                    if chrom == "chr1"
            ),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_writer_rejection_message_is_the_line_then_the_writers_message() {
        let result = bedgraph_with("chr1\t50\t60\t1\nchr1\t0\t10\t2\n", &WriterOptions::new());
        assert_eq!(
            result.unwrap_err().to_string(),
            "line 2: chr1:0 starts before the previous item ends at 60 (unsorted or overlapping)"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_an_unknown_chromosome_at_its_line() {
        let (line, error) = writer_rejection(bedgraph_with(
            "chr1\t0\t10\t1\nchrX\t0\t10\t2\n",
            &WriterOptions::new(),
        ));
        assert_eq!(line, 2);
        assert!(matches!(error, Error::UnknownChrom(ref chrom) if chrom == "chrX"), "{error:?}");
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_nan_value_at_its_line() {
        let (line, error) = writer_rejection(bedgraph_with(
            "chr1\t0\t10\t1\nchr1\t10\t20\tnan\n",
            &WriterOptions::new(),
        ));
        assert_eq!(line, 2);
        assert!(
            matches!(error, Error::NonFiniteValue { position: 10, value, .. } if value.is_nan()),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_returning_chromosome_at_its_line() {
        let text = "chr1\t0\t10\t1\nchr2\t0\t10\t2\nchr1\t20\t30\t3\n";
        let (line, error) = writer_rejection(bedgraph_with(text, &WriterOptions::new()));
        assert_eq!(line, 3);
        assert!(
            matches!(error, Error::ChromAlreadyFinished(ref chrom) if chrom == "chr1"),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_chromosome_out_of_listed_order_at_its_line() {
        let options = WriterOptions::new().chrom_id_order(ChromIdOrder::Listed);
        let (line, error) =
            writer_rejection(bedgraph_with("chr2\t0\t10\t1\nchr1\t0\t10\t2\n", &options));
        assert_eq!(line, 2);
        assert!(
            matches!(
                error,
                Error::ChromOutOfOrder { ref chrom, ref current } if chrom == "chr1" && current == "chr2"
            ),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_an_interval_past_the_chromosome_end_at_its_line() {
        let (line, error) =
            writer_rejection(bedgraph_with("chr2\t490\t510\t1\n", &WriterOptions::new()));
        assert_eq!(line, 1);
        assert!(
            matches!(error, Error::OutOfBounds { start: 490, end: 510, length: 500, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_reports_a_start_after_the_end_at_its_line() {
        let (line, error) =
            writer_rejection(bedgraph_with("chr1\t20\t10\t1\n", &WriterOptions::new()));
        assert_eq!(line, 1);
        assert!(
            matches!(error, Error::InvalidItem { start: 20, reason, .. }
                if reason.contains("start must be below its end")),
            "{error:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_counts_comment_and_blank_lines_in_line_numbers() {
        let text = "# header\n\nchr1\t0\t10\t1\n  # indented\nchr1\t5\t15\t2\n";
        assert_eq!(writer_rejection(bedgraph_with(text, &WriterOptions::new())).0, 5);
    }

    #[test]
    fn bedgraph_to_bigwig_passes_reader_parse_errors_through() {
        let (line, message) =
            parse_failure(bedgraph_with("chr1\t0\t10\t1\nchr1\t10\t20\n", &WriterOptions::new()));
        assert_eq!(line, 2);
        assert_eq!(message, "expected 4 fields, found 3");
    }

    #[test]
    fn bedgraph_to_bigwig_rejects_an_invalid_chromosome_list_before_reading() {
        let mut sink = Cursor::new(Vec::new());
        let result =
            bedgraph_to_bigwig(&b"chr1\t0\t10\t1\n"[..], &[], &WriterOptions::new(), &mut sink);
        assert!(matches!(result, Err(Error::InvalidChromList(_))), "{result:?}");
    }

    #[test]
    fn bedgraph_to_bigwig_rejects_a_sink_that_is_not_empty_before_reading() {
        let sink = Cursor::new(vec![0; 4]);
        let result =
            bedgraph_to_bigwig(&b"chr1\t0\t10\t1\n"[..], &chroms(), &WriterOptions::new(), sink);
        assert!(
            matches!(result, Err(Error::Io(ref error)) if error.kind() == io::ErrorKind::InvalidInput),
            "{result:?}"
        );
    }

    #[test]
    fn bedgraph_to_bigwig_returns_sink_failures_as_io_errors() {
        let input = &b"chr1\t0\t10\t1\n"[..];
        let result = bedgraph_to_bigwig(input, &chroms(), &WriterOptions::new(), FailingSink);
        assert!(matches!(result, Err(Error::Io(_))), "{result:?}");
    }

    #[test]
    fn bedgraph_to_bigwig_returns_a_sink_failure_while_pushing_unwrapped() {
        // The unparseable last line would fail the conversion if the sink failure came only at
        // the end.
        let text = format!("{}chr1\t0\t10\n", one_base_intervals(20_000));
        let options = options_writing_to_the_sink_while_pushing();
        let result = bedgraph_to_bigwig(text.as_bytes(), &chroms(), &options, FailingSink);
        assert!(matches!(result, Err(Error::Io(_))), "{result:?}");
    }

    // ---- wig_to_bigwig -----------------------------------------------------------------------

    #[test]
    fn wig_to_bigwig_writes_fixed_step_values_at_their_positions() {
        let (_, bytes) = wig("fixedStep chrom=chr1 start=101 step=10 span=5\n1\n2\n3\n");
        assert_eq!(
            items(&bytes),
            [item("chr1", 100, 105, 1.0), item("chr1", 110, 115, 2.0), item("chr1", 120, 125, 3.0)]
        );
    }

    #[test]
    fn wig_to_bigwig_keeps_fixed_step_blocks_as_fixed_step_sections() {
        let (_, bytes) = wig("fixedStep chrom=chr1 start=101 step=10 span=5\n1\n2\n3\n");
        assert_eq!(section_kinds(&bytes), [SectionKind::FixedStep { step: 10, span: 5 }]);
    }

    #[test]
    fn wig_to_bigwig_writes_variable_step_items_at_their_positions() {
        let (_, bytes) = wig("variableStep chrom=chr1 span=5\n11\t1\n101\t2\n");
        assert_eq!(items(&bytes), [item("chr1", 10, 15, 1.0), item("chr1", 100, 105, 2.0)]);
    }

    #[test]
    fn wig_to_bigwig_keeps_variable_step_blocks_as_variable_step_sections() {
        let (_, bytes) = wig("variableStep chrom=chr1 span=5\n11\t1\n101\t2\n");
        assert_eq!(section_kinds(&bytes), [SectionKind::VariableStep { span: 5 }]);
    }

    #[test]
    fn wig_to_bigwig_writes_bedgraph_lines_as_bedgraph_sections() {
        let (_, bytes) = wig("chr1\t0\t10\t1\nchr1\t10\t20\t2\n");
        assert_eq!(section_kinds(&bytes), [SectionKind::BedGraph]);
        assert_eq!(items(&bytes), [item("chr1", 0, 10, 1.0), item("chr1", 10, 20, 2.0)]);
    }

    #[test]
    fn wig_to_bigwig_keeps_the_section_type_of_each_block_in_a_mixed_file() {
        let text = "fixedStep chrom=chr1 start=1 step=2 span=1\n1\n2\n\
                    variableStep chrom=chr1 span=3\n11\t3\n\
                    chr1\t20\t30\t4\n\
                    fixedStep chrom=chr1 start=41 step=1\n5\n";
        let (_, bytes) = wig(text);
        assert_eq!(
            section_kinds(&bytes),
            [
                SectionKind::FixedStep { step: 2, span: 1 },
                SectionKind::VariableStep { span: 3 },
                SectionKind::BedGraph,
                SectionKind::FixedStep { step: 1, span: 1 },
            ]
        );
    }

    #[test]
    fn wig_to_bigwig_reports_the_items_written() {
        let text = "fixedStep chrom=chr1 start=1 step=1\n1\n2\nvariableStep chrom=chr1\n11 3\n\
                    chr1\t20\t30\t4\n";
        assert_eq!(wig(text).0, Report { items: 4, ..Report::default() });
    }

    #[test]
    fn wig_to_bigwig_writes_nothing_for_a_declaration_without_values() {
        let text = "fixedStep chrom=chr2 start=1 step=1\nvariableStep chrom=chr1\n5 1\n";
        let (report, bytes) = wig(text);
        assert_eq!(report.items, 1);
        assert_eq!(items(&bytes), [item("chr1", 4, 5, 1.0)]);
    }

    #[test]
    fn wig_to_bigwig_writes_a_fixed_step_block_longer_than_a_batch_at_the_right_positions() {
        let values = BATCH_LEN * 2 + 10;
        let text = format!("fixedStep chrom=chr1 start=1 step=2 span=1\n{}", "7\n".repeat(values));
        let (report, bytes) = wig(&text);
        assert_eq!(report.items, u64::try_from(values).unwrap());
        let starts: Vec<u32> = items(&bytes).iter().map(|&(_, start, _, _)| start).collect();
        let expected: Vec<u32> = (0..).step_by(2).take(values).collect();
        assert_eq!(starts, expected);
    }

    #[test]
    fn wig_to_bigwig_writes_a_variable_step_block_longer_than_a_batch() {
        let values = BATCH_LEN + 10;
        let mut text = String::from("variableStep chrom=chr1\n");
        for position in 1..=values {
            writeln!(text, "{position}\t1").unwrap();
        }
        let (report, bytes) = wig(&text);
        assert_eq!(report.items, u64::try_from(values).unwrap());
        assert_eq!(items(&bytes).len(), values);
    }

    #[test]
    fn wig_to_bigwig_skips_and_counts_zero_length_bedgraph_lines() {
        let (report, bytes) = wig("chr1\t0\t10\t1\nchr1\t10\t10\t2\n");
        assert_eq!(report, Report { items: 1, skipped_empty: 1, ..Report::default() });
        assert_eq!(items(&bytes), [item("chr1", 0, 10, 1.0)]);
    }

    #[test]
    fn wig_to_bigwig_counts_fixed_step_values_truncated_or_dropped_by_clip() {
        // Values start at 495, 497, 499 and 501 on the 500-base chr2.
        let text = "fixedStep chrom=chr2 start=496 step=2 span=2\n1\n2\n3\n4\n";
        let (report, bytes) = wig_with(text, &WriterOptions::new().clip(true)).unwrap();
        assert_eq!(report, Report { items: 4, clipped: 1, dropped: 1, ..Report::default() });
        assert_eq!(
            items(&bytes),
            [item("chr2", 495, 497, 1.0), item("chr2", 497, 499, 2.0), item("chr2", 499, 500, 3.0)]
        );
    }

    #[test]
    fn wig_to_bigwig_counts_variable_step_items_truncated_or_dropped_by_clip() {
        // Items cover 480-490, 494-504 and 504-514 on the 500-base chr2.
        let text = "variableStep chrom=chr2 span=10\n481\t1\n495\t2\n505\t3\n";
        let (report, bytes) = wig_with(text, &WriterOptions::new().clip(true)).unwrap();
        assert_eq!(report, Report { items: 3, clipped: 1, dropped: 1, ..Report::default() });
        assert_eq!(items(&bytes), [item("chr2", 480, 490, 1.0), item("chr2", 494, 500, 2.0)]);
    }

    /// A wig item running past its chromosome's end must be reported, not silently written
    /// (<https://github.com/igvteam/igv/issues/425>).
    #[test]
    fn wig_item_past_chromosome_end_is_an_error_with_line_number_unless_clipped() {
        // On the 500-base chr2, the last item of each block covers 495-505.
        let fixed_step = "fixedStep chrom=chr2 start=476 step=10 span=10\n1\n2\n# last\n3\n";
        let variable_step = "variableStep chrom=chr2 span=10\n476\t1\n\n496\t2\n";

        let (line, error) = writer_rejection(wig_with(fixed_step, &WriterOptions::new()));
        assert_eq!(line, 5);
        assert!(
            matches!(error, Error::OutOfBounds { start: 495, end: 505, length: 500, .. }),
            "{error:?}"
        );
        let (line, error) = writer_rejection(wig_with(variable_step, &WriterOptions::new()));
        assert_eq!(line, 4);
        assert!(
            matches!(error, Error::OutOfBounds { start: 495, end: 505, length: 500, .. }),
            "{error:?}"
        );

        let clip = WriterOptions::new().clip(true);
        let (report, bytes) = wig_with(fixed_step, &clip).unwrap();
        assert_eq!(report, Report { items: 3, clipped: 1, ..Report::default() });
        assert_eq!(items(&bytes).last(), Some(&item("chr2", 495, 500, 3.0)));
        let (report, bytes) = wig_with(variable_step, &clip).unwrap();
        assert_eq!(report, Report { items: 2, clipped: 1, ..Report::default() });
        assert_eq!(items(&bytes).last(), Some(&item("chr2", 495, 500, 2.0)));
    }

    #[test]
    fn wig_to_bigwig_reports_an_item_overlapping_a_truncated_one_at_its_line() {
        // On the 300-base chr3 items cover 250-260, 296-306 (truncated) and 303-313.
        let text = "variableStep chrom=chr3 span=10\n251\t1\n297\t2\n304\t3\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new().clip(true)));
        assert_eq!(line, 4);
        assert!(
            matches!(
                error,
                Error::UnsortedOrOverlapping { ref chrom, start: 303, previous_end: 306 }
                    if chrom == "chr3"
            ),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_a_nan_fixed_step_value_at_its_line() {
        let text = "fixedStep chrom=chr1 start=1 step=1\n1\n2\nnan\n4\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 4);
        assert!(
            matches!(
                error,
                Error::NonFiniteValue { ref chrom, position: 2, value }
                    if chrom == "chr1" && value.is_nan()
            ),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_a_fixed_step_value_past_the_chromosome_end_at_its_line() {
        // Values start at 498, 499 and 500 on the 500-base chr2.
        let text = "fixedStep chrom=chr2 start=499 step=1\n1\n2\n3\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 4);
        assert!(
            matches!(error, Error::OutOfBounds { start: 500, end: 501, length: 500, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_the_first_rejected_value_of_a_batch() {
        // The NaN on line 3 comes before the value past the end of chr2 on line 12.
        let text = format!("fixedStep chrom=chr2 start=491 step=1\n1\nnan\n{}", "1\n".repeat(9));
        let (line, error) = writer_rejection(wig_with(&text, &WriterOptions::new()));
        assert_eq!(line, 3);
        assert!(matches!(error, Error::NonFiniteValue { position: 491, .. }), "{error:?}");
    }

    #[test]
    fn wig_to_bigwig_reports_a_rejected_value_after_comments_at_its_line() {
        let text = "fixedStep chrom=chr1 start=1 step=1\n1\n# note\n\n2\nnan\n";
        assert_eq!(writer_rejection(wig_with(text, &WriterOptions::new())).0, 6);
    }

    #[test]
    fn wig_to_bigwig_reports_a_rejected_value_in_a_later_batch_at_its_line() {
        let bad_index = BATCH_LEN + 100;
        let values: String = (0..BATCH_LEN * 2)
            .map(|index| if index == bad_index { "inf\n" } else { "1\n" })
            .collect();
        let text = format!("fixedStep chrom=chr1 start=1 step=1\n{values}");
        let (line, error) = writer_rejection(wig_with(&text, &WriterOptions::new()));
        // The header is line 1, so value `index` is on line `index + 2`.
        assert_eq!(line, u64::try_from(bad_index + 2).unwrap());
        let bad_position = u64::try_from(bad_index).unwrap();
        assert!(
            matches!(
                error,
                Error::NonFiniteValue { position, value, .. }
                    if position == bad_position && value == f32::INFINITY
            ),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_a_repeated_variable_step_position_at_its_line() {
        let text = "variableStep chrom=chr1\n10\t1\n10\t2\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 3);
        assert!(
            matches!(error, Error::UnsortedOrOverlapping { start: 9, previous_end: 10, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_an_unknown_chromosome_at_its_first_value() {
        let text = "fixedStep chrom=chrX start=1 step=1\n# values follow\n1\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 3);
        assert!(matches!(error, Error::UnknownChrom(ref chrom) if chrom == "chrX"), "{error:?}");
    }

    #[test]
    fn wig_to_bigwig_reports_a_returning_chromosome_at_its_line() {
        let text = "variableStep chrom=chr1\n1\t1\nvariableStep chrom=chr2\n1\t1\n\
                    variableStep chrom=chr1\n5\t1\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 6);
        assert!(
            matches!(error, Error::ChromAlreadyFinished(ref chrom) if chrom == "chr1"),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_a_span_wider_than_the_step_at_the_first_value() {
        let text = "fixedStep chrom=chr1 start=1 step=1 span=2\n1\n2\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 2);
        assert!(
            matches!(error, Error::InvalidItem { reason, .. }
                if reason.contains("span must not exceed step")),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_overlap_between_blocks_at_the_line_of_the_later_value() {
        let text =
            "fixedStep chrom=chr1 start=1 step=1 span=1\n1\n2\nvariableStep chrom=chr1\n2\t3\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 5);
        assert!(
            matches!(error, Error::UnsortedOrOverlapping { start: 1, previous_end: 2, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_reports_a_bedgraph_line_rejected_by_the_writer_at_its_line() {
        let text = "variableStep chrom=chr1\n10\t1\nchr1\t5\t20\t2\n";
        let (line, error) = writer_rejection(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 3);
        assert!(
            matches!(error, Error::UnsortedOrOverlapping { start: 5, previous_end: 10, .. }),
            "{error:?}"
        );
    }

    #[test]
    fn wig_to_bigwig_returns_a_sink_failure_while_pushing_unwrapped() {
        // The unparseable last line would fail the conversion if the sink failure came only at
        // the end.
        let text = format!("fixedStep chrom=chr1 start=1 step=1\n{}x\n", "1\n".repeat(20_000));
        let options = options_writing_to_the_sink_while_pushing();
        let result = wig_to_bigwig(text.as_bytes(), &chroms(), &options, FailingSink);
        assert!(matches!(result, Err(Error::Io(_))), "{result:?}");
    }

    #[test]
    fn wig_to_bigwig_passes_reader_parse_errors_through() {
        let text = "fixedStep chrom=chr1 start=1 step=1\n1\nchr1\t5\t6\t1\n2\n";
        let (line, message) = parse_failure(wig_with(text, &WriterOptions::new()));
        assert_eq!(line, 4);
        assert!(message.contains("before any fixedStep or variableStep declaration"), "{message}");
    }

    // ---- bigwig_to_bedgraph ------------------------------------------------------------------

    #[test]
    fn bigwig_to_bedgraph_writes_each_item_as_a_line() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 10, 5, &[1.0, 1.0])?;
            writer.push_variable_step("chr1", 2, &[(30, 2.5)])?;
            writer.push_interval("chr1", 40, 50, -1.0)
        });
        assert_eq!(
            to_bedgraph(&bytes, &Selection::All, false),
            "chr1\t0\t5\t1\nchr1\t10\t15\t1\nchr1\t30\t32\t2.5\nchr1\t40\t50\t-1\n"
        );
    }

    #[test]
    fn bigwig_to_bedgraph_merges_adjacent_equal_values_when_asked() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &[(0, 10, 1.0), (10, 20, 1.0), (20, 30, 2.0)])
        });
        assert_eq!(to_bedgraph(&bytes, &Selection::All, true), "chr1\t0\t20\t1\nchr1\t20\t30\t2\n");
    }

    #[test]
    fn bigwig_to_bedgraph_keeps_adjacent_equal_values_apart_without_merging() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &[(0, 10, 1.0), (10, 20, 1.0)])
        });
        assert_eq!(
            to_bedgraph(&bytes, &Selection::All, false),
            "chr1\t0\t10\t1\nchr1\t10\t20\t1\n"
        );
    }

    #[test]
    fn bigwig_to_bedgraph_does_not_merge_equal_values_across_a_gap() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &[(0, 10, 1.0), (11, 20, 1.0)])
        });
        assert_eq!(to_bedgraph(&bytes, &Selection::All, true), "chr1\t0\t10\t1\nchr1\t11\t20\t1\n");
    }

    #[test]
    fn bigwig_to_bedgraph_does_not_merge_equal_values_across_chromosomes() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_interval("chr1", 5, 10, 1.0)?;
            writer.push_interval("chr2", 10, 20, 1.0)
        });
        assert_eq!(to_bedgraph(&bytes, &Selection::All, true), "chr1\t5\t10\t1\nchr2\t10\t20\t1\n");
    }

    #[test]
    fn bigwig_to_bedgraph_merges_equal_values_across_sections_and_section_types() {
        let bytes = build(&WriterOptions::new().items_per_section(2), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &[1.0, 1.0, 1.0])?;
            writer.push_interval("chr1", 3, 10, 1.0)?;
            writer.push_variable_step("chr1", 5, &[(10, 2.0)])
        });
        assert_eq!(to_bedgraph(&bytes, &Selection::All, true), "chr1\t0\t10\t1\nchr1\t10\t15\t2\n");
    }

    #[test]
    fn bigwig_to_bedgraph_writes_chromosomes_in_chrom_id_order() {
        // With arrival order, chr2 receives its chromId first.
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_interval("chr2", 0, 10, 2.0)?;
            writer.push_interval("chr1", 0, 10, 1.0)
        });
        assert_eq!(to_bedgraph(&bytes, &Selection::All, false), "chr2\t0\t10\t2\nchr1\t0\t10\t1\n");
    }

    #[test]
    fn bigwig_to_bedgraph_selection_chrom_writes_only_that_chromosome() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_interval("chr1", 0, 10, 1.0)?;
            writer.push_interval("chr2", 0, 10, 2.0)
        });
        let selection = Selection::Chrom("chr2".to_string());
        assert_eq!(to_bedgraph(&bytes, &selection, false), "chr2\t0\t10\t2\n");
    }

    #[test]
    fn bigwig_to_bedgraph_region_clips_items_at_both_ends() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_intervals(
                "chr1",
                &[(0, 10, 1.0), (10, 20, 2.0), (20, 30, 3.0), (40, 50, 4.0)],
            )
        });
        assert_eq!(
            to_bedgraph(&bytes, &region("chr1", 5, 25), false),
            "chr1\t5\t10\t1\nchr1\t10\t20\t2\nchr1\t20\t25\t3\n"
        );
    }

    #[test]
    fn bigwig_to_bedgraph_writes_nothing_for_a_chromosome_without_data() {
        let bytes = build(&WriterOptions::new(), |writer| writer.push_interval("chr1", 0, 10, 1.0));
        assert_eq!(to_bedgraph(&bytes, &Selection::Chrom("chr3".to_string()), true), "");
    }

    #[test]
    fn bigwig_to_bedgraph_fails_for_an_unknown_chromosome() {
        let bytes = build(&WriterOptions::new(), |writer| writer.push_interval("chr1", 0, 10, 1.0));
        let result = try_to_bedgraph(&bytes, &Selection::Chrom("chrX".to_string()), true);
        assert!(
            matches!(result, Err(Error::UnknownChrom(ref chrom)) if chrom == "chrX"),
            "{result:?}"
        );
    }

    #[test]
    fn bigwig_to_bedgraph_fails_for_a_region_past_the_chromosome_end() {
        let bytes = build(&WriterOptions::new(), |writer| writer.push_interval("chr2", 0, 10, 1.0));
        let result = try_to_bedgraph(&bytes, &region("chr2", 0, 501), true);
        assert!(matches!(result, Err(Error::InvalidRegion { .. })), "{result:?}");
    }

    // ---- bigwig_to_wig -----------------------------------------------------------------------

    #[test]
    fn bigwig_to_wig_writes_a_fixed_step_section_as_a_fixed_step_block() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 100, 10, 5, &[1.0, 2.0, 3.0])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "fixedStep chrom=chr1 start=101 step=10 span=5\n1\n2\n3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_writes_a_variable_step_section_as_a_variable_step_block() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_variable_step("chr1", 5, &[(10, 1.5), (100, 2.0)])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "variableStep chrom=chr1 span=5\n11\t1.5\n101\t2\n"
        );
    }

    #[test]
    fn bigwig_to_wig_writes_a_bedgraph_section_after_a_comment_giving_its_extent() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &[(0, 10, 1.0), (20, 30, 2.0)])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "#bedGraph section chr1:0-30\nchr1\t0\t10\t1\nchr1\t20\t30\t2\n"
        );
    }

    #[test]
    fn bigwig_to_wig_comments_each_stored_bedgraph_section() {
        let bytes = build(&WriterOptions::new().items_per_section(2), |writer| {
            writer.push_intervals("chr1", &[(0, 10, 1.0), (10, 20, 2.0), (30, 40, 3.0)])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "#bedGraph section chr1:0-20\nchr1\t0\t10\t1\nchr1\t10\t20\t2\n\
             #bedGraph section chr1:30-40\nchr1\t30\t40\t3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_writes_a_fixed_step_run_of_several_sections_under_one_header() {
        let bytes = build(&WriterOptions::new().items_per_section(2), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &[1.0, 2.0, 3.0, 4.0, 5.0])
        });
        assert_eq!(section_kinds(&bytes).len(), 3);
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "fixedStep chrom=chr1 start=1 step=1 span=1\n1\n2\n3\n4\n5\n"
        );
    }

    #[test]
    fn bigwig_to_wig_writes_a_variable_step_run_of_several_sections_under_one_header() {
        let bytes = build(&WriterOptions::new().items_per_section(2), |writer| {
            writer.push_variable_step("chr1", 1, &[(0, 1.0), (5, 2.0), (9, 3.0)])
        });
        assert_eq!(section_kinds(&bytes).len(), 2);
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "variableStep chrom=chr1 span=1\n1\t1\n6\t2\n10\t3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_starts_a_new_fixed_step_header_after_a_gap() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &[1.0, 2.0])?;
            writer.push_fixed_step("chr1", 10, 1, 1, &[3.0])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "fixedStep chrom=chr1 start=1 step=1 span=1\n1\n2\n\
             fixedStep chrom=chr1 start=11 step=1 span=1\n3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_starts_a_new_variable_step_header_when_the_span_changes() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_variable_step("chr1", 1, &[(0, 1.0)])?;
            writer.push_variable_step("chr1", 2, &[(10, 2.0)])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "variableStep chrom=chr1 span=1\n1\t1\nvariableStep chrom=chr1 span=2\n11\t2\n"
        );
    }

    #[test]
    fn bigwig_to_wig_starts_a_new_header_for_each_chromosome() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 1, 1, &[1.0, 2.0])?;
            writer.push_fixed_step("chr2", 2, 1, 1, &[3.0])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "fixedStep chrom=chr1 start=1 step=1 span=1\n1\n2\n\
             fixedStep chrom=chr2 start=3 step=1 span=1\n3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_starts_a_new_header_after_bedgraph_lines() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_variable_step("chr1", 1, &[(0, 1.0)])?;
            writer.push_interval("chr1", 5, 6, 2.0)?;
            writer.push_variable_step("chr1", 1, &[(10, 3.0)])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "variableStep chrom=chr1 span=1\n1\t1\n\
             #bedGraph section chr1:5-6\nchr1\t5\t6\t2\n\
             variableStep chrom=chr1 span=1\n11\t3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_writes_chromosomes_in_chrom_id_order() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_variable_step("chr2", 1, &[(0, 2.0)])?;
            writer.push_variable_step("chr1", 1, &[(0, 1.0)])
        });
        assert_eq!(
            to_wig(&bytes, &Selection::All),
            "variableStep chrom=chr2 span=1\n1\t2\nvariableStep chrom=chr1 span=1\n1\t1\n"
        );
    }

    #[test]
    fn bigwig_to_wig_selection_chrom_writes_only_that_chromosome() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_variable_step("chr1", 1, &[(0, 1.0)])?;
            writer.push_variable_step("chr2", 1, &[(0, 2.0)])
        });
        let selection = Selection::Chrom("chr2".to_string());
        assert_eq!(to_wig(&bytes, &selection), "variableStep chrom=chr2 span=1\n1\t2\n");
    }

    #[test]
    fn bigwig_to_wig_region_gives_fixed_step_items_cut_at_both_ends_headers_of_their_own() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 10, 10, &[1.0, 2.0, 3.0, 4.0])
        });
        assert_eq!(
            to_wig(&bytes, &region("chr1", 5, 35)),
            "fixedStep chrom=chr1 start=6 step=10 span=5\n1\n\
             fixedStep chrom=chr1 start=11 step=10 span=10\n2\n3\n\
             fixedStep chrom=chr1 start=31 step=10 span=5\n4\n"
        );
    }

    #[test]
    fn bigwig_to_wig_region_gives_variable_step_items_cut_at_both_ends_headers_of_their_own() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_variable_step("chr1", 10, &[(0, 1.0), (20, 2.0), (40, 3.0)])
        });
        assert_eq!(
            to_wig(&bytes, &region("chr1", 5, 45)),
            "variableStep chrom=chr1 span=5\n6\t1\n\
             variableStep chrom=chr1 span=10\n21\t2\n\
             variableStep chrom=chr1 span=5\n41\t3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_region_clips_bedgraph_items_at_both_ends_and_comments_the_clipped_extent() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_intervals(
                "chr1",
                &[(0, 10, 1.0), (10, 20, 2.0), (20, 30, 3.0), (40, 50, 4.0)],
            )
        });
        assert_eq!(
            to_wig(&bytes, &region("chr1", 5, 25)),
            "#bedGraph section chr1:5-25\nchr1\t5\t10\t1\nchr1\t10\t20\t2\nchr1\t20\t25\t3\n"
        );
    }

    #[test]
    fn bigwig_to_wig_writes_no_header_for_a_region_between_sparse_fixed_step_items() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_fixed_step("chr1", 0, 10, 1, &[1.0, 2.0, 3.0])
        });
        assert_eq!(to_wig(&bytes, &region("chr1", 2, 8)), "");
    }

    #[test]
    fn bigwig_to_wig_writes_no_comment_for_a_region_between_bedgraph_items() {
        let bytes = build(&WriterOptions::new(), |writer| {
            writer.push_intervals("chr1", &[(0, 10, 1.0), (50, 60, 2.0)])
        });
        assert_eq!(to_wig(&bytes, &region("chr1", 20, 30)), "");
    }

    #[test]
    fn bigwig_to_wig_writes_nothing_for_a_chromosome_without_data() {
        let bytes = build(&WriterOptions::new(), |writer| writer.push_interval("chr1", 0, 10, 1.0));
        assert_eq!(to_wig(&bytes, &Selection::Chrom("chr3".to_string())), "");
    }

    #[test]
    fn bigwig_to_wig_fails_for_an_unknown_chromosome() {
        let bytes = build(&WriterOptions::new(), |writer| writer.push_interval("chr1", 0, 10, 1.0));
        let result = try_to_wig(&bytes, &region("chrX", 0, 10));
        assert!(
            matches!(result, Err(Error::UnknownChrom(ref chrom)) if chrom == "chrX"),
            "{result:?}"
        );
    }

    // ---- Round trips -------------------------------------------------------------------------

    #[test]
    fn bedgraph_round_trips_through_bigwig_unchanged() {
        let text = "chr1\t0\t10\t1.5\nchr1\t10\t20\t1.5\nchr1\t30\t31\t1.2345679\nchr2\t5\t6\t-2\n";
        let (_, bytes) = bedgraph(text);
        assert_eq!(to_bedgraph(&bytes, &Selection::All, false), text);
    }

    #[test]
    fn fixed_step_wig_round_trips_through_bigwig_unchanged() {
        let text = "fixedStep chrom=chr1 start=101 step=10 span=5\n1.5\n2\n-3\n\
                    fixedStep chrom=chr2 start=1 step=1 span=1\n0.25\n";
        let (_, bytes) = wig(text);
        assert_eq!(to_wig(&bytes, &Selection::All), text);
    }

    #[test]
    fn variable_step_wig_round_trips_through_bigwig_unchanged() {
        let text = "variableStep chrom=chr1 span=5\n11\t1.5\n101\t2\n\
                    variableStep chrom=chr1 span=1\n201\t3\n";
        let (_, bytes) = wig(text);
        assert_eq!(to_wig(&bytes, &Selection::All), text);
    }

    #[test]
    fn wig_with_bedgraph_lines_round_trips_through_bigwig_unchanged() {
        let text = "variableStep chrom=chr1 span=5\n1\t1\n\
                    #bedGraph section chr1:10-30\nchr1\t10\t20\t2\nchr1\t20\t30\t3\n\
                    fixedStep chrom=chr1 start=101 step=1 span=1\n4\n5\n";
        let (_, bytes) = wig(text);
        assert_eq!(to_wig(&bytes, &Selection::All), text);
    }

    #[test]
    fn long_fixed_step_wig_round_trips_through_bigwig_under_one_header() {
        let mut text = String::from("fixedStep chrom=chr1 start=1 step=1 span=1\n");
        for index in 0..BATCH_LEN * 2 {
            writeln!(text, "{}", index % 7).unwrap();
        }
        let (_, bytes) = wig_with(&text, &WriterOptions::new().items_per_section(1_000)).unwrap();
        assert_eq!(to_wig(&bytes, &Selection::All), text);
    }
}
