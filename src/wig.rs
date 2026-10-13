//! Wig text (`fixedStep` and `variableStep`): a streaming reader and a writer.
//!
//! Wig positions are 1-based; every coordinate this module hands out or accepts is 0-based, like
//! the rest of the crate. The reader follows what kent's `wigToBigWig` does rather than what the
//! UCSC documentation says where the two differ (notably, a `fixedStep` block without `span`
//! uses `span = step`). It does not check ordering or overlaps; the bigWig writer does.

use std::io::{self, BufRead, Write};
use std::mem;

use crate::bedgraph::{self, parse_record};
use crate::text::{
    LineReader, check_chrom, is_separator, parse_chrom, parse_error, parse_u32, parse_unsigned,
    parse_value, split_fields, write_bedgraph_line, write_value,
};
use crate::{Error, Result};

/// One parsed wig event. Coordinates are 0-based.
#[derive(Clone, Debug, PartialEq)]
pub enum Event {
    /// A `fixedStep` declaration. The `FixedValue`s that follow sit at `start`, `start + step`,
    /// `start + 2 * step`, ..., each covering `span` bases.
    FixedStep {
        /// Chromosome name.
        chrom: String,
        /// 0-based position of the first value.
        start: u32,
        /// Distance between consecutive values.
        step: u32,
        /// Bases covered by each value; equals `step` when the file gave no `span`.
        span: u32,
    },
    /// A `variableStep` declaration. The `VariableValue`s that follow each cover `span` bases.
    VariableStep {
        /// Chromosome name.
        chrom: String,
        /// Bases covered by each value; 1 when the file gave no `span`.
        span: u32,
    },
    /// A value inside a `fixedStep` block; its position follows from the declaration and the
    /// number of values already read.
    FixedValue(f32),
    /// A position and value inside a `variableStep` block.
    VariableValue {
        /// 0-based position.
        start: u32,
        /// The value.
        value: f32,
    },
    /// A four-column bedGraph line inside a wig file, which kent accepts and so do we. It ends
    /// any open `fixedStep` or `variableStep` block.
    BedGraph(bedgraph::Record),
}

/// The wig block the reader is currently inside.
enum Section {
    /// No declaration yet, or a bedGraph line closed the previous block.
    None,
    Fixed {
        /// 0-based position of the next value; `u64` so advancing cannot overflow before the
        /// range check against `u32`.
        next_start: u64,
        step: u32,
        span: u32,
    },
    Variable {
        span: u32,
    },
}

/// Streaming wig reader over any [`BufRead`].
///
/// Understands `fixedStep` and `variableStep` blocks and four-column bedGraph lines mixed freely
/// among them. Blank lines, `#` comments and leading `track`/`browser` lines follow the same
/// rules as [`bedgraph::Reader`]. Declarations must carry `chrom=` and, for `fixedStep`, `start=`
/// (at least 1) and `step=`; `variableStep` rejects `start=` and `step=`. Unknown or repeated
/// keys, non-integer positions, zero `step` or `span`, and positions that overflow `u32` are
/// parse errors, as is a value line before any declaration. A parse error does not make the
/// reader unusable.
///
/// ```
/// use bigwig::wig::{Event, Reader};
///
/// let text = "fixedStep chrom=chr1 start=101 step=10\n1.5\n2.5\n";
/// let events: Vec<Event> = Reader::new(text.as_bytes()).collect::<Result<_, _>>()?;
/// assert_eq!(
///     events,
///     [
///         // `span` defaults to `step`, and `start` is now 0-based.
///         Event::FixedStep { chrom: "chr1".into(), start: 100, step: 10, span: 10 },
///         Event::FixedValue(1.5),
///         Event::FixedValue(2.5),
///     ]
/// );
/// # Ok::<(), bigwig::Error>(())
/// ```
pub struct Reader<R> {
    lines: LineReader<R>,
    section: Section,
    /// Set after an I/O error so that iterating does not spin on a persistent failure.
    finished: bool,
}

impl<R: BufRead> Reader<R> {
    /// Creates a reader over `inner`.
    pub fn new(inner: R) -> Self {
        Self { lines: LineReader::new(inner), section: Section::None, finished: false }
    }

    /// Reads the next event into `event`, reusing its `String` where possible. Returns
    /// `Ok(false)` at end of input. After an error `event` holds an unspecified value.
    pub fn read_event(&mut self, event: &mut Event) -> Result<bool> {
        if !self.lines.next_line()? {
            return Ok(false);
        }
        let line = self.lines.line();
        let line_number = self.lines.line_number();
        match line.split(is_separator).next().unwrap_or("") {
            "fixedStep" => {
                // Closed first, so that after a bad declaration the values that follow are
                // rejected rather than attributed to the previous block.
                self.section = Section::None;
                self.section = read_fixed_step_header(line, line_number, event)?;
            }
            "variableStep" => {
                self.section = Section::None;
                self.section = read_variable_step_header(line, line_number, event)?;
            }
            _ => read_data_line(&mut self.section, line, line_number, event)?,
        }
        Ok(true)
    }

    /// The leading `track` lines seen so far, verbatim and trimmed.
    pub fn track_lines(&self) -> &[String] {
        &self.lines.track_lines
    }

    /// The 1-based number of the line the last event was read from.
    pub(crate) fn line_number(&self) -> u64 {
        self.lines.line_number()
    }
}

impl<R: BufRead> Iterator for Reader<R> {
    type Item = Result<Event>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let mut event = Event::FixedValue(0.0);
        match self.read_event(&mut event) {
            Ok(true) => Some(Ok(event)),
            Ok(false) => None,
            Err(error) => {
                self.finished = matches!(error, Error::Io(_));
                Some(Err(error))
            }
        }
    }
}

/// Which section kind a writer is inside, so values can be checked against their header.
#[derive(Clone, Copy, PartialEq, Eq)]
enum WriterSection {
    None,
    Fixed,
    Variable,
}

/// Writes wig text: `fixedStep` and `variableStep` blocks, bedGraph lines, and comments.
///
/// Headers always carry an explicit `span=` and 1-based positions, so the output means the same
/// under the UCSC documentation's `span` default and kent's. Values use the shortest form that
/// parses back to the same `f32`. A value call must match the most recent header
/// (`fixed_value` after `fixed_step`, `variable_value` after `variable_step`); a bedGraph line
/// closes the current block, so a new header is needed before more values. Each line reaches the
/// underlying writer in a single `write_all`; wrap unbuffered sinks in a
/// [`std::io::BufWriter`].
///
/// ```
/// use bigwig::wig::Writer;
///
/// let mut writer = Writer::new(Vec::new());
/// writer.fixed_step("chr1", 100, 10, 10)?;
/// writer.fixed_value(1.5)?;
/// writer.fixed_value(2.0)?;
/// writer.variable_step("chr2", 5)?;
/// writer.variable_value(7, 3.0)?;
/// let text = String::from_utf8(writer.into_inner()?).unwrap();
/// assert_eq!(
///     text,
///     "fixedStep chrom=chr1 start=101 step=10 span=10\n1.5\n2\n\
///      variableStep chrom=chr2 span=5\n8\t3\n"
/// );
/// # Ok::<(), bigwig::Error>(())
/// ```
pub struct Writer<W> {
    inner: W,
    line: Vec<u8>,
    section: WriterSection,
}

impl<W: Write> Writer<W> {
    /// Creates a writer over `inner`.
    pub fn new(inner: W) -> Self {
        Self { inner, line: Vec::new(), section: WriterSection::None }
    }

    /// Starts a `fixedStep` block whose first value is at 0-based `start`; the following
    /// [`Self::fixed_value`] calls are at `start`, `start + step`, ..., each covering `span`
    /// bases. Fails with [`Error::InvalidItem`] for an unwritable `chrom` or a zero `step` or
    /// `span`.
    pub fn fixed_step(&mut self, chrom: &str, start: u32, step: u32, span: u32) -> Result<()> {
        check_step_and_span(chrom, start, step, span)?;
        self.line.clear();
        writeln!(
            self.line,
            "fixedStep chrom={chrom} start={} step={step} span={span}",
            u64::from(start) + 1
        )?;
        self.write_line(WriterSection::Fixed)
    }

    /// Starts a `variableStep` block whose values each cover `span` bases. Fails with
    /// [`Error::InvalidItem`] for an unwritable `chrom` or a zero `span`.
    pub fn variable_step(&mut self, chrom: &str, span: u32) -> Result<()> {
        check_step_and_span(chrom, 0, 1, span)?;
        self.line.clear();
        writeln!(self.line, "variableStep chrom={chrom} span={span}")?;
        self.write_line(WriterSection::Variable)
    }

    /// Writes one value of the current `fixedStep` block.
    pub fn fixed_value(&mut self, value: f32) -> Result<()> {
        self.require_section(WriterSection::Fixed, "fixed_value", "fixed_step")?;
        self.line.clear();
        write_value(&mut self.line, value)?;
        self.line.push(b'\n');
        self.write_line(WriterSection::Fixed)
    }

    /// Writes one `(start, value)` item of the current `variableStep` block; `start` is 0-based.
    pub fn variable_value(&mut self, start: u32, value: f32) -> Result<()> {
        self.require_section(WriterSection::Variable, "variable_value", "variable_step")?;
        self.line.clear();
        write!(self.line, "{}\t", u64::from(start) + 1)?;
        write_value(&mut self.line, value)?;
        self.line.push(b'\n');
        self.write_line(WriterSection::Variable)
    }

    /// Writes a bedGraph line, ending any open block. Fails with [`Error::InvalidItem`] for an
    /// unwritable `chrom`.
    pub fn bedgraph_record(&mut self, chrom: &str, start: u32, end: u32, value: f32) -> Result<()> {
        check_chrom(chrom, u64::from(start))?;
        self.line.clear();
        write_bedgraph_line(&mut self.line, chrom, start, end, value)?;
        self.write_line(WriterSection::None)
    }

    /// Writes `#` followed by `text` as a comment line, for example
    /// `bedGraph section chr1:1-100` ahead of bedGraph data, as kent does. Readers skip it
    /// wherever it appears and it does not end the current block. Fails if `text` contains a
    /// line break.
    pub fn comment(&mut self, text: &str) -> Result<()> {
        if text.contains(['\n', '\r']) {
            return Err(invalid_use("comment text contains a line break"));
        }
        self.line.clear();
        writeln!(self.line, "#{text}")?;
        self.write_line(self.section)
    }

    /// Flushes and returns the underlying writer.
    pub fn into_inner(mut self) -> Result<W> {
        self.inner.flush()?;
        Ok(self.inner)
    }

    /// Writes the buffered line and, only if that succeeds, moves to `next_section`.
    fn write_line(&mut self, next_section: WriterSection) -> Result<()> {
        self.inner.write_all(&self.line)?;
        self.section = next_section;
        Ok(())
    }

    fn require_section(
        &self,
        needed: WriterSection,
        method: &str,
        header_method: &str,
    ) -> Result<()> {
        if self.section == needed {
            Ok(())
        } else {
            Err(invalid_use(format!("{method} called outside a block started by {header_method}")))
        }
    }
}

/// Interprets a line that is not a declaration: a bedGraph line, or a value in the open block.
fn read_data_line(
    section: &mut Section,
    line: &str,
    line_number: u64,
    event: &mut Event,
) -> Result<()> {
    let (count, fields) = split_fields(line);
    if count == 4 {
        // A bedGraph line ends the block even when it turns out to be malformed.
        *section = Section::None;
        let mut record = match mem::replace(event, Event::FixedValue(0.0)) {
            Event::BedGraph(record) => record,
            other => bedgraph::Record { chrom: chrom_buffer(other), ..Default::default() },
        };
        let parsed = parse_record(line, line_number, &mut record);
        *event = Event::BedGraph(record);
        return parsed;
    }
    match section {
        Section::None => Err(parse_error(
            line_number,
            "data line before any fixedStep or variableStep declaration",
        )),
        Section::Fixed { next_start, step, span } => {
            // Advance before validating: a rejected line still occupies its grid slot, so the
            // values after it keep their true positions if the caller skips the error.
            let position = *next_start;
            *next_start = position + u64::from(*step);
            if count != 1 {
                return Err(parse_error(
                    line_number,
                    format!("expected 1 field (a value) in fixedStep data, found {count}"),
                ));
            }
            let value = parse_value(fields[0], line_number)?;
            checked_end(position, *span, line_number)?;
            *event = Event::FixedValue(value);
            Ok(())
        }
        Section::Variable { span } => {
            if count != 2 {
                return Err(parse_error(
                    line_number,
                    format!(
                        "expected 2 fields (position and value) in variableStep data, found {count}"
                    ),
                ));
            }
            let start = one_based_to_zero_based(fields[0], line_number, "position")?;
            checked_end(u64::from(start), *span, line_number)?;
            let value = parse_value(fields[1], line_number)?;
            *event = Event::VariableValue { start, value };
            Ok(())
        }
    }
}

/// Parses a `fixedStep` declaration into `event` and returns the section it opens.
fn read_fixed_step_header(line: &str, line_number: u64, event: &mut Event) -> Result<Section> {
    let header = parse_header_keys(line, line_number, true)?;
    let chrom =
        header.chrom.ok_or_else(|| parse_error(line_number, "fixedStep requires chrom="))?;
    let chrom = parse_chrom(chrom, line_number)?;
    let start_field =
        header.start.ok_or_else(|| parse_error(line_number, "fixedStep requires start="))?;
    let step = header.step.ok_or_else(|| parse_error(line_number, "fixedStep requires step="))?;
    let start = one_based_to_zero_based(start_field, line_number, "start")?;
    let step = nonzero_u32(step, line_number, "step")?;
    let span = match header.span {
        Some(span) => nonzero_u32(span, line_number, "span")?,
        None => step,
    };
    let mut buffer = chrom_buffer(mem::replace(event, Event::FixedValue(0.0)));
    buffer.push_str(chrom);
    *event = Event::FixedStep { chrom: buffer, start, step, span };
    Ok(Section::Fixed { next_start: u64::from(start), step, span })
}

/// Parses a `variableStep` declaration into `event` and returns the section it opens.
fn read_variable_step_header(line: &str, line_number: u64, event: &mut Event) -> Result<Section> {
    let header = parse_header_keys(line, line_number, false)?;
    let chrom =
        header.chrom.ok_or_else(|| parse_error(line_number, "variableStep requires chrom="))?;
    let chrom = parse_chrom(chrom, line_number)?;
    let span = match header.span {
        Some(span) => nonzero_u32(span, line_number, "span")?,
        None => 1,
    };
    let mut buffer = chrom_buffer(mem::replace(event, Event::FixedValue(0.0)));
    buffer.push_str(chrom);
    *event = Event::VariableStep { chrom: buffer, span };
    Ok(Section::Variable { span })
}

/// The `key=value` pairs of a declaration line; values are still text.
#[derive(Default)]
struct HeaderKeys<'a> {
    chrom: Option<&'a str>,
    start: Option<&'a str>,
    step: Option<&'a str>,
    span: Option<&'a str>,
}

/// Splits a declaration line into its keys. `is_fixed` selects which keys are allowed.
fn parse_header_keys(line: &str, line_number: u64, is_fixed: bool) -> Result<HeaderKeys<'_>> {
    let kind = if is_fixed { "fixedStep" } else { "variableStep" };
    let mut keys = HeaderKeys::default();
    for word in line.split(is_separator).filter(|word| !word.is_empty()).skip(1) {
        let (key, value) = word.split_once('=').ok_or_else(|| {
            parse_error(
                line_number,
                format!("expected key=value in {kind} declaration, found {word:?}"),
            )
        })?;
        let slot = match key {
            "chrom" => &mut keys.chrom,
            "span" => &mut keys.span,
            "start" | "step" if is_fixed => {
                if key == "start" {
                    &mut keys.start
                } else {
                    &mut keys.step
                }
            }
            "start" | "step" => {
                return Err(parse_error(
                    line_number,
                    format!("{key}= is not allowed in variableStep"),
                ));
            }
            _ => return Err(parse_error(line_number, format!("unknown key {key:?} in {kind}"))),
        };
        if value.is_empty() {
            return Err(parse_error(line_number, format!("{key}= has no value")));
        }
        if slot.replace(value).is_some() {
            return Err(parse_error(line_number, format!("repeated key {key}=")));
        }
    }
    Ok(keys)
}

/// Takes the chromosome `String` out of `event` (cleared) so it can be reused.
fn chrom_buffer(event: Event) -> String {
    let mut chrom = match event {
        Event::FixedStep { chrom, .. } | Event::VariableStep { chrom, .. } => chrom,
        Event::BedGraph(record) => record.chrom,
        Event::FixedValue(_) | Event::VariableValue { .. } => String::new(),
    };
    chrom.clear();
    chrom
}

/// Parses a 1-based position (at least 1) and converts it to 0-based.
fn one_based_to_zero_based(field: &str, line_number: u64, what: &str) -> Result<u32> {
    let position = parse_unsigned(field, line_number, what)?;
    position.checked_sub(1).and_then(|zero_based| u32::try_from(zero_based).ok()).ok_or_else(|| {
        parse_error(line_number, format!("{what} {field:?} is out of range; positions are 1-based"))
    })
}

/// Parses a `step` or `span` value, which must be a positive `u32`.
fn nonzero_u32(field: &str, line_number: u64, what: &str) -> Result<u32> {
    match parse_u32(field, line_number, what)? {
        0 => Err(parse_error(line_number, format!("{what} must be at least 1"))),
        value => Ok(value),
    }
}

/// Checks that an item at 0-based `start` covering `span` bases ends at or before `u32::MAX`.
fn checked_end(start: u64, span: u32, line_number: u64) -> Result<u32> {
    start.checked_add(u64::from(span)).and_then(|end| u32::try_from(end).ok()).ok_or_else(|| {
        parse_error(line_number, "position plus span exceeds the maximum coordinate")
    })
}

/// Validates the arguments of a block header.
fn check_step_and_span(chrom: &str, start: u32, step: u32, span: u32) -> Result<()> {
    check_chrom(chrom, u64::from(start))?;
    let reason = if step == 0 {
        "step must be at least 1"
    } else if span == 0 {
        "span must be at least 1"
    } else {
        return Ok(());
    };
    Err(Error::InvalidItem { chrom: chrom.to_string(), start: u64::from(start), reason })
}

/// An error for a call that cannot produce valid wig text.
fn invalid_use(message: impl Into<String>) -> Error {
    Error::Io(io::Error::new(io::ErrorKind::InvalidInput, message.into()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bedgraph::Record;

    fn read_all(text: &str) -> Result<Vec<Event>> {
        Reader::new(text.as_bytes()).collect()
    }

    fn parse_failure(text: &str) -> (u64, String) {
        match read_all(text) {
            Err(Error::Parse { line, message }) => (line, message),
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    fn fixed(chrom: &str, start: u32, step: u32, span: u32) -> Event {
        Event::FixedStep { chrom: chrom.to_string(), start, step, span }
    }

    fn variable(chrom: &str, span: u32) -> Event {
        Event::VariableStep { chrom: chrom.to_string(), span }
    }

    fn bedgraph_event(chrom: &str, start: u32, end: u32, value: f32) -> Event {
        Event::BedGraph(Record { chrom: chrom.to_string(), start, end, value })
    }

    fn written(build: impl FnOnce(&mut Writer<Vec<u8>>) -> Result<()>) -> String {
        let mut writer = Writer::new(Vec::new());
        build(&mut writer).unwrap();
        String::from_utf8(writer.into_inner().unwrap()).unwrap()
    }

    // ---- reader: fixedStep -------------------------------------------------------------

    #[test]
    fn fixed_step_converts_start_to_zero_based() {
        let events = read_all("fixedStep chrom=chr1 start=1 step=1 span=1\n5\n").unwrap();
        assert_eq!(events, vec![fixed("chr1", 0, 1, 1), Event::FixedValue(5.0)]);
    }

    #[test]
    fn fixed_step_span_defaults_to_step() {
        let events = read_all("fixedStep chrom=chr1 start=11 step=5\n1\n").unwrap();
        assert_eq!(events[0], fixed("chr1", 10, 5, 5));
    }

    #[test]
    fn fixed_step_explicit_span_is_kept() {
        let events = read_all("fixedStep chrom=chr1 start=11 step=5 span=2\n1\n").unwrap();
        assert_eq!(events[0], fixed("chr1", 10, 5, 2));
    }

    #[test]
    fn fixed_step_keys_may_come_in_any_order() {
        let events = read_all("fixedStep span=2 step=5 start=11 chrom=chr1\n").unwrap();
        assert_eq!(events, vec![fixed("chr1", 10, 5, 2)]);
    }

    #[test]
    fn fixed_step_span_greater_than_step_is_not_rejected() {
        let events = read_all("fixedStep chrom=c start=1 step=2 span=5\n1\n").unwrap();
        assert_eq!(events[0], fixed("c", 0, 2, 5));
    }

    #[test]
    fn fixed_step_requires_step() {
        let (line, message) = parse_failure("fixedStep chrom=c start=1\n");
        assert_eq!(line, 1);
        assert!(message.contains("step"), "{message}");
    }

    #[test]
    fn fixed_step_requires_start() {
        assert!(parse_failure("fixedStep chrom=c step=1\n").1.contains("start"));
    }

    #[test]
    fn control_character_in_fixed_step_chromosome_name_is_a_parse_error() {
        assert!(parse_failure("fixedStep chrom=c\u{8} start=1 step=1\n").1.contains("control"));
    }

    #[test]
    fn unicode_whitespace_in_variable_step_chromosome_name_is_a_parse_error() {
        assert!(parse_failure("variableStep chrom=c\u{a0}\n").1.contains("whitespace"));
    }

    #[test]
    fn fixed_step_requires_chrom() {
        assert!(parse_failure("fixedStep start=1 step=1\n").1.contains("chrom"));
    }

    #[test]
    fn fixed_step_start_zero_is_an_error() {
        assert_eq!(parse_failure("\nfixedStep chrom=c start=0 step=1\n").0, 2);
    }

    #[test]
    fn fixed_step_zero_step_is_an_error() {
        assert!(parse_failure("fixedStep chrom=c start=1 step=0\n").1.contains("step"));
    }

    #[test]
    fn fixed_step_zero_span_is_an_error() {
        assert!(parse_failure("fixedStep chrom=c start=1 step=1 span=0\n").1.contains("span"));
    }

    #[test]
    fn fixed_step_non_integer_start_is_an_error() {
        assert_eq!(parse_failure("fixedStep chrom=c start=5.7 step=1\n").0, 1);
    }

    #[test]
    fn fixed_step_non_integer_step_is_an_error() {
        assert_eq!(parse_failure("fixedStep chrom=c start=1 step=2.5\n").0, 1);
    }

    #[test]
    fn fixed_step_start_one_past_u32_max_converts_to_u32_max() {
        let events = read_all("fixedStep chrom=c start=4294967296 step=1\n").unwrap();
        assert_eq!(events[0], fixed("c", u32::MAX, 1, 1));
    }

    #[test]
    fn fixed_step_start_beyond_u32_is_an_error() {
        assert_eq!(parse_failure("fixedStep chrom=c start=4294967297 step=1\n").0, 1);
    }

    #[test]
    fn fixed_values_are_validated_against_position_overflow() {
        let (line, message) =
            parse_failure("fixedStep chrom=c start=4294967290 step=3 span=1\n1\n1\n1\n");
        assert_eq!(line, 4);
        assert!(message.contains("exceeds"), "{message}");
    }

    #[test]
    fn fixed_values_whose_span_overflows_are_an_error() {
        assert_eq!(parse_failure("fixedStep chrom=c start=4294967295 step=1 span=2\n1\n").0, 2);
    }

    #[test]
    fn fixed_step_huge_step_does_not_panic() {
        let _ = read_all("fixedStep chrom=c start=1 step=4294967295\n1\n1\n1\n");
    }

    #[test]
    fn fixed_value_line_with_two_fields_is_an_error() {
        assert_eq!(parse_failure("fixedStep chrom=c start=1 step=1\n1 2\n").0, 2);
    }

    #[test]
    fn fixed_value_that_is_not_a_number_is_an_error() {
        assert_eq!(parse_failure("fixedStep chrom=c start=1 step=1\n1\nabc\n").0, 3);
    }

    // ---- reader: variableStep ----------------------------------------------------------

    #[test]
    fn variable_step_converts_positions_to_zero_based() {
        let events = read_all("variableStep chrom=chr2 span=3\n10 1.5\n20\t2\n").unwrap();
        assert_eq!(
            events,
            vec![
                variable("chr2", 3),
                Event::VariableValue { start: 9, value: 1.5 },
                Event::VariableValue { start: 19, value: 2.0 },
            ]
        );
    }

    #[test]
    fn variable_step_span_defaults_to_one() {
        assert_eq!(read_all("variableStep chrom=c\n").unwrap(), vec![variable("c", 1)]);
    }

    #[test]
    fn variable_step_requires_chrom() {
        assert!(parse_failure("variableStep span=2\n").1.contains("chrom"));
    }

    #[test]
    fn variable_step_rejects_start() {
        let (line, message) = parse_failure("variableStep chrom=c start=1\n");
        assert_eq!(line, 1);
        assert!(message.contains("start"), "{message}");
    }

    #[test]
    fn variable_step_rejects_step() {
        assert!(parse_failure("variableStep chrom=c step=1\n").1.contains("step"));
    }

    #[test]
    fn variable_step_position_zero_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=c\n0 1\n").0, 2);
    }

    #[test]
    fn variable_step_non_integer_position_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=c\n5.7 1\n").0, 2);
    }

    #[test]
    fn variable_step_negative_position_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=c\n-5 1\n").0, 2);
    }

    #[test]
    fn variable_step_position_plus_span_overflow_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=c span=2\n4294967295 1\n").0, 2);
    }

    #[test]
    fn variable_step_value_line_with_one_field_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=c\n5\n").0, 2);
    }

    #[test]
    fn variable_step_value_line_with_three_fields_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=c\n5 1 2\n").0, 2);
    }

    // ---- reader: headers ---------------------------------------------------------------

    #[test]
    fn unknown_header_key_is_an_error() {
        let (line, message) = parse_failure("variableStep chrom=c color=red\n");
        assert_eq!(line, 1);
        assert!(message.contains("color"), "{message}");
    }

    #[test]
    fn repeated_header_key_is_an_error() {
        let (_, message) = parse_failure("variableStep chrom=c span=1 span=2\n");
        assert!(message.contains("repeated"), "{message}");
    }

    #[test]
    fn repeated_chrom_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=a chrom=b\n").0, 1);
    }

    #[test]
    fn header_word_without_equals_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=c span\n").0, 1);
    }

    #[test]
    fn header_key_with_empty_value_is_an_error() {
        assert_eq!(parse_failure("variableStep chrom=\n").0, 1);
    }

    #[test]
    fn header_fields_may_be_separated_by_tabs_and_runs_of_spaces() {
        let events = read_all("variableStep\tchrom=c   span=2\n").unwrap();
        assert_eq!(events, vec![variable("c", 2)]);
    }

    // ---- reader: structure -------------------------------------------------------------

    #[test]
    fn data_before_any_declaration_is_an_error() {
        let (line, message) = parse_failure("# c\n\n5\n");
        assert_eq!(line, 3);
        assert!(message.contains("declaration"), "{message}");
    }

    #[test]
    fn two_field_data_before_any_declaration_is_an_error() {
        assert_eq!(parse_failure("10 1\n").0, 1);
    }

    #[test]
    fn bedgraph_line_inside_wig_is_a_bedgraph_event() {
        let events = read_all("chr1 0 10 2.5\n").unwrap();
        assert_eq!(events, vec![bedgraph_event("chr1", 0, 10, 2.5)]);
    }

    #[test]
    fn bedgraph_lines_mix_freely_with_stepped_blocks() {
        let text = "fixedStep chrom=a start=1 step=1\n1\nb 0 5 2\nvariableStep chrom=c\n3 4\n";
        assert_eq!(
            read_all(text).unwrap(),
            vec![
                fixed("a", 0, 1, 1),
                Event::FixedValue(1.0),
                bedgraph_event("b", 0, 5, 2.0),
                variable("c", 1),
                Event::VariableValue { start: 2, value: 4.0 },
            ]
        );
    }

    #[test]
    fn bedgraph_line_ends_the_open_block() {
        let (line, _) = parse_failure("fixedStep chrom=a start=1 step=1\nb 0 5 2\n7\n");
        assert_eq!(line, 3);
    }

    #[test]
    fn malformed_bedgraph_line_in_wig_is_an_error_with_line_number() {
        assert_eq!(parse_failure("variableStep chrom=c\nb 1.5 5 2\n").0, 2);
    }

    #[test]
    fn a_new_declaration_restarts_fixed_positions() {
        let text =
            "fixedStep chrom=a start=1 step=10\n1\n2\nfixedStep chrom=a start=101 step=10\n3\n";
        let events = read_all(text).unwrap();
        assert_eq!(events[3], fixed("a", 100, 10, 10));
    }

    #[test]
    fn blank_and_comment_lines_are_skipped_everywhere() {
        let text = "\n# top\nfixedStep chrom=a start=1 step=1\n\n  # mid\n1\n\t\n2\n#end\n";
        let events = read_all(text).unwrap();
        assert_eq!(events.len(), 3);
    }

    #[test]
    fn leading_track_and_browser_lines_are_skipped() {
        let text = "browser position chr1:1-9\ntrack type=wiggle_0\nvariableStep chrom=c\n1 1\n";
        let mut reader = Reader::new(text.as_bytes());
        let events: Vec<Event> = reader.by_ref().collect::<Result<_>>().unwrap();
        assert_eq!(events.len(), 2);
        assert_eq!(reader.track_lines(), ["track type=wiggle_0"]);
    }

    #[test]
    fn track_line_after_a_declaration_is_an_error() {
        let (line, message) = parse_failure("variableStep chrom=c\ntrack name=x\n");
        assert_eq!(line, 2);
        assert!(message.contains("multiple tracks"), "{message}");
    }

    #[test]
    fn chromosome_named_trackx_in_a_bedgraph_line_is_data() {
        assert_eq!(read_all("trackX 0 5 1\n").unwrap(), vec![bedgraph_event("trackX", 0, 5, 1.0)]);
    }

    #[test]
    fn unsorted_data_is_not_rejected() {
        let text = "variableStep chrom=c\n50 1\n10 1\n";
        assert_eq!(read_all(text).unwrap().len(), 3);
    }

    #[test]
    fn crlf_line_endings_are_accepted() {
        let events = read_all("fixedStep chrom=c start=1 step=1\r\n1\r\n").unwrap();
        assert_eq!(events, vec![fixed("c", 0, 1, 1), Event::FixedValue(1.0)]);
    }

    #[test]
    fn nan_and_infinite_values_are_accepted() {
        let events = read_all("variableStep chrom=c\n1 nan\n2 inf\n").unwrap();
        assert!(matches!(events[1], Event::VariableValue { value, .. } if value.is_nan()));
        assert!(matches!(events[2], Event::VariableValue { value, .. } if value == f32::INFINITY));
    }

    #[test]
    fn read_event_reuses_the_event_and_returns_false_at_eof() {
        let mut reader = Reader::new("variableStep chrom=c\n1 2\n".as_bytes());
        let mut event = Event::FixedValue(0.0);
        assert!(reader.read_event(&mut event).unwrap());
        assert_eq!(event, variable("c", 1));
        assert!(reader.read_event(&mut event).unwrap());
        assert_eq!(event, Event::VariableValue { start: 0, value: 2.0 });
        assert!(!reader.read_event(&mut event).unwrap());
    }

    #[test]
    fn reader_continues_after_a_parse_error() {
        let reader = Reader::new("variableStep chrom=c\nbad line here\n3 4\n".as_bytes());
        let results: Vec<_> = reader.collect();
        assert!(results[0].is_ok() && results[1].is_err() && results[2].is_ok());
    }

    #[test]
    fn values_after_a_bad_fixed_step_declaration_are_rejected() {
        let text =
            "fixedStep chrom=chr1 start=1 step=1\n1\nfixedStep chrom=chr2 start=0 step=1\n5\n6\n";
        let results: Vec<_> = Reader::new(text.as_bytes()).collect();
        assert!(results[2].is_err(), "start=0 is invalid");
        for result in &results[3..] {
            match result {
                Err(Error::Parse { message, .. }) => assert!(message.contains("before any")),
                other => panic!("expected a data-before-declaration error, got {other:?}"),
            }
        }
    }

    #[test]
    fn values_after_a_bad_variable_step_declaration_are_rejected() {
        let text = "fixedStep chrom=chr1 start=1 step=1\n1\nvariableStep span=1\n5\n";
        let results: Vec<_> = Reader::new(text.as_bytes()).collect();
        assert!(results[2].is_err(), "chrom= is required");
        assert!(results[3].is_err());
    }

    #[test]
    fn values_after_a_malformed_bedgraph_line_are_rejected() {
        let text = "fixedStep chrom=chr1 start=1 step=1\n1\nchr1 x 10 1.0\n7\n";
        let results: Vec<_> = Reader::new(text.as_bytes()).collect();
        assert!(results[2].is_err());
        assert!(
            matches!(&results[3], Err(Error::Parse { message, .. }) if message.contains("before any"))
        );
    }

    #[test]
    fn a_bad_fixed_step_value_still_advances_the_position() {
        let mut reader = Reader::new("fixedStep chrom=chr1 start=1 step=2\n1\nabc\n3\n".as_bytes());
        let mut event = Event::FixedValue(0.0);
        assert!(reader.read_event(&mut event).unwrap());
        assert!(reader.read_event(&mut event).unwrap());
        assert!(reader.read_event(&mut event).is_err());
        assert!(matches!(reader.section, Section::Fixed { next_start: 4, .. }));
    }

    #[test]
    fn a_wrong_field_count_in_fixed_step_still_advances_the_position() {
        let mut reader = Reader::new("fixedStep chrom=chr1 start=1 step=1\n1 2\n".as_bytes());
        let mut event = Event::FixedValue(0.0);
        assert!(reader.read_event(&mut event).unwrap());
        assert!(reader.read_event(&mut event).is_err());
        assert!(matches!(reader.section, Section::Fixed { next_start: 1, .. }));
    }

    #[test]
    fn invalid_utf8_is_a_parse_error_with_line_number() {
        let mut reader = Reader::new(&b"variableStep chrom=c\n1\xff 2\n"[..]);
        let mut event = Event::FixedValue(0.0);
        assert!(reader.read_event(&mut event).unwrap());
        match reader.read_event(&mut event) {
            Err(Error::Parse { line, .. }) => assert_eq!(line, 2),
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    #[test]
    fn arbitrary_garbage_does_not_panic() {
        for text in [
            "\0",
            "fixedStep",
            "variableStep=",
            "fixedStep chrom= start= step=",
            "fixedStep chrom=c start=18446744073709551616 step=1",
            "fixedStep chrom=c start=1 step=1\n1 2 3 4 5 6\n",
            "variableStep chrom=c span=4294967295\n1 1\n",
            "a b c d e",
            "= = =",
        ] {
            let _ = read_all(text);
        }
    }

    // ---- writer ------------------------------------------------------------------------

    #[test]
    fn fixed_step_header_has_one_based_start_and_explicit_span() {
        let text = written(|w| w.fixed_step("chr1", 0, 10, 5));
        assert_eq!(text, "fixedStep chrom=chr1 start=1 step=10 span=5\n");
    }

    #[test]
    fn variable_step_header_has_explicit_span() {
        let text = written(|w| w.variable_step("chr1", 1));
        assert_eq!(text, "variableStep chrom=chr1 span=1\n");
    }

    #[test]
    fn variable_values_are_written_one_based() {
        let text = written(|w| {
            w.variable_step("c", 2)?;
            w.variable_value(0, 1.5)?;
            w.variable_value(9, 2.0)
        });
        assert_eq!(text, "variableStep chrom=c span=2\n1\t1.5\n10\t2\n");
    }

    #[test]
    fn comment_line_is_hash_followed_by_text() {
        let text = written(|w| w.comment("bedGraph section chr1:1-100"));
        assert_eq!(text, "#bedGraph section chr1:1-100\n");
    }

    #[test]
    fn bedgraph_record_is_written_zero_based() {
        let text = written(|w| w.bedgraph_record("chr1", 10, 20, 3.0));
        assert_eq!(text, "chr1\t10\t20\t3\n");
    }

    #[test]
    fn writer_rejects_bad_chromosome_names() {
        let mut writer = Writer::new(Vec::new());
        assert!(matches!(writer.fixed_step("a b", 0, 1, 1), Err(Error::InvalidItem { .. })));
        assert!(matches!(writer.variable_step("", 1), Err(Error::InvalidItem { .. })));
        assert!(matches!(writer.bedgraph_record("#x", 0, 1, 1.0), Err(Error::InvalidItem { .. })));
        assert!(writer.into_inner().unwrap().is_empty());
    }

    #[test]
    fn writer_rejects_zero_step_and_span() {
        let mut writer = Writer::new(Vec::new());
        assert!(matches!(writer.fixed_step("c", 0, 0, 1), Err(Error::InvalidItem { .. })));
        assert!(matches!(writer.fixed_step("c", 0, 1, 0), Err(Error::InvalidItem { .. })));
        assert!(matches!(writer.variable_step("c", 0), Err(Error::InvalidItem { .. })));
    }

    #[test]
    fn writer_rejects_values_outside_a_matching_block() {
        let mut writer = Writer::new(Vec::new());
        assert!(writer.fixed_value(1.0).is_err());
        assert!(writer.variable_value(0, 1.0).is_err());
        writer.fixed_step("c", 0, 1, 1).unwrap();
        assert!(writer.variable_value(0, 1.0).is_err());
        writer.bedgraph_record("c", 0, 1, 1.0).unwrap();
        assert!(writer.fixed_value(1.0).is_err());
    }

    #[test]
    fn comment_does_not_end_the_current_block() {
        let text = written(|w| {
            w.fixed_step("c", 0, 1, 1)?;
            w.fixed_value(1.0)?;
            w.comment(" note")?;
            w.fixed_value(2.0)
        });
        assert!(text.ends_with("1\n# note\n2\n"), "{text}");
    }

    #[test]
    fn comment_with_a_line_break_is_rejected() {
        let mut writer = Writer::new(Vec::new());
        assert!(writer.comment("a\nb").is_err());
        assert!(writer.comment("a\rb").is_err());
        assert!(writer.into_inner().unwrap().is_empty());
    }

    #[test]
    fn fixed_step_round_trips() {
        let text = written(|w| {
            w.fixed_step("chr1", 99, 10, 4)?;
            w.fixed_value(1.5)?;
            w.fixed_value(-2.0)?;
            w.fixed_value(1e-7)
        });
        assert_eq!(
            read_all(&text).unwrap(),
            vec![
                fixed("chr1", 99, 10, 4),
                Event::FixedValue(1.5),
                Event::FixedValue(-2.0),
                Event::FixedValue(1e-7),
            ]
        );
    }

    #[test]
    fn variable_step_round_trips() {
        let text = written(|w| {
            w.variable_step("chr2", 3)?;
            w.variable_value(0, 1.0)?;
            w.variable_value(1000, 0.25)?;
            w.variable_value(u32::MAX - 3, 7.0)
        });
        assert_eq!(
            read_all(&text).unwrap(),
            vec![
                variable("chr2", 3),
                Event::VariableValue { start: 0, value: 1.0 },
                Event::VariableValue { start: 1000, value: 0.25 },
                Event::VariableValue { start: u32::MAX - 3, value: 7.0 },
            ]
        );
    }

    #[test]
    fn mixed_blocks_comments_and_bedgraph_round_trip() {
        let text = written(|w| {
            w.fixed_step("chr1", 0, 1, 1)?;
            w.fixed_value(1.0)?;
            w.comment("bedGraph section chr1:11-30")?;
            w.bedgraph_record("chr1", 10, 20, 2.0)?;
            w.bedgraph_record("chr1", 20, 30, 3.0)?;
            w.variable_step("chr2", 2)?;
            w.variable_value(4, 5.0)?;
            w.fixed_step("chr3", u32::MAX - 1, 1, 1)?;
            w.fixed_value(6.0)
        });
        assert_eq!(
            read_all(&text).unwrap(),
            vec![
                fixed("chr1", 0, 1, 1),
                Event::FixedValue(1.0),
                bedgraph_event("chr1", 10, 20, 2.0),
                bedgraph_event("chr1", 20, 30, 3.0),
                variable("chr2", 2),
                Event::VariableValue { start: 4, value: 5.0 },
                fixed("chr3", u32::MAX - 1, 1, 1),
                Event::FixedValue(6.0),
            ]
        );
    }
}
