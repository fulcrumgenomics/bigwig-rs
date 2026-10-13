//! bedGraph text: a streaming reader and a writer.
//!
//! bedGraph lines are `chrom start end value`, 0-based half-open, separated by tabs or spaces.
//! The reader is deliberately forgiving about layout (comments, blank lines, leading `track`
//! lines) and strict about numbers; it does not check ordering, which is the bigWig writer's job.

use std::io::{BufRead, Write};

use crate::text::{
    LineReader, check_chrom, parse_chrom, parse_error, parse_u32, parse_value, split_fields,
    write_bedgraph_line,
};
use crate::{Error, Result};

/// One bedGraph line: a value over the 0-based half-open range `[start, end)` of `chrom`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Record {
    /// Chromosome name.
    pub chrom: String,
    /// 0-based inclusive start.
    pub start: u32,
    /// 0-based exclusive end.
    pub end: u32,
    /// The value over the range; may be NaN or infinite, which the bigWig writer rejects.
    pub value: f32,
}

/// Streaming bedGraph reader over any [`BufRead`] (plain text; decompress gzip/bgzf first).
///
/// Blank lines and lines whose first non-blank character is `#` are skipped anywhere. `track` and
/// `browser` lines are accepted only before the first data line, matched on the exact first
/// word (so a chromosome named `trackX` is data); a later one is an [`Error::Parse`]. Data
/// lines must have exactly four fields separated by runs of tabs or spaces, with strict
/// unsigned-integer coordinates. A parse error does not make the reader unusable.
///
/// ```
/// use bigwig::bedgraph::{Reader, Record};
///
/// let text = "track type=bedGraph\nchr1\t0\t10\t1.5\nchr1 10 20 2\n";
/// let records: Vec<Record> = Reader::new(text.as_bytes()).collect::<Result<_, _>>()?;
/// assert_eq!(records.len(), 2);
/// assert_eq!(records[1], Record { chrom: "chr1".into(), start: 10, end: 20, value: 2.0 });
/// # Ok::<(), bigwig::Error>(())
/// ```
pub struct Reader<R> {
    lines: LineReader<R>,
    /// Set after an I/O error so that iterating does not spin on a persistent failure.
    finished: bool,
}

impl<R: BufRead> Reader<R> {
    /// Creates a reader over `inner`.
    pub fn new(inner: R) -> Self {
        Self { lines: LineReader::new(inner), finished: false }
    }

    /// Reads the next record into `record`, reusing its `String`. Returns `Ok(false)` at end of
    /// input. On error the contents of `record` are unchanged.
    pub fn read_record(&mut self, record: &mut Record) -> Result<bool> {
        if !self.lines.next_line()? {
            return Ok(false);
        }
        parse_record(self.lines.line(), self.lines.line_number(), record)?;
        Ok(true)
    }

    /// The leading `track` lines seen so far, verbatim and trimmed, for callers that want the
    /// display hints they carry.
    pub fn track_lines(&self) -> &[String] {
        &self.lines.track_lines
    }
}

impl<R: BufRead> Iterator for Reader<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let mut record = Record::default();
        match self.read_record(&mut record) {
            Ok(true) => Some(Ok(record)),
            Ok(false) => None,
            Err(error) => {
                self.finished = matches!(error, Error::Io(_));
                Some(Err(error))
            }
        }
    }
}

/// Writes `chrom\tstart\tend\tvalue\n` lines, with values in the shortest form that parses back
/// to the same `f32`.
///
/// Each line reaches the underlying writer in a single `write_all`; wrap unbuffered sinks in a
/// [`std::io::BufWriter`].
///
/// ```
/// use bigwig::bedgraph::Writer;
///
/// let mut writer = Writer::new(Vec::new());
/// writer.write_record("chr1", 0, 10, 1.5)?;
/// writer.write_record("chr1", 10, 20, 2.0)?;
/// assert_eq!(writer.into_inner()?, b"chr1\t0\t10\t1.5\nchr1\t10\t20\t2\n");
/// # Ok::<(), bigwig::Error>(())
/// ```
pub struct Writer<W> {
    inner: W,
    line: Vec<u8>,
}

impl<W: Write> Writer<W> {
    /// Creates a writer over `inner`.
    pub fn new(inner: W) -> Self {
        Self { inner, line: Vec::new() }
    }

    /// Writes one line. Fails with [`Error::InvalidItem`] if `chrom` could not be read back as a
    /// chromosome name (empty, containing whitespace or control characters, starting with `#`,
    /// or exactly `track` or `browser`); nothing is written in that case. Neither ordering nor
    /// `start < end` is checked.
    pub fn write_record(&mut self, chrom: &str, start: u32, end: u32, value: f32) -> Result<()> {
        check_chrom(chrom, u64::from(start))?;
        self.line.clear();
        write_bedgraph_line(&mut self.line, chrom, start, end, value)?;
        self.inner.write_all(&self.line)?;
        Ok(())
    }

    /// Flushes and returns the underlying writer.
    pub fn into_inner(mut self) -> Result<W> {
        self.inner.flush()?;
        Ok(self.inner)
    }
}

/// Parses a four-field bedGraph line into `record`; `record` is untouched on error.
pub(crate) fn parse_record(line: &str, line_number: u64, record: &mut Record) -> Result<()> {
    let (count, fields) = split_fields(line);
    if count != 4 {
        let found = if count > 4 { "more".to_string() } else { count.to_string() };
        return Err(parse_error(line_number, format!("expected 4 fields, found {found}")));
    }
    let chrom = parse_chrom(fields[0], line_number)?;
    let start = parse_u32(fields[1], line_number, "start")?;
    let end = parse_u32(fields[2], line_number, "end")?;
    let value = parse_value(fields[3], line_number)?;
    record.chrom.clear();
    record.chrom.push_str(chrom);
    record.start = start;
    record.end = end;
    record.value = value;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_all(text: &str) -> Result<Vec<Record>> {
        Reader::new(text.as_bytes()).collect()
    }

    fn parse_failure(text: &str) -> (u64, String) {
        match read_all(text) {
            Err(Error::Parse { line, message }) => (line, message),
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    fn record(chrom: &str, start: u32, end: u32, value: f32) -> Record {
        Record { chrom: chrom.to_string(), start, end, value }
    }

    // ---- reader ----------------------------------------------------------------------

    #[test]
    fn control_character_in_chromosome_name_is_a_parse_error() {
        assert!(parse_failure("chr\u{8}1\t0\t10\t1\n").1.contains("control character"));
    }

    #[test]
    fn unicode_whitespace_in_chromosome_name_is_a_parse_error() {
        assert!(parse_failure("chr\u{a0}1\t0\t10\t1\n").1.contains("whitespace"));
    }

    #[test]
    fn reads_tab_separated_records() {
        let records = read_all("chr1\t0\t10\t1.5\nchr2\t5\t6\t-2\n").unwrap();
        assert_eq!(records, vec![record("chr1", 0, 10, 1.5), record("chr2", 5, 6, -2.0)]);
    }

    #[test]
    fn reads_runs_of_mixed_tabs_and_spaces() {
        let records = read_all("chr1 \t 0\t\t10   1.5\n").unwrap();
        assert_eq!(records, vec![record("chr1", 0, 10, 1.5)]);
    }

    #[test]
    fn reads_crlf_line_endings() {
        let records = read_all("chr1\t0\t10\t1\r\nchr1\t10\t20\t2\r\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[1].value, 2.0);
    }

    #[test]
    fn reads_final_line_without_newline() {
        assert_eq!(read_all("chr1\t0\t10\t1").unwrap().len(), 1);
    }

    #[test]
    fn empty_input_yields_no_records() {
        assert!(read_all("").unwrap().is_empty());
    }

    #[test]
    fn skips_blank_and_comment_lines_even_when_indented() {
        let text = "\n# comment\n   # indented comment\n\t\nchr1\t0\t10\t1\n  \n#end\n";
        assert_eq!(read_all(text).unwrap(), vec![record("chr1", 0, 10, 1.0)]);
    }

    #[test]
    fn skips_comments_between_data_lines() {
        let records = read_all("chr1 0 10 1\n# mid\nchr1 10 20 2\n").unwrap();
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn leading_track_and_browser_lines_are_skipped() {
        let text = "browser position chr1:1-100\ntrack type=bedGraph name=x\nchr1 0 10 1\n";
        assert_eq!(read_all(text).unwrap(), vec![record("chr1", 0, 10, 1.0)]);
    }

    #[test]
    fn track_lines_are_kept_verbatim() {
        let text = "track type=bedGraph name=x\nbrowser hide all\nchr1 0 1 1\n";
        let mut reader = Reader::new(text.as_bytes());
        assert!(reader.read_record(&mut Record::default()).unwrap());
        assert_eq!(reader.track_lines(), ["track type=bedGraph name=x"]);
    }

    #[test]
    fn track_line_after_comments_but_before_data_is_accepted() {
        assert_eq!(read_all("# hi\ntrack name=x\nchr1 0 1 1\n").unwrap().len(), 1);
    }

    #[test]
    fn chromosome_named_trackx_is_data() {
        assert_eq!(read_all("trackX 0 10 1\n").unwrap(), vec![record("trackX", 0, 10, 1.0)]);
    }

    #[test]
    fn track_line_after_data_is_an_error_on_that_line() {
        let (line, message) = parse_failure("chr1 0 10 1\n\ntrack name=x\n");
        assert_eq!(line, 3);
        assert!(message.contains("multiple tracks"), "{message}");
    }

    #[test]
    fn browser_line_after_data_is_an_error() {
        assert_eq!(parse_failure("chr1 0 10 1\nbrowser hide all\n").0, 2);
    }

    #[test]
    fn too_few_fields_is_an_error_with_line_number() {
        let (line, message) = parse_failure("chr1 0 10 1\nchr1 10 20\n");
        assert_eq!(line, 2);
        assert!(message.contains("found 3"), "{message}");
    }

    #[test]
    fn too_many_fields_is_an_error() {
        assert_eq!(parse_failure("chr1 0 10 1 extra\n").0, 1);
    }

    #[test]
    fn fractional_coordinate_is_an_error() {
        assert_eq!(parse_failure("chr1 10.5 20 1\n").0, 1);
    }

    #[test]
    fn negative_coordinate_is_an_error() {
        assert_eq!(parse_failure("\nchr1 -5 20 1\n").0, 2);
    }

    #[test]
    fn exponent_coordinate_is_an_error() {
        assert_eq!(parse_failure("chr1 0 1e3 1\n").0, 1);
    }

    #[test]
    fn plus_signed_coordinate_is_an_error() {
        assert_eq!(parse_failure("chr1 +5 20 1\n").0, 1);
    }

    #[test]
    fn coordinate_above_u32_is_an_error() {
        assert_eq!(parse_failure("chr1 0 4294967296 1\n").0, 1);
    }

    #[test]
    fn coordinate_u32_max_is_accepted() {
        assert_eq!(read_all("chr1 0 4294967295 1\n").unwrap()[0].end, u32::MAX);
    }

    #[test]
    fn non_numeric_value_is_an_error() {
        let (line, message) = parse_failure("chr1 0 10 abc\n");
        assert_eq!(line, 1);
        assert!(message.contains("abc"), "{message}");
    }

    #[test]
    fn nan_and_infinite_values_are_accepted() {
        let records = read_all("c 0 1 nan\nc 1 2 inf\nc 2 3 -inf\n").unwrap();
        assert!(records[0].value.is_nan());
        assert_eq!(records[1].value, f32::INFINITY);
        assert_eq!(records[2].value, f32::NEG_INFINITY);
    }

    #[test]
    fn scientific_values_are_accepted() {
        assert_eq!(read_all("c 0 1 1e3\n").unwrap()[0].value, 1000.0);
    }

    #[test]
    fn unsorted_and_overlapping_records_are_not_rejected() {
        assert_eq!(read_all("b 50 60 1\na 10 20 1\na 5 30 1\n").unwrap().len(), 3);
    }

    #[test]
    fn start_after_end_is_not_rejected() {
        assert_eq!(read_all("c 20 10 1\n").unwrap().len(), 1);
    }

    #[test]
    fn invalid_utf8_is_a_parse_error_with_line_number() {
        let mut reader = Reader::new(&b"chr1 0 1 1\nch\xff 0 1 1\n"[..]);
        let mut record = Record::default();
        assert!(reader.read_record(&mut record).unwrap());
        match reader.read_record(&mut record) {
            Err(Error::Parse { line, .. }) => assert_eq!(line, 2),
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    #[test]
    fn reader_continues_after_a_parse_error() {
        let reader = Reader::new("chr1 0 1 1\nbad\nchr1 1 2 2\n".as_bytes());
        let results: Vec<_> = reader.collect();
        assert!(results[0].is_ok() && results[1].is_err() && results[2].is_ok());
    }

    #[test]
    fn read_record_leaves_record_unchanged_on_error() {
        let mut reader = Reader::new("chr2 1 x 3\n".as_bytes());
        let mut target = record("chr1", 0, 10, 1.0);
        assert!(reader.read_record(&mut target).is_err());
        assert_eq!(target, record("chr1", 0, 10, 1.0));
    }

    #[test]
    fn read_record_reuses_the_record_and_returns_false_at_eof() {
        let mut reader = Reader::new("chr1 0 1 1\nchr2 1 2 2\n".as_bytes());
        let mut target = Record::default();
        assert!(reader.read_record(&mut target).unwrap());
        assert_eq!(target.chrom, "chr1");
        assert!(reader.read_record(&mut target).unwrap());
        assert_eq!(target.chrom, "chr2");
        assert!(!reader.read_record(&mut target).unwrap());
    }

    #[test]
    fn arbitrary_garbage_does_not_panic() {
        for text in ["\0\0\0", "\t\t\t", "track", "chr1 \u{a0} 1 1", "####", "a b c d\ne f g h\n"] {
            let _ = read_all(text);
        }
    }

    // ---- writer ----------------------------------------------------------------------

    #[test]
    fn writes_tab_separated_lines() {
        let mut writer = Writer::new(Vec::new());
        writer.write_record("chr1", 0, 10, 1.5).unwrap();
        writer.write_record("chr2", 3, 4, -5.0).unwrap();
        assert_eq!(writer.into_inner().unwrap(), b"chr1\t0\t10\t1.5\nchr2\t3\t4\t-5\n");
    }

    #[test]
    fn writer_rejects_chromosome_names_that_cannot_be_read_back() {
        let mut writer = Writer::new(Vec::new());
        for chrom in ["", "a b", "a\tb", "#x", "track", "browser", "a\nb"] {
            assert!(
                matches!(writer.write_record(chrom, 0, 1, 1.0), Err(Error::InvalidItem { .. })),
                "{chrom:?}"
            );
        }
        assert!(writer.into_inner().unwrap().is_empty());
    }

    #[test]
    fn writer_accepts_chromosome_names_that_merely_resemble_keywords() {
        let mut writer = Writer::new(Vec::new());
        writer.write_record("trackX", 0, 1, 1.0).unwrap();
        writer.write_record("chr#1", 0, 1, 1.0).unwrap();
    }

    #[test]
    fn write_then_read_round_trips_records() {
        let written = vec![
            record("chr1", 0, 10, 1.5),
            record("chr1", 10, 4_294_967_295, -0.000_123_45),
            record("chr2", 7, 8, 3.402_823_5e38),
            record("chr2", 8, 9, 1.0e-30),
            record("chrM", 0, 1, 0.0),
        ];
        let mut writer = Writer::new(Vec::new());
        for r in &written {
            writer.write_record(&r.chrom, r.start, r.end, r.value).unwrap();
        }
        let bytes = writer.into_inner().unwrap();
        let read: Vec<Record> = Reader::new(&bytes[..]).collect::<Result<_>>().unwrap();
        assert_eq!(read, written);
    }
}
