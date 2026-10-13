//! Line-oriented text helpers shared by the bedGraph and wig readers and writers.

use std::io::{self, BufRead, Write};
use std::str::FromStr;

use crate::{Error, Result};

/// skipping, and the rules for leading `track`/`browser` lines.
pub(crate) struct LineReader<R> {
    inner: R,
    line: String,
    line_number: u64,
    seen_data: bool,
    pub(crate) track_lines: Vec<String>,
}

impl<R: BufRead> LineReader<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            line: String::new(),
            line_number: 0,
            seen_data: false,
            track_lines: Vec::new(),
        }
    }

    /// Advances to the next data (or wig declaration) line, returning `false` at end of input.
    /// The line is then available from [`Self::line`].
    pub(crate) fn next_line(&mut self) -> Result<bool> {
        loop {
            self.line.clear();
            self.line_number += 1;
            match self.inner.read_line(&mut self.line) {
                Ok(0) => {
                    self.line_number -= 1;
                    return Ok(false);
                }
                Ok(_) => {}
                Err(error) if error.kind() == io::ErrorKind::InvalidData => {
                    return Err(parse_error(self.line_number, "line is not valid UTF-8"));
                }
                Err(error) => return Err(error.into()),
            }
            let trimmed = self.line.trim_matches(is_separator);
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let first_word = trimmed.split(is_separator).next().unwrap_or("");
            if first_word == "track" || first_word == "browser" {
                if self.seen_data {
                    return Err(parse_error(
                        self.line_number,
                        format!(
                            "`{first_word}` line after data; multiple tracks per file are not supported"
                        ),
                    ));
                }
                if first_word == "track" {
                    self.track_lines.push(trimmed.to_string());
                }
                continue;
            }
            self.seen_data = true;
            return Ok(true);
        }
    }

    /// The current line with surrounding separators removed.
    pub(crate) fn line(&self) -> &str {
        self.line.trim_matches(is_separator)
    }

    /// The 1-based number of the current line.
    pub(crate) fn line_number(&self) -> u64 {
        self.line_number
    }
}

/// Builds an [`Error::Parse`].
pub(crate) fn parse_error(line: u64, message: impl Into<String>) -> Error {
    Error::Parse { line, message: message.into() }
}

/// True for the characters that separate fields. `\r` is included so CRLF files parse.
pub(crate) fn is_separator(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\r' | '\n')
}

/// Splits `line` into its fields, returning the total field count and the first five of them.
/// Five is one more than any line type needs, enough to tell "too many" from "exactly right".
pub(crate) fn split_fields(line: &str) -> (usize, [&str; 5]) {
    let mut fields = [""; 5];
    let mut count = 0;
    for field in line.split(is_separator).filter(|field| !field.is_empty()) {
        if count < fields.len() {
            fields[count] = field;
        }
        count += 1;
    }
    (count, fields)
}

/// Parses a non-negative integer made only of ASCII digits. Rust's own parser would also accept
/// a leading `+`, which neither format allows.
pub(crate) fn parse_unsigned(field: &str, line: u64, what: &str) -> Result<u64> {
    if field.is_empty() || !field.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(parse_error(
            line,
            format!("invalid {what} {field:?}: expected an unsigned integer"),
        ));
    }
    field.parse::<u64>().map_err(|_| parse_error(line, format!("{what} {field:?} is out of range")))
}

/// Parses an unsigned integer that must fit in a `u32`.
pub(crate) fn parse_u32(field: &str, line: u64, what: &str) -> Result<u32> {
    u32::try_from(parse_unsigned(field, line, what)?)
        .map_err(|_| parse_error(line, format!("{what} {field:?} is out of range")))
}

/// Parses a data value with `f32::from_str`, which accepts `nan` and `inf`.
pub(crate) fn parse_value(field: &str, line: u64) -> Result<f32> {
    f32::from_str(field).map_err(|_| parse_error(line, format!("invalid value {field:?}")))
}

/// Validates a chromosome name read from text. Fields are split on tabs and spaces only, so a
/// name can still carry other control or whitespace characters, which [`bedgraph::Writer`] refuses to write.
/// Rejecting them on read keeps every name the readers return writable.
pub(crate) fn parse_chrom(field: &str, line: u64) -> Result<&str> {
    if field.chars().any(|c| c.is_whitespace() || c.is_control()) {
        return Err(parse_error(
            line,
            format!("chromosome name {field:?} contains whitespace or a control character"),
        ));
    }
    Ok(field)
}

/// Rejects chromosome names the text readers could not read back: empty, containing whitespace or
/// control characters, starting with `#`, or exactly `track` or `browser`. `position` is only
/// reported in the error.
pub(crate) fn check_chrom(chrom: &str, position: u64) -> Result<()> {
    let reason = if chrom.is_empty() {
        "chromosome name is empty"
    } else if chrom.chars().any(|c| c.is_whitespace() || c.is_control()) {
        "chromosome name contains whitespace or a control character"
    } else if chrom.starts_with('#') {
        "chromosome name starts with '#'"
    } else if chrom == "track" || chrom == "browser" {
        "chromosome name would be read as a track or browser line"
    } else {
        return Ok(());
    };
    Err(Error::InvalidItem { chrom: chrom.to_string(), start: position, reason })
}

/// Appends `chrom\tstart\tend\tvalue\n` to `out`.
pub(crate) fn write_bedgraph_line(
    out: &mut Vec<u8>,
    chrom: &str,
    start: u32,
    end: u32,
    value: f32,
) -> io::Result<()> {
    write!(out, "{chrom}\t{start}\t{end}\t")?;
    write_value(out, value)?;
    out.write_all(b"\n")
}

/// Writes `value` in the shortest form that parses back to the same `f32`: no fractional part
/// for integral values (`5`, not `5.0`), and scientific notation only when the magnitude is
/// outside `[1e-5, 1e16)`. Zero, NaN and the infinities use their plain forms.
pub(crate) fn write_value(out: &mut impl Write, value: f32) -> io::Result<()> {
    let magnitude = value.abs();
    if !value.is_finite() || value == 0.0 || (1e-5..1e16).contains(&magnitude) {
        write!(out, "{value}")
    } else {
        write!(out, "{value:e}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn formatted(value: f32) -> String {
        let mut out = Vec::new();
        write_value(&mut out, value).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn integral_values_have_no_fractional_part() {
        assert_eq!(formatted(5.0), "5");
        assert_eq!(formatted(-5.0), "-5");
        assert_eq!(formatted(100.0), "100");
        assert_eq!(formatted(0.0), "0");
    }

    #[test]
    fn negative_zero_keeps_its_sign() {
        assert_eq!(formatted(-0.0), "-0");
    }

    #[test]
    fn fractional_values_use_the_shortest_round_trip_form() {
        assert_eq!(formatted(0.1), "0.1");
        assert_eq!(formatted(1.5), "1.5");
        assert_eq!(formatted(-2.25), "-2.25");
        assert_eq!(formatted(1.234_567_9), "1.2345679");
    }

    #[test]
    fn one_e_minus_five_is_plain_decimal() {
        assert_eq!(formatted(1e-5), "0.00001");
        assert_eq!(formatted(-1e-5), "-0.00001");
    }

    #[test]
    fn values_just_below_one_e_minus_five_use_scientific_notation() {
        let below = f32::from_bits(1e-5f32.to_bits() - 1);
        let text = formatted(below);
        assert!(text.contains('e'), "{text}");
        assert_eq!(text.parse::<f32>().unwrap(), below);
        assert_eq!(formatted(1e-6), "1e-6");
        assert_eq!(formatted(-1e-6), "-1e-6");
    }

    #[test]
    fn one_e_sixteen_uses_scientific_notation() {
        assert_eq!(formatted(1e16), "1e16");
        assert_eq!(formatted(-1e16), "-1e16");
        assert_eq!(formatted(1.5e20), "1.5e20");
    }

    #[test]
    fn values_just_below_one_e_sixteen_are_plain_decimal() {
        let below = f32::from_bits(1e16f32.to_bits() - 1);
        let text = formatted(below);
        assert!(!text.contains('e'), "{text}");
        assert_eq!(text.parse::<f32>().unwrap(), below);
        assert_eq!(formatted(1e15), "1000000000000000");
    }

    #[test]
    fn extreme_magnitudes_round_trip() {
        for value in [
            f32::MIN_POSITIVE,
            -f32::MIN_POSITIVE,
            f32::MIN_POSITIVE / 2.0,
            f32::from_bits(1),
            f32::MAX,
            f32::MIN,
            f32::EPSILON,
        ] {
            assert_eq!(formatted(value).parse::<f32>().unwrap(), value, "{value:e}");
        }
    }

    #[test]
    fn min_positive_and_smallest_subnormal_are_scientific() {
        assert_eq!(formatted(f32::MIN_POSITIVE), "1.1754944e-38");
        assert_eq!(formatted(f32::from_bits(1)), "1e-45");
    }

    #[test]
    fn f32_max_is_scientific() {
        assert_eq!(formatted(f32::MAX), "3.4028235e38");
    }

    #[test]
    fn non_finite_values_use_forms_the_reader_accepts() {
        for value in [f32::INFINITY, f32::NEG_INFINITY] {
            assert_eq!(formatted(value).parse::<f32>().unwrap(), value);
        }
        assert!(formatted(f32::NAN).parse::<f32>().unwrap().is_nan());
    }

    #[test]
    fn widely_spread_values_round_trip_exactly() {
        for bits in (1..0x7f80_0000u32).step_by(99_991) {
            let value = f32::from_bits(bits);
            assert_eq!(formatted(value).parse::<f32>().unwrap(), value, "{bits:#x}");
        }
    }
}
