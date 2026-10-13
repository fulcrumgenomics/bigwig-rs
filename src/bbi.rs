//! The BBI container layer shared by bigWig and bigBed: magic numbers, byte order, the fixed
//! header, zoom headers and total summary, chromosome metadata, and the bounds-checked binary codec
//! every on-disk structure is read and written with.

use std::io::{self, BufRead, Read, Seek, SeekFrom};

use crate::error::{Error, Result};
use crate::text::{is_separator, parse_error, parse_u32, split_fields};
use crate::zoom::Summary;

/// Magic number at offset 0 (and, since version 4, at the end) of a bigWig file.
pub(crate) const BIGWIG_MAGIC: u32 = 0x888F_FC26;
/// Magic number of a bigBed file, which shares the container but is not read by this crate.
pub(crate) const BIGBED_MAGIC: u32 = 0x8789_F2EB;
/// Magic number of the chromosome B+ tree header.
pub(crate) const BPTREE_MAGIC: u32 = 0x78CA_8C91;
/// Magic number of an R-tree ("cirTree") index header.
pub(crate) const RTREE_MAGIC: u32 = 0x2468_ACE0;

/// The format version this crate writes. Versions 1 up to this one are readable.
pub(crate) const BIGWIG_VERSION: u16 = 4;
/// Size of the fixed header at offset 0.
pub(crate) const FIXED_HEADER_SIZE: usize = 64;
/// Size of one zoom header; `zoomLevels` of them follow the fixed header.
pub(crate) const ZOOM_HEADER_SIZE: usize = 24;
/// Size of the total summary at `totalSummaryOffset`.
pub(crate) const TOTAL_SUMMARY_SIZE: usize = 40;
/// Most node levels the B+ tree and R-tree decoders follow before declaring a file corrupt.
///
/// A tree with fan-out 2 (the smallest that terminates) holding 2^64 items is 64 levels deep, so
/// a deeper chain of child offsets can only come from a corrupt or hostile file.
pub(crate) const MAX_TREE_DEPTH: usize = 64;

/// The 64-byte fixed header at offset 0.
///
/// The magic is not stored: [`FixedHeader::decode`] accepts only the bigWig magic and reports the
/// byte order it was found in, and [`FixedHeader::encode`] always writes it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct FixedHeader {
    /// Format version: 1 to [`BIGWIG_VERSION`].
    pub(crate) version: u16,
    /// Number of zoom headers following the fixed header.
    pub(crate) zoom_levels: u16,
    /// Offset of the chromosome B+ tree.
    pub(crate) chrom_tree_offset: u64,
    /// Offset of `sectionCount` (u64), which the data sections follow.
    pub(crate) full_data_offset: u64,
    /// Offset of the data R-tree index.
    pub(crate) full_index_offset: u64,
    /// bigBed only: number of BED fields; 0 for bigWig.
    pub(crate) field_count: u16,
    /// bigBed only: number of standard BED fields; 0 for bigWig.
    pub(crate) defined_field_count: u16,
    /// bigBed only: offset of the autoSql text; 0 for bigWig.
    pub(crate) auto_sql_offset: u64,
    /// Offset of the total summary, or 0 when the file has none.
    pub(crate) total_summary_offset: u64,
    /// Largest decompressed block size, or 0 when blocks are stored uncompressed.
    pub(crate) uncompress_buf_size: u32,
    /// bigBed only: offset of the extension header; 0 for bigWig.
    pub(crate) extension_offset: u64,
}

impl FixedHeader {
    /// Decodes the header from the start of a file and returns it with the file's byte order,
    /// which is detected from the magic.
    ///
    /// Fails with [`Error::NotBigWig`] when the magic matches in neither byte order, including
    /// when `bytes` is too short to hold one, [`Error::IsBigBed`] for the bigBed magic,
    /// [`Error::UnsupportedVersion`] for a version outside 1 to [`BIGWIG_VERSION`], and
    /// [`Error::Corrupt`] when `bytes` ends after the magic but before the header does.
    pub(crate) fn decode(bytes: &[u8]) -> Result<(Self, ByteOrder)> {
        // Missing magic bytes read as zeros. Neither magic has a zero byte, so input too short to
        // hold one is `NotBigWig`, as other input that isn't a bigWig is.
        let mut magic = [0; 4];
        let present = bytes.len().min(magic.len());
        magic[..present].copy_from_slice(&bytes[..present]);
        let found = u32::from_le_bytes(magic);
        let byte_order = match ByteOrder::detect(found, BIGWIG_MAGIC) {
            Some(order) => order,
            None if ByteOrder::detect(found, BIGBED_MAGIC).is_some() => {
                return Err(Error::IsBigBed);
            }
            None => return Err(Error::NotBigWig { found }),
        };

        let mut reader = ByteReader::new(bytes, byte_order, 0);
        reader.skip(4)?;
        let version = reader.u16()?;
        if !(1..=BIGWIG_VERSION).contains(&version) {
            return Err(Error::UnsupportedVersion(version));
        }
        let header = Self {
            version,
            zoom_levels: reader.u16()?,
            chrom_tree_offset: reader.u64()?,
            full_data_offset: reader.u64()?,
            full_index_offset: reader.u64()?,
            field_count: reader.u16()?,
            defined_field_count: reader.u16()?,
            auto_sql_offset: reader.u64()?,
            total_summary_offset: reader.u64()?,
            uncompress_buf_size: reader.u32()?,
            extension_offset: reader.u64()?,
        };
        Ok((header, byte_order))
    }

    /// Appends the 64-byte header, magic first, to `buf` in byte order `order`.
    pub(crate) fn encode(&self, order: ByteOrder, buf: &mut Vec<u8>) {
        order.put_u32(buf, BIGWIG_MAGIC);
        order.put_u16(buf, self.version);
        order.put_u16(buf, self.zoom_levels);
        order.put_u64(buf, self.chrom_tree_offset);
        order.put_u64(buf, self.full_data_offset);
        order.put_u64(buf, self.full_index_offset);
        order.put_u16(buf, self.field_count);
        order.put_u16(buf, self.defined_field_count);
        order.put_u64(buf, self.auto_sql_offset);
        order.put_u64(buf, self.total_summary_offset);
        order.put_u32(buf, self.uncompress_buf_size);
        order.put_u64(buf, self.extension_offset);
    }
}

/// One zoom level as decoded from the file header: a summary of the data at a coarser resolution.
///
/// A level comes only from a reader and can't be constructed, so handing one back to a query
/// identifies the level exactly, without searching by reduction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ZoomLevel {
    /// Nominal bases per zoom record (`reductionLevel`).
    reduction: u32,
    /// Offset of this level's `zoomCount` and zoom blocks.
    data_offset: u64,
    /// Offset of this level's R-tree index.
    index_offset: u64,
}

impl ZoomLevel {
    /// Creates a level from its zoom-header fields.
    pub(crate) fn new(reduction: u32, data_offset: u64, index_offset: u64) -> Self {
        Self { reduction, data_offset, index_offset }
    }

    /// Decodes one 24-byte zoom header. The reserved field is ignored, as kent does.
    pub(crate) fn decode(reader: &mut ByteReader<'_>) -> Result<Self> {
        let reduction = reader.u32()?;
        reader.skip(4)?;
        Ok(Self { reduction, data_offset: reader.u64()?, index_offset: reader.u64()? })
    }

    /// Nominal bases summarised by each zoom record of this level.
    pub fn reduction(&self) -> u32 {
        self.reduction
    }

    /// Offset of this level's `zoomCount` and zoom blocks.
    pub(crate) fn data_offset(&self) -> u64 {
        self.data_offset
    }

    /// Offset of this level's R-tree index.
    pub(crate) fn index_offset(&self) -> u64 {
        self.index_offset
    }

    /// Appends the 24-byte zoom header, with a zero reserved field, to `buf` in byte order `order`.
    pub(crate) fn encode(&self, order: ByteOrder, buf: &mut Vec<u8>) {
        order.put_u32(buf, self.reduction);
        order.put_u32(buf, 0);
        order.put_u64(buf, self.data_offset);
        order.put_u64(buf, self.index_offset);
    }
}

/// Byte order of every multi-byte integer and float in a bigWig file.
///
/// Readers detect it from the magic number. Writers produce little-endian files unless told
/// otherwise; big-endian files are readable only by kent's tools and this crate.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ByteOrder {
    /// Least significant byte first, as written on every common host.
    #[default]
    Little,
    /// Most significant byte first.
    Big,
}

impl ByteOrder {
    /// Detects the byte order of a structure from its first four bytes, given as `found` (those
    /// bytes read as a little-endian u32): the order in which they decode as `magic`, or `None` if
    /// neither does.
    pub(crate) fn detect(found: u32, magic: u32) -> Option<Self> {
        if found == magic {
            Some(Self::Little)
        } else if found.swap_bytes() == magic {
            Some(Self::Big)
        } else {
            None
        }
    }

    /// Appends `value` to `buf` in this byte order.
    pub(crate) fn put_u16(self, buf: &mut Vec<u8>, value: u16) {
        buf.extend_from_slice(&self.u16_bytes(value));
    }

    /// Appends `value` to `buf` in this byte order.
    pub(crate) fn put_u32(self, buf: &mut Vec<u8>, value: u32) {
        buf.extend_from_slice(&self.u32_bytes(value));
    }

    /// Appends `value` to `buf` in this byte order.
    pub(crate) fn put_u64(self, buf: &mut Vec<u8>, value: u64) {
        buf.extend_from_slice(&self.u64_bytes(value));
    }

    /// Appends `value` to `buf` in this byte order.
    pub(crate) fn put_f32(self, buf: &mut Vec<u8>, value: f32) {
        self.put_u32(buf, value.to_bits());
    }

    /// Appends `value` to `buf` in this byte order.
    pub(crate) fn put_f64(self, buf: &mut Vec<u8>, value: f64) {
        self.put_u64(buf, value.to_bits());
    }

    /// Encodes `value` in this byte order, for patching a field in place.
    pub(crate) fn u16_bytes(self, value: u16) -> [u8; 2] {
        match self {
            Self::Little => value.to_le_bytes(),
            Self::Big => value.to_be_bytes(),
        }
    }

    /// Encodes `value` in this byte order, for patching a field in place.
    pub(crate) fn u32_bytes(self, value: u32) -> [u8; 4] {
        match self {
            Self::Little => value.to_le_bytes(),
            Self::Big => value.to_be_bytes(),
        }
    }

    /// Encodes `value` in this byte order, for patching a field in place.
    pub(crate) fn u64_bytes(self, value: u64) -> [u8; 8] {
        match self {
            Self::Little => value.to_le_bytes(),
            Self::Big => value.to_be_bytes(),
        }
    }

    /// Decodes a u16 stored in this byte order.
    pub(crate) fn u16_from(self, bytes: [u8; 2]) -> u16 {
        match self {
            Self::Little => u16::from_le_bytes(bytes),
            Self::Big => u16::from_be_bytes(bytes),
        }
    }

    /// Decodes a u32 stored in this byte order.
    pub(crate) fn u32_from(self, bytes: [u8; 4]) -> u32 {
        match self {
            Self::Little => u32::from_le_bytes(bytes),
            Self::Big => u32::from_be_bytes(bytes),
        }
    }

    /// Decodes a u64 stored in this byte order.
    pub(crate) fn u64_from(self, bytes: [u8; 8]) -> u64 {
        match self {
            Self::Little => u64::from_le_bytes(bytes),
            Self::Big => u64::from_be_bytes(bytes),
        }
    }
}

/// One chromosome (sequence): its name and length in bases.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ChromInfo {
    /// Sequence name, e.g. `chr1`.
    pub name: String,
    /// Length in bases.
    pub length: u32,
}

impl ChromInfo {
    /// Creates a chromosome entry. Names and lengths are validated by the writer, not here.
    pub fn new(name: impl Into<String>, length: u32) -> Self {
        Self { name: name.into(), length }
    }

    /// Parses a UCSC `chrom.sizes` file: one `name length` pair per line, separated by tabs or
    /// spaces, returned in file order. Blank lines and lines whose first non-blank character is
    /// `#` are skipped. Lengths must be plain unsigned integers that fit in a `u32`.
    ///
    /// Names, duplicates and zero lengths are not checked here; the writer validates the list
    /// when it is created. Fails with [`Error::Parse`] carrying the 1-based line number, or with
    /// [`Error::Io`].
    ///
    /// ```
    /// use bigwig::ChromInfo;
    ///
    /// let chroms = ChromInfo::read_chrom_sizes("chr1\t248956422\nchrM 16569\n".as_bytes())?;
    /// assert_eq!(chroms, [ChromInfo::new("chr1", 248_956_422), ChromInfo::new("chrM", 16_569)]);
    /// # Ok::<(), bigwig::Error>(())
    /// ```
    pub fn read_chrom_sizes(mut input: impl BufRead) -> Result<Vec<ChromInfo>> {
        let mut chroms = Vec::new();
        let mut line = String::new();
        let mut line_number = 0;
        loop {
            line.clear();
            line_number += 1;
            match input.read_line(&mut line) {
                Ok(0) => return Ok(chroms),
                Ok(_) => {}
                Err(error) if error.kind() == io::ErrorKind::InvalidData => {
                    return Err(parse_error(line_number, "line is not valid UTF-8"));
                }
                Err(error) => return Err(error.into()),
            }
            let trimmed = line.trim_matches(is_separator);
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let (count, fields) = split_fields(trimmed);
            if count != 2 {
                let found = if count > 2 { "more".to_string() } else { count.to_string() };
                return Err(parse_error(
                    line_number,
                    format!("expected 2 fields (name and length), found {found}"),
                ));
            }
            let length = parse_u32(fields[1], line_number, "length")?;
            chroms.push(Self::new(fields[0], length));
        }
    }
}

/// A bounds-checked cursor over bytes read from a file, decoding fields in the file's byte order.
///
/// It knows the absolute file offset of its first byte, so a truncated structure is reported as
/// [`Error::Corrupt`] at the offset where the missing field starts.
#[derive(Clone, Debug)]
pub(crate) struct ByteReader<'a> {
    /// The bytes being decoded.
    bytes: &'a [u8],
    /// Index of the next byte to decode.
    position: usize,
    /// Byte order of multi-byte fields.
    byte_order: ByteOrder,
    /// Absolute file offset of `bytes[0]`.
    base_offset: u64,
}

impl<'a> ByteReader<'a> {
    /// Creates a reader over `bytes`, which were read from absolute file offset `base_offset`.
    pub(crate) fn new(bytes: &'a [u8], byte_order: ByteOrder, base_offset: u64) -> Self {
        Self { bytes, position: 0, byte_order, base_offset }
    }

    /// Decodes the next byte.
    pub(crate) fn u8(&mut self) -> Result<u8> {
        Ok(self.array::<1>()?[0])
    }

    /// Decodes the next u16.
    pub(crate) fn u16(&mut self) -> Result<u16> {
        let bytes = self.array()?;
        Ok(self.byte_order.u16_from(bytes))
    }

    /// Decodes the next u32.
    pub(crate) fn u32(&mut self) -> Result<u32> {
        let bytes = self.array()?;
        Ok(self.byte_order.u32_from(bytes))
    }

    /// Decodes the next u64.
    pub(crate) fn u64(&mut self) -> Result<u64> {
        let bytes = self.array()?;
        Ok(self.byte_order.u64_from(bytes))
    }

    /// Decodes the next f32.
    pub(crate) fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_bits(self.u32()?))
    }

    /// Decodes the next f64.
    pub(crate) fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_bits(self.u64()?))
    }

    /// Returns the next `len` bytes undecoded.
    pub(crate) fn bytes(&mut self, len: usize) -> Result<&'a [u8]> {
        let start = self.position;
        let end = self.checked_end(len)?;
        self.position = end;
        Ok(&self.bytes[start..end])
    }

    /// Skips the next `len` bytes.
    pub(crate) fn skip(&mut self, len: usize) -> Result<()> {
        self.position = self.checked_end(len)?;
        Ok(())
    }

    /// Absolute file offset of the next byte.
    pub(crate) fn file_offset(&self) -> u64 {
        // A slice index always fits in u64 on supported (≤ 64-bit) targets.
        self.base_offset + self.position as u64
    }

    /// Number of bytes not yet decoded.
    pub(crate) fn remaining(&self) -> usize {
        self.bytes.len() - self.position
    }

    /// Returns the next `N` bytes as an array.
    fn array<const N: usize>(&mut self) -> Result<[u8; N]> {
        let bytes = self.bytes(N)?;
        let mut array = [0u8; N];
        array.copy_from_slice(bytes);
        Ok(array)
    }

    /// Returns the index just past the next `len` bytes, or `Corrupt` if fewer remain.
    fn checked_end(&self, len: usize) -> Result<usize> {
        if len > self.remaining() {
            return Err(Error::corrupt(
                self.file_offset(),
                format!("truncated: needed {len} bytes but only {} remain", self.remaining()),
            ));
        }
        Ok(self.position + len)
    }
}

/// Decodes the 40-byte total summary.
#[allow(clippy::cast_precision_loss, reason = "base counts are far below 2^53, so exact")]
pub(crate) fn decode_total_summary(reader: &mut ByteReader<'_>) -> Result<Summary> {
    let bases_covered = reader.u64()?;
    Ok(Summary {
        covered_bases: bases_covered as f64,
        min: reader.f64()?,
        max: reader.f64()?,
        sum: reader.f64()?,
        sum_squares: reader.f64()?,
    })
}

/// Appends the 40-byte total summary to `buf` in byte order `order`.
///
/// The file stores `basesCovered` as a u64, so `summary.covered_bases` must be a whole number of
/// bases, as it is when accumulated from raw data.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "covered_bases is a whole, non-negative count below 2^53, so the cast is exact"
)]
pub(crate) fn encode_total_summary(summary: &Summary, order: ByteOrder, buf: &mut Vec<u8>) {
    order.put_u64(buf, summary.covered_bases as u64);
    order.put_f64(buf, summary.min);
    order.put_f64(buf, summary.max);
    order.put_f64(buf, summary.sum);
    order.put_f64(buf, summary.sum_squares);
}

/// Reads exactly `len` bytes at absolute `offset` from `source` into a new buffer.
///
/// `file_len` is the length of the source. A range extending past it is rejected as
/// [`Error::Corrupt`] before anything is allocated, so a corrupt length can't cause a huge
/// allocation.
pub(crate) fn read_exact_at<R: Read + Seek>(
    source: &mut R,
    offset: u64,
    len: usize,
    file_len: u64,
) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    read_exact_at_into(source, offset, len, file_len, &mut buf)?;
    Ok(buf)
}

/// As [`read_exact_at`], reusing `buf`, which is resized to `len`.
pub(crate) fn read_exact_at_into<R: Read + Seek>(
    source: &mut R,
    offset: u64,
    len: usize,
    file_len: u64,
    buf: &mut Vec<u8>,
) -> Result<()> {
    let end = u64::try_from(len).ok().and_then(|len| offset.checked_add(len));
    if end.is_none_or(|end| end > file_len) {
        return Err(Error::corrupt(
            offset,
            format!("{len}-byte read extends past the end of the {file_len}-byte file"),
        ));
    }
    buf.clear();
    buf.resize(len, 0);
    source.seek(SeekFrom::Start(offset))?;
    source.read_exact(buf)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn little_endian_puts_least_significant_byte_first() {
        let mut buf = Vec::new();
        ByteOrder::Little.put_u32(&mut buf, 0x0102_0304);
        assert_eq!(buf, [0x04, 0x03, 0x02, 0x01]);
    }

    #[test]
    fn big_endian_puts_most_significant_byte_first() {
        let mut buf = Vec::new();
        ByteOrder::Big.put_u32(&mut buf, 0x0102_0304);
        assert_eq!(buf, [0x01, 0x02, 0x03, 0x04]);
    }

    #[test]
    fn default_byte_order_is_little() {
        assert_eq!(ByteOrder::default(), ByteOrder::Little);
    }

    #[test]
    fn every_field_width_round_trips_in_both_byte_orders() {
        for order in [ByteOrder::Little, ByteOrder::Big] {
            let mut buf = vec![0xAB];
            order.put_u16(&mut buf, 0xBEEF);
            order.put_u32(&mut buf, 0xDEAD_BEEF);
            order.put_u64(&mut buf, 0x0123_4567_89AB_CDEF);
            order.put_f32(&mut buf, -1.5);
            order.put_f64(&mut buf, 1e300);

            let mut reader = ByteReader::new(&buf, order, 0);
            assert_eq!(reader.u8().unwrap(), 0xAB);
            assert_eq!(reader.u16().unwrap(), 0xBEEF);
            assert_eq!(reader.u32().unwrap(), 0xDEAD_BEEF);
            assert_eq!(reader.u64().unwrap(), 0x0123_4567_89AB_CDEF);
            assert_eq!(reader.f32().unwrap(), -1.5);
            assert_eq!(reader.f64().unwrap(), 1e300);
            assert_eq!(reader.remaining(), 0);
        }
    }

    #[test]
    fn byte_arrays_match_appended_bytes() {
        for order in [ByteOrder::Little, ByteOrder::Big] {
            let mut buf = Vec::new();
            order.put_u16(&mut buf, 7);
            order.put_u32(&mut buf, 8);
            order.put_u64(&mut buf, 9);
            let expected: Vec<u8> =
                [&order.u16_bytes(7)[..], &order.u32_bytes(8)[..], &order.u64_bytes(9)[..]]
                    .concat();
            assert_eq!(buf, expected);
        }
    }

    #[test]
    fn same_bytes_decode_differently_by_byte_order() {
        let bytes = [0x00, 0x00, 0x00, 0x01];
        assert_eq!(ByteReader::new(&bytes, ByteOrder::Big, 0).u32().unwrap(), 1);
        assert_eq!(ByteReader::new(&bytes, ByteOrder::Little, 0).u32().unwrap(), 0x0100_0000);
    }

    #[test]
    fn truncated_field_is_corrupt_at_its_absolute_offset() {
        let bytes = [1, 2, 3, 4, 5, 6];
        let mut reader = ByteReader::new(&bytes, ByteOrder::Little, 1000);
        reader.u32().unwrap();
        match reader.u32() {
            Err(Error::Corrupt { offset, .. }) => assert_eq!(offset, 1004),
            other => panic!("expected Corrupt, got {other:?}"),
        }
    }

    #[test]
    fn failed_read_does_not_advance_the_reader() {
        let bytes = [1, 2];
        let mut reader = ByteReader::new(&bytes, ByteOrder::Little, 0);
        assert!(reader.u64().is_err());
        assert_eq!(reader.file_offset(), 0);
        assert_eq!(reader.u16().unwrap(), 0x0201);
    }

    #[test]
    fn bytes_returns_the_next_slice_and_advances() {
        let bytes = [1, 2, 3, 4, 5];
        let mut reader = ByteReader::new(&bytes, ByteOrder::Little, 10);
        assert_eq!(reader.bytes(2).unwrap(), &[1, 2]);
        assert_eq!(reader.file_offset(), 12);
        assert_eq!(reader.bytes(3).unwrap(), &[3, 4, 5]);
        assert!(reader.bytes(1).is_err());
    }

    #[test]
    fn skip_advances_and_rejects_skipping_past_the_end() {
        let bytes = [1, 2, 3, 4];
        let mut reader = ByteReader::new(&bytes, ByteOrder::Big, 0);
        reader.skip(3).unwrap();
        assert_eq!(reader.u8().unwrap(), 4);
        assert!(matches!(reader.skip(1), Err(Error::Corrupt { offset: 4, .. })));
    }

    #[test]
    fn read_exact_at_reads_the_requested_range() {
        let mut source = Cursor::new((0u8..100).collect::<Vec<_>>());
        assert_eq!(read_exact_at(&mut source, 10, 4, 100).unwrap(), [10, 11, 12, 13]);
    }

    #[test]
    fn read_exact_at_accepts_a_range_ending_exactly_at_eof() {
        let mut source = Cursor::new(vec![7u8; 8]);
        assert_eq!(read_exact_at(&mut source, 4, 4, 8).unwrap(), [7, 7, 7, 7]);
    }

    #[test]
    fn read_exact_at_rejects_a_range_past_eof_as_corrupt() {
        let mut source = Cursor::new(vec![0u8; 8]);
        assert!(matches!(
            read_exact_at(&mut source, 6, 4, 8),
            Err(Error::Corrupt { offset: 6, .. })
        ));
    }

    #[test]
    fn read_exact_at_rejects_an_overflowing_range_without_allocating() {
        let mut source = Cursor::new(vec![0u8; 8]);
        let result = read_exact_at(&mut source, u64::MAX - 1, usize::MAX, 8);
        assert!(matches!(result, Err(Error::Corrupt { .. })));
    }

    #[test]
    fn read_exact_at_into_reuses_and_resizes_the_buffer() {
        let mut source = Cursor::new((0u8..10).collect::<Vec<_>>());
        let mut buf = vec![99; 32];
        read_exact_at_into(&mut source, 2, 3, 10, &mut buf).unwrap();
        assert_eq!(buf, [2, 3, 4]);
    }

    #[test]
    fn magic_numbers_match_the_format() {
        assert_eq!(BIGWIG_MAGIC, 0x888F_FC26);
        assert_eq!(BIGBED_MAGIC, 0x8789_F2EB);
        assert_eq!(BPTREE_MAGIC, 0x78CA_8C91);
        assert_eq!(RTREE_MAGIC, 0x2468_ACE0);
    }

    /// A header with every field distinct, so a misplaced field can't go unnoticed.
    fn distinct_header() -> FixedHeader {
        FixedHeader {
            version: BIGWIG_VERSION,
            zoom_levels: 3,
            chrom_tree_offset: 0x0101_0101_0101_0101,
            full_data_offset: 0x0202_0202_0202_0202,
            full_index_offset: 0x0303_0303_0303_0303,
            field_count: 0x0404,
            defined_field_count: 0x0505,
            auto_sql_offset: 0x0606_0606_0606_0606,
            total_summary_offset: 0x0707_0707_0707_0707,
            uncompress_buf_size: 0x0808_0808,
            extension_offset: 0x0909_0909_0909_0909,
        }
    }

    /// Encodes `header` in byte order `order`.
    fn encoded_header(header: &FixedHeader, order: ByteOrder) -> Vec<u8> {
        let mut buf = Vec::new();
        header.encode(order, &mut buf);
        buf
    }

    #[test]
    fn detect_finds_little_endian_magic() {
        let found = u32::from_le_bytes(BIGWIG_MAGIC.to_le_bytes());
        assert_eq!(ByteOrder::detect(found, BIGWIG_MAGIC), Some(ByteOrder::Little));
    }

    #[test]
    fn detect_finds_big_endian_magic() {
        let found = u32::from_le_bytes(BIGWIG_MAGIC.to_be_bytes());
        assert_eq!(ByteOrder::detect(found, BIGWIG_MAGIC), Some(ByteOrder::Big));
    }

    #[test]
    fn detect_rejects_other_magic() {
        assert_eq!(ByteOrder::detect(BIGBED_MAGIC, BIGWIG_MAGIC), None);
    }

    #[test]
    fn fixed_header_round_trips_in_both_byte_orders() {
        for order in [ByteOrder::Little, ByteOrder::Big] {
            let bytes = encoded_header(&distinct_header(), order);
            assert_eq!(FixedHeader::decode(&bytes).unwrap(), (distinct_header(), order));
        }
    }

    #[test]
    fn fixed_header_encodes_to_its_fixed_size() {
        assert_eq!(encoded_header(&distinct_header(), ByteOrder::Little).len(), FIXED_HEADER_SIZE);
    }

    #[test]
    fn fixed_header_starts_with_the_magic_in_the_chosen_byte_order() {
        let little = encoded_header(&FixedHeader::default(), ByteOrder::Little);
        let big = encoded_header(&FixedHeader::default(), ByteOrder::Big);
        assert_eq!(little[..4], BIGWIG_MAGIC.to_le_bytes());
        assert_eq!(big[..4], BIGWIG_MAGIC.to_be_bytes());
    }

    #[test]
    fn fixed_header_fields_sit_at_their_format_offsets() {
        let bytes = encoded_header(&distinct_header(), ByteOrder::Little);
        let u16_at = |at: usize| u16::from_le_bytes(bytes[at..at + 2].try_into().unwrap());
        let u32_at = |at: usize| u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
        let u64_at = |at: usize| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
        assert_eq!(u16_at(4), BIGWIG_VERSION);
        assert_eq!(u16_at(6), 3);
        assert_eq!(u64_at(8), 0x0101_0101_0101_0101);
        assert_eq!(u64_at(16), 0x0202_0202_0202_0202);
        assert_eq!(u64_at(24), 0x0303_0303_0303_0303);
        assert_eq!(u16_at(32), 0x0404);
        assert_eq!(u16_at(34), 0x0505);
        assert_eq!(u64_at(36), 0x0606_0606_0606_0606);
        assert_eq!(u64_at(44), 0x0707_0707_0707_0707);
        assert_eq!(u32_at(52), 0x0808_0808);
        assert_eq!(u64_at(56), 0x0909_0909_0909_0909);
    }

    #[test]
    fn fixed_header_accepts_versions_one_to_four() {
        for version in 1..=4 {
            let header = FixedHeader { version, ..FixedHeader::default() };
            let bytes = encoded_header(&header, ByteOrder::Little);
            assert_eq!(FixedHeader::decode(&bytes).unwrap().0.version, version);
        }
    }

    #[test]
    fn fixed_header_rejects_version_zero() {
        let bytes = encoded_header(&FixedHeader::default(), ByteOrder::Little);
        assert!(matches!(FixedHeader::decode(&bytes), Err(Error::UnsupportedVersion(0))));
    }

    #[test]
    fn fixed_header_rejects_version_five() {
        let header = FixedHeader { version: 5, ..FixedHeader::default() };
        let bytes = encoded_header(&header, ByteOrder::Big);
        assert!(matches!(FixedHeader::decode(&bytes), Err(Error::UnsupportedVersion(5))));
    }

    #[test]
    fn fixed_header_rejects_unknown_magic_reporting_it_little_endian() {
        let mut bytes = encoded_header(&distinct_header(), ByteOrder::Little);
        bytes[..4].copy_from_slice(&[0x01, 0x02, 0x03, 0x04]);
        assert!(matches!(
            FixedHeader::decode(&bytes),
            Err(Error::NotBigWig { found: 0x0403_0201 })
        ));
    }

    #[test]
    fn fixed_header_rejects_bigbed_magic_in_either_byte_order() {
        for magic in [BIGBED_MAGIC.to_le_bytes(), BIGBED_MAGIC.to_be_bytes()] {
            let mut bytes = encoded_header(&distinct_header(), ByteOrder::Little);
            bytes[..4].copy_from_slice(&magic);
            assert!(matches!(FixedHeader::decode(&bytes), Err(Error::IsBigBed)));
        }
    }

    #[test]
    fn fixed_header_truncated_after_the_magic_is_corrupt_at_the_missing_field() {
        let bytes = encoded_header(&distinct_header(), ByteOrder::Little);
        // autoSqlOffset spans bytes 36..44, so 40 bytes cut it short.
        assert!(matches!(
            FixedHeader::decode(&bytes[..40]),
            Err(Error::Corrupt { offset: 36, .. })
        ));
    }

    #[test]
    fn fixed_header_shorter_than_the_magic_is_not_bigwig_with_the_missing_bytes_as_zero() {
        assert!(matches!(
            FixedHeader::decode(&[0x26, 0xFC]),
            Err(Error::NotBigWig { found: 0x0000_FC26 })
        ));
    }

    #[test]
    fn zoom_level_round_trips_in_both_byte_orders() {
        let level = ZoomLevel::new(256, 0x1234_5678_9ABC, 0x0FED_CBA9_8765);
        for order in [ByteOrder::Little, ByteOrder::Big] {
            let mut buf = Vec::new();
            level.encode(order, &mut buf);
            let mut reader = ByteReader::new(&buf, order, 64);
            assert_eq!(ZoomLevel::decode(&mut reader).unwrap(), level);
            assert_eq!(reader.remaining(), 0);
        }
    }

    #[test]
    fn zoom_level_encodes_to_its_fixed_size_with_zero_reserved() {
        let mut buf = Vec::new();
        ZoomLevel::new(u32::MAX, u64::MAX, u64::MAX).encode(ByteOrder::Little, &mut buf);
        assert_eq!(buf.len(), ZOOM_HEADER_SIZE);
        assert_eq!(buf[4..8], [0, 0, 0, 0]);
    }

    #[test]
    fn zoom_level_decode_ignores_the_reserved_field() {
        let mut buf = Vec::new();
        ZoomLevel::new(40, 1000, 2000).encode(ByteOrder::Big, &mut buf);
        buf[4..8].copy_from_slice(&[0xFF; 4]);
        let level = ZoomLevel::decode(&mut ByteReader::new(&buf, ByteOrder::Big, 0)).unwrap();
        assert_eq!(level, ZoomLevel::new(40, 1000, 2000));
    }

    #[test]
    fn zoom_level_reports_its_reduction_and_offsets() {
        let level = ZoomLevel::new(1024, 5000, 9000);
        assert_eq!(level.reduction(), 1024);
        assert_eq!(level.data_offset(), 5000);
        assert_eq!(level.index_offset(), 9000);
    }

    #[test]
    fn truncated_zoom_header_is_corrupt_at_the_missing_field() {
        let mut buf = Vec::new();
        ZoomLevel::new(40, 1000, 2000).encode(ByteOrder::Little, &mut buf);
        let mut reader = ByteReader::new(&buf[..20], ByteOrder::Little, 64);
        assert!(matches!(ZoomLevel::decode(&mut reader), Err(Error::Corrupt { offset: 80, .. })));
    }

    #[test]
    fn total_summary_round_trips_in_both_byte_orders() {
        let summary =
            Summary { covered_bases: 1234.0, min: -2.5, max: 99.25, sum: 1e6, sum_squares: 3e9 };
        for order in [ByteOrder::Little, ByteOrder::Big] {
            let mut buf = Vec::new();
            encode_total_summary(&summary, order, &mut buf);
            assert_eq!(buf.len(), TOTAL_SUMMARY_SIZE);
            let mut reader = ByteReader::new(&buf, order, 0);
            assert_eq!(decode_total_summary(&mut reader).unwrap(), summary);
        }
    }

    #[test]
    fn total_summary_stores_bases_covered_as_an_integer() {
        let summary =
            Summary { covered_bases: 3e9, min: 0.0, max: 0.0, sum: 0.0, sum_squares: 0.0 };
        let mut buf = Vec::new();
        encode_total_summary(&summary, ByteOrder::Little, &mut buf);
        assert_eq!(buf[..8], 3_000_000_000u64.to_le_bytes());
    }

    #[test]
    fn truncated_total_summary_is_corrupt() {
        let mut reader = ByteReader::new(&[0; 39], ByteOrder::Little, 100);
        assert!(matches!(
            decode_total_summary(&mut reader),
            Err(Error::Corrupt { offset: 132, .. })
        ));
    }

    #[test]
    fn chrom_info_new_stores_name_and_length() {
        assert_eq!(
            ChromInfo::new("chr1", 248_956_422),
            ChromInfo { name: "chr1".to_string(), length: 248_956_422 }
        );
    }

    /// The line and message of the parse error `read_chrom_sizes` returns for `text`.
    fn chrom_sizes_failure(text: &[u8]) -> (u64, String) {
        match ChromInfo::read_chrom_sizes(text) {
            Err(Error::Parse { line, message }) => (line, message),
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    #[test]
    fn read_chrom_sizes_returns_name_and_length_pairs_in_file_order() {
        let chroms = ChromInfo::read_chrom_sizes(&b"chr2\t500\nchr1\t1000\n"[..]).unwrap();
        assert_eq!(chroms, [ChromInfo::new("chr2", 500), ChromInfo::new("chr1", 1000)]);
    }

    #[test]
    fn read_chrom_sizes_accepts_runs_of_spaces_and_tabs_between_fields() {
        let chroms = ChromInfo::read_chrom_sizes(&b"chr1   10\nchr2\t \t20\n"[..]).unwrap();
        assert_eq!(chroms, [ChromInfo::new("chr1", 10), ChromInfo::new("chr2", 20)]);
    }

    #[test]
    fn read_chrom_sizes_skips_comments_and_blank_lines() {
        let text = b"# assembly sizes\n\n  # indented comment\n   \nchr1\t10\n";
        let chroms = ChromInfo::read_chrom_sizes(&text[..]).unwrap();
        assert_eq!(chroms, [ChromInfo::new("chr1", 10)]);
    }

    #[test]
    fn read_chrom_sizes_accepts_crlf_line_endings_and_a_missing_final_newline() {
        let chroms = ChromInfo::read_chrom_sizes(&b"chr1\t10\r\nchr2\t20"[..]).unwrap();
        assert_eq!(chroms, [ChromInfo::new("chr1", 10), ChromInfo::new("chr2", 20)]);
    }

    #[test]
    fn read_chrom_sizes_of_empty_input_is_empty() {
        assert_eq!(ChromInfo::read_chrom_sizes(&b""[..]).unwrap(), []);
    }

    #[test]
    fn read_chrom_sizes_leaves_duplicates_and_zero_lengths_to_the_writer() {
        let chroms = ChromInfo::read_chrom_sizes(&b"chr1\t10\nchr1\t0\n"[..]).unwrap();
        assert_eq!(chroms, [ChromInfo::new("chr1", 10), ChromInfo::new("chr1", 0)]);
    }

    #[test]
    fn read_chrom_sizes_rejects_a_line_with_one_field_at_its_line_number() {
        let (line, message) = chrom_sizes_failure(b"# sizes\nchr1\t10\nchr2\n");
        assert_eq!(line, 3);
        assert!(message.contains("expected 2 fields"), "{message}");
    }

    #[test]
    fn read_chrom_sizes_rejects_a_line_with_three_fields() {
        let (line, message) = chrom_sizes_failure(b"chr1\t10\tchr1.fa\n");
        assert_eq!(line, 1);
        assert!(message.contains("found more"), "{message}");
    }

    #[test]
    fn read_chrom_sizes_rejects_a_fractional_length() {
        assert_eq!(chrom_sizes_failure(b"chr1\t10.5\n").0, 1);
    }

    #[test]
    fn read_chrom_sizes_rejects_a_negative_length() {
        assert_eq!(chrom_sizes_failure(b"chr1\t-5\n").0, 1);
    }

    #[test]
    fn read_chrom_sizes_rejects_a_length_with_a_plus_sign() {
        assert_eq!(chrom_sizes_failure(b"chr1\t+5\n").0, 1);
    }

    #[test]
    fn read_chrom_sizes_rejects_a_length_in_scientific_notation() {
        assert_eq!(chrom_sizes_failure(b"chr1\t1e3\n").0, 1);
    }

    #[test]
    fn read_chrom_sizes_rejects_a_length_past_u32_max() {
        let (line, message) = chrom_sizes_failure(b"chr1\t4294967295\nchr2\t4294967296\n");
        assert_eq!(line, 2);
        assert!(message.contains("out of range"), "{message}");
    }

    #[test]
    fn read_chrom_sizes_reports_invalid_utf8_at_its_line_number() {
        assert_eq!(chrom_sizes_failure(b"chr1\t10\nchr\xff\t20\n").0, 2);
    }
}
