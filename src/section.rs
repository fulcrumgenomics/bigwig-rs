//! bigWig data sections: the [`DataItem`] type, and encoding and decoding of the three section
//! types (bedGraph, variableStep, fixedStep).
//!
//! A decompressed section is a 24-byte header followed by `itemCount` items whose layout depends
//! on the type. [`SectionBuilder`] produces that form for the writer, one section at a time;
//! [`Section`] decodes it for the reader.

use crate::bbi::{ByteOrder, ByteReader};
use crate::error::{Error, Result};

/// Length of the header that starts every section.
pub(crate) const SECTION_HEADER_LEN: usize = 24;
/// Most items one section can hold: `itemCount` is a u16.
pub(crate) const MAX_ITEMS_PER_SECTION: u16 = u16::MAX;

/// One stored data value: `value` applies to every base in `[start, end)` on chromosome
/// `chrom_id`.
///
/// An item is raw data, not a summary. Reading expands fixedStep and variableStep items into
/// explicit ranges, so items from all three section types have this one form.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DataItem {
    /// Index into the reader's chromosome list (the file's chromId).
    pub chrom_id: u32,
    /// First base covered (0-based).
    pub start: u32,
    /// One past the last base covered.
    pub end: u32,
    /// The value of every base in `[start, end)`.
    pub value: f32,
}

/// How a section stores its items: the header's `type`, with the step and span that go with it.
///
/// A section never mixes types, steps or spans, so comparing kinds with `==` tells the writer
/// whether the next item may join the open section.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SectionKind {
    /// Type 1: each item stores `start`, `end` and `value` (12 bytes).
    BedGraph,
    /// Type 2: each item stores `start` and `value` (8 bytes) and covers `span` bases.
    VariableStep {
        /// Bases covered by every item.
        span: u32,
    },
    /// Type 3: each item stores only `value` (4 bytes); item *i* starts `i × step` bases after
    /// the section start and covers `span` bases.
    FixedStep {
        /// Distance between consecutive item starts.
        step: u32,
        /// Bases covered by every item.
        span: u32,
    },
}

impl SectionKind {
    /// Rebuilds the kind from a decoded header, or `None` for an unknown type. A bedGraph
    /// header's `itemStep` and `itemSpan` are ignored because its items are explicit.
    fn from_header(type_code: u8, step: u32, span: u32) -> Option<Self> {
        match type_code {
            1 => Some(Self::BedGraph),
            2 => Some(Self::VariableStep { span }),
            3 => Some(Self::FixedStep { step, span }),
            _ => None,
        }
    }

    /// The header's `type` byte.
    fn type_code(self) -> u8 {
        match self {
            Self::BedGraph => 1,
            Self::VariableStep { .. } => 2,
            Self::FixedStep { .. } => 3,
        }
    }

    /// The header's `itemStep` field, which only fixedStep uses.
    fn item_step(self) -> u32 {
        match self {
            Self::BedGraph | Self::VariableStep { .. } => 0,
            Self::FixedStep { step, .. } => step,
        }
    }

    /// The header's `itemSpan` field, which bedGraph leaves at 0.
    fn item_span(self) -> u32 {
        match self {
            Self::BedGraph => 0,
            Self::VariableStep { span } | Self::FixedStep { span, .. } => span,
        }
    }

    /// Encoded size of one item.
    fn item_len(self) -> usize {
        match self {
            Self::BedGraph => 12,
            Self::VariableStep { .. } => 8,
            Self::FixedStep { .. } => 4,
        }
    }
}

/// Accumulates one section's items, already in file byte order, for the writer to flush as one
/// block.
///
/// Items are appended as they arrive and the header is patched in by [`Self::finish`]. The
/// builder validates nothing about ordering, overlap or bounds; the writer rejects bad items
/// before pushing them. One builder is reused for every section of a file.
#[derive(Debug)]
pub(crate) struct SectionBuilder {
    byte_order: ByteOrder,
    /// The header placeholder followed by the encoded items.
    buf: Vec<u8>,
    chrom_id: u32,
    kind: SectionKind,
    /// Start of the first item.
    start: u32,
    /// End of the last item by kent's formulas; equal to `start` while empty.
    end: u32,
    /// Start of the next fixedStep item.
    next_start: u32,
    item_count: u16,
}

impl SectionBuilder {
    /// Creates an empty builder. Call one of the `begin_*` methods before pushing items.
    pub(crate) fn new(byte_order: ByteOrder) -> Self {
        let mut builder = Self {
            byte_order,
            buf: Vec::new(),
            chrom_id: 0,
            kind: SectionKind::BedGraph,
            start: 0,
            end: 0,
            next_start: 0,
            item_count: 0,
        };
        builder.reset(0, SectionKind::BedGraph, 0);
        builder
    }

    /// Starts a bedGraph section on `chrom_id`, discarding anything not yet finished.
    pub(crate) fn begin_bed_graph(&mut self, chrom_id: u32) {
        self.reset(chrom_id, SectionKind::BedGraph, 0);
    }

    /// Starts a variableStep section on `chrom_id` whose items all cover `span` bases.
    pub(crate) fn begin_variable_step(&mut self, chrom_id: u32, span: u32) {
        self.reset(chrom_id, SectionKind::VariableStep { span }, 0);
    }

    /// Starts a fixedStep section on `chrom_id` whose first item starts at `start`; item *i*
    /// covers `[start + i × step, start + i × step + span)`.
    pub(crate) fn begin_fixed_step(&mut self, chrom_id: u32, start: u32, step: u32, span: u32) {
        self.reset(chrom_id, SectionKind::FixedStep { step, span }, start);
    }

    /// Appends a bedGraph item covering `[start, end)`.
    pub(crate) fn push_bed_graph(&mut self, start: u32, end: u32, value: f32) {
        assert_eq!(self.kind, SectionKind::BedGraph, "push_bed_graph on a {:?} section", self.kind);
        self.count_item(start);
        self.byte_order.put_u32(&mut self.buf, start);
        self.byte_order.put_u32(&mut self.buf, end);
        self.byte_order.put_f32(&mut self.buf, value);
        self.end = end;
    }

    /// Appends a variableStep item starting at `start`.
    pub(crate) fn push_variable_step(&mut self, start: u32, value: f32) {
        let SectionKind::VariableStep { span } = self.kind else {
            panic!("push_variable_step on a {:?} section", self.kind);
        };
        self.count_item(start);
        self.byte_order.put_u32(&mut self.buf, start);
        self.byte_order.put_f32(&mut self.buf, value);
        // The writer has already rejected items that end past u32::MAX.
        self.end = start + span;
    }

    /// Appends the next fixedStep item. The writer pushes in bulk; tests use this one-at-a-time
    /// form as a reference for [`Self::push_fixed_step_values`].
    #[cfg(test)]
    pub(crate) fn push_fixed_step(&mut self, value: f32) {
        let SectionKind::FixedStep { step, span } = self.kind else {
            panic!("push_fixed_step on a {:?} section", self.kind);
        };
        self.count_item(self.next_start);
        self.byte_order.put_f32(&mut self.buf, value);
        // The writer has already rejected items that end past u32::MAX. The start after the last
        // item may still lie past it when span < step, and is then never used.
        self.end = self.next_start + span;
        self.next_start = self.next_start.wrapping_add(step);
    }

    /// Appends the next fixedStep items, one per value: each starts `step` bases after the one
    /// before. Dense per-base data comes through here, so the values are encoded in one pass.
    pub(crate) fn push_fixed_step_values(&mut self, values: &[f32]) {
        let SectionKind::FixedStep { step, span } = self.kind else {
            panic!("push_fixed_step_values on a {:?} section", self.kind);
        };
        let Some(last_index) = values.len().checked_sub(1) else { return };
        let item_count = u16::try_from(usize::from(self.item_count) + values.len())
            .expect("a section holds at most 65,535 items");
        if self.item_count == 0 {
            self.start = self.next_start;
        }
        self.item_count = item_count;

        let first_byte = self.buf.len();
        self.buf.resize(first_byte + values.len() * 4, 0);
        let encoded = self.buf[first_byte..].chunks_exact_mut(4).zip(values);
        // One loop per byte order, so neither re-checks the order for every value.
        match self.byte_order {
            ByteOrder::Little => {
                encoded.for_each(|(bytes, value)| bytes.copy_from_slice(&value.to_le_bytes()));
            }
            ByteOrder::Big => {
                encoded.for_each(|(bytes, value)| bytes.copy_from_slice(&value.to_be_bytes()));
            }
        }

        // The writer has already rejected items that end past u32::MAX, so the last item's start
        // fits; the start after it may not when span < step, and is then never used.
        let last_index = u32::try_from(last_index).expect("a section holds at most 65,535 items");
        let last_start = self.next_start + last_index * step;
        self.end = last_start + span;
        self.next_start = last_start.wrapping_add(step);
    }

    /// Fills in the header and returns the whole encoded section, header first. The builder keeps
    /// the section (and its accessors keep answering) until the next `begin_*`.
    pub(crate) fn finish(&mut self) -> &[u8] {
        let order = self.byte_order;
        let header = &mut self.buf[..SECTION_HEADER_LEN];
        header[0..4].copy_from_slice(&order.u32_bytes(self.chrom_id));
        header[4..8].copy_from_slice(&order.u32_bytes(self.start));
        header[8..12].copy_from_slice(&order.u32_bytes(self.end));
        header[12..16].copy_from_slice(&order.u32_bytes(self.kind.item_step()));
        header[16..20].copy_from_slice(&order.u32_bytes(self.kind.item_span()));
        header[20] = self.kind.type_code();
        header[21] = 0;
        header[22..24].copy_from_slice(&order.u16_bytes(self.item_count));
        &self.buf
    }

    /// Items pushed since the last `begin_*`.
    pub(crate) fn item_count(&self) -> u16 {
        self.item_count
    }

    /// Chromosome of the section.
    pub(crate) fn chrom_id(&self) -> u32 {
        self.chrom_id
    }

    /// Start of the first item (the fixedStep start given to `begin_fixed_step`, otherwise 0,
    /// until an item is pushed).
    pub(crate) fn start(&self) -> u32 {
        self.start
    }

    /// End of the last item, as kent computes it: fixedStep `start + (n − 1) × step + span`,
    /// variableStep `last start + span`, bedGraph the last item's end. Equal to [`Self::start`]
    /// while the section is empty.
    pub(crate) fn end(&self) -> u32 {
        self.end
    }

    /// Kind of the section being built.
    pub(crate) fn kind(&self) -> SectionKind {
        self.kind
    }

    /// Counts one more item, taking `start` as the section start when it is the first.
    fn count_item(&mut self, start: u32) {
        assert!(
            self.item_count < MAX_ITEMS_PER_SECTION,
            "a section holds at most {MAX_ITEMS_PER_SECTION} items"
        );
        if self.item_count == 0 {
            self.start = start;
        }
        self.item_count += 1;
    }

    /// Empties the builder and reserves the header for a new section.
    fn reset(&mut self, chrom_id: u32, kind: SectionKind, start: u32) {
        self.buf.clear();
        self.buf.resize(SECTION_HEADER_LEN, 0);
        self.chrom_id = chrom_id;
        self.kind = kind;
        self.start = start;
        self.end = start;
        self.next_start = start;
        self.item_count = 0;
    }
}

/// A decoded view of one decompressed section: its header fields and an allocation-free
/// iterator over its items.
///
/// [`Self::parse`] checks what the format fixes: a known type, and `itemCount` items that fit in
/// the block (trailing bytes are tolerated). Item extents come from the header's step and span
/// and are checked as the items are visited. The header's `end` plays no part in them, because
/// libBigWig wrote it too large before pyBigWig #166 was fixed, so only tests read it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Section<'a> {
    chrom_id: u32,
    start: u32,
    /// The header's `end`, which tests compare with what the writer meant to store.
    #[cfg(test)]
    end_hint: u32,
    kind: SectionKind,
    item_count: u16,
    /// Exactly `item_count` encoded items.
    items: &'a [u8],
    byte_order: ByteOrder,
    /// File offset of the block, for error reports.
    offset: u64,
}

impl<'a> Section<'a> {
    /// Decodes the header of `block`, the decompressed bytes of a section read from file offset
    /// `offset`. Every error is [`Error::Corrupt`] at that offset.
    pub(crate) fn parse(block: &'a [u8], byte_order: ByteOrder, offset: u64) -> Result<Self> {
        if block.len() < SECTION_HEADER_LEN {
            return Err(Error::corrupt(
                offset,
                format!(
                    "section header needs {SECTION_HEADER_LEN} bytes but the block holds {}",
                    block.len()
                ),
            ));
        }
        let (header, body) = block.split_at(SECTION_HEADER_LEN);
        let mut fields = ByteReader::new(header, byte_order, offset);
        let chrom_id = fields.u32()?;
        let start = fields.u32()?;
        #[cfg(test)]
        let end_hint = fields.u32()?;
        #[cfg(not(test))]
        fields.skip(4)?;
        let step = fields.u32()?;
        let span = fields.u32()?;
        let type_code = fields.u8()?;
        fields.skip(1)?;
        let item_count = fields.u16()?;

        let kind = SectionKind::from_header(type_code, step, span)
            .ok_or_else(|| Error::corrupt(offset, format!("unknown section type {type_code}")))?;
        let items_len = usize::from(item_count) * kind.item_len();
        if items_len > body.len() {
            return Err(Error::corrupt(
                offset,
                format!(
                    "itemCount {item_count} needs {items_len} bytes but only {} follow the header",
                    body.len()
                ),
            ));
        }
        Ok(Self {
            chrom_id,
            start,
            #[cfg(test)]
            end_hint,
            kind,
            item_count,
            items: &body[..items_len],
            byte_order,
            offset,
        })
    }

    /// Iterates over the items in stored order, deriving each extent from the header's step
    /// and span.
    pub(crate) fn items(self) -> SectionItems<'a> {
        SectionItems { section: self, index: 0 }
    }

    /// Chromosome of every item.
    pub(crate) fn chrom_id(self) -> u32 {
        self.chrom_id
    }

    /// Start of the first item, as the header states it.
    #[cfg(test)]
    pub(crate) fn start(self) -> u32 {
        self.start
    }

    /// The header's `end`. It is not trusted for item extents, and some writers overstate it.
    #[cfg(test)]
    pub(crate) fn end_hint(self) -> u32 {
        self.end_hint
    }

    /// Kind of the section.
    pub(crate) fn kind(self) -> SectionKind {
        self.kind
    }

    /// Number of items in the section.
    #[cfg(test)]
    pub(crate) fn item_count(self) -> u16 {
        self.item_count
    }

    /// Decodes item `index`, which must be below `item_count`. An item whose extent overflows
    /// u32 or has `start > end` is [`Error::Corrupt`]. A zero-length item (`start == end`) is
    /// returned as stored: it is harmless, and the reader decides what to do with it.
    fn item(self, index: usize) -> Result<DataItem> {
        let item_len = self.kind.item_len();
        let bytes = &self.items[index * item_len..][..item_len];
        let mut fields = ByteReader::new(bytes, self.byte_order, self.offset);
        let extent = match self.kind {
            SectionKind::BedGraph => Some((fields.u32()?, fields.u32()?)),
            SectionKind::VariableStep { span } => {
                let start = fields.u32()?;
                start.checked_add(span).map(|end| (start, end))
            }
            SectionKind::FixedStep { step, span } => u32::try_from(index)
                .ok()
                .and_then(|i| i.checked_mul(step))
                .and_then(|from_section_start| self.start.checked_add(from_section_start))
                .and_then(|start| start.checked_add(span).map(|end| (start, end))),
        };
        let value = fields.f32()?;
        let Some((start, end)) = extent else {
            return Err(Error::corrupt(self.offset, format!("item {index} extends past u32::MAX")));
        };
        if start > end {
            return Err(Error::corrupt(
                self.offset,
                format!("item {index} has start {start} > end {end}"),
            ));
        }
        Ok(DataItem { chrom_id: self.chrom_id, start, end, value })
    }
}

/// Iterator over a section's items, decoding each on demand. After an error it yields nothing
/// more.
#[derive(Clone, Debug)]
pub(crate) struct SectionItems<'a> {
    section: Section<'a>,
    /// Index of the next item to decode.
    index: usize,
}

impl Iterator for SectionItems<'_> {
    type Item = Result<DataItem>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= usize::from(self.section.item_count) {
            return None;
        }
        let item = self.section.item(self.index);
        self.index =
            if item.is_ok() { self.index + 1 } else { usize::from(self.section.item_count) };
        Some(item)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = usize::from(self.section.item_count).saturating_sub(self.index);
        (remaining, Some(remaining))
    }

    /// Skips `n` items without decoding them; items are fixed-size, so this is O(1).
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.index = self.index.saturating_add(n);
        self.next()
    }
}

impl ExactSizeIterator for SectionItems<'_> {}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    const BOTH_ORDERS: [ByteOrder; 2] = [ByteOrder::Little, ByteOrder::Big];

    /// A distinct, exactly representable value per item index.
    fn value(i: usize) -> f32 {
        f32::from(u16::try_from(i % 65_536).unwrap()) * 0.25 - 100.0
    }

    /// Encodes `items` (start, end, value) as one bedGraph section on chromosome 3.
    fn encode_bed_graph(order: ByteOrder, items: &[(u32, u32, f32)]) -> Vec<u8> {
        let mut builder = SectionBuilder::new(order);
        builder.begin_bed_graph(3);
        for &(start, end, value) in items {
            builder.push_bed_graph(start, end, value);
        }
        builder.finish().to_vec()
    }

    /// Encodes `items` (start, value) as one variableStep section on chromosome 3.
    fn encode_variable_step(order: ByteOrder, span: u32, items: &[(u32, f32)]) -> Vec<u8> {
        let mut builder = SectionBuilder::new(order);
        builder.begin_variable_step(3, span);
        for &(start, value) in items {
            builder.push_variable_step(start, value);
        }
        builder.finish().to_vec()
    }

    /// Encodes `values` as one fixedStep section on chromosome 3.
    fn encode_fixed_step(
        order: ByteOrder,
        start: u32,
        step: u32,
        span: u32,
        values: &[f32],
    ) -> Vec<u8> {
        let mut builder = SectionBuilder::new(order);
        builder.begin_fixed_step(3, start, step, span);
        for &value in values {
            builder.push_fixed_step(value);
        }
        builder.finish().to_vec()
    }

    /// Parses and fully decodes a section read from offset 1000.
    fn decode(block: &[u8], order: ByteOrder) -> Vec<DataItem> {
        Section::parse(block, order, 1000).unwrap().items().collect::<Result<Vec<_>>>().unwrap()
    }

    fn item(start: u32, end: u32, value: f32) -> DataItem {
        DataItem { chrom_id: 3, start, end, value }
    }

    /// Reads the u32 header field at `offset`.
    fn header_u32(block: &[u8], order: ByteOrder, offset: usize) -> u32 {
        order.u32_from(block[offset..offset + 4].try_into().unwrap())
    }

    fn bed_graph_items(n: usize) -> Vec<(u32, u32, f32)> {
        (0..n)
            .map(|i| (u32::try_from(i).unwrap() * 10, u32::try_from(i).unwrap() * 10 + 7, value(i)))
            .collect()
    }

    fn variable_step_items(n: usize) -> Vec<(u32, f32)> {
        (0..n).map(|i| (u32::try_from(i).unwrap() * 50 + 5, value(i))).collect()
    }

    fn fixed_step_values(n: usize) -> Vec<f32> {
        (0..n).map(value).collect()
    }

    fn assert_bed_graph_round_trips(n: usize) {
        let items = bed_graph_items(n);
        for order in BOTH_ORDERS {
            let block = encode_bed_graph(order, &items);
            assert_eq!(block.len(), SECTION_HEADER_LEN + 12 * n);
            let decoded = decode(&block, order);
            let expected: Vec<_> = items.iter().map(|&(s, e, v)| item(s, e, v)).collect();
            assert_eq!(decoded, expected);
        }
    }

    fn assert_variable_step_round_trips(n: usize) {
        let items = variable_step_items(n);
        for order in BOTH_ORDERS {
            let block = encode_variable_step(order, 20, &items);
            assert_eq!(block.len(), SECTION_HEADER_LEN + 8 * n);
            let decoded = decode(&block, order);
            let expected: Vec<_> = items.iter().map(|&(s, v)| item(s, s + 20, v)).collect();
            assert_eq!(decoded, expected);
        }
    }

    fn assert_fixed_step_round_trips(n: usize) {
        let values = fixed_step_values(n);
        for order in BOTH_ORDERS {
            let block = encode_fixed_step(order, 1_000_000, 10, 10, &values);
            assert_eq!(block.len(), SECTION_HEADER_LEN + 4 * n);
            let decoded = decode(&block, order);
            let expected: Vec<_> = values
                .iter()
                .enumerate()
                .map(|(i, &v)| {
                    let start = 1_000_000 + u32::try_from(i).unwrap() * 10;
                    item(start, start + 10, v)
                })
                .collect();
            assert_eq!(decoded, expected);
        }
    }

    // ---- round trips ------------------------------------------------------------------------

    #[test]
    fn bed_graph_section_with_one_item_round_trips() {
        assert_bed_graph_round_trips(1);
    }

    #[test]
    fn bed_graph_section_with_two_items_round_trips() {
        assert_bed_graph_round_trips(2);
    }

    #[test]
    fn bed_graph_section_with_65535_items_round_trips() {
        assert_bed_graph_round_trips(65_535);
    }

    #[test]
    fn variable_step_section_with_one_item_round_trips() {
        assert_variable_step_round_trips(1);
    }

    #[test]
    fn variable_step_section_with_two_items_round_trips() {
        assert_variable_step_round_trips(2);
    }

    #[test]
    fn variable_step_section_with_65535_items_round_trips() {
        assert_variable_step_round_trips(65_535);
    }

    #[test]
    fn fixed_step_section_with_one_item_round_trips() {
        assert_fixed_step_round_trips(1);
    }

    #[test]
    fn fixed_step_section_with_two_items_round_trips() {
        assert_fixed_step_round_trips(2);
    }

    #[test]
    fn fixed_step_section_with_65535_items_round_trips() {
        assert_fixed_step_round_trips(65_535);
    }

    #[test]
    fn fixed_step_values_pushed_in_bulk_encode_as_pushed_one_by_one() {
        let values = fixed_step_values(1_000);
        for order in BOTH_ORDERS {
            let mut builder = SectionBuilder::new(order);
            builder.begin_fixed_step(3, 500, 7, 5);
            builder.push_fixed_step(values[0]);
            builder.push_fixed_step_values(&values[1..400]);
            builder.push_fixed_step_values(&[]);
            builder.push_fixed_step_values(&values[400..]);
            assert_eq!(builder.item_count(), 1_000);
            assert_eq!(builder.end(), 500 + 999 * 7 + 5);
            assert_eq!(builder.finish(), encode_fixed_step(order, 500, 7, 5, &values));
        }
    }

    #[test]
    fn fixed_step_values_pushed_in_bulk_into_an_empty_section_start_at_its_start() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_fixed_step(3, 40, 10, 10);
        builder.push_fixed_step_values(&[1.0, 2.0]);
        assert_eq!((builder.start(), builder.end()), (40, 60));
        let items = decode(builder.finish(), ByteOrder::Little);
        assert_eq!(items, vec![item(40, 50, 1.0), item(50, 60, 2.0)]);
    }

    #[test]
    fn fixed_step_with_span_below_step_leaves_gaps_between_items() {
        for order in BOTH_ORDERS {
            let block = encode_fixed_step(order, 100, 10, 3, &[1.0, 2.0, 3.0]);
            assert_eq!(
                decode(&block, order),
                [item(100, 103, 1.0), item(110, 113, 2.0), item(120, 123, 3.0)]
            );
            assert_eq!(header_u32(&block, order, 8), 123, "end = start + (n-1)*step + span");
        }
    }

    #[test]
    fn values_round_trip_bit_for_bit() {
        let specials = [0.0, -0.0, f32::MIN_POSITIVE, f32::MAX, -f32::MAX, f32::NAN, f32::INFINITY];
        let block = encode_fixed_step(ByteOrder::Little, 0, 1, 1, &specials);
        let decoded = decode(&block, ByteOrder::Little);
        for (decoded, expected) in decoded.iter().zip(specials) {
            assert_eq!(decoded.value.to_bits(), expected.to_bits());
        }
    }

    #[test]
    fn same_items_encode_to_different_bytes_per_byte_order() {
        let items = bed_graph_items(2);
        let little = encode_bed_graph(ByteOrder::Little, &items);
        let big = encode_bed_graph(ByteOrder::Big, &items);
        assert_ne!(little, big);
        assert_eq!(&little[0..4], &[3, 0, 0, 0], "chromId little-endian");
        assert_eq!(&big[0..4], &[0, 0, 0, 3], "chromId big-endian");
    }

    // ---- header fields --------------------------------------------------------------------

    #[test]
    fn bed_graph_header_has_type_one_zero_step_and_span_and_last_end() {
        for order in BOTH_ORDERS {
            let block = encode_bed_graph(order, &[(10, 20, 1.0), (30, 45, 2.0)]);
            assert_eq!(header_u32(&block, order, 0), 3, "chromId");
            assert_eq!(header_u32(&block, order, 4), 10, "start");
            assert_eq!(header_u32(&block, order, 8), 45, "end");
            assert_eq!(header_u32(&block, order, 12), 0, "itemStep");
            assert_eq!(header_u32(&block, order, 16), 0, "itemSpan");
            assert_eq!(block[20], 1, "type");
            assert_eq!(block[21], 0, "reserved");
            assert_eq!(order.u16_from([block[22], block[23]]), 2, "itemCount");
        }
    }

    #[test]
    fn variable_step_header_has_type_two_span_and_last_start_plus_span() {
        for order in BOTH_ORDERS {
            let block = encode_variable_step(order, 25, &[(100, 1.0), (500, 2.0), (900, 3.0)]);
            assert_eq!(header_u32(&block, order, 4), 100, "start");
            assert_eq!(header_u32(&block, order, 8), 925, "end");
            assert_eq!(header_u32(&block, order, 12), 0, "itemStep");
            assert_eq!(header_u32(&block, order, 16), 25, "itemSpan");
            assert_eq!(block[20], 2, "type");
            assert_eq!(order.u16_from([block[22], block[23]]), 3, "itemCount");
        }
    }

    #[test]
    fn fixed_step_header_has_type_three_step_span_and_kent_end() {
        for order in BOTH_ORDERS {
            let block = encode_fixed_step(order, 1000, 50, 50, &[1.0, 2.0, 3.0, 4.0]);
            assert_eq!(header_u32(&block, order, 4), 1000, "start");
            assert_eq!(header_u32(&block, order, 8), 1200, "end = 1000 + 3*50 + 50");
            assert_eq!(header_u32(&block, order, 12), 50, "itemStep");
            assert_eq!(header_u32(&block, order, 16), 50, "itemSpan");
            assert_eq!(block[20], 3, "type");
            assert_eq!(order.u16_from([block[22], block[23]]), 4, "itemCount");
        }
    }

    #[test]
    fn item_count_is_a_u16_at_offset_22() {
        let block = encode_fixed_step(ByteOrder::Little, 0, 1, 1, &fixed_step_values(65_535));
        assert_eq!(&block[22..24], &[0xFF, 0xFF]);
    }

    // ---- builder state --------------------------------------------------------------------

    #[test]
    fn builder_accessors_track_the_open_section() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_variable_step(7, 5);
        assert_eq!(builder.item_count(), 0);
        assert_eq!(builder.chrom_id(), 7);
        assert_eq!(builder.kind(), SectionKind::VariableStep { span: 5 });
        builder.push_variable_step(100, 1.0);
        builder.push_variable_step(200, 2.0);
        assert_eq!(builder.item_count(), 2);
        assert_eq!(builder.start(), 100);
        assert_eq!(builder.end(), 205);
    }

    #[test]
    fn fixed_step_builder_reports_start_before_the_first_push() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_fixed_step(1, 5000, 10, 10);
        assert_eq!(builder.start(), 5000);
        assert_eq!(builder.end(), 5000);
        builder.push_fixed_step(1.0);
        assert_eq!(builder.start(), 5000);
        assert_eq!(builder.end(), 5010);
    }

    #[test]
    fn fixed_step_builder_end_tracks_each_push() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_fixed_step(1, 0, 100, 1);
        for i in 1..=5 {
            builder.push_fixed_step(0.0);
            assert_eq!(builder.end(), (i - 1) * 100 + 1);
        }
    }

    #[test]
    fn fixed_step_last_item_may_end_within_step_of_u32_max() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_fixed_step(3, u32::MAX - 15, 10, 5);
        builder.push_fixed_step(1.0);
        builder.push_fixed_step(2.0);
        assert_eq!(builder.end(), u32::MAX);
        let block = builder.finish().to_vec();
        assert_eq!(
            decode(&block, ByteOrder::Little),
            [item(u32::MAX - 15, u32::MAX - 10, 1.0), item(u32::MAX - 5, u32::MAX, 2.0)]
        );
    }

    #[test]
    fn empty_section_is_a_bare_header_with_end_equal_to_start() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_fixed_step(2, 300, 10, 10);
        let block = builder.finish().to_vec();
        assert_eq!(block.len(), SECTION_HEADER_LEN);
        assert_eq!(header_u32(&block, ByteOrder::Little, 4), 300);
        assert_eq!(header_u32(&block, ByteOrder::Little, 8), 300);
        let section = Section::parse(&block, ByteOrder::Little, 0).unwrap();
        assert_eq!(section.item_count(), 0);
        assert_eq!(section.items().count(), 0);
    }

    #[test]
    fn begin_discards_the_previous_section() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_bed_graph(1);
        builder.push_bed_graph(0, 10, 1.0);
        builder.push_bed_graph(10, 20, 2.0);
        builder.finish();
        builder.begin_fixed_step(2, 50, 5, 5);
        builder.push_fixed_step(9.0);
        let block = builder.finish().to_vec();
        assert_eq!(block.len(), SECTION_HEADER_LEN + 4);
        assert_eq!(header_u32(&block, ByteOrder::Little, 0), 2);
        assert_eq!(
            decode(&block, ByteOrder::Little),
            [DataItem { chrom_id: 2, start: 50, end: 55, value: 9.0 }]
        );
    }

    #[test]
    fn finish_can_be_called_again_without_changing_the_section() {
        let mut builder = SectionBuilder::new(ByteOrder::Big);
        builder.begin_bed_graph(1);
        builder.push_bed_graph(0, 10, 1.0);
        let first = builder.finish().to_vec();
        assert_eq!(builder.finish(), first);
        assert_eq!(builder.item_count(), 1);
    }

    #[test]
    #[should_panic(expected = "at most 65535 items")]
    fn pushing_a_65536th_item_panics() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_fixed_step(1, 0, 1, 1);
        for _ in 0..=65_535 {
            builder.push_fixed_step(0.0);
        }
    }

    #[test]
    #[should_panic(expected = "push_variable_step on a BedGraph section")]
    fn pushing_an_item_of_the_wrong_kind_panics() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_bed_graph(1);
        builder.push_variable_step(0, 0.0);
    }

    // ---- decoder: header `end` ------------------------------------------------------------

    #[test]
    fn fixed_step_decode_ignores_an_overstated_header_end() {
        // pyBigWig #166: libBigWig wrote fixedStep `end` six steps too far.
        for order in BOTH_ORDERS {
            let values = fixed_step_values(100);
            let correct = encode_fixed_step(order, 2000, 10, 10, &values);
            let mut overstated = correct.clone();
            let bad_end = header_u32(&correct, order, 8) + 6 * 10;
            overstated[8..12].copy_from_slice(&order.u32_bytes(bad_end));

            let section = Section::parse(&overstated, order, 0).unwrap();
            assert_eq!(section.end_hint(), bad_end);
            assert_eq!(decode(&overstated, order), decode(&correct, order));
            assert_eq!(decode(&overstated, order).last().unwrap().end, 2000 + 99 * 10 + 10);
        }
    }

    #[test]
    fn decode_ignores_an_understated_header_end() {
        let correct = encode_bed_graph(ByteOrder::Little, &bed_graph_items(10));
        let mut understated = correct.clone();
        understated[8..12].copy_from_slice(&0u32.to_le_bytes());
        assert_eq!(decode(&understated, ByteOrder::Little), decode(&correct, ByteOrder::Little));
    }

    // ---- decoder: header checks -----------------------------------------------------------

    #[test]
    fn parsed_header_fields_are_exposed() {
        let block = encode_fixed_step(ByteOrder::Big, 500, 20, 10, &[1.0, 2.0]);
        let section = Section::parse(&block, ByteOrder::Big, 0).unwrap();
        assert_eq!(section.chrom_id(), 3);
        assert_eq!(section.start(), 500);
        assert_eq!(section.end_hint(), 530);
        assert_eq!(section.kind(), SectionKind::FixedStep { step: 20, span: 10 });
        assert_eq!(section.item_count(), 2);
    }

    #[test]
    fn unknown_section_type_is_corrupt() {
        let mut block = encode_bed_graph(ByteOrder::Little, &bed_graph_items(1));
        for type_code in [0, 4, 0x80, 0xFF] {
            block[20] = type_code;
            match Section::parse(&block, ByteOrder::Little, 555) {
                Err(Error::Corrupt { offset, reason }) => {
                    assert_eq!(offset, 555);
                    assert!(reason.contains("type"), "{reason}");
                }
                other => panic!("expected Corrupt, got {other:?}"),
            }
        }
    }

    #[test]
    fn item_count_past_the_block_is_corrupt() {
        let mut block = encode_bed_graph(ByteOrder::Little, &bed_graph_items(2));
        block[22..24].copy_from_slice(&3u16.to_le_bytes());
        match Section::parse(&block, ByteOrder::Little, 555) {
            Err(Error::Corrupt { offset, reason }) => {
                assert_eq!(offset, 555);
                assert!(reason.contains("itemCount"), "{reason}");
            }
            other => panic!("expected Corrupt, got {other:?}"),
        }
    }

    #[test]
    fn item_count_filling_the_block_exactly_is_accepted() {
        let block = encode_fixed_step(ByteOrder::Little, 0, 1, 1, &fixed_step_values(5));
        assert_eq!(Section::parse(&block, ByteOrder::Little, 0).unwrap().item_count(), 5);
    }

    #[test]
    fn trailing_bytes_after_the_items_are_tolerated() {
        let mut block = encode_variable_step(ByteOrder::Little, 5, &variable_step_items(3));
        block.extend_from_slice(&[0xAB; 11]);
        assert_eq!(decode(&block, ByteOrder::Little).len(), 3);
    }

    #[test]
    fn truncated_header_is_corrupt() {
        let block = encode_bed_graph(ByteOrder::Little, &bed_graph_items(1));
        for len in [0, 1, 23] {
            assert!(matches!(
                Section::parse(&block[..len], ByteOrder::Little, 9),
                Err(Error::Corrupt { offset: 9, .. })
            ));
        }
    }

    #[test]
    fn bed_graph_header_step_and_span_are_ignored() {
        let mut block = encode_bed_graph(ByteOrder::Little, &[(10, 20, 1.0)]);
        block[12..16].copy_from_slice(&99u32.to_le_bytes());
        block[16..20].copy_from_slice(&77u32.to_le_bytes());
        let section = Section::parse(&block, ByteOrder::Little, 0).unwrap();
        assert_eq!(section.kind(), SectionKind::BedGraph);
        assert_eq!(decode(&block, ByteOrder::Little), [item(10, 20, 1.0)]);
    }

    // ---- decoder: item checks -------------------------------------------------------------

    #[test]
    fn bed_graph_item_with_start_after_end_is_corrupt() {
        for (start, end) in [(11, 10), (u32::MAX, 0)] {
            let block = encode_bed_graph(ByteOrder::Little, &[(0, 5, 1.0), (start, end, 2.0)]);
            let section = Section::parse(&block, ByteOrder::Little, 321).unwrap();
            let mut items = section.items();
            assert_eq!(items.next().unwrap().unwrap(), item(0, 5, 1.0));
            assert!(matches!(items.next(), Some(Err(Error::Corrupt { offset: 321, .. }))));
        }
    }

    #[test]
    fn zero_length_bed_graph_item_decodes_without_stopping_its_neighbours() {
        let block =
            encode_bed_graph(ByteOrder::Little, &[(0, 5, 1.0), (10, 10, 2.0), (20, 30, 3.0)]);
        assert_eq!(
            decode(&block, ByteOrder::Little),
            [item(0, 5, 1.0), item(10, 10, 2.0), item(20, 30, 3.0)]
        );
    }

    #[test]
    fn variable_step_item_ending_past_u32_max_is_corrupt() {
        let block = encode_variable_step(ByteOrder::Little, 10, &[(u32::MAX - 10, 1.0)]);
        assert_eq!(decode(&block, ByteOrder::Little), [item(u32::MAX - 10, u32::MAX, 1.0)]);
        // The builder never encodes an overflowing item, so widen the header's span by hand.
        let mut block = encode_variable_step(ByteOrder::Little, 1, &[(u32::MAX - 8, 1.0)]);
        block[16..20].copy_from_slice(&10u32.to_le_bytes());
        let mut items = Section::parse(&block, ByteOrder::Little, 5).unwrap().items();
        assert!(matches!(items.next(), Some(Err(Error::Corrupt { offset: 5, .. }))));
    }

    #[test]
    fn variable_step_with_zero_span_yields_zero_length_items() {
        let block = encode_variable_step(ByteOrder::Big, 0, &[(100, 1.0), (200, 2.0)]);
        assert_eq!(decode(&block, ByteOrder::Big), [item(100, 100, 1.0), item(200, 200, 2.0)]);
    }

    /// Encodes a valid two-item fixedStep section, then moves its header `start` to
    /// `header_start` so the second item runs past `u32::MAX` (the builder never encodes that).
    fn fixed_step_near_u32_max(header_start: u32, step: u32, span: u32) -> Vec<u8> {
        let mut block = encode_fixed_step(ByteOrder::Little, 0, step, span, &[1.0, 2.0]);
        block[4..8].copy_from_slice(&header_start.to_le_bytes());
        block
    }

    #[test]
    fn fixed_step_item_starting_past_u32_max_is_corrupt() {
        let block = fixed_step_near_u32_max(u32::MAX - 5, 10, 5);
        let mut items = Section::parse(&block, ByteOrder::Little, 5).unwrap().items();
        assert_eq!(items.next().unwrap().unwrap(), item(u32::MAX - 5, u32::MAX, 1.0));
        assert!(matches!(items.next(), Some(Err(Error::Corrupt { offset: 5, .. }))));
    }

    #[test]
    fn fixed_step_item_ending_past_u32_max_is_corrupt() {
        let block = fixed_step_near_u32_max(u32::MAX - 15, 10, 10);
        let mut items = Section::parse(&block, ByteOrder::Little, 5).unwrap().items();
        assert_eq!(items.next().unwrap().unwrap(), item(u32::MAX - 15, u32::MAX - 5, 1.0));
        assert!(matches!(items.next(), Some(Err(Error::Corrupt { offset: 5, .. }))));
    }

    #[test]
    fn fixed_step_with_zero_span_yields_zero_length_items() {
        let block = encode_fixed_step(ByteOrder::Little, 100, 10, 0, &[1.0, 2.0]);
        assert_eq!(decode(&block, ByteOrder::Little), [item(100, 100, 1.0), item(110, 110, 2.0)]);
    }

    #[test]
    fn fixed_step_with_zero_step_stacks_items_at_the_start() {
        let block = encode_fixed_step(ByteOrder::Little, 100, 0, 10, &[1.0, 2.0]);
        assert_eq!(decode(&block, ByteOrder::Little), [item(100, 110, 1.0), item(100, 110, 2.0)]);
    }

    #[test]
    fn iterator_yields_nothing_after_an_error() {
        let block = encode_bed_graph(ByteOrder::Little, &[(6, 5, 1.0), (10, 20, 2.0)]);
        let mut items = Section::parse(&block, ByteOrder::Little, 0).unwrap().items();
        assert!(items.next().unwrap().is_err());
        assert!(items.next().is_none());
        assert_eq!(items.len(), 0);
    }

    // ---- decoder: iterator ----------------------------------------------------------------

    #[test]
    fn iterator_reports_its_exact_length() {
        let block = encode_fixed_step(ByteOrder::Little, 0, 1, 1, &fixed_step_values(4));
        let mut items = Section::parse(&block, ByteOrder::Little, 0).unwrap().items();
        assert_eq!(items.len(), 4);
        items.next();
        assert_eq!(items.len(), 3);
        assert_eq!(items.size_hint(), (3, Some(3)));
    }

    #[test]
    fn nth_skips_straight_to_the_requested_item() {
        let block = encode_fixed_step(ByteOrder::Little, 1000, 10, 10, &fixed_step_values(100));
        let mut items = Section::parse(&block, ByteOrder::Little, 0).unwrap().items();
        assert_eq!(items.nth(42).unwrap().unwrap(), item(1420, 1430, value(42)));
        assert_eq!(items.next().unwrap().unwrap(), item(1430, 1440, value(43)));
        assert!(items.nth(1000).is_none());
        assert_eq!(items.len(), 0);
    }

    #[test]
    fn items_carry_the_section_chrom_id() {
        let mut builder = SectionBuilder::new(ByteOrder::Little);
        builder.begin_bed_graph(41);
        builder.push_bed_graph(1, 2, 0.5);
        let block = builder.finish().to_vec();
        assert_eq!(decode(&block, ByteOrder::Little)[0].chrom_id, 41);
    }

    proptest! {
        #[test]
        fn arbitrary_bytes_never_panic_the_decoder(
            block in prop::collection::vec(any::<u8>(), 0..1024),
            big_endian: bool,
        ) {
            let order = if big_endian { ByteOrder::Big } else { ByteOrder::Little };
            if let Ok(section) = Section::parse(&block, order, 0) {
                for item in section.items() {
                    let _ = item;
                }
            }
        }

        #[test]
        fn arbitrary_well_formed_items_round_trip(
            items in prop::collection::vec((0u32..1_000_000, 1u32..10_000, any::<f32>()), 0..200),
            big_endian: bool,
        ) {
            let order = if big_endian { ByteOrder::Big } else { ByteOrder::Little };
            let intervals: Vec<_> = items.iter().map(|&(s, len, v)| (s, s + len, v)).collect();
            let block = encode_bed_graph(order, &intervals);
            let decoded = decode(&block, order);
            prop_assert_eq!(decoded.len(), intervals.len());
            for (decoded, &(s, e, v)) in decoded.iter().zip(&intervals) {
                prop_assert_eq!((decoded.start, decoded.end, decoded.value.to_bits()), (s, e, v.to_bits()));
            }
        }
    }
}
