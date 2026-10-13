//! The R-tree ("cirTree") index over data sections and zoom blocks: bulk writing and overlap
//! queries.

use std::collections::HashSet;
use std::io::{Read, Seek};

use crate::bbi::{
    ByteOrder, ByteReader, MAX_TREE_DEPTH, RTREE_MAGIC, read_exact_at, read_exact_at_into,
};
use crate::error::{Error, Result};

/// Size of the index header.
const HEADER_SIZE: usize = 48;
/// Size of a node header: `isLeaf` (u8), reserved (u8), `count` (u16).
const NODE_HEADER_SIZE: usize = 4;
/// Size of a leaf slot: four u32 bounds, then `dataOffset` and `dataSize` (u64 each).
const LEAF_SLOT_SIZE: usize = 32;
/// Size of a non-leaf slot: four u32 bounds, then `childOffset` (u64).
const NON_LEAF_SLOT_SIZE: usize = 24;
/// Most slots a node can hold: `count` is a u16.
const MAX_SLOTS: usize = u16::MAX as usize;

/// An R-tree index in a file, from its header: enough to walk its nodes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RTree {
    /// Byte order of the index, detected from its own magic.
    byte_order: ByteOrder,
    /// Absolute offset of the root node, or `None` when the index is empty and has nothing to
    /// read.
    root_offset: Option<u64>,
    /// Bytes read for each node: a node of `blockSize` leaf slots, the larger kind.
    node_read_len: usize,
    /// Length of the file the index is in; every node and block must lie inside it.
    file_len: u64,
}

impl RTree {
    /// Reads the header of the index at absolute `offset` in a source of `file_len` bytes.
    ///
    /// An `offset` of 0 (libBigWig's empty files) or an `itemCount` of 0 (kent writes no root
    /// node then) gives an empty index without reading further. Fails with [`Error::Corrupt`] on
    /// a truncated header or a bad magic.
    pub(crate) fn read<R: Read + Seek>(source: &mut R, offset: u64, file_len: u64) -> Result<Self> {
        // The fixed header occupies offset 0, so no index can be there.
        if offset == 0 {
            return Ok(Self {
                byte_order: ByteOrder::default(),
                root_offset: None,
                node_read_len: NODE_HEADER_SIZE,
                file_len,
            });
        }
        let header = read_exact_at(source, offset, HEADER_SIZE, file_len)?;
        let found = ByteReader::new(&header, ByteOrder::Little, offset).u32()?;
        let byte_order = ByteOrder::detect(found, RTREE_MAGIC)
            .ok_or_else(|| Error::corrupt(offset, format!("bad R-tree magic {found:#010x}")))?;
        let mut reader = ByteReader::new(&header, byte_order, offset);
        reader.skip(4)?; // magic
        let block_size = usize::try_from(reader.u32()?).unwrap_or(MAX_SLOTS).min(MAX_SLOTS);
        let item_count = reader.u64()?;
        let root_offset = (item_count > 0).then_some(offset + HEADER_SIZE as u64);
        let node_read_len = NODE_HEADER_SIZE + block_size * LEAF_SLOT_SIZE;
        Ok(Self { byte_order, root_offset, node_read_len, file_len })
    }

    /// Returns the leaf entries overlapping `[start, end)` on `chrom_id`, in index order.
    ///
    /// Overlap is kent's lexicographic test on (chromIx, base) pairs, so an entry spanning several
    /// chromosomes overlaps every query between its bounds, and an empty query `[p, p)` matches
    /// entries that strictly contain `p`. Every overlapping node is searched, with no early exit,
    /// so the result is exact even when leaves are not in `chromId` order.
    ///
    /// Fails with [`Error::Corrupt`] as [`RTree::leaves`] does.
    pub(crate) fn overlapping<R: Read + Seek>(
        &self,
        source: &mut R,
        chrom_id: u32,
        start: u32,
        end: u32,
    ) -> Result<Vec<LeafEntry>> {
        let mut entries = Vec::new();
        let mut buffers = WalkBuffers::default();
        self.overlapping_into(source, chrom_id, start, end, &mut buffers, &mut entries)?;
        Ok(entries)
    }

    /// As [`RTree::overlapping`], replacing the contents of `entries` and reusing `buffers`, so
    /// repeated queries allocate nothing once the buffers have grown.
    pub(crate) fn overlapping_into<R: Read + Seek>(
        &self,
        source: &mut R,
        chrom_id: u32,
        start: u32,
        end: u32,
        buffers: &mut WalkBuffers,
        entries: &mut Vec<LeafEntry>,
    ) -> Result<()> {
        self.walk(source, |bounds| bounds.overlaps(chrom_id, start, end), buffers, entries)
    }

    /// Returns every leaf entry, in index order.
    ///
    /// Fails with [`Error::Corrupt`] at the offending offset on a node or child offset past the
    /// end of the file, an `isLeaf` other than 0 or 1, a node reached twice (a cyclic or shared
    /// child offset), a tree deeper than [`MAX_TREE_DEPTH`], or a returned entry whose block
    /// extends past the end of the file.
    pub(crate) fn leaves<R: Read + Seek>(&self, source: &mut R) -> Result<Vec<LeafEntry>> {
        let mut entries = Vec::new();
        self.walk(source, |_| true, &mut WalkBuffers::default(), &mut entries)?;
        Ok(entries)
    }

    /// Walks the tree depth first, descending into the children whose bounds `keep` accepts, and
    /// replaces the contents of `entries` with the accepted leaf entries in index order.
    ///
    /// Each node is fetched with one read of `blockSize` leaf slots, stopping at the end of the
    /// file, and only its `count` slots are decoded, so padding is never relied on (kent pads its
    /// last leaf with 24-byte slots, leaving it shorter than a full leaf). Bounds are never
    /// checked for order, so zero-width and inverted entries are tolerated.
    fn walk<R: Read + Seek>(
        &self,
        source: &mut R,
        keep: impl Fn(&Bounds) -> bool,
        buffers: &mut WalkBuffers,
        entries: &mut Vec<LeafEntry>,
    ) -> Result<()> {
        entries.clear();
        let Some(root_offset) = self.root_offset else {
            return Ok(());
        };
        let WalkBuffers { pending, visited, node: buf } = buffers;
        pending.clear();
        visited.clear();
        pending.push((root_offset, 0));
        while let Some((node_offset, depth)) = pending.pop() {
            if depth >= MAX_TREE_DEPTH {
                return Err(Error::corrupt(
                    node_offset,
                    format!("R-tree is deeper than {MAX_TREE_DEPTH} levels"),
                ));
            }
            if !visited.insert(node_offset) {
                return Err(Error::corrupt(
                    node_offset,
                    "R-tree node reached twice (cyclic or shared child offset)",
                ));
            }
            // A node's header and slots come in one read: two would cost two seeks, and a
            // `BufReader` discards its buffer on every seek.
            let to_end = usize::try_from(self.file_len.saturating_sub(node_offset));
            let node_len = self.node_read_len.min(to_end.unwrap_or(usize::MAX));
            let node_len = node_len.max(NODE_HEADER_SIZE);
            read_exact_at_into(source, node_offset, node_len, self.file_len, buf)?;
            let mut header = ByteReader::new(buf, self.byte_order, node_offset);
            let is_leaf = header.u8()?;
            header.skip(1)?;
            let count = usize::from(header.u16()?);
            let slot_size = match is_leaf {
                0 => NON_LEAF_SLOT_SIZE,
                1 => LEAF_SLOT_SIZE,
                _ => {
                    return Err(Error::corrupt(node_offset, format!("R-tree isLeaf is {is_leaf}")));
                }
            };

            // A node with more slots than `blockSize`, or one cut short by the end of the file,
            // is read again from its slots, where a truncation is then reported.
            let slots_offset = node_offset + NODE_HEADER_SIZE as u64;
            let slots_len = count * slot_size;
            let slots_start = if NODE_HEADER_SIZE + slots_len <= buf.len() {
                NODE_HEADER_SIZE
            } else {
                read_exact_at_into(source, slots_offset, slots_len, self.file_len, buf)?;
                0
            };
            let mut slots = ByteReader::new(&buf[slots_start..], self.byte_order, slots_offset);
            if is_leaf == 1 {
                for _ in 0..count {
                    let slot_offset = slots.file_offset();
                    let entry = LeafEntry::decode(&mut slots)?;
                    if keep(&entry.bounds()) {
                        entry.check_block_in_file(slot_offset, self.file_len)?;
                        entries.push(entry);
                    }
                }
            } else {
                let first_child = pending.len();
                for _ in 0..count {
                    let bounds = Bounds::decode(&mut slots)?;
                    let child_offset = slots.u64()?;
                    if keep(&bounds) {
                        pending.push((child_offset, depth + 1));
                    }
                }
                // The stack pops the last push first; reversing visits children in slot order.
                pending[first_child..].reverse();
            }
        }
        Ok(())
    }
}

/// Scratch space for walking an index, kept between queries so they don't allocate.
#[derive(Debug, Default)]
pub(crate) struct WalkBuffers {
    /// Nodes still to visit, as (offset, depth).
    pending: Vec<(u64, usize)>,
    /// Offsets of the nodes visited so far, to detect cycles and shared children.
    visited: HashSet<u64>,
    /// The bytes of the node being decoded.
    node: Vec<u8>,
}

/// One leaf slot: the extent of one block (a data section or a zoom block) and where the block
/// is in the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct LeafEntry {
    /// `chromId` of the block's first item.
    pub(crate) start_chrom_ix: u32,
    /// Start of the block's first item.
    pub(crate) start_base: u32,
    /// `chromId` of the block's last item.
    pub(crate) end_chrom_ix: u32,
    /// End (exclusive) of the block's last item.
    pub(crate) end_base: u32,
    /// Absolute file offset of the (compressed) block.
    pub(crate) data_offset: u64,
    /// Size in bytes of the (compressed) block.
    pub(crate) data_size: u64,
}

impl LeafEntry {
    /// Decodes one 32-byte leaf slot.
    fn decode(reader: &mut ByteReader<'_>) -> Result<Self> {
        let bounds = Bounds::decode(reader)?;
        Ok(Self {
            start_chrom_ix: bounds.start.0,
            start_base: bounds.start.1,
            end_chrom_ix: bounds.end.0,
            end_base: bounds.end.1,
            data_offset: reader.u64()?,
            data_size: reader.u64()?,
        })
    }

    /// The entry's extent as lexicographic (chromIx, base) bounds.
    fn bounds(&self) -> Bounds {
        Bounds {
            start: (self.start_chrom_ix, self.start_base),
            end: (self.end_chrom_ix, self.end_base),
        }
    }

    /// Fails with [`Error::Corrupt`] at the slot's `dataOffset` field unless the block lies
    /// inside a file of `file_len` bytes. Some real files carry garbage sizes, so this is checked
    /// before any caller sizes a read from them.
    fn check_block_in_file(&self, slot_offset: u64, file_len: u64) -> Result<()> {
        if self.data_offset.checked_add(self.data_size).is_none_or(|end| end > file_len) {
            return Err(Error::corrupt(
                slot_offset + 16,
                format!(
                    "R-tree block of {} bytes at offset {} extends past the end of the \
                     {file_len}-byte file",
                    self.data_size, self.data_offset
                ),
            ));
        }
        Ok(())
    }

    /// Appends the 32-byte leaf slot to `buf` in byte order `order`.
    fn encode(&self, order: ByteOrder, buf: &mut Vec<u8>) {
        self.bounds().encode(order, buf);
        order.put_u64(buf, self.data_offset);
        order.put_u64(buf, self.data_size);
    }
}

/// The extent of an entry or a node, from `start` to `end`, each a (chromIx, base) pair compared
/// lexicographically.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Bounds {
    /// (chromIx, base) of the start.
    start: (u32, u32),
    /// (chromIx, base) of the (exclusive) end.
    end: (u32, u32),
}

impl Bounds {
    /// Decodes four u32 bounds: `startChromIx`, `startBase`, `endChromIx`, `endBase`.
    fn decode(reader: &mut ByteReader<'_>) -> Result<Self> {
        Ok(Self { start: (reader.u32()?, reader.u32()?), end: (reader.u32()?, reader.u32()?) })
    }

    /// The smallest bounds enclosing all of `bounds`: the lexicographic minimum start and maximum
    /// end, as kent computes them. All zero when `bounds` is empty.
    fn enclosing(bounds: impl IntoIterator<Item = Bounds>) -> Self {
        bounds
            .into_iter()
            .reduce(|a, b| Self { start: a.start.min(b.start), end: a.end.max(b.end) })
            .unwrap_or_default()
    }

    /// kent's `cirTreeOverlaps`: whether `[start, end)` on `chrom_id` overlaps these bounds,
    /// comparing (chromIx, base) pairs lexicographically and treating both ranges as half-open.
    fn overlaps(&self, chrom_id: u32, start: u32, end: u32) -> bool {
        (chrom_id, start) < self.end && (chrom_id, end) > self.start
    }

    /// Appends the four u32 bounds to `buf` in byte order `order`.
    fn encode(&self, order: ByteOrder, buf: &mut Vec<u8>) {
        order.put_u32(buf, self.start.0);
        order.put_u32(buf, self.start.1);
        order.put_u32(buf, self.end.0);
        order.put_u32(buf, self.end.1);
    }
}

/// Encodes `entries` as an R-tree to be written at absolute file offset `offset`, in kent's
/// bulk-loaded layout: the header, then the root, then each level down to the leaves.
///
/// Leaves hold up to `block_size` entries in the order given. Each node's bounds are the
/// lexicographic minimum and maximum of its children's, and every node is padded to `block_size`
/// zero-filled slots of its own slot size. The header records `itemCount` as the number of
/// entries and `items_per_slot` as given. Its `endFileOffset` is `offset`, since the data the
/// index covers ends where the index starts. An empty index is the header and one empty root
/// leaf.
///
/// # Panics
///
/// If `block_size` is outside 2..=65,535: smaller fan-outs never reduce to a root, and a node's
/// `count` is a u16.
pub(crate) fn encode(
    entries: &[LeafEntry],
    offset: u64,
    block_size: u32,
    items_per_slot: u32,
    order: ByteOrder,
) -> Vec<u8> {
    assert!(
        (2..=u32::from(u16::MAX)).contains(&block_size),
        "R-tree block size {block_size} is outside 2..=65535"
    );
    let fan_out = block_size as usize;
    let leaf_node_size = NODE_HEADER_SIZE + fan_out * LEAF_SLOT_SIZE;
    let non_leaf_node_size = NODE_HEADER_SIZE + fan_out * NON_LEAF_SLOT_SIZE;

    // Bounds of every node, level by level from the leaves (level 0) up to the single root.
    let mut levels: Vec<Vec<Bounds>> = vec![
        entries
            .chunks(fan_out)
            .map(|leaf| Bounds::enclosing(leaf.iter().map(LeafEntry::bounds)))
            .collect(),
    ];
    if entries.is_empty() {
        levels[0].push(Bounds::default());
    }
    while levels[levels.len() - 1].len() > 1 {
        let parents = levels[levels.len() - 1]
            .chunks(fan_out)
            .map(|children| Bounds::enclosing(children.iter().copied()))
            .collect();
        levels.push(parents);
    }
    let node_size = |level: usize| if level == 0 { leaf_node_size } else { non_leaf_node_size };

    // Absolute offset of each level's first node; levels are written root first.
    let mut level_offsets = vec![0; levels.len()];
    let mut next_offset = offset + HEADER_SIZE as u64;
    for level in (0..levels.len()).rev() {
        level_offsets[level] = next_offset;
        next_offset += (levels[level].len() * node_size(level)) as u64;
    }
    let mut bytes = Vec::with_capacity(usize::try_from(next_offset - offset).unwrap_or(0));

    order.put_u32(&mut bytes, RTREE_MAGIC);
    order.put_u32(&mut bytes, block_size);
    order.put_u64(&mut bytes, entries.len() as u64);
    levels[levels.len() - 1][0].encode(order, &mut bytes);
    order.put_u64(&mut bytes, offset);
    order.put_u32(&mut bytes, items_per_slot);
    order.put_u32(&mut bytes, 0);

    // The i-th slot of node n on level L points at node n * fan_out + i on level L - 1.
    for level in (1..levels.len()).rev() {
        for (node, children) in levels[level - 1].chunks(fan_out).enumerate() {
            put_node_header(&mut bytes, order, false, children.len());
            for (slot, child_bounds) in children.iter().enumerate() {
                let child = node * fan_out + slot;
                child_bounds.encode(order, &mut bytes);
                let child_offset = level_offsets[level - 1] + (child * node_size(level - 1)) as u64;
                order.put_u64(&mut bytes, child_offset);
            }
            bytes.resize(bytes.len() + (fan_out - children.len()) * NON_LEAF_SLOT_SIZE, 0);
        }
    }

    if entries.is_empty() {
        put_node_header(&mut bytes, order, true, 0);
        bytes.resize(bytes.len() + fan_out * LEAF_SLOT_SIZE, 0);
    }
    for leaf in entries.chunks(fan_out) {
        put_node_header(&mut bytes, order, true, leaf.len());
        for entry in leaf {
            entry.encode(order, &mut bytes);
        }
        bytes.resize(bytes.len() + (fan_out - leaf.len()) * LEAF_SLOT_SIZE, 0);
    }
    debug_assert_eq!(offset + bytes.len() as u64, next_offset);
    bytes
}

/// Appends a node header: `isLeaf`, a zero reserved byte and the slot count.
fn put_node_header(bytes: &mut Vec<u8>, order: ByteOrder, is_leaf: bool, count: usize) {
    bytes.push(u8::from(is_leaf));
    bytes.push(0);
    order.put_u16(bytes, u16::try_from(count).expect("a node holds at most 65,535 slots"));
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, SeekFrom};

    use proptest::prelude::*;

    use super::*;

    /// Where test indexes sit in their files, as a slice index: after some filler, since offset 0
    /// means "no index".
    const INDEX_START: usize = 100;
    /// [`INDEX_START`] as a file offset.
    const INDEX_OFFSET: u64 = INDEX_START as u64;

    /// Leaf entries for `(chromIx, start, end)` blocks, each on one chromosome. Block `i` has data
    /// offset `i` and size 1, so every block is distinct and lies inside any test file.
    fn numbered(blocks: &[(u32, u32, u32)]) -> Vec<LeafEntry> {
        blocks
            .iter()
            .enumerate()
            .map(|(index, &(chrom, start, end))| LeafEntry {
                start_chrom_ix: chrom,
                start_base: start,
                end_chrom_ix: chrom,
                end_base: end,
                data_offset: index as u64,
                data_size: 1,
            })
            .collect()
    }

    /// `per_chrom` adjacent 100-base blocks on each of `chroms`, one chromosome after another in
    /// the order given.
    fn tiled(chroms: &[u32], per_chrom: u32) -> Vec<LeafEntry> {
        let blocks: Vec<(u32, u32, u32)> = chroms
            .iter()
            .flat_map(|&chrom| (0..per_chrom).map(move |i| (chrom, i * 100, (i + 1) * 100)))
            .collect();
        numbered(&blocks)
    }

    /// A file holding the index of `entries` at [`INDEX_OFFSET`], after filler.
    fn file_with_index(entries: &[LeafEntry], block_size: u32, order: ByteOrder) -> Vec<u8> {
        let mut file = vec![0xAA; INDEX_START];
        file.extend(encode(entries, INDEX_OFFSET, block_size, 1, order));
        file
    }

    /// Every leaf entry of the index at [`INDEX_OFFSET`] in `file`.
    fn leaves_of(file: &[u8]) -> Result<Vec<LeafEntry>> {
        let mut source = Cursor::new(file);
        RTree::read(&mut source, INDEX_OFFSET, file.len() as u64)?.leaves(&mut source)
    }

    /// The leaf entries of the index at [`INDEX_OFFSET`] in `file` that overlap `[start, end)` on
    /// `chrom`.
    fn overlapping_in(file: &[u8], chrom: u32, start: u32, end: u32) -> Result<Vec<LeafEntry>> {
        let mut source = Cursor::new(file);
        RTree::read(&mut source, INDEX_OFFSET, file.len() as u64)?.overlapping(
            &mut source,
            chrom,
            start,
            end,
        )
    }

    /// The single-chromosome `entries` that overlap `[start, end)` on `chrom`, in order, by a
    /// direct interval test.
    fn brute_force(entries: &[LeafEntry], chrom: u32, start: u32, end: u32) -> Vec<LeafEntry> {
        entries
            .iter()
            .filter(|e| e.start_chrom_ix == chrom && e.start_base < end && e.end_base > start)
            .copied()
            .collect()
    }

    /// Asserts that queries over many windows of every chromosome in `0..chroms`, plus one past
    /// them, return exactly the brute-force overlaps.
    fn assert_queries_match_brute_force(file: &[u8], entries: &[LeafEntry], chroms: u32) {
        for chrom in 0..=chroms {
            for start in (0..1200).step_by(50) {
                for len in [0, 1, 49, 100, 250, 1000] {
                    let got = overlapping_in(file, chrom, start, start + len).unwrap();
                    let want = brute_force(entries, chrom, start, start + len);
                    assert_eq!(got, want, "chrom {chrom} [{start}, {})", start + len);
                }
            }
        }
    }

    /// The little-endian u32 at `at` in `bytes`.
    fn u32_at(bytes: &[u8], at: usize) -> u32 {
        u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap())
    }

    /// The little-endian u64 at `at` in `bytes`.
    fn u64_at(bytes: &[u8], at: usize) -> u64 {
        u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap())
    }

    /// The `(isLeaf, count)` of the node at file offset `at`.
    fn node_header(file: &[u8], at: usize) -> (u8, u16) {
        (file[at], u16::from_le_bytes(file[at + 2..at + 4].try_into().unwrap()))
    }

    /// File offset of the root node of the index at [`INDEX_OFFSET`].
    const ROOT: usize = INDEX_START + HEADER_SIZE;

    /// File offset of the child offset in the root's `slot`-th (non-leaf) slot.
    fn root_child_offset_field(slot: usize) -> usize {
        ROOT + NODE_HEADER_SIZE + slot * NON_LEAF_SLOT_SIZE + 16
    }

    /// A little-endian file holding, at [`INDEX_OFFSET`], an index of `levels` levels:
    /// single-child non-leaf nodes, without padding, chained down to a leaf holding `block`.
    fn chain_index(levels: usize, block: LeafEntry) -> Vec<u8> {
        let order = ByteOrder::Little;
        let mut file = file_with_index(&[block], 2, order);
        file.truncate(ROOT);
        let node_size = (NODE_HEADER_SIZE + NON_LEAF_SLOT_SIZE) as u64;
        for level in 1..levels as u64 {
            put_node_header(&mut file, order, false, 1);
            block.bounds().encode(order, &mut file);
            order.put_u64(&mut file, ROOT as u64 + level * node_size);
        }
        put_node_header(&mut file, order, true, 1);
        block.encode(order, &mut file);
        file
    }

    /// A `Read + Seek` source that counts the reads made from it.
    struct CountingReads<R> {
        inner: R,
        reads: usize,
    }

    impl<R: Read> Read for CountingReads<R> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            self.reads += 1;
            self.inner.read(buf)
        }
    }

    impl<R: Seek> Seek for CountingReads<R> {
        fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
            self.inner.seek(position)
        }
    }

    /// The number of reads walking every node of the index at [`INDEX_OFFSET`] in `file` takes,
    /// not counting the read of the index header.
    fn reads_to_list_leaves(file: &[u8]) -> usize {
        let mut source = CountingReads { inner: Cursor::new(file), reads: 0 };
        let index = RTree::read(&mut source, INDEX_OFFSET, file.len() as u64).unwrap();
        source.reads = 0;
        index.leaves(&mut source).unwrap();
        source.reads
    }

    #[test]
    fn empty_index_is_the_header_and_one_padded_empty_root_leaf() {
        let index = encode(&[], 5000, 256, 1, ByteOrder::Little);
        assert_eq!(index.len(), HEADER_SIZE + NODE_HEADER_SIZE + 256 * LEAF_SLOT_SIZE);
        assert_eq!(u64_at(&index, 8), 0);
        assert_eq!(node_header(&index, HEADER_SIZE), (1, 0));
        assert!(index[HEADER_SIZE + 1..].iter().all(|&byte| byte == 0));
    }

    #[test]
    fn empty_index_has_no_leaves() {
        assert_eq!(leaves_of(&file_with_index(&[], 256, ByteOrder::Little)).unwrap(), []);
    }

    #[test]
    fn empty_index_overlaps_nothing() {
        let file = file_with_index(&[], 256, ByteOrder::Big);
        assert_eq!(overlapping_in(&file, 0, 0, u32::MAX).unwrap(), []);
    }

    #[test]
    fn kent_empty_index_without_a_root_reads_as_empty() {
        let mut file = file_with_index(&[], 256, ByteOrder::Little);
        file.truncate(ROOT);
        assert_eq!(leaves_of(&file).unwrap(), []);
        assert_eq!(overlapping_in(&file, 0, 0, 100).unwrap(), []);
    }

    #[test]
    fn index_offset_zero_reads_as_empty() {
        let file = vec![0xAA; 64];
        let mut source = Cursor::new(&file);
        let index = RTree::read(&mut source, 0, 64).unwrap();
        assert_eq!(index.leaves(&mut source).unwrap(), []);
    }

    #[test]
    fn header_records_the_parameters_and_the_entry_count() {
        let index = encode(&tiled(&[0], 5), 7000, 64, 1024, ByteOrder::Little);
        assert_eq!(u32_at(&index, 0), RTREE_MAGIC);
        assert_eq!(u32_at(&index, 4), 64);
        assert_eq!(u64_at(&index, 8), 5);
        assert_eq!(u64_at(&index, 32), 7000, "endFileOffset is the index's own offset");
        assert_eq!(u32_at(&index, 40), 1024);
        assert_eq!(u32_at(&index, 44), 0);
    }

    #[test]
    fn header_bounds_are_the_lexicographic_extent_of_all_entries() {
        let entries = [
            LeafEntry {
                start_chrom_ix: 2,
                start_base: 50,
                end_chrom_ix: 2,
                end_base: 60,
                ..tiled(&[0], 1)[0]
            },
            LeafEntry {
                start_chrom_ix: 0,
                start_base: 500,
                end_chrom_ix: 0,
                end_base: 600,
                ..tiled(&[0], 1)[0]
            },
            LeafEntry {
                start_chrom_ix: 1,
                start_base: 0,
                end_chrom_ix: 3,
                end_base: 10,
                ..tiled(&[0], 1)[0]
            },
        ];
        let index = encode(&entries, 7000, 2, 1, ByteOrder::Little);
        let bounds: Vec<u32> = (16..32).step_by(4).map(|at| u32_at(&index, at)).collect();
        assert_eq!(bounds, [0, 500, 3, 10]);
    }

    #[test]
    fn up_to_block_size_entries_fit_a_single_leaf_root() {
        let file = file_with_index(&tiled(&[0], 4), 4, ByteOrder::Little);
        assert_eq!(node_header(&file, ROOT), (1, 4));
        assert_eq!(file.len(), ROOT + NODE_HEADER_SIZE + 4 * LEAF_SLOT_SIZE);
    }

    #[test]
    fn one_entry_past_block_size_adds_a_level() {
        let file = file_with_index(&tiled(&[0], 5), 4, ByteOrder::Little);
        assert_eq!(node_header(&file, ROOT), (0, 2));
    }

    #[test]
    fn nodes_are_padded_with_slots_of_their_own_size() {
        // A root using 2 of 2 non-leaf slots, a full leaf, and a leaf using 1 of 2 leaf slots.
        let file = file_with_index(&tiled(&[0], 3), 2, ByteOrder::Little);
        let non_leaf_node = NODE_HEADER_SIZE + 2 * NON_LEAF_SLOT_SIZE;
        let leaf_node = NODE_HEADER_SIZE + 2 * LEAF_SLOT_SIZE;
        assert_eq!(file.len(), ROOT + non_leaf_node + 2 * leaf_node);
        assert!(file[file.len() - LEAF_SLOT_SIZE..].iter().all(|&byte| byte == 0));
    }

    #[test]
    fn non_leaf_bounds_enclose_their_children() {
        let entries = numbered(&[(0, 500, 600), (0, 100, 200), (0, 900, 950), (1, 0, 50)]);
        let file = file_with_index(&entries, 2, ByteOrder::Little);
        let slot_bounds = |slot: usize| -> Vec<u32> {
            let at = ROOT + NODE_HEADER_SIZE + slot * NON_LEAF_SLOT_SIZE;
            (0..4).map(|field| u32_at(&file, at + 4 * field)).collect()
        };
        assert_eq!(slot_bounds(0), [0, 100, 0, 600]);
        assert_eq!(slot_bounds(1), [0, 900, 1, 50]);
    }

    #[test]
    fn three_level_tree_lists_every_entry_in_order() {
        let entries = tiled(&[0, 1], 10);
        let file = file_with_index(&entries, 3, ByteOrder::Little);
        assert_eq!(node_header(&file, ROOT), (0, 3));
        assert_eq!(leaves_of(&file).unwrap(), entries);
    }

    #[test]
    fn three_level_tree_overlaps_match_brute_force() {
        let entries = tiled(&[0, 1], 10);
        let file = file_with_index(&entries, 3, ByteOrder::Little);
        assert_queries_match_brute_force(&file, &entries, 2);
    }

    #[test]
    fn deep_tree_of_fan_out_two_lists_and_finds_every_entry() {
        let entries = tiled(&[0, 1, 2], 11);
        let file = file_with_index(&entries, 2, ByteOrder::Little);
        assert_eq!(leaves_of(&file).unwrap(), entries);
        assert_queries_match_brute_force(&file, &entries, 3);
    }

    #[test]
    fn big_endian_index_reads_back() {
        let entries = tiled(&[0, 1], 10);
        let file = file_with_index(&entries, 3, ByteOrder::Big);
        assert_eq!(file[INDEX_START..INDEX_START + 4], RTREE_MAGIC.to_be_bytes());
        assert_eq!(leaves_of(&file).unwrap(), entries);
        assert_eq!(overlapping_in(&file, 1, 250, 450).unwrap(), brute_force(&entries, 1, 250, 450));
    }

    #[test]
    fn little_endian_index_reads_back() {
        let entries = tiled(&[0, 1], 10);
        let file = file_with_index(&entries, 3, ByteOrder::Little);
        assert_eq!(file[INDEX_START..INDEX_START + 4], RTREE_MAGIC.to_le_bytes());
        assert_eq!(leaves_of(&file).unwrap(), entries);
    }

    #[test]
    fn leaf_entries_keep_their_data_offsets_and_sizes() {
        let entry = LeafEntry {
            start_chrom_ix: 1,
            start_base: 2,
            end_chrom_ix: 3,
            end_base: 4,
            data_offset: 10,
            data_size: 90,
        };
        assert_eq!(leaves_of(&file_with_index(&[entry], 256, ByteOrder::Big)).unwrap(), [entry]);
    }

    #[test]
    fn overlap_is_half_open_at_both_ends() {
        let entries = numbered(&[(0, 100, 200)]);
        let file = file_with_index(&entries, 256, ByteOrder::Little);
        assert_eq!(overlapping_in(&file, 0, 200, 300).unwrap(), []);
        assert_eq!(overlapping_in(&file, 0, 0, 100).unwrap(), []);
        assert_eq!(overlapping_in(&file, 0, 199, 200).unwrap(), entries);
        assert_eq!(overlapping_in(&file, 0, 100, 101).unwrap(), entries);
    }

    #[test]
    fn entry_spanning_chromosomes_overlaps_queries_between_its_bounds() {
        let entry = LeafEntry { end_chrom_ix: 2, end_base: 50, ..numbered(&[(0, 100, 0)])[0] };
        let file = file_with_index(&[entry], 256, ByteOrder::Little);
        assert_eq!(overlapping_in(&file, 1, 0, 10).unwrap(), [entry]);
        assert_eq!(overlapping_in(&file, 0, 99, 101).unwrap(), [entry]);
        assert_eq!(overlapping_in(&file, 2, 49, 50).unwrap(), [entry]);
        assert_eq!(overlapping_in(&file, 0, 0, 100).unwrap(), []);
        assert_eq!(overlapping_in(&file, 2, 50, 60).unwrap(), []);
        assert_eq!(overlapping_in(&file, 3, 0, 10).unwrap(), []);
    }

    #[test]
    fn empty_query_matches_only_entries_strictly_containing_its_position() {
        let entries = numbered(&[(0, 100, 200)]);
        let file = file_with_index(&entries, 256, ByteOrder::Little);
        assert_eq!(overlapping_in(&file, 0, 150, 150).unwrap(), entries);
        assert_eq!(overlapping_in(&file, 0, 100, 100).unwrap(), []);
        assert_eq!(overlapping_in(&file, 0, 200, 200).unwrap(), []);
    }

    #[test]
    fn query_on_an_absent_chromosome_returns_nothing() {
        let file = file_with_index(&tiled(&[0, 2], 10), 3, ByteOrder::Little);
        assert_eq!(overlapping_in(&file, 1, 0, u32::MAX).unwrap(), []);
        assert_eq!(overlapping_in(&file, 7, 0, u32::MAX).unwrap(), []);
    }

    #[test]
    fn unordered_chrom_ids_still_return_exactly_the_overlapping_blocks() {
        // kent >= 2.10 writes data in input order with name-sorted IDs, so leaves can run 2, 0, 1.
        let entries = tiled(&[2, 0, 1], 10);
        let file = file_with_index(&entries, 2, ByteOrder::Little);
        assert_queries_match_brute_force(&file, &entries, 3);
    }

    #[test]
    fn interleaved_chrom_ids_still_return_exactly_the_overlapping_blocks() {
        let blocks: Vec<(u32, u32, u32)> =
            (0..30).map(|i| (i % 3, (i / 3) * 100, (i / 3) * 100 + 80)).collect();
        let entries = numbered(&blocks);
        let file = file_with_index(&entries, 2, ByteOrder::Little);
        assert_queries_match_brute_force(&file, &entries, 3);
    }

    #[test]
    fn query_spanning_multiple_rtree_leaves_returns_each_block_once() {
        let entries = tiled(&[0], 9);
        let file = file_with_index(&entries, 2, ByteOrder::Little);
        assert_eq!(overlapping_in(&file, 0, 0, 900).unwrap(), entries);
        assert_eq!(overlapping_in(&file, 0, 150, 650).unwrap(), entries[1..7]);
    }

    #[test]
    fn query_spanning_leaves_of_a_non_monotonic_tree_returns_each_block_once() {
        // chr1, chr10 and chr2 take name-sorted IDs 0, 1 and 2, but data come in numeric order:
        // chr1, chr2, chr10, so the leaves run through IDs 0, 2, 1.
        let entries = tiled(&[0, 2, 1], 5);
        let file = file_with_index(&entries, 2, ByteOrder::Little);
        for chrom in 0..3 {
            let got = overlapping_in(&file, chrom, 0, 500).unwrap();
            assert_eq!(got, brute_force(&entries, chrom, 0, 500));
            assert_eq!(got.len(), 5);
        }
    }

    #[test]
    fn reader_accepts_kent_short_padding_in_last_rtree_leaf() {
        // kent pads the last leaf with 24-byte slots, so with the index at the end of the file
        // that leaf is 8 bytes shorter than a full one.
        let entries = tiled(&[0], 3);
        let mut file = file_with_index(&entries, 2, ByteOrder::Little);
        file.truncate(file.len() - (LEAF_SLOT_SIZE - NON_LEAF_SLOT_SIZE));
        assert_eq!(leaves_of(&file).unwrap(), entries);
        assert_eq!(overlapping_in(&file, 0, 250, 260).unwrap(), entries[2..]);
    }

    #[test]
    fn each_node_is_fetched_with_one_read() {
        // Fan-out 2 over 4 entries: a root and two leaves.
        let file = file_with_index(&tiled(&[0], 4), 2, ByteOrder::Little);
        assert_eq!(reads_to_list_leaves(&file), 3);
    }

    #[test]
    fn node_with_more_slots_than_the_block_size_is_still_read_whole() {
        let entries = tiled(&[0], 4);
        let mut file = file_with_index(&entries, 4, ByteOrder::Little);
        file[INDEX_START + 4..INDEX_START + 8].copy_from_slice(&2u32.to_le_bytes());
        assert_eq!(leaves_of(&file).unwrap(), entries);
    }

    #[test]
    fn huge_block_size_still_reads_nodes_up_to_the_end_of_the_file() {
        let entries = tiled(&[0], 3);
        let mut file = file_with_index(&entries, 2, ByteOrder::Little);
        file[INDEX_START + 4..INDEX_START + 8].copy_from_slice(&u32::MAX.to_le_bytes());
        assert_eq!(leaves_of(&file).unwrap(), entries);
    }

    #[test]
    fn reader_tolerates_rtree_entries_with_zero_width_bounds() {
        let entries = numbered(&[(0, 100, 100), (0, 200, 300), (0, 500, 400), (1, 0, 0)]);
        let file = file_with_index(&entries, 2, ByteOrder::Little);
        assert_eq!(leaves_of(&file).unwrap(), entries);
        assert_queries_match_brute_force(&file, &entries, 2);
    }

    #[test]
    fn returned_block_past_the_end_of_the_file_is_corrupt_at_its_data_offset_field() {
        let mut entries = tiled(&[0], 2);
        entries[1].data_size = 0xFFFF_FFFF_8855_EE2F;
        let file = file_with_index(&entries, 256, ByteOrder::Little);
        let data_offset_field = (ROOT + NODE_HEADER_SIZE + LEAF_SLOT_SIZE + 16) as u64;
        assert!(matches!(
            leaves_of(&file),
            Err(Error::Corrupt { offset, .. }) if offset == data_offset_field
        ));
    }

    #[test]
    fn block_past_the_end_of_the_file_does_not_fail_queries_that_miss_it() {
        let mut entries = tiled(&[0], 2);
        entries[1].data_size = 0xFFFF_FFFF_8855_EE2F;
        let file = file_with_index(&entries, 256, ByteOrder::Little);
        assert_eq!(overlapping_in(&file, 0, 0, 100).unwrap(), entries[..1]);
    }

    #[test]
    fn block_whose_end_overflows_is_corrupt() {
        let mut entries = tiled(&[0], 1);
        entries[0].data_offset = u64::MAX;
        let file = file_with_index(&entries, 256, ByteOrder::Little);
        assert!(matches!(overlapping_in(&file, 0, 0, 100), Err(Error::Corrupt { .. })));
    }

    #[test]
    fn bad_magic_is_corrupt_at_the_index_offset() {
        let mut file = file_with_index(&tiled(&[0], 3), 256, ByteOrder::Little);
        file[INDEX_START] ^= 0xFF;
        let mut source = Cursor::new(&file);
        assert!(matches!(
            RTree::read(&mut source, INDEX_OFFSET, file.len() as u64),
            Err(Error::Corrupt { offset: INDEX_OFFSET, .. })
        ));
    }

    #[test]
    fn truncated_header_is_corrupt() {
        let file = file_with_index(&tiled(&[0], 3), 256, ByteOrder::Little);
        assert!(matches!(
            leaves_of(&file[..ROOT - 1]),
            Err(Error::Corrupt { offset: INDEX_OFFSET, .. })
        ));
    }

    #[test]
    fn truncated_node_is_corrupt_at_its_slots() {
        let file = file_with_index(&tiled(&[0], 3), 3, ByteOrder::Little);
        let slots_offset = (ROOT + NODE_HEADER_SIZE) as u64;
        assert!(matches!(
            leaves_of(&file[..ROOT + NODE_HEADER_SIZE + 2 * LEAF_SLOT_SIZE]),
            Err(Error::Corrupt { offset, .. }) if offset == slots_offset
        ));
    }

    #[test]
    fn items_without_a_root_node_are_corrupt() {
        let mut file = file_with_index(&tiled(&[0], 3), 256, ByteOrder::Little);
        file.truncate(ROOT);
        assert!(matches!(
            leaves_of(&file),
            Err(Error::Corrupt { offset, .. }) if offset == ROOT as u64
        ));
    }

    #[test]
    fn is_leaf_other_than_zero_or_one_is_corrupt_at_the_node() {
        let mut file = file_with_index(&tiled(&[0], 3), 256, ByteOrder::Little);
        file[ROOT] = 7;
        assert!(matches!(
            leaves_of(&file),
            Err(Error::Corrupt { offset, .. }) if offset == ROOT as u64
        ));
    }

    #[test]
    fn cyclic_child_offset_is_corrupt() {
        let mut file = file_with_index(&tiled(&[0], 4), 2, ByteOrder::Little);
        let field = root_child_offset_field(1);
        file[field..field + 8].copy_from_slice(&(ROOT as u64).to_le_bytes());
        assert!(matches!(
            leaves_of(&file),
            Err(Error::Corrupt { offset, .. }) if offset == ROOT as u64
        ));
    }

    #[test]
    fn shared_child_offset_is_corrupt() {
        let mut file = file_with_index(&tiled(&[0], 4), 2, ByteOrder::Little);
        let (first, second) = (root_child_offset_field(0), root_child_offset_field(1));
        file.copy_within(first..first + 8, second);
        assert!(matches!(overlapping_in(&file, 0, 0, 400), Err(Error::Corrupt { .. })));
    }

    #[test]
    fn child_offset_past_the_end_of_the_file_is_corrupt() {
        let mut file = file_with_index(&tiled(&[0], 4), 2, ByteOrder::Little);
        let field = root_child_offset_field(1);
        let past_end = file.len() as u64 + 100;
        file[field..field + 8].copy_from_slice(&past_end.to_le_bytes());
        assert!(matches!(
            leaves_of(&file),
            Err(Error::Corrupt { offset, .. }) if offset == past_end
        ));
    }

    #[test]
    fn index_at_the_depth_limit_is_read() {
        let block = tiled(&[0], 1)[0];
        assert_eq!(leaves_of(&chain_index(MAX_TREE_DEPTH, block)).unwrap(), [block]);
    }

    #[test]
    fn index_deeper_than_the_limit_is_corrupt() {
        let block = tiled(&[0], 1)[0];
        let file = chain_index(MAX_TREE_DEPTH + 1, block);
        assert!(matches!(overlapping_in(&file, 0, 0, 100), Err(Error::Corrupt { .. })));
    }

    #[test]
    #[should_panic(expected = "outside 2..=65535")]
    fn encode_rejects_a_block_size_below_two() {
        encode(&tiled(&[0], 3), INDEX_OFFSET, 1, 1, ByteOrder::Little);
    }

    #[test]
    #[should_panic(expected = "outside 2..=65535")]
    fn encode_rejects_a_block_size_above_a_u16_count() {
        encode(&tiled(&[0], 3), INDEX_OFFSET, 65_536, 1, ByteOrder::Little);
    }

    proptest! {
        #[test]
        fn overlaps_match_brute_force_for_any_entries(
            blocks in prop::collection::vec((0u32..4, 0u32..1000, 0u32..300), 0..120),
            block_size in 2u32..6,
            chrom in 0u32..5,
            start in 0u32..1300,
            len in 0u32..400,
        ) {
            let blocks: Vec<(u32, u32, u32)> =
                blocks.into_iter().map(|(chrom, start, len)| (chrom, start, start + len)).collect();
            let entries = numbered(&blocks);
            let file = file_with_index(&entries, block_size, ByteOrder::Little);
            prop_assert_eq!(leaves_of(&file).unwrap(), entries.clone());
            prop_assert_eq!(
                overlapping_in(&file, chrom, start, start + len).unwrap(),
                brute_force(&entries, chrom, start, start + len)
            );
        }

        #[test]
        fn decoding_mutated_indexes_never_panics(
            per_chrom in 0u32..20,
            block_size in 2u32..5,
            mutations in prop::collection::vec((any::<prop::sample::Index>(), any::<u8>()), 1..8),
            keep in any::<prop::sample::Index>(),
            chrom in 0u32..3,
            start in 0u32..2000,
            len in 0u32..2000,
        ) {
            let mut file = file_with_index(&tiled(&[0, 1], per_chrom), block_size, ByteOrder::Big);
            let index_len = file.len() - INDEX_START;
            for (index, byte) in mutations {
                file[INDEX_START + index.index(index_len)] = byte;
            }
            file.truncate(INDEX_START + keep.index(index_len) + 1);
            let _ = leaves_of(&file);
            let _ = overlapping_in(&file, chrom, start, start + len);
        }
    }
}
