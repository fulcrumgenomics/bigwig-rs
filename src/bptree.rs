//! The chromosome B+ tree: writing a byte-sorted tree (a single leaf when possible) and reading every
//! leaf into the chromosome list.

use std::collections::HashSet;
use std::io::{Read, Seek};

use crate::bbi::{
    BPTREE_MAGIC, ByteOrder, ByteReader, ChromInfo, MAX_TREE_DEPTH, read_exact_at,
    read_exact_at_into,
};
use crate::error::{Error, Result};

/// Size of the tree header.
const HEADER_SIZE: usize = 32;
/// Size of a node header: `isLeaf` (u8), reserved (u8), `count` (u16).
const NODE_HEADER_SIZE: usize = 4;
/// Size of a leaf slot's value (`chromId` u32 then `chromSize` u32). A non-leaf slot's child offset
/// (u64) is the same size, so leaf and non-leaf nodes have the same layout and size.
const VALUE_SIZE: u32 = 8;
/// Most chromosomes in one node, so up to this many fit a single leaf. It stays below 32,768
/// because some Java readers treat the u16 `count` as signed.
const MAX_BLOCK_SIZE: usize = 32_767;

/// One chromosome to store in the tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ChromEntry<'a> {
    /// Sequence name: non-empty, free of NUL bytes and unique within the tree.
    pub(crate) name: &'a str,
    /// The `chromId` to write; a placeholder when the writer patches it later.
    pub(crate) chrom_id: u32,
    /// Length in bases.
    pub(crate) length: u32,
}

/// A tree encoded for a particular file offset.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct EncodedChromTree {
    /// The tree, header first, to be written at the offset it was encoded for.
    pub(crate) bytes: Vec<u8>,
    /// For each input entry, in input order, the absolute file offset of its 4-byte `chromId`
    /// field, so IDs can be overwritten in place once they are known.
    pub(crate) chrom_id_offsets: Vec<u64>,
}

/// Encodes `entries` as a B+ tree to be written at absolute file offset `offset`.
///
/// Keys are sorted bytewise and NUL-padded to the longest name, and every node is padded to
/// `blockSize` slots with zeros. Up to 32,767 chromosomes fit a single leaf with `blockSize` equal
/// to the count; above that the tree has several levels with fan-out 32,767. An empty list gives
/// just the header, as kent writes.
///
/// # Panics
///
/// If a name is longer than `u32::MAX` bytes.
pub(crate) fn encode(
    entries: &[ChromEntry<'_>],
    offset: u64,
    order: ByteOrder,
) -> EncodedChromTree {
    encode_with_max_block_size(entries, offset, order, MAX_BLOCK_SIZE)
}

/// As [`encode`], with `blockSize = min(max_block_size, entries.len())`. kent's tools use 256.
///
/// # Panics
///
/// If `max_block_size` is outside 2..=65,535 (smaller fan-outs never reduce to a root, and a
/// node's `count` is a u16), or a name is longer than `u32::MAX` bytes.
fn encode_with_max_block_size(
    entries: &[ChromEntry<'_>],
    offset: u64,
    order: ByteOrder,
    max_block_size: usize,
) -> EncodedChromTree {
    assert!(
        (2..=usize::from(u16::MAX)).contains(&max_block_size),
        "B+ tree block size {max_block_size} is outside 2..=65535"
    );
    let block_size = entries.len().min(max_block_size);
    let key_size = entries.iter().map(|entry| entry.name.len()).max().unwrap_or(0);
    let node_size = NODE_HEADER_SIZE + block_size * (key_size + VALUE_SIZE as usize);

    let mut bytes = Vec::new();
    order.put_u32(&mut bytes, BPTREE_MAGIC);
    order.put_u32(&mut bytes, u32::try_from(block_size).expect("block size is at most 65,535"));
    order.put_u32(&mut bytes, u32::try_from(key_size).expect("name longer than u32::MAX bytes"));
    order.put_u32(&mut bytes, VALUE_SIZE);
    order.put_u64(&mut bytes, entries.len() as u64);
    order.put_u64(&mut bytes, 0);
    if entries.is_empty() {
        return EncodedChromTree { bytes, chrom_id_offsets: Vec::new() };
    }

    // Byte order of NUL-padded keys equals byte order of the names, since names have no NULs.
    let mut sorted: Vec<usize> = (0..entries.len()).collect();
    sorted.sort_unstable_by_key(|&index| entries[index].name.as_bytes());
    debug_assert!(
        sorted.windows(2).all(|pair| entries[pair[0]].name != entries[pair[1]].name),
        "duplicate chromosome names"
    );

    // Node count of each level, leaves (level 0) first, up to the single root.
    let mut level_sizes = vec![entries.len().div_ceil(block_size)];
    while let Some(&nodes) = level_sizes.last().filter(|&&nodes| nodes > 1) {
        level_sizes.push(nodes.div_ceil(block_size));
    }
    // Absolute offset of each level's first node; levels are written root first.
    let mut level_offsets = vec![0; level_sizes.len()];
    let mut next_offset = offset + HEADER_SIZE as u64;
    for level in (0..level_sizes.len()).rev() {
        level_offsets[level] = next_offset;
        next_offset += (level_sizes[level] * node_size) as u64;
    }
    bytes.reserve(level_sizes.iter().sum::<usize>() * node_size);

    // A node on level L covers block_size^(L+1) chromosomes, and its i-th slot points at the
    // node with the same global index on the level below, keyed by that node's first name.
    for level in (1..level_sizes.len()).rev() {
        let level_exponent = u32::try_from(level).expect("a tree has far fewer than 2^32 levels");
        let names_per_child = block_size.pow(level_exponent);
        let children = level_sizes[level - 1];
        for first_child in (0..children).step_by(block_size) {
            let last_child = children.min(first_child + block_size);
            put_node_header(&mut bytes, order, false, last_child - first_child);
            for child in first_child..last_child {
                put_key(&mut bytes, entries[sorted[child * names_per_child]].name, key_size);
                order.put_u64(&mut bytes, level_offsets[level - 1] + (child * node_size) as u64);
            }
            put_padding(&mut bytes, block_size - (last_child - first_child), key_size);
        }
    }

    let mut chrom_id_offsets = vec![0; entries.len()];
    for leaf in sorted.chunks(block_size) {
        put_node_header(&mut bytes, order, true, leaf.len());
        for &index in leaf {
            let entry = &entries[index];
            put_key(&mut bytes, entry.name, key_size);
            chrom_id_offsets[index] = offset + bytes.len() as u64;
            order.put_u32(&mut bytes, entry.chrom_id);
            order.put_u32(&mut bytes, entry.length);
        }
        put_padding(&mut bytes, block_size - leaf.len(), key_size);
    }
    debug_assert_eq!(offset + bytes.len() as u64, next_offset);
    EncodedChromTree { bytes, chrom_id_offsets }
}

/// Appends a node header: `isLeaf`, a zero reserved byte and the slot count.
fn put_node_header(bytes: &mut Vec<u8>, order: ByteOrder, is_leaf: bool, count: usize) {
    bytes.push(u8::from(is_leaf));
    bytes.push(0);
    order.put_u16(bytes, u16::try_from(count).expect("a node holds at most 65,535 slots"));
}

/// Appends `name` NUL-padded to `key_size` bytes.
fn put_key(bytes: &mut Vec<u8>, name: &str, key_size: usize) {
    bytes.extend_from_slice(name.as_bytes());
    bytes.resize(bytes.len() + key_size - name.len(), 0);
}

/// Appends `slots` unused, zero-filled slots.
fn put_padding(bytes: &mut Vec<u8>, slots: usize, key_size: usize) {
    bytes.resize(bytes.len() + slots * (key_size + VALUE_SIZE as usize), 0);
}

/// Reads the B+ tree at absolute `offset` in a source of `file_len` bytes and returns its
/// chromosomes ordered by `chromId`, so a chromosome's index in the result is its `chromId`.
///
/// Every node is visited rather than searched, so keys need not be sorted. The byte order comes
/// from the tree's own magic. Keys are trimmed at their first NUL and must be UTF-8. An
/// `itemCount` of 0 gives an empty list without reading any node, since kent writes none.
///
/// Fails with [`Error::Corrupt`] at the offending offset on a bad magic, a `valSize` other than 8,
/// an `itemCount` whose slots can't fit in the file, a node or child offset past the end of the
/// file, an `isLeaf` other than 0 or 1, a node reached twice (a cyclic or shared child offset), a
/// tree deeper than [`MAX_TREE_DEPTH`], a key that is not UTF-8, or `chromId`s that are not
/// exactly `0..itemCount`.
pub(crate) fn read<R: Read + Seek>(
    source: &mut R,
    offset: u64,
    file_len: u64,
) -> Result<Vec<ChromInfo>> {
    let header = read_exact_at(source, offset, HEADER_SIZE, file_len)?;
    let found = ByteReader::new(&header, ByteOrder::Little, offset).u32()?;
    let byte_order = ByteOrder::detect(found, BPTREE_MAGIC)
        .ok_or_else(|| Error::corrupt(offset, format!("bad B+ tree magic {found:#010x}")))?;
    let mut reader = ByteReader::new(&header, byte_order, offset);
    reader.skip(8)?; // magic, blockSize
    let key_size = reader.u32()?;
    let value_size_offset = reader.file_offset();
    let value_size = reader.u32()?;
    let item_count_offset = reader.file_offset();
    let item_count = reader.u64()?;

    if value_size != VALUE_SIZE {
        return Err(Error::corrupt(
            value_size_offset,
            format!("B+ tree valSize is {value_size}, but bigWig needs {VALUE_SIZE}"),
        ));
    }
    if item_count == 0 {
        return Ok(Vec::new());
    }
    // Every chromosome has a leaf slot of its own, so slots that can't fit in the file mean a
    // corrupt count; this bounds everything allocated below by the file's size.
    let slot_size = u64::from(key_size) + u64::from(VALUE_SIZE);
    if item_count.checked_mul(slot_size).is_none_or(|bytes| bytes > file_len) {
        return Err(Error::corrupt(
            item_count_offset,
            format!(
                "B+ tree itemCount {item_count} of {slot_size}-byte slots can't fit in the \
                 {file_len}-byte file"
            ),
        ));
    }
    let slot_size = usize::try_from(slot_size)
        .map_err(|_| Error::corrupt(offset, format!("B+ tree keySize {key_size} is too large")))?;
    let key_len = slot_size - VALUE_SIZE as usize;

    // (chromId, file offset of the chromId, chromosome) for every leaf slot, in visit order.
    let mut found: Vec<(u32, u64, ChromInfo)> = Vec::new();
    let mut pending = vec![(offset + HEADER_SIZE as u64, 0)];
    let mut visited = HashSet::new();
    let mut buf = Vec::new();
    while let Some((node_offset, depth)) = pending.pop() {
        if depth >= MAX_TREE_DEPTH {
            return Err(Error::corrupt(
                node_offset,
                format!("B+ tree is deeper than {MAX_TREE_DEPTH} levels"),
            ));
        }
        if !visited.insert(node_offset) {
            return Err(Error::corrupt(
                node_offset,
                "B+ tree node reached twice (cyclic or shared child offset)",
            ));
        }
        read_exact_at_into(source, node_offset, NODE_HEADER_SIZE, file_len, &mut buf)?;
        let mut node = ByteReader::new(&buf, byte_order, node_offset);
        let is_leaf = node.u8()?;
        node.skip(1)?;
        let count = node.u16()?;
        if is_leaf > 1 {
            return Err(Error::corrupt(node_offset, format!("B+ tree isLeaf is {is_leaf}")));
        }

        // Only `count` slots are read: a node's unused slots are padding, which isn't needed.
        let slots_offset = node_offset + NODE_HEADER_SIZE as u64;
        let slots_len = usize::from(count).checked_mul(slot_size).ok_or_else(|| {
            Error::corrupt(node_offset, format!("{count} B+ tree slots are too large to read"))
        })?;
        read_exact_at_into(source, slots_offset, slots_len, file_len, &mut buf)?;
        let mut slots = ByteReader::new(&buf, byte_order, slots_offset);
        for _ in 0..count {
            if is_leaf == 0 {
                slots.skip(key_len)?;
                pending.push((slots.u64()?, depth + 1));
                continue;
            }
            let key_offset = slots.file_offset();
            let name = decode_key(slots.bytes(key_len)?, key_offset)?;
            let chrom_id_offset = slots.file_offset();
            let chrom_id = slots.u32()?;
            let length = slots.u32()?;
            if u64::from(chrom_id) >= item_count {
                return Err(Error::corrupt(
                    chrom_id_offset,
                    format!("chromId {chrom_id} is not below itemCount {item_count}"),
                ));
            }
            if found.len() as u64 == item_count {
                return Err(Error::corrupt(
                    key_offset,
                    format!("B+ tree leaves hold more than itemCount {item_count} chromosomes"),
                ));
            }
            found.push((chrom_id, chrom_id_offset, ChromInfo::new(name, length)));
        }
    }

    order_by_chrom_id(found, item_count, item_count_offset)
}

/// Orders the leaf slots `found` (chromId, offset of the chromId, chromosome) by `chromId`,
/// checking that their IDs are exactly `0..item_count`. All IDs are already known to be below
/// `item_count`, and there are at most `item_count` of them.
fn order_by_chrom_id(
    mut found: Vec<(u32, u64, ChromInfo)>,
    item_count: u64,
    item_count_offset: u64,
) -> Result<Vec<ChromInfo>> {
    if (found.len() as u64) < item_count {
        return Err(Error::corrupt(
            item_count_offset,
            format!(
                "B+ tree itemCount is {item_count} but its leaves hold {} chromosomes",
                found.len()
            ),
        ));
    }
    // With exactly itemCount IDs, all below itemCount, no duplicates means exactly 0..itemCount.
    // The sort is stable, so of two equal IDs the second is the one visited later.
    found.sort_by_key(|&(chrom_id, ..)| chrom_id);
    if let Some(pair) = found.windows(2).find(|pair| pair[0].0 == pair[1].0) {
        return Err(Error::corrupt(pair[1].1, format!("duplicate chromId {}", pair[1].0)));
    }
    Ok(found.into_iter().map(|(.., chrom)| chrom).collect())
}

/// Decodes a key: the bytes before its first NUL, which must be UTF-8.
fn decode_key(key: &[u8], key_offset: u64) -> Result<String> {
    let name = key.iter().position(|&byte| byte == 0).map_or(key, |end| &key[..end]);
    match std::str::from_utf8(name) {
        Ok(name) => Ok(name.to_owned()),
        Err(error) => {
            Err(Error::corrupt(key_offset, format!("B+ tree key is not valid UTF-8: {error}")))
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::io::Cursor;

    use proptest::prelude::*;

    use super::*;

    /// kent's `bptFileFind` (`kent:lib/bPlusTree.c:119-165`) over an encoded tree whose header is
    /// `tree[0]`, at absolute offset `base`: it descends non-leaf nodes by
    /// `memcmp(key, slotKey) < 0`, then looks for the NUL-padded key in one leaf. Each leaf is
    /// indexed on its first visit, which answers as kent's linear scan does while keeping tens of
    /// thousands of lookups fast.
    struct KentLookup<'a> {
        /// The tree, header first.
        tree: &'a [u8],
        /// Absolute offset of `tree[0]`.
        base: u64,
        /// Byte order of the tree.
        order: ByteOrder,
        /// The header's `keySize`.
        key_size: usize,
        /// Leaf contents by node offset: padded key to (chromId, chromSize).
        leaves: HashMap<u64, HashMap<Vec<u8>, (u32, u32)>>,
    }

    impl<'a> KentLookup<'a> {
        /// Prepares lookups in `tree`, which sits at absolute offset `base`.
        fn new(tree: &'a [u8], base: u64, order: ByteOrder) -> Self {
            let key_size = ByteReader::new(&tree[8..12], order, 8).u32().unwrap() as usize;
            Self { tree, base, order, key_size, leaves: HashMap::new() }
        }

        /// Returns (chromId, chromSize) for `name`, or `None` where kent finds nothing.
        fn find(&mut self, name: &str) -> Option<(u32, u32)> {
            if name.len() > self.key_size {
                return None;
            }
            let key_size = self.key_size;
            let mut key = name.as_bytes().to_vec();
            key.resize(key_size, 0);
            let mut node_offset = self.base + HEADER_SIZE as u64;
            loop {
                let mut node = self.node(node_offset);
                let is_leaf = node.u8().unwrap();
                node.skip(1).unwrap();
                let count = node.u16().unwrap();
                if is_leaf == 1 {
                    let leaf = self.leaves.entry(node_offset).or_insert_with(|| {
                        (0..count)
                            .map(|_| {
                                let slot_key = node.bytes(key_size).unwrap().to_vec();
                                (slot_key, (node.u32().unwrap(), node.u32().unwrap()))
                            })
                            .collect()
                    });
                    return leaf.get(&key).copied();
                }
                node.skip(key_size).unwrap();
                let mut child_offset = node.u64().unwrap();
                for _ in 1..count {
                    if key.as_slice() < node.bytes(key_size).unwrap() {
                        break;
                    }
                    child_offset = node.u64().unwrap();
                }
                node_offset = child_offset;
            }
        }

        /// A reader positioned at the node at absolute `offset`.
        fn node(&self, offset: u64) -> ByteReader<'a> {
            let tree: &'a [u8] = self.tree;
            ByteReader::new(
                &tree[usize::try_from(offset - self.base).unwrap()..],
                self.order,
                offset,
            )
        }
    }

    /// `count` distinct names in numeric rather than byte order, many of them prefixes of others:
    /// `contig1`, `contig2`, ..., `contig10`, ...
    fn contig_names(count: usize) -> Vec<String> {
        (1..=count).map(|n| format!("contig{n}")).collect()
    }

    /// Entries for `names` with `chromId` equal to the name's position and a distinct length.
    fn entries(names: &[String]) -> Vec<ChromEntry<'_>> {
        names
            .iter()
            .enumerate()
            .map(|(position, name)| {
                let chrom_id = u32::try_from(position).unwrap();
                ChromEntry { name, chrom_id, length: 1000 + chrom_id }
            })
            .collect()
    }

    /// The chromosomes that reading the tree of `entries(names)` must return.
    fn expected_chroms(names: &[String]) -> Vec<ChromInfo> {
        entries(names).iter().map(|entry| ChromInfo::new(entry.name, entry.length)).collect()
    }

    /// The tree of `entries(names)` at offset 0 in little-endian order, kent-style with
    /// `blockSize = min(256, count)`.
    fn kent_style_tree(names: &[String]) -> Vec<u8> {
        encode_with_max_block_size(&entries(names), 0, ByteOrder::Little, 256).bytes
    }

    /// The default-shaped little-endian tree of `entries(names)` at offset 0.
    fn our_tree(names: &[String]) -> Vec<u8> {
        encode(&entries(names), 0, ByteOrder::Little).bytes
    }

    /// Reads the tree at `offset` in `file`.
    fn read_file(file: &[u8], offset: u64) -> Result<Vec<ChromInfo>> {
        read(&mut Cursor::new(file), offset, file.len() as u64)
    }

    /// The header's `blockSize` of a little-endian tree starting at `tree[0]`.
    fn block_size(tree: &[u8]) -> u32 {
        u32::from_le_bytes(tree[4..8].try_into().unwrap())
    }

    /// Whether the root of a tree starting at `tree[0]` is a leaf.
    fn root_is_leaf(tree: &[u8]) -> bool {
        tree[HEADER_SIZE] == 1
    }

    /// The slot count of the little-endian root node of a tree starting at `tree[0]`.
    fn root_count(tree: &[u8]) -> u16 {
        u16::from_le_bytes(tree[HEADER_SIZE + 2..HEADER_SIZE + 4].try_into().unwrap())
    }

    /// Asserts that kent's lookup finds every one of `names` in `tree` with its ID and length.
    fn assert_kent_finds_every_name(tree: &[u8], names: &[String]) {
        let mut lookup = KentLookup::new(tree, 0, ByteOrder::Little);
        for entry in entries(names) {
            assert_eq!(lookup.find(entry.name), Some((entry.chrom_id, entry.length)), "{entry:?}");
        }
    }

    /// Names absent from `contig_names(_)` trees, several sharing a prefix with names in them.
    const ABSENT_NAMES: [&str; 8] =
        ["", "c", "contig", "contig0", "contig1_", "contig99999", "Contig1", "zzz"];

    /// Asserts that kent's lookup finds none of [`ABSENT_NAMES`] in `tree`.
    fn assert_kent_misses_absent_names(tree: &[u8]) {
        let mut lookup = KentLookup::new(tree, 0, ByteOrder::Little);
        for name in ABSENT_NAMES {
            assert_eq!(lookup.find(name), None, "{name:?}");
        }
    }

    /// A little-endian tree of `levels` levels: single-child non-leaf nodes chained down to one
    /// leaf holding chromosome `a`.
    fn chain_tree(levels: usize) -> Vec<u8> {
        let order = ByteOrder::Little;
        let node_size = (NODE_HEADER_SIZE + 1 + VALUE_SIZE as usize) as u64;
        let mut bytes = Vec::new();
        for field in [BPTREE_MAGIC, 1, 1, VALUE_SIZE] {
            order.put_u32(&mut bytes, field);
        }
        order.put_u64(&mut bytes, 1);
        order.put_u64(&mut bytes, 0);
        for level in 1..levels as u64 {
            put_node_header(&mut bytes, order, false, 1);
            bytes.push(b'a');
            order.put_u64(&mut bytes, HEADER_SIZE as u64 + level * node_size);
        }
        put_node_header(&mut bytes, order, true, 1);
        bytes.push(b'a');
        order.put_u32(&mut bytes, 0);
        order.put_u32(&mut bytes, 100);
        bytes
    }

    /// Offset within `kent_style_tree(&contig_names(300))` of the root's `slot`-th child offset;
    /// its keys are 9 bytes.
    fn root_child_offset_field(slot: usize) -> usize {
        HEADER_SIZE + NODE_HEADER_SIZE + slot * (9 + VALUE_SIZE as usize) + 9
    }

    /// Owned names from string literals.
    fn owned(names: &[&str]) -> Vec<String> {
        names.iter().map(|&name| name.to_owned()).collect()
    }

    #[test]
    fn one_name_is_a_single_leaf() {
        let tree = our_tree(&contig_names(1));
        assert_eq!(block_size(&tree), 1);
        assert!(root_is_leaf(&tree));
        assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&contig_names(1)));
    }

    #[test]
    fn three_thousand_names_are_a_single_leaf() {
        let tree = our_tree(&contig_names(3000));
        assert_eq!(block_size(&tree), 3000);
        assert!(root_is_leaf(&tree));
        assert_eq!(root_count(&tree), 3000);
    }

    #[test]
    fn names_up_to_the_single_leaf_limit_are_a_single_leaf() {
        let tree = our_tree(&contig_names(32_767));
        assert_eq!(block_size(&tree), 32_767);
        assert!(root_is_leaf(&tree));
        assert_eq!(root_count(&tree), 32_767);
    }

    #[test]
    fn one_name_past_the_single_leaf_limit_is_multi_level() {
        let tree = our_tree(&contig_names(32_768));
        assert_eq!(block_size(&tree), 32_767);
        assert!(!root_is_leaf(&tree));
        assert_eq!(root_count(&tree), 2);
    }

    #[test]
    fn forty_thousand_names_are_multi_level() {
        let tree = our_tree(&contig_names(40_000));
        assert_eq!(block_size(&tree), 32_767);
        assert!(!root_is_leaf(&tree));
        assert_eq!(root_count(&tree), 2);
    }

    #[test]
    fn kent_style_tree_of_three_thousand_names_has_two_levels() {
        let tree = kent_style_tree(&contig_names(3000));
        assert_eq!(block_size(&tree), 256);
        assert!(!root_is_leaf(&tree));
        assert_eq!(root_count(&tree), 12);
    }

    #[test]
    fn kent_style_block_size_is_the_count_below_256() {
        assert_eq!(block_size(&kent_style_tree(&contig_names(100))), 100);
    }

    #[test]
    fn single_leaf_of_three_thousand_reads_back_by_chrom_id() {
        let names = contig_names(3000);
        assert_eq!(read_file(&our_tree(&names), 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn single_leaf_at_the_limit_reads_back_by_chrom_id() {
        let names = contig_names(32_767);
        assert_eq!(read_file(&our_tree(&names), 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn multi_level_tree_of_32768_reads_back_by_chrom_id() {
        let names = contig_names(32_768);
        assert_eq!(read_file(&our_tree(&names), 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn multi_level_tree_of_40000_reads_back_by_chrom_id() {
        let names = contig_names(40_000);
        assert_eq!(read_file(&our_tree(&names), 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn kent_style_tree_of_three_thousand_reads_back_by_chrom_id() {
        let names = contig_names(3000);
        assert_eq!(read_file(&kent_style_tree(&names), 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn read_chromosomes_omit_absent_names() {
        let chroms = read_file(&kent_style_tree(&contig_names(3000)), 0).unwrap();
        let names: HashSet<&str> = chroms.iter().map(|chrom| chrom.name.as_str()).collect();
        for name in ABSENT_NAMES {
            assert!(!names.contains(name), "{name:?}");
        }
    }

    #[test]
    fn kent_lookup_finds_every_name_in_a_single_leaf_of_three_thousand() {
        let names = contig_names(3000);
        assert_kent_finds_every_name(&our_tree(&names), &names);
    }

    #[test]
    fn kent_lookup_finds_every_name_in_a_single_leaf_at_the_limit() {
        let names = contig_names(32_767);
        assert_kent_finds_every_name(&our_tree(&names), &names);
    }

    #[test]
    fn kent_lookup_finds_every_name_in_a_multi_level_tree_of_32768() {
        let names = contig_names(32_768);
        assert_kent_finds_every_name(&our_tree(&names), &names);
    }

    #[test]
    fn kent_lookup_finds_every_name_in_a_multi_level_tree_of_40000() {
        let names = contig_names(40_000);
        assert_kent_finds_every_name(&our_tree(&names), &names);
    }

    #[test]
    fn kent_lookup_finds_every_name_in_a_kent_style_tree_of_three_thousand() {
        let names = contig_names(3000);
        assert_kent_finds_every_name(&kent_style_tree(&names), &names);
    }

    #[test]
    fn three_level_kent_style_tree_reads_back_and_is_searchable() {
        let names = contig_names(70_000);
        let tree = kent_style_tree(&names);
        assert_kent_finds_every_name(&tree, &names);
        assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn kent_lookup_misses_absent_names_in_a_single_leaf() {
        assert_kent_misses_absent_names(&our_tree(&contig_names(3000)));
    }

    #[test]
    fn kent_lookup_misses_absent_names_in_a_multi_level_tree() {
        assert_kent_misses_absent_names(&our_tree(&contig_names(40_000)));
    }

    #[test]
    fn kent_lookup_misses_absent_names_in_a_kent_style_tree() {
        assert_kent_misses_absent_names(&kent_style_tree(&contig_names(3000)));
    }

    #[test]
    fn prefix_names_are_stored_in_byte_order() {
        let tree = our_tree(&owned(&["chr2", "chr1_alt", "chr10", "chr1"]));
        // keySize is 8, so each slot is 16 bytes and the first starts after the node header.
        let keys: Vec<&[u8]> = (0..4)
            .map(|slot| {
                let start = HEADER_SIZE + NODE_HEADER_SIZE + slot * 16;
                &tree[start..start + 8]
            })
            .collect();
        assert_eq!(keys, [b"chr1\0\0\0\0", b"chr10\0\0\0", b"chr1_alt", b"chr2\0\0\0\0"]);
    }

    #[test]
    fn prefix_names_read_back_and_are_found_by_kent_lookup() {
        let names = owned(&["chr2", "chr1_alt", "chr10", "chr1"]);
        let tree = our_tree(&names);
        assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&names));
        assert_kent_finds_every_name(&tree, &names);
        let mut lookup = KentLookup::new(&tree, 0, ByteOrder::Little);
        assert_eq!(lookup.find("chr1_"), None);
        assert_eq!(lookup.find("chr100"), None);
    }

    #[test]
    fn prefix_names_in_different_leaves_are_found_by_kent_lookup() {
        // A fan-out of 2 spreads the prefix names over several leaves, three levels deep.
        let names = owned(&["chr1", "chr10", "chr1_alt", "chr2", "chr20", "chr2_alt", "chrX"]);
        let tree = encode_with_max_block_size(&entries(&names), 0, ByteOrder::Little, 2).bytes;
        assert_kent_finds_every_name(&tree, &names);
        assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn keys_are_nul_padded_to_the_longest_name() {
        let names = owned(&["chrUn_long", "chr1"]);
        let tree = our_tree(&names);
        assert_eq!(u32::from_le_bytes(tree[8..12].try_into().unwrap()), 10);
        let first_key = HEADER_SIZE + NODE_HEADER_SIZE;
        assert_eq!(&tree[first_key..first_key + 10], b"chr1\0\0\0\0\0\0");
        assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn val_size_is_eight() {
        let tree = our_tree(&contig_names(3));
        assert_eq!(u32::from_le_bytes(tree[12..16].try_into().unwrap()), 8);
    }

    #[test]
    fn unused_slots_are_zero_filled() {
        // 300 names in a kent-style tree: a root using 2 of 256 slots, a full leaf, and a leaf
        // using 44 of 256 slots. Keys are 9 bytes.
        let tree = kent_style_tree(&contig_names(300));
        let slot_size = 9 + VALUE_SIZE as usize;
        let node_size = NODE_HEADER_SIZE + 256 * slot_size;
        assert_eq!(tree.len(), HEADER_SIZE + 3 * node_size);
        let root_unused = HEADER_SIZE + NODE_HEADER_SIZE + 2 * slot_size..HEADER_SIZE + node_size;
        let last_leaf_unused = HEADER_SIZE + 2 * node_size + NODE_HEADER_SIZE + 44 * slot_size..;
        assert!(tree[root_unused].iter().all(|&byte| byte == 0));
        assert!(tree[last_leaf_unused].iter().all(|&byte| byte == 0));
    }

    #[test]
    fn chrom_id_offsets_are_absolute_and_point_at_each_chrom_id() {
        let names = contig_names(300);
        let encoded = encode_with_max_block_size(&entries(&names), 1000, ByteOrder::Big, 256);
        for (entry, &field) in entries(&names).iter().zip(&encoded.chrom_id_offsets) {
            let at = usize::try_from(field - 1000).unwrap();
            let stored = u32::from_be_bytes(encoded.bytes[at..at + 4].try_into().unwrap());
            assert_eq!(stored, entry.chrom_id);
        }
    }

    #[test]
    fn chrom_ids_patched_through_their_offsets_are_read_back() {
        let names = contig_names(300);
        let encoded = encode_with_max_block_size(&entries(&names), 1000, ByteOrder::Little, 256);
        let mut file = vec![0xAA; 1000];
        file.extend_from_slice(&encoded.bytes);
        // Reverse the IDs, so the last name becomes chromId 0.
        for (position, &field) in encoded.chrom_id_offsets.iter().enumerate() {
            let new_id = u32::try_from(names.len() - 1 - position).unwrap();
            let at = usize::try_from(field).unwrap();
            file[at..at + 4].copy_from_slice(&ByteOrder::Little.u32_bytes(new_id));
        }
        let mut expected = expected_chroms(&names);
        expected.reverse();
        assert_eq!(read_file(&file, 1000).unwrap(), expected);
    }

    #[test]
    fn big_endian_tree_reads_back() {
        let names = contig_names(3000);
        let tree = encode_with_max_block_size(&entries(&names), 0, ByteOrder::Big, 256).bytes;
        assert_eq!(tree[..4], BPTREE_MAGIC.to_be_bytes());
        assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn little_endian_tree_reads_back_at_a_nonzero_offset() {
        let names = contig_names(3000);
        let mut file = vec![0xAA; 77];
        file.extend(encode(&entries(&names), 77, ByteOrder::Little).bytes);
        assert_eq!(file[77..81], BPTREE_MAGIC.to_le_bytes());
        assert_eq!(read_file(&file, 77).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn empty_tree_is_just_a_header_with_block_size_zero() {
        let encoded = encode(&[], 0, ByteOrder::Little);
        assert_eq!(encoded.bytes.len(), HEADER_SIZE);
        assert_eq!(block_size(&encoded.bytes), 0);
        assert!(encoded.chrom_id_offsets.is_empty());
    }

    #[test]
    fn empty_tree_without_a_root_reads_as_no_chromosomes() {
        assert_eq!(read_file(&encode(&[], 0, ByteOrder::Little).bytes, 0).unwrap(), []);
    }

    #[test]
    fn item_count_zero_ignores_whatever_follows_the_header() {
        let mut file = encode(&[], 0, ByteOrder::Big).bytes;
        file.extend_from_slice(&[0xFF; 16]);
        assert_eq!(read_file(&file, 0).unwrap(), []);
    }

    #[test]
    fn unsorted_keys_are_read() {
        let names = contig_names(3);
        let mut tree = our_tree(&names);
        // Swap the first and last leaf slots (7-byte keys), leaving the keys out of order.
        let slot =
            |index: usize| HEADER_SIZE + NODE_HEADER_SIZE + index * (7 + VALUE_SIZE as usize);
        let first = tree[slot(0)..slot(1)].to_vec();
        tree.copy_within(slot(2)..slot(3), slot(0));
        tree[slot(2)..slot(3)].copy_from_slice(&first);
        assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&names));
    }

    #[test]
    fn chromosomes_are_ordered_by_chrom_id_not_by_name() {
        let entries = [
            ChromEntry { name: "b", chrom_id: 0, length: 1 },
            ChromEntry { name: "c", chrom_id: 1, length: 2 },
            ChromEntry { name: "a", chrom_id: 2, length: 3 },
        ];
        let tree = encode(&entries, 0, ByteOrder::Little).bytes;
        let names: Vec<String> =
            read_file(&tree, 0).unwrap().into_iter().map(|chrom| chrom.name).collect();
        assert_eq!(names, ["b", "c", "a"]);
    }

    #[test]
    fn keys_are_trimmed_at_the_first_nul() {
        let mut tree = our_tree(&owned(&["chr1", "chr100"]));
        // "chr1" is the first key, NUL-padded to 6 bytes; put a byte after its first NUL.
        tree[HEADER_SIZE + NODE_HEADER_SIZE + 5] = b'X';
        assert_eq!(read_file(&tree, 0).unwrap()[0].name, "chr1");
    }

    #[test]
    fn tree_at_the_depth_limit_is_read() {
        assert_eq!(read_file(&chain_tree(MAX_TREE_DEPTH), 0).unwrap(), [ChromInfo::new("a", 100)]);
    }

    #[test]
    fn tree_deeper_than_the_limit_is_corrupt() {
        assert!(matches!(
            read_file(&chain_tree(MAX_TREE_DEPTH + 1), 0),
            Err(Error::Corrupt { .. })
        ));
    }

    #[test]
    fn bad_magic_is_corrupt_at_the_tree_offset() {
        let mut file = vec![0; 10];
        file.extend(encode(&entries(&contig_names(3)), 10, ByteOrder::Little).bytes);
        file[10] ^= 0xFF;
        assert!(matches!(read_file(&file, 10), Err(Error::Corrupt { offset: 10, .. })));
    }

    #[test]
    fn truncated_header_is_corrupt() {
        let tree = our_tree(&contig_names(3));
        assert!(matches!(read_file(&tree[..20], 0), Err(Error::Corrupt { offset: 0, .. })));
    }

    #[test]
    fn truncated_leaf_is_corrupt_at_its_slots() {
        let tree = our_tree(&contig_names(3));
        let slots_offset = (HEADER_SIZE + NODE_HEADER_SIZE) as u64;
        assert!(matches!(
            read_file(&tree[..tree.len() - 1], 0),
            Err(Error::Corrupt { offset, .. }) if offset == slots_offset
        ));
    }

    #[test]
    fn val_size_other_than_eight_is_corrupt_at_the_field() {
        let mut tree = our_tree(&contig_names(3));
        tree[12..16].copy_from_slice(&12u32.to_le_bytes());
        assert!(matches!(read_file(&tree, 0), Err(Error::Corrupt { offset: 12, .. })));
    }

    #[test]
    fn item_count_larger_than_the_file_is_corrupt_at_the_field() {
        let mut tree = our_tree(&contig_names(3));
        tree[16..24].copy_from_slice(&(1u64 << 40).to_le_bytes());
        assert!(matches!(read_file(&tree, 0), Err(Error::Corrupt { offset: 16, .. })));
    }

    #[test]
    fn item_count_overflowing_the_slot_arithmetic_is_corrupt() {
        let mut tree = our_tree(&contig_names(3));
        tree[16..24].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(matches!(read_file(&tree, 0), Err(Error::Corrupt { offset: 16, .. })));
    }

    #[test]
    fn fewer_leaf_entries_than_item_count_is_corrupt() {
        let mut tree = kent_style_tree(&contig_names(300));
        tree[16..24].copy_from_slice(&301u64.to_le_bytes());
        assert!(matches!(read_file(&tree, 0), Err(Error::Corrupt { offset: 16, .. })));
    }

    #[test]
    fn more_leaf_entries_than_item_count_is_corrupt() {
        let mut tree = kent_style_tree(&contig_names(300));
        tree[16..24].copy_from_slice(&299u64.to_le_bytes());
        assert!(matches!(read_file(&tree, 0), Err(Error::Corrupt { .. })));
    }

    #[test]
    fn duplicate_chrom_id_is_corrupt_at_the_later_id() {
        let entries = [
            ChromEntry { name: "a", chrom_id: 0, length: 1 },
            ChromEntry { name: "b", chrom_id: 1, length: 2 },
            ChromEntry { name: "c", chrom_id: 1, length: 3 },
        ];
        let encoded = encode(&entries, 0, ByteOrder::Little);
        let later = encoded.chrom_id_offsets[2];
        assert!(matches!(
            read_file(&encoded.bytes, 0),
            Err(Error::Corrupt { offset, .. }) if offset == later
        ));
    }

    #[test]
    fn chrom_id_not_below_item_count_is_corrupt_at_the_id() {
        let entries = [
            ChromEntry { name: "a", chrom_id: 0, length: 1 },
            ChromEntry { name: "b", chrom_id: 2, length: 2 },
        ];
        let encoded = encode(&entries, 0, ByteOrder::Little);
        let bad = encoded.chrom_id_offsets[1];
        assert!(matches!(
            read_file(&encoded.bytes, 0),
            Err(Error::Corrupt { offset, .. }) if offset == bad
        ));
    }

    #[test]
    fn invalid_utf8_key_is_corrupt_at_the_key() {
        let mut tree = our_tree(&contig_names(3));
        let first_key = HEADER_SIZE + NODE_HEADER_SIZE;
        tree[first_key] = 0xFF;
        assert!(matches!(
            read_file(&tree, 0),
            Err(Error::Corrupt { offset, .. }) if offset == first_key as u64
        ));
    }

    #[test]
    fn is_leaf_other_than_zero_or_one_is_corrupt_at_the_node() {
        let mut tree = our_tree(&contig_names(3));
        tree[HEADER_SIZE] = 2;
        assert!(matches!(
            read_file(&tree, 0),
            Err(Error::Corrupt { offset, .. }) if offset == HEADER_SIZE as u64
        ));
    }

    #[test]
    fn cyclic_child_offset_is_corrupt() {
        let mut tree = kent_style_tree(&contig_names(300));
        let field = root_child_offset_field(1);
        tree[field..field + 8].copy_from_slice(&(HEADER_SIZE as u64).to_le_bytes());
        assert!(matches!(
            read_file(&tree, 0),
            Err(Error::Corrupt { offset, .. }) if offset == HEADER_SIZE as u64
        ));
    }

    #[test]
    fn shared_child_offset_is_corrupt() {
        let mut tree = kent_style_tree(&contig_names(300));
        let (first, second) = (root_child_offset_field(0), root_child_offset_field(1));
        tree.copy_within(first..first + 8, second);
        assert!(matches!(read_file(&tree, 0), Err(Error::Corrupt { .. })));
    }

    #[test]
    fn child_offset_past_the_end_of_the_file_is_corrupt() {
        let mut tree = kent_style_tree(&contig_names(300));
        let field = root_child_offset_field(1);
        let past_end = tree.len() as u64 + 100;
        tree[field..field + 8].copy_from_slice(&past_end.to_le_bytes());
        assert!(matches!(
            read_file(&tree, 0),
            Err(Error::Corrupt { offset, .. }) if offset == past_end
        ));
    }

    #[test]
    #[should_panic(expected = "outside 2..=65535")]
    fn encode_rejects_a_fan_out_that_never_reduces_to_a_root() {
        encode_with_max_block_size(&entries(&contig_names(3)), 0, ByteOrder::Little, 1);
    }

    proptest! {
        #[test]
        fn decoding_mutated_trees_never_panics(
            count in 1usize..40,
            max_block_size in 2usize..6,
            mutations in prop::collection::vec((any::<prop::sample::Index>(), any::<u8>()), 1..8),
            keep in any::<prop::sample::Index>(),
        ) {
            let names = contig_names(count);
            let mut tree =
                encode_with_max_block_size(&entries(&names), 0, ByteOrder::Little, max_block_size)
                    .bytes;
            for (index, byte) in mutations {
                let at = index.index(tree.len());
                tree[at] = byte;
            }
            tree.truncate(keep.index(tree.len()) + 1);
            let _ = read_file(&tree, 0);
        }

        #[test]
        fn any_printable_names_read_back_and_are_found_by_kent_lookup(
            names in prop::collection::hash_set("[!-~]{1,12}", 1..200),
            max_block_size in 2usize..8,
        ) {
            let names: Vec<String> = names.into_iter().collect();
            let tree =
                encode_with_max_block_size(&entries(&names), 0, ByteOrder::Little, max_block_size)
                    .bytes;
            prop_assert_eq!(read_file(&tree, 0).unwrap(), expected_chroms(&names));
            assert_kent_finds_every_name(&tree, &names);
        }
    }
}
