//! Block compression: the [`Compression`] setting, and the backend dispatch that deflates and
//! inflates data sections and zoom blocks.
//!
//! Exactly one backend is compiled in: libdeflate whenever the `libdeflate` feature is on,
//! otherwise zlib-rs through flate2. The two live in a `backend` module with one surface, so no
//! backend type reaches a signature and swapping backends changes no API.

use crate::error::{Error, Result};

/// Hard cap on the size a block may inflate to. It is far above the format's real maximum (a
/// data section is at most 24 + 65,535 × 12 = 786,444 bytes, and zoom blocks are typically
/// 32 KiB) and bounds the allocation a corrupt or hostile block can cause.
pub(crate) const MAX_INFLATED_SIZE: usize = 64 << 20;

/// How data sections and zoom blocks are stored.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Compression {
    /// Blocks stored raw; header `uncompressBufSize` = 0 (kent `-unc`). igv.js mis-decodes such
    /// files when a query touches more than one block, so this is for testing and kent-only
    /// pipelines.
    None,
    /// zlib-wrapped deflate at `level` (0..=9 portable; 10..=12 accepted only with libdeflate).
    /// Level 0 still writes zlib streams, made of stored deflate blocks.
    Deflate {
        /// Compression level.
        level: u8,
    },
}

impl Compression {
    /// Highest `Deflate` level the compiled backend accepts: 12 with libdeflate, 9 with zlib-rs.
    pub const MAX_LEVEL: u8 = backend::MAX_LEVEL;

    /// Checks the level against the compiled backend. An out-of-range level is
    /// [`Error::InvalidOption`]; it is never clamped.
    pub(crate) fn validate(self) -> Result<()> {
        match self {
            Self::None => Ok(()),
            Self::Deflate { level } if level <= Self::MAX_LEVEL => Ok(()),
            Self::Deflate { level } => Err(invalid_level(level)),
        }
    }
}

impl Default for Compression {
    /// Level 1: about 11% larger than level 6 on dense depth, for about 5.4× less CPU. Use
    /// `Deflate { level: 6 }` for output comparable in size to kent's zlib level 6.
    fn default() -> Self {
        Self::Deflate { level: 1 }
    }
}

/// Turns raw blocks into their stored form: one zlib (RFC 1950) stream each, or a plain copy with
/// [`Compression::None`]. It holds the backend's state, so one instance serves a whole file.
pub(crate) struct BlockCompressor {
    /// `None` when blocks are stored raw.
    deflater: Option<backend::Deflater>,
}

impl BlockCompressor {
    /// Creates a compressor, rejecting a level the compiled backend does not support.
    pub(crate) fn new(compression: Compression) -> Result<Self> {
        compression.validate()?;
        let deflater = match compression {
            Compression::None => None,
            Compression::Deflate { level } => Some(backend::Deflater::new(level)?),
        };
        Ok(Self { deflater })
    }

    /// Runs every job on the calling thread, in order, stopping at the first failure.
    ///
    /// This is the seam a parallel executor plugs into later: each job owns its buffers, so a
    /// batch can be split across threads and the outputs still written in order.
    pub(crate) fn compress_batch(&mut self, jobs: &mut [CompressionJob]) -> Result<()> {
        for job in jobs {
            self.compress(&job.raw, &mut job.compressed)?;
        }
        Ok(())
    }

    /// Encodes `raw` into `out`, replacing its contents. Level 0 still produces a zlib stream
    /// (of stored deflate blocks); only [`Compression::None`] copies the bytes through.
    pub(crate) fn compress(&mut self, raw: &[u8], out: &mut Vec<u8>) -> Result<()> {
        match &mut self.deflater {
            None => {
                out.clear();
                out.extend_from_slice(raw);
                Ok(())
            }
            Some(deflater) => deflater.compress(raw, out),
        }
    }
}

/// One block to compress: the encoded bytes in, the stored bytes out.
///
/// A job owns both buffers so it can be moved to another thread and back (`std::mem::take`),
/// and so the allocations can be recycled from one section to the next.
#[derive(Debug, Default)]
pub(crate) struct CompressionJob {
    /// The block as encoded, before compression.
    pub(crate) raw: Vec<u8>,
    /// The bytes to write to the file, filled in by [`BlockCompressor::compress_batch`].
    pub(crate) compressed: Vec<u8>,
}

/// Recovers blocks' decompressed bytes, reusing one inflate buffer across blocks.
///
/// With `uncompressBufSize == 0` blocks are stored raw and come back untouched. Otherwise the
/// buffer starts at `uncompressBufSize` bytes and doubles whenever a block does not fit, up to a
/// hard cap, because some writers understate the field (libBigWig writes a constant 32,768).
pub(crate) struct BlockDecompressor {
    /// `None` when blocks are stored raw.
    inflater: Option<backend::Inflater>,
    /// The inflate buffer; its length is the current size limit.
    buf: Vec<u8>,
    /// Largest size `buf` may grow to.
    max_inflated_size: usize,
    /// Whether the buffer has had to grow, so the warning is logged once per instance.
    warned: bool,
}

impl BlockDecompressor {
    /// Creates a decompressor for a file whose header carries `uncompress_buf_size`.
    pub(crate) fn new(uncompress_buf_size: u32) -> Self {
        Self::with_cap(uncompress_buf_size, MAX_INFLATED_SIZE)
    }

    /// As [`Self::new`], with the cap on the inflated size (at least 1 byte) made explicit so
    /// tests can use a small one.
    fn with_cap(uncompress_buf_size: u32, max_inflated_size: usize) -> Self {
        if uncompress_buf_size == 0 {
            return Self { inflater: None, buf: Vec::new(), max_inflated_size, warned: false };
        }
        // A hostile header can claim any size; the cap bounds the allocation.
        let initial = usize::try_from(uncompress_buf_size)
            .map_or(max_inflated_size, |size| size.min(max_inflated_size));
        Self {
            inflater: Some(backend::Inflater::new()),
            buf: vec![0; initial],
            max_inflated_size,
            warned: false,
        }
    }

    /// Returns the decompressed form of `block`, which was read from file offset `offset`.
    ///
    /// Raw blocks are returned as given; inflated ones borrow the reusable buffer. Every
    /// failure, including a block that inflates past the cap, is [`Error::Decompress`] at
    /// `offset`.
    pub(crate) fn decompress<'a>(&'a mut self, block: &'a [u8], offset: u64) -> Result<&'a [u8]> {
        let Self { inflater, buf, max_inflated_size, warned } = self;
        let Some(inflater) = inflater else { return Ok(block) };
        loop {
            match inflater.inflate(block, buf) {
                Ok(len) => return Ok(&buf[..len]),
                Err(InflateError::Invalid(reason)) => {
                    return Err(Error::Decompress { offset, reason });
                }
                Err(InflateError::InsufficientSpace) => {
                    if buf.len() >= *max_inflated_size {
                        return Err(Error::Decompress {
                            offset,
                            reason: format!("block inflates past the {max_inflated_size}-byte cap"),
                        });
                    }
                    if !*warned {
                        log::warn!(
                            "block at offset {offset} inflates past the {}-byte uncompressBufSize \
                             the file declares; growing the buffer (reported once per file)",
                            buf.len()
                        );
                        *warned = true;
                    }
                    let doubled = buf.len().saturating_mul(2).min(*max_inflated_size);
                    buf.resize(doubled, 0);
                }
            }
        }
    }
}

/// Why a backend could not inflate a block.
enum InflateError {
    /// The output buffer is too small; the caller grows it and retries.
    InsufficientSpace,
    /// The bytes are not a valid zlib stream; the backend's description of the problem.
    Invalid(String),
}

/// Builds the error for a level the compiled backend does not support.
fn invalid_level(level: u8) -> Error {
    Error::InvalidOption(format!(
        "compression: level {level} is outside 0..={} for the {} backend",
        backend::MAX_LEVEL,
        backend::NAME
    ))
}

/// Builds the error for a compression failure. It is reported as I/O so that no backend's error
/// type reaches [`Error`], and switching backends changes no API.
fn compression_failed(reason: impl std::fmt::Display) -> Error {
    Error::Io(std::io::Error::other(format!("zlib compression failed: {reason}")))
}

/// The libdeflate backend: whole-buffer deflate and inflate, with levels 0..=12.
#[cfg(feature = "libdeflate")]
mod backend {
    use libdeflater::{CompressionLvl, Compressor, DecompressionError, Decompressor};

    use super::{InflateError, compression_failed, invalid_level};
    use crate::error::Result;

    /// Name used in error messages.
    pub(super) const NAME: &str = "libdeflate";
    /// Highest level the backend accepts.
    pub(super) const MAX_LEVEL: u8 = 12;

    /// A zlib compressor fixed at one level.
    pub(super) struct Deflater {
        inner: Compressor,
    }

    impl Deflater {
        /// Creates a compressor at `level`, which must already be validated.
        pub(super) fn new(level: u8) -> Result<Self> {
            let lvl = CompressionLvl::new(i32::from(level)).map_err(|_| invalid_level(level))?;
            Ok(Self { inner: Compressor::new(lvl) })
        }

        /// Encodes `raw` as one zlib stream in `out`, replacing its contents.
        pub(super) fn compress(&mut self, raw: &[u8], out: &mut Vec<u8>) -> Result<()> {
            // The bound makes `InsufficientSpace` impossible; whatever `out` held is overwritten.
            out.resize(self.inner.zlib_compress_bound(raw.len()), 0);
            let len = self.inner.zlib_compress(raw, out).map_err(compression_failed)?;
            out.truncate(len);
            Ok(())
        }
    }

    /// A zlib decompressor.
    pub(super) struct Inflater {
        inner: Decompressor,
    }

    impl Inflater {
        /// Creates a decompressor.
        pub(super) fn new() -> Self {
            Self { inner: Decompressor::new() }
        }

        /// Inflates `block` into `buf`, returning the inflated length.
        pub(super) fn inflate(
            &mut self,
            block: &[u8],
            buf: &mut [u8],
        ) -> std::result::Result<usize, InflateError> {
            self.inner.zlib_decompress(block, buf).map_err(|err| match err {
                DecompressionError::InsufficientSpace => InflateError::InsufficientSpace,
                DecompressionError::BadData => InflateError::Invalid(
                    "not a valid zlib stream (bad header, data or checksum)".to_string(),
                ),
            })
        }
    }
}

/// The zlib-rs backend through flate2: streaming deflate and inflate, with levels 0..=9.
#[cfg(all(feature = "zlib-rs", not(feature = "libdeflate")))]
mod backend {
    use flate2::{Compress, Decompress, FlushCompress, FlushDecompress, Status};

    use super::{InflateError, compression_failed};
    use crate::error::Result;

    /// Name used in error messages.
    pub(super) const NAME: &str = "zlib-rs";
    /// Highest level the backend accepts.
    pub(super) const MAX_LEVEL: u8 = 9;

    /// A zlib compressor fixed at one level.
    pub(super) struct Deflater {
        inner: Compress,
    }

    impl Deflater {
        /// Creates a compressor at `level`, which must already be validated.
        #[expect(
            clippy::unnecessary_wraps,
            reason = "shares the libdeflate constructor's signature, which can fail"
        )]
        pub(super) fn new(level: u8) -> Result<Self> {
            Ok(Self { inner: Compress::new(flate2::Compression::new(u32::from(level)), true) })
        }

        /// Encodes `raw` as one zlib stream in `out`, replacing its contents.
        pub(super) fn compress(&mut self, raw: &[u8], out: &mut Vec<u8>) -> Result<()> {
            self.inner.reset();
            out.clear();
            // Enough for incompressible input: stored blocks add 5 bytes per 64 KiB, plus the
            // 6-byte zlib wrapper. Anything more is grown below.
            out.reserve(raw.len() + raw.len() / 1000 + 64);
            let mut input = raw;
            loop {
                let (in_before, out_before) = (self.inner.total_in(), self.inner.total_out());
                let status = self
                    .inner
                    .compress_vec(input, out, FlushCompress::Finish)
                    .map_err(compression_failed)?;
                let consumed = progress(in_before, self.inner.total_in());
                let produced = progress(out_before, self.inner.total_out());
                input = &input[consumed..];
                match status {
                    Status::StreamEnd => return Ok(()),
                    Status::Ok | Status::BufError if consumed > 0 || produced > 0 => {
                        if out.len() == out.capacity() {
                            out.reserve(out.len().max(1024));
                        }
                    }
                    Status::Ok | Status::BufError => {
                        return Err(compression_failed("the stream stalled"));
                    }
                }
            }
        }
    }

    /// A zlib decompressor.
    pub(super) struct Inflater {
        inner: Decompress,
    }

    impl Inflater {
        /// Creates a decompressor.
        pub(super) fn new() -> Self {
            Self { inner: Decompress::new(true) }
        }

        /// Inflates `block` into `buf`, returning the inflated length.
        ///
        /// The stream is inflated from the start on every call, so a retry with a larger buffer
        /// repeats the work; the total stays within twice the final pass because the buffer
        /// doubles each time.
        pub(super) fn inflate(
            &mut self,
            block: &[u8],
            buf: &mut [u8],
        ) -> std::result::Result<usize, InflateError> {
            self.inner.reset(true);
            let mut input = block;
            let mut filled = 0;
            loop {
                let (in_before, out_before) = (self.inner.total_in(), self.inner.total_out());
                let status = self
                    .inner
                    .decompress(input, &mut buf[filled..], FlushDecompress::None)
                    .map_err(|err| InflateError::Invalid(err.to_string()))?;
                let consumed = progress(in_before, self.inner.total_in());
                let produced = progress(out_before, self.inner.total_out());
                input = &input[consumed..];
                filled += produced;
                match status {
                    Status::StreamEnd => return Ok(filled),
                    // A full buffer is only insufficient once the trailer cannot be consumed
                    // either, so an exactly-sized buffer needs no retry.
                    Status::Ok | Status::BufError if consumed > 0 || produced > 0 => {}
                    Status::Ok | Status::BufError if filled == buf.len() => {
                        return Err(InflateError::InsufficientSpace);
                    }
                    Status::Ok | Status::BufError => {
                        return Err(InflateError::Invalid("truncated zlib stream".to_string()));
                    }
                }
            }
        }
    }

    /// Bytes consumed or produced since `before`. Each is bounded by a slice length, so the
    /// conversion cannot fail.
    fn progress(before: u64, after: u64) -> usize {
        usize::try_from(after - before).unwrap_or(usize::MAX)
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    /// A compressible block: a slow ramp with a short repeating ripple, so every level shrinks it.
    fn sample_block(len: usize) -> Vec<u8> {
        (0..len).map(|i| u8::try_from((i / 3 + i % 7) % 251).unwrap()).collect()
    }

    /// Compresses `raw` at `compression` with a fresh compressor.
    fn compress(compression: Compression, raw: &[u8]) -> Vec<u8> {
        let mut compressor = BlockCompressor::new(compression).unwrap();
        let mut out = Vec::new();
        compressor.compress(raw, &mut out).unwrap();
        out
    }

    /// A zlib stream of one stored block holding `abc`, built by hand from RFC 1950 and 1951.
    fn hand_made_zlib_abc() -> Vec<u8> {
        vec![0x78, 0x01, 0x01, 0x03, 0x00, 0xFC, 0xFF, b'a', b'b', b'c', 0x02, 0x4D, 0x01, 0x27]
    }

    /// Asserts that `stream` starts with a zlib header: deflate method, window 32 KiB, and a
    /// check value that makes `CMF × 256 + FLG` a multiple of 31.
    fn assert_zlib_header(stream: &[u8]) {
        assert_eq!(stream[0], 0x78, "CMF should say deflate with a 32 KiB window");
        assert_eq!((u32::from(stream[0]) << 8 | u32::from(stream[1])) % 31, 0);
    }

    fn assert_send<T: Send>() {}

    #[test]
    fn default_is_deflate_level_one() {
        assert_eq!(Compression::default(), Compression::Deflate { level: 1 });
    }

    #[test]
    fn max_level_matches_the_compiled_backend() {
        let expected = if cfg!(feature = "libdeflate") { 12 } else { 9 };
        assert_eq!(Compression::MAX_LEVEL, expected);
    }

    #[test]
    fn none_is_always_valid() {
        assert!(Compression::None.validate().is_ok());
    }

    #[test]
    fn every_level_up_to_the_backend_maximum_is_valid() {
        for level in 0..=Compression::MAX_LEVEL {
            assert!(Compression::Deflate { level }.validate().is_ok(), "level {level}");
        }
    }

    #[test]
    fn level_above_the_backend_maximum_is_invalid_option_naming_the_option() {
        for level in [Compression::MAX_LEVEL + 1, u8::MAX] {
            match (Compression::Deflate { level }).validate() {
                Err(Error::InvalidOption(message)) => {
                    assert!(message.contains("compression"), "{message}");
                    assert!(message.contains(&level.to_string()), "{message}");
                }
                other => panic!("expected InvalidOption, got {other:?}"),
            }
        }
    }

    #[cfg(not(feature = "libdeflate"))]
    #[test]
    fn level_ten_is_rejected_without_libdeflate() {
        assert!(matches!(
            Compression::Deflate { level: 10 }.validate(),
            Err(Error::InvalidOption(_))
        ));
    }

    #[test]
    fn new_compressor_rejects_an_out_of_range_level() {
        let result =
            BlockCompressor::new(Compression::Deflate { level: Compression::MAX_LEVEL + 1 });
        assert!(matches!(result, Err(Error::InvalidOption(_))));
    }

    #[test]
    fn none_copies_the_block_through() {
        let raw = sample_block(1000);
        assert_eq!(compress(Compression::None, &raw), raw);
    }

    #[test]
    fn level_zero_emits_a_zlib_stream_of_stored_blocks() {
        let raw = sample_block(100_000);
        let stream = compress(Compression::Deflate { level: 0 }, &raw);
        assert_zlib_header(&stream);
        assert!(stream.len() > raw.len(), "stored blocks add framing, never remove bytes");
        assert!(stream.len() < raw.len() + raw.len() / 1000 + 64);
        let mut decompressor = BlockDecompressor::new(100_000);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
    }

    #[test]
    fn level_one_round_trips() {
        let raw = sample_block(50_000);
        let stream = compress(Compression::Deflate { level: 1 }, &raw);
        assert_zlib_header(&stream);
        assert!(stream.len() < raw.len() / 2, "the ramp is highly compressible");
        let mut decompressor = BlockDecompressor::new(50_000);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
    }

    #[test]
    fn level_six_round_trips() {
        let raw = sample_block(50_000);
        let stream = compress(Compression::Deflate { level: 6 }, &raw);
        assert_zlib_header(&stream);
        let mut decompressor = BlockDecompressor::new(50_000);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
    }

    #[test]
    fn level_nine_round_trips() {
        let raw = sample_block(50_000);
        let stream = compress(Compression::Deflate { level: 9 }, &raw);
        assert_zlib_header(&stream);
        let mut decompressor = BlockDecompressor::new(50_000);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
    }

    #[cfg(feature = "libdeflate")]
    #[test]
    fn level_twelve_round_trips() {
        let raw = sample_block(50_000);
        let stream = compress(Compression::Deflate { level: 12 }, &raw);
        assert_zlib_header(&stream);
        let mut decompressor = BlockDecompressor::new(50_000);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
    }

    #[test]
    fn level_six_output_is_smaller_than_level_zero_output() {
        let raw = sample_block(50_000);
        let stored = compress(Compression::Deflate { level: 0 }, &raw);
        let packed = compress(Compression::Deflate { level: 6 }, &raw);
        assert!(packed.len() < stored.len());
    }

    #[test]
    fn empty_block_round_trips() {
        let stream = compress(Compression::Deflate { level: 6 }, &[]);
        assert_zlib_header(&stream);
        let mut decompressor = BlockDecompressor::new(1);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), &[] as &[u8]);
    }

    #[test]
    fn compress_replaces_whatever_the_output_buffer_held() {
        let mut compressor = BlockCompressor::new(Compression::Deflate { level: 1 }).unwrap();
        let raw = sample_block(2000);
        let mut out = vec![0xEE; 5000];
        compressor.compress(&raw, &mut out).unwrap();
        let mut fresh = Vec::new();
        compressor.compress(&raw, &mut fresh).unwrap();
        assert_eq!(out, fresh);
    }

    #[test]
    fn compressor_is_reusable_across_blocks() {
        let mut compressor = BlockCompressor::new(Compression::Deflate { level: 1 }).unwrap();
        let mut decompressor = BlockDecompressor::new(4000);
        let mut out = Vec::new();
        for len in [10, 4000, 1] {
            let raw = sample_block(len);
            compressor.compress(&raw, &mut out).unwrap();
            assert_eq!(decompressor.decompress(&out, 0).unwrap(), raw);
        }
    }

    #[test]
    fn compress_batch_fills_every_job_in_order() {
        let mut compressor = BlockCompressor::new(Compression::Deflate { level: 1 }).unwrap();
        let mut jobs: Vec<CompressionJob> = (1..=3)
            .map(|i| CompressionJob { raw: sample_block(i * 1000), compressed: vec![0xEE; 7] })
            .collect();
        compressor.compress_batch(&mut jobs).unwrap();
        let mut decompressor = BlockDecompressor::new(3000);
        for (i, job) in jobs.iter().enumerate() {
            assert_eq!(
                decompressor.decompress(&job.compressed, 0).unwrap(),
                sample_block((i + 1) * 1000)
            );
        }
    }

    #[test]
    fn compress_batch_with_no_jobs_succeeds() {
        let mut compressor = BlockCompressor::new(Compression::None).unwrap();
        compressor.compress_batch(&mut []).unwrap();
    }

    #[test]
    fn compress_batch_with_none_copies_each_job() {
        let mut compressor = BlockCompressor::new(Compression::None).unwrap();
        let mut jobs = [CompressionJob { raw: sample_block(10), ..CompressionJob::default() }];
        compressor.compress_batch(&mut jobs).unwrap();
        assert_eq!(jobs[0].compressed, sample_block(10));
    }

    #[test]
    fn jobs_and_codecs_can_move_between_threads() {
        assert_send::<CompressionJob>();
        assert_send::<BlockCompressor>();
        assert_send::<BlockDecompressor>();
    }

    #[test]
    fn raw_mode_returns_any_block_unchanged() {
        let mut decompressor = BlockDecompressor::new(0);
        let block = sample_block(100);
        assert_eq!(decompressor.decompress(&block, 0).unwrap(), block);
        assert_eq!(decompressor.decompress(&[0xFF, 0x00], 9).unwrap(), &[0xFF, 0x00]);
        assert!(decompressor.buf.is_empty(), "raw mode allocates nothing");
    }

    #[test]
    fn hand_made_zlib_stream_inflates() {
        let mut decompressor = BlockDecompressor::new(3);
        assert_eq!(decompressor.decompress(&hand_made_zlib_abc(), 0).unwrap(), b"abc");
    }

    #[test]
    fn exact_uncompress_buf_size_needs_no_retry() {
        let raw = sample_block(10_000);
        let stream = compress(Compression::default(), &raw);
        let mut decompressor = BlockDecompressor::new(10_000);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
        assert!(!decompressor.warned);
        assert_eq!(decompressor.buf.len(), 10_000);
    }

    #[test]
    fn undersized_uncompress_buf_size_grows_and_succeeds() {
        let raw = sample_block(100_000);
        let stream = compress(Compression::default(), &raw);
        let mut decompressor = BlockDecompressor::new(1);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
        assert!(decompressor.warned);
        assert!(decompressor.buf.len() >= 100_000);
    }

    #[test]
    fn buffer_doubles_from_uncompress_buf_size_until_the_block_fits() {
        let raw = sample_block(10_000);
        let stream = compress(Compression::default(), &raw);
        let mut decompressor = BlockDecompressor::new(1000);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
        assert_eq!(decompressor.buf.len(), 16_000);
    }

    #[test]
    fn one_byte_short_retries_once_and_succeeds() {
        let raw = sample_block(10_000);
        let stream = compress(Compression::default(), &raw);
        let mut decompressor = BlockDecompressor::new(9_999);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
        assert_eq!(decompressor.buf.len(), 19_998);
    }

    #[test]
    fn block_inflating_past_the_cap_is_a_decompress_error() {
        let raw = vec![0u8; 4096];
        let stream = compress(Compression::default(), &raw);
        let mut decompressor = BlockDecompressor::with_cap(16, 1024);
        match decompressor.decompress(&stream, 77) {
            Err(Error::Decompress { offset, reason }) => {
                assert_eq!(offset, 77);
                assert!(reason.contains("1024"), "{reason}");
            }
            other => panic!("expected Decompress, got {other:?}"),
        }
        assert!(decompressor.buf.len() <= 1024);
    }

    #[test]
    fn block_inflating_to_exactly_the_cap_succeeds() {
        let raw = vec![0u8; 1024];
        let stream = compress(Compression::default(), &raw);
        let mut decompressor = BlockDecompressor::with_cap(16, 1024);
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
    }

    #[test]
    fn initial_buffer_is_clamped_to_the_cap() {
        let decompressor = BlockDecompressor::with_cap(u32::MAX, 4096);
        assert_eq!(decompressor.buf.len(), 4096);
    }

    #[test]
    fn default_cap_is_64_mib() {
        assert_eq!(MAX_INFLATED_SIZE, 64 * 1024 * 1024);
        assert_eq!(BlockDecompressor::new(7).max_inflated_size, MAX_INFLATED_SIZE);
    }

    #[test]
    fn corrupt_stream_is_a_decompress_error_at_the_block_offset() {
        let raw = sample_block(10_000);
        let mut stream = compress(Compression::Deflate { level: 6 }, &raw);
        let middle = stream.len() / 2;
        stream[middle] ^= 0xFF;
        let mut decompressor = BlockDecompressor::new(10_000);
        assert!(matches!(
            decompressor.decompress(&stream, 4242),
            Err(Error::Decompress { offset: 4242, .. })
        ));
    }

    #[test]
    fn bad_checksum_is_a_decompress_error() {
        let mut stream = hand_made_zlib_abc();
        let last = stream.len() - 1;
        stream[last] ^= 0x01;
        let mut decompressor = BlockDecompressor::new(3);
        assert!(matches!(decompressor.decompress(&stream, 0), Err(Error::Decompress { .. })));
    }

    #[test]
    fn truncated_stream_is_a_decompress_error() {
        let raw = sample_block(10_000);
        let stream = compress(Compression::Deflate { level: 6 }, &raw);
        let truncated = &stream[..stream.len() / 2];
        let mut decompressor = BlockDecompressor::new(10_000);
        assert!(matches!(decompressor.decompress(truncated, 0), Err(Error::Decompress { .. })));
    }

    #[test]
    fn truncated_stream_with_an_undersized_buffer_is_a_decompress_error() {
        let raw = sample_block(10_000);
        let stream = compress(Compression::Deflate { level: 6 }, &raw);
        let truncated = &stream[..stream.len() / 2];
        let mut decompressor = BlockDecompressor::new(16);
        assert!(matches!(decompressor.decompress(truncated, 0), Err(Error::Decompress { .. })));
    }

    #[test]
    fn empty_block_is_a_decompress_error_when_compressed() {
        let mut decompressor = BlockDecompressor::new(16);
        assert!(matches!(decompressor.decompress(&[], 0), Err(Error::Decompress { .. })));
    }

    #[test]
    fn decompressor_stays_usable_after_an_error() {
        let raw = sample_block(1000);
        let stream = compress(Compression::default(), &raw);
        let mut decompressor = BlockDecompressor::new(1000);
        assert!(decompressor.decompress(&[1, 2, 3], 0).is_err());
        assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw);
    }

    proptest! {
        #[test]
        fn arbitrary_bytes_round_trip_at_level_one(raw in prop::collection::vec(any::<u8>(), 0..4096)) {
            let stream = compress(Compression::Deflate { level: 1 }, &raw);
            let mut decompressor = BlockDecompressor::new(1);
            prop_assert_eq!(decompressor.decompress(&stream, 0).unwrap(), raw.as_slice());
        }

        #[test]
        fn arbitrary_bytes_never_panic_the_decompressor(block in prop::collection::vec(any::<u8>(), 0..512)) {
            let mut decompressor = BlockDecompressor::with_cap(8, 4096);
            let _ = decompressor.decompress(&block, 0);
        }
    }
}
