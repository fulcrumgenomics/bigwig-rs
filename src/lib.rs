//! Read and write [bigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) files, and the
//! [wig](https://genome.ucsc.edu/goldenPath/help/wiggle.html) and
//! [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) text formats.
//!
//! bigWig stores one `f32` value per base range, indexed for random access and summarised at
//! several zoom levels so genome browsers can draw any scale quickly. This crate writes files in
//! the layout UCSC's own tools produce, so they open in UCSC, IGV, igv.js, JBrowse 2,
//! libBigWig, pyBigWig and bigtools, and reads files written by any of those tools. It is
//! synchronous, has no `unsafe` code, and writes with bounded memory.
//!
//! # Writing and reading
//!
//! [`WriterOptions`] configures a [`BigWigWriter`], which takes data in sorted order as
//! fixedStep runs ([`push_fixed_step`](BigWigWriter::push_fixed_step)), variableStep items
//! ([`push_variable_step`](BigWigWriter::push_variable_step)) or bedGraph intervals
//! ([`push_interval`](BigWigWriter::push_interval)), and builds the zoom levels in the same
//! pass. [`BigWigReader`] reads the result back from a path or any `Read + Seek` source.
//!
//! ```
//! use std::io::Cursor;
//!
//! use bigwig::{BigWigReader, ChromInfo, SummarySource, WriterOptions};
//!
//! // Write to memory; `WriterOptions::create` writes to a path instead.
//! let chroms = [ChromInfo::new("chr1", 1_000), ChromInfo::new("chr2", 500)];
//! let mut writer = WriterOptions::new().create_with(Cursor::new(Vec::new()), &chroms)?;
//! writer.push_fixed_step("chr1", 100, 10, 10, &[1.0, 2.0, 3.0])?;
//! writer.push_interval("chr2", 0, 50, 0.5)?;
//! let bytes = writer.finish()?.into_inner();
//!
//! // Read it back.
//! let mut reader = BigWigReader::new(Cursor::new(bytes))?;
//! assert_eq!(reader.chroms().len(), 2);
//! assert_eq!(reader.values("chr1", 108, 112)?[..2], [1.0, 1.0]);
//!
//! let items = reader.items("chr1", 0, 1_000)?.collect::<Result<Vec<_>, _>>()?;
//! assert_eq!((items[1].start, items[1].end, items[1].value), (110, 120, 2.0));
//!
//! let summary = reader.summary("chr1", 100, 130, SummarySource::Raw)?.expect("data in region");
//! assert_eq!(summary.mean(), 2.0);
//! # Ok::<(), bigwig::Error>(())
//! ```
//!
//! # Calling `finish`
//!
//! [`BigWigWriter::finish`] is mandatory. It writes the indexes, the zoom levels and the real
//! header, and reports any error doing so. Until then the file carries a zero magic number, so
//! an unfinished file is rejected by every reader rather than looking valid but truncated.
//! Dropping a writer without calling `finish` logs a warning through the [`log`] crate and does
//! not finish the file, because `Drop` cannot report errors.
//!
//! # Text formats and conversions
//!
//! The [`bedgraph`] and [`wig`] modules hold streaming text readers and writers. The [`convert`]
//! module goes between text and bigWig in either direction, reading text once from any
//! [`BufRead`](std::io::BufRead). A chromosome list comes from a UCSC `chrom.sizes` file through
//! [`ChromInfo::read_chrom_sizes`].
//!
//! ```
//! use std::io::Cursor;
//!
//! use bigwig::convert::{Selection, bedgraph_to_bigwig, bigwig_to_bedgraph};
//! use bigwig::{BigWigReader, ChromInfo, WriterOptions};
//!
//! let chroms = ChromInfo::read_chrom_sizes("chr1\t1000\n".as_bytes())?;
//! let text = "chr1\t0\t100\t1.5\nchr1\t150\t200\t2\n";
//!
//! let mut bigwig = Cursor::new(Vec::new());
//! let report = bedgraph_to_bigwig(text.as_bytes(), &chroms, &WriterOptions::new(), &mut bigwig)?;
//! assert_eq!(report.items, 2);
//!
//! let mut reader = BigWigReader::new(bigwig)?;
//! let mut round_trip = Vec::new();
//! bigwig_to_bedgraph(&mut reader, &Selection::All, true, &mut round_trip)?;
//! assert_eq!(round_trip, text.as_bytes());
//! # Ok::<(), bigwig::Error>(())
//! ```
//!
//! # Main types
//!
//! - Writing: [`WriterOptions`], [`BigWigWriter`], [`ChromInfo`], [`ChromIdOrder`], [`ZoomPolicy`]
//!   and [`Compression`].
//! - Reading: [`BigWigReader`], [`Info`], [`DataItem`], [`Items`], [`Summary`],
//!   [`SummarySource`], [`ZoomLevel`], [`ZoomRecord`] and [`ZoomRecords`].
//! - Text and conversion: [`bedgraph`], [`wig`] and [`convert`].
//! - Errors: [`Error`] and [`Result`].
//!
//! # Compression backends
//!
//! Blocks are zlib streams. The `libdeflate` feature (default) uses libdeflate through
//! `libdeflater`; the `zlib-rs` feature uses the pure-Rust zlib-rs through `flate2`. At least one
//! must be enabled; with both, libdeflate is used.
#![forbid(unsafe_code)]
#![cfg_attr(docsrs, feature(doc_cfg))]
#![allow(clippy::must_use_candidate, clippy::missing_errors_doc, clippy::missing_panics_doc)]
// Tests check exact round trips of stored values, where bit-for-bit equality is the requirement.
#![cfg_attr(test, allow(clippy::float_cmp))]

#[cfg(not(any(feature = "libdeflate", feature = "zlib-rs")))]
compile_error!("bigwig needs a compression backend: enable the `libdeflate` or `zlib-rs` feature");

mod bbi;
mod bptree;
mod compress;
mod error;
mod reader;
mod rtree;
mod section;
mod text;
mod writer;
mod zoom;

pub mod bedgraph;
pub mod convert;
pub mod wig;

pub use bbi::{ByteOrder, ChromInfo, ZoomLevel};
pub use compress::Compression;
pub use error::{Error, Result};
pub use reader::{BigWigReader, Info, Items, SummarySource, ZoomRecords};
pub use section::DataItem;
pub use writer::{BigWigWriter, ChromIdOrder, WriterOptions};
pub use zoom::{Summary, ZoomPolicy, ZoomRecord};
