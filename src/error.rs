//! The crate-wide error type.

/// Result alias used throughout the crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Every error the crate can return.
///
/// Data-validation errors from the writer (`UnknownChrom`, `ChromOutOfOrder`,
/// `ChromAlreadyFinished`, `UnsortedOrOverlapping`, `OutOfBounds`, `InvalidItem`,
/// `NonFiniteValue`) are detected before the writer changes any state, so the caller may skip the
/// offending item and continue. I/O and compression errors leave a writer unusable: later calls
/// return [`Error::WriterFailed`]. Reader errors never make a reader unusable.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    // ---- I/O -------------------------------------------------------------------------------
    /// An underlying read, write or seek failed.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    // ---- Reading: the input is not a (supported) bigWig -----------------------------------
    /// The first four bytes are neither the bigWig magic nor its byte-swapped form.
    #[error("not a bigWig file: magic {found:#010x}")]
    NotBigWig {
        /// The first four bytes, read as a little-endian u32. Input shorter than four bytes is
        /// padded with zeros.
        found: u32,
    },
    /// The file is a bigBed (or another BBI type) rather than a bigWig.
    #[error("file is a bigBed, not a bigWig")]
    IsBigBed,
    /// The header's version is outside 1..=4.
    #[error("unsupported bigWig version {0}")]
    UnsupportedVersion(u16),
    /// A structure is internally inconsistent: bad sub-structure magic, an offset past EOF, a
    /// count that cannot fit in the file, an unknown section type, a truncated block, ...
    #[error("corrupt bigWig at offset {offset}: {reason}")]
    Corrupt {
        /// Absolute file offset at which the problem was found.
        offset: u64,
        /// What was wrong.
        reason: String,
    },
    /// A block failed to inflate, or inflated past the hard cap.
    #[error("failed to decompress block at offset {offset}: {reason}")]
    Decompress {
        /// Absolute file offset of the block.
        offset: u64,
        /// The backend's description of the failure.
        reason: String,
    },

    // ---- Queries --------------------------------------------------------------------------
    /// The named chromosome is not in the file (or not in the writer's chromosome list).
    #[error("unknown chromosome {0:?}")]
    UnknownChrom(String),
    /// start > end, end past the chromosome length, an output buffer of the wrong length, or
    /// zero bins requested.
    #[error("invalid region {chrom}:{start}-{end}: {reason}")]
    InvalidRegion {
        /// Chromosome of the region.
        chrom: String,
        /// Start of the region (0-based).
        start: u64,
        /// End of the region (exclusive).
        end: u64,
        /// What was wrong.
        reason: &'static str,
    },
    /// A `ZoomLevel` that did not come from this file's header.
    #[error("zoom level with reduction {reduction} does not belong to this file")]
    ForeignZoomLevel {
        /// The foreign level's reduction.
        reduction: u32,
    },

    // ---- Writing: configuration -----------------------------------------------------------
    /// An option failed validation when the writer was created (message names the option).
    #[error("invalid writer option: {0}")]
    InvalidOption(String),
    /// The chromosome list is empty, or has a duplicate name, a name that is empty or not
    /// printable ASCII other than space, or a zero length.
    #[error("invalid chromosome list: {0}")]
    InvalidChromList(String),

    // ---- Writing: data validation (the writer is left unchanged) ---------------------------
    /// With `ChromIdOrder::Listed`, data for a chromosome that comes before the current one in
    /// the chromosome list.
    #[error("data for {chrom} after {current}, which comes later in the chromosome list")]
    ChromOutOfOrder {
        /// The chromosome the rejected data was for.
        chrom: String,
        /// The chromosome currently being written.
        current: String,
    },
    /// Data for a chromosome arrived after that chromosome was finished.
    #[error("data for {0} after it was finished; all data for a chromosome must be contiguous")]
    ChromAlreadyFinished(String),
    /// An item starts before the end of the previous item on the same chromosome.
    #[error(
        "{chrom}:{start} starts before the previous item ends at {previous_end} (unsorted or overlapping)"
    )]
    UnsortedOrOverlapping {
        /// Chromosome of the rejected item.
        chrom: String,
        /// Start of the rejected item.
        start: u64,
        /// End of the previous item on the chromosome, before any clipping.
        previous_end: u64,
    },
    /// An item extends past the end of its chromosome (and clipping is off).
    #[error("{chrom}:{start}-{end} extends past the chromosome length {length}")]
    OutOfBounds {
        /// Chromosome of the rejected item.
        chrom: String,
        /// Start of the rejected item.
        start: u64,
        /// End of the rejected item (exclusive).
        end: u64,
        /// Length of the chromosome.
        length: u32,
    },
    /// An interval with start >= end, a zero step or span, span > step, or an empty slice.
    #[error("invalid item on {chrom} at {start}: {reason}")]
    InvalidItem {
        /// Chromosome of the rejected item.
        chrom: String,
        /// Start of the rejected item.
        start: u64,
        /// What was wrong.
        reason: &'static str,
    },
    /// A NaN or infinite value (UCSC: "NaN values are not supported").
    #[error("non-finite value {value} at {chrom}:{position}")]
    NonFiniteValue {
        /// Chromosome of the rejected value.
        chrom: String,
        /// Position (0-based) of the rejected value.
        position: u64,
        /// The rejected value.
        value: f32,
    },
    /// A previous I/O or compression error left the writer in an unknown state.
    #[error("writer is unusable after an earlier error")]
    WriterFailed,

    // ---- Text formats ---------------------------------------------------------------------
    /// A bedGraph, wig or `chrom.sizes` line could not be parsed (1-based line number). An item
    /// that parses but that the writer rejects during a conversion is an [`Error::AtLine`].
    #[error("line {line}: {message}")]
    Parse {
        /// 1-based line number.
        line: u64,
        /// What was wrong.
        message: String,
    },
    /// The writer rejected the item on an input line of a conversion (1-based line number).
    ///
    /// `error` is one of the data-validation variants. Its message is part of this one, so it is
    /// not also given as the [`source`](std::error::Error::source).
    #[error("line {line}: {error}")]
    AtLine {
        /// 1-based line number.
        line: u64,
        /// The writer's error.
        error: Box<Error>,
    },
}

impl Error {
    /// Builds an [`Error::Corrupt`] for a structure found to be inconsistent at `offset`.
    pub(crate) fn corrupt(offset: u64, reason: impl Into<String>) -> Self {
        Self::Corrupt { offset, reason: reason.into() }
    }

    /// Whether this is a data-validation error, with which the writer rejects a push before
    /// changing any state, leaving it usable.
    pub(crate) fn is_data_validation(&self) -> bool {
        matches!(
            self,
            Self::UnknownChrom(_)
                | Self::ChromOutOfOrder { .. }
                | Self::ChromAlreadyFinished(_)
                | Self::UnsortedOrOverlapping { .. }
                | Self::OutOfBounds { .. }
                | Self::InvalidItem { .. }
                | Self::NonFiniteValue { .. }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An unsorted item rejected on line 1823.
    fn unsorted_at_line() -> Error {
        let unsorted = Error::UnsortedOrOverlapping {
            chrom: "chr1".to_string(),
            start: 500,
            previous_end: 510,
        };
        Error::AtLine { line: 1823, error: Box::new(unsorted) }
    }

    #[test]
    fn at_line_displays_the_line_then_the_wrapped_error() {
        assert_eq!(
            unsorted_at_line().to_string(),
            "line 1823: chr1:500 starts before the previous item ends at 510 (unsorted or overlapping)"
        );
    }

    #[test]
    fn at_line_has_no_source() {
        assert!(std::error::Error::source(&unsorted_at_line()).is_none());
    }
}
