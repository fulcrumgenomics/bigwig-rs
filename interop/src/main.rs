//! Cross-implementation interop harness for the `bigwig` crate.
//!
//! `generate` writes every fixture deterministically from a seed; `dump` prints this crate's
//! reading of a bigWig as TSV, for `compare.py` to check against kent, pyBigWig, bigtools and
//! `@gmod/bbi`. `run.sh` drives the whole matrix.

mod dump;
mod generate;
mod scan;

use std::error::Error;
use std::path::PathBuf;
use std::process::ExitCode;

use bigwig::convert::Selection;
use clap::{Args, Parser, Subcommand, ValueEnum};

/// Fixture generator and reader dump for the bigWig interop matrix.
#[derive(Debug, Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Writes every fixture: our bigWigs under `<OUTDIR>/ours`, and text inputs for the other
    /// tools' writers, with the items each must store, under `<OUTDIR>/text`.
    Generate {
        /// Seed for the fixture data; the same seed always gives byte-identical files.
        seed: u64,
        /// Directory to write into; created if missing.
        outdir: PathBuf,
    },
    /// Prints this crate's reading of a bigWig as TSV with a header line. Floats are printed as
    /// the shortest decimal that parses back to the same f64, and stored f32 values are widened
    /// to f64 first, so a parser that reads doubles recovers every value exactly.
    Dump {
        #[command(subcommand)]
        what: DumpCommand,
    },
    /// Reads every item of a bigWig, recomputes its total summary and, with `--round-trip`,
    /// rewrites the items with our writer and compares the copy; prints one TSV row. For the
    /// public-file corpus (corpus.sh).
    Scan {
        bigwig: PathBuf,
        /// Scratch path for the rewritten copy, removed afterwards.
        #[arg(long)]
        round_trip: Option<PathBuf>,
    },
    /// Prints `count` items of a bigWig chosen uniformly at random, as item TSV rows, so query
    /// regions can be placed where there is data (corpus.sh).
    Sample { bigwig: PathBuf, count: usize, seed: u64 },
}

#[derive(Debug, Subcommand)]
enum DumpCommand {
    /// One row of header fields and the stored total summary.
    Info { bigwig: PathBuf },
    /// The chromosomes in `chromId` order.
    Chroms { bigwig: PathBuf },
    /// The stored items, unclipped: the whole file, one chromosome, or those overlapping a
    /// region.
    Items {
        bigwig: PathBuf,
        #[command(flatten)]
        selection: SelectionArgs,
    },
    /// bedGraph text from `convert::bigwig_to_bedgraph` with abutting equal values merged, as
    /// kent `bigWigToBedGraph` prints them; items are clipped to a region.
    Bedgraph {
        bigwig: PathBuf,
        #[command(flatten)]
        selection: SelectionArgs,
    },
    /// Wig text from `convert::bigwig_to_wig`; items are clipped to a region.
    Wig {
        bigwig: PathBuf,
        #[command(flatten)]
        selection: SelectionArgs,
    },
    /// The value of every base of a region, NaN where there is no data.
    Values { bigwig: PathBuf, chrom: String, start: u32, end: u32 },
    /// Summaries of `bins` equal bins of a region; every statistic is empty for a bin with no
    /// data.
    Summaries {
        bigwig: PathBuf,
        chrom: String,
        start: u32,
        end: u32,
        bins: u32,
        /// Raw data, or the zoom level kent's `bigWigSummary` would choose.
        #[arg(long, value_enum, default_value_t = Source::BestZoom)]
        source: Source,
    },
    /// The records of one zoom level (by header index) overlapping a region.
    Zoom { bigwig: PathBuf, level: usize, chrom: String, start: u32, end: u32 },
}

/// The whole file, one chromosome (`--chrom`), or a region of one (`--chrom --start --end`).
#[derive(Debug, Args)]
struct SelectionArgs {
    /// Restrict to this chromosome.
    #[arg(long)]
    chrom: Option<String>,
    /// Region start (0-based); needs `--chrom` and `--end`.
    #[arg(long, requires_all = ["chrom", "end"])]
    start: Option<u32>,
    /// Region end (exclusive); needs `--chrom` and `--start`.
    #[arg(long, requires_all = ["chrom", "start"])]
    end: Option<u32>,
}

impl SelectionArgs {
    /// The selection these arguments name.
    fn selection(self) -> Selection {
        match (self.chrom, self.start.zip(self.end)) {
            (None, _) => Selection::All,
            (Some(chrom), None) => Selection::Chrom(chrom),
            (Some(chrom), Some((start, end))) => Selection::Region { chrom, start, end },
        }
    }
}

/// Where `dump summaries` takes its data from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum Source {
    Raw,
    BestZoom,
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("bigwig-interop: {error}");
            ExitCode::FAILURE
        }
    }
}

/// Dispatches one subcommand.
fn run(cli: Cli) -> Result<(), Box<dyn Error>> {
    match cli.command {
        Command::Generate { seed, outdir } => generate::generate(seed, &outdir),
        Command::Scan { bigwig, round_trip } => scan::scan(&bigwig, round_trip.as_deref()),
        Command::Sample { bigwig, count, seed } => scan::sample(&bigwig, count, seed),
        Command::Dump { what } => match what {
            DumpCommand::Info { bigwig } => dump::info(&bigwig),
            DumpCommand::Chroms { bigwig } => dump::chroms(&bigwig),
            DumpCommand::Items { bigwig, selection } => {
                dump::items(&bigwig, &selection.selection())
            }
            DumpCommand::Bedgraph { bigwig, selection } => {
                dump::bedgraph(&bigwig, &selection.selection())
            }
            DumpCommand::Wig { bigwig, selection } => dump::wig(&bigwig, &selection.selection()),
            DumpCommand::Values { bigwig, chrom, start, end } => {
                dump::values(&bigwig, &chrom, start, end)
            }
            DumpCommand::Summaries { bigwig, chrom, start, end, bins, source } => {
                let source = match source {
                    Source::Raw => bigwig::SummarySource::Raw,
                    Source::BestZoom => bigwig::SummarySource::BestZoom,
                };
                dump::summaries(&bigwig, &chrom, start, end, bins, source)
            }
            DumpCommand::Zoom { bigwig, level, chrom, start, end } => {
                dump::zoom_records(&bigwig, level, &chrom, start, end)
            }
        },
    }
}
