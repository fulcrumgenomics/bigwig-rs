#!/usr/bin/env python3
"""The cross-implementation interop checks, run by run.sh.

Usage: compare.py <fixture-dir> [fixture ...]

<fixture-dir> is what `bigwig-interop generate` wrote. Naming fixtures (as in its manifest.tsv,
plus `theirs` for the other tools' writers) restricts the run to them. Tool locations come from
environment variables (see run.sh), each defaulting to the tool's name on PATH; run.sh also sets
IGV_CLASSPATH, for IgvDump.java compiled against IGV's jars, and IGVJS_DIR, the igv.js source.
Every check prints one PASS/FAIL line (NOTE for a recorded, not required, failure), a summary
follows, and the exit status is 1 if any required check failed.

Both directions are covered:

- our files, read by them: kent bigWigInfo, bigWigToBedGraph, bigWigToWig and bigWigSummary;
  pyBigWig; bigtools bigwiginfo and bigwigtobedgraph; @gmod/bbi; IGV desktop; igv.js;
- their files, read by us: kent bedGraphToBigWig and wigToBigWig, pyBigWig addEntries and
  bigtools bedgraphtobigwig, fed the generated text, which our reader must reproduce exactly.

Tolerances: raw values, intervals, bases covered, min and max are exact; raw means and standard
deviations agree to a relative 1e-6; zoom-based summaries to f32 precision. Where a tool prints
fewer digits than that (kent's %g and %f), the comparison allows half a unit in the last printed
digit on top. Big-endian files must be read by kent and by us; the other readers' results on
them are recorded, not required, as are igv.js's on uncompressed files.
"""

import csv
import io
import json
import math
import os
import struct
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field

# Relative tolerance for raw-data means and standard deviations.
RAW_REL_TOL = 1e-6
# Relative tolerance for zoom-based statistics: f32 has a 24-bit significand, and both sides
# read the same f32 record fields, so agreement to a few f32 ulps is expected.
ZOOM_REL_TOL = 1e-6
# Fixtures whose data must not be required of readers that are only little-endian.
BIG_ENDIAN_FIXTURES = {"big_endian"}
# Fixtures with no data at all.
EMPTY_FIXTURES = {"empty"}
# Fixtures written with Compression::None.
UNCOMPRESSED_FIXTURES = {"uncompressed"}
# Most sampled chromosomes with data per fixture, and the widest region queried for values.
MAX_SAMPLED_CHROMS = 3
WINDOW_LENGTH = 30_000
# Bins per raw-summary query, and most bins per zoom-summary query.
RAW_BINS = 10
MAX_ZOOM_BINS = 8
# kent bigWigSummary statistics compared, as its -type names them.
KENT_SUMMARY_TYPES = ("mean", "min", "max", "std", "coverage")
# pyBigWig stats types compared. Zoom-based "sum" is not: libBigWig truncates each record's
# overlap fraction to an integer (bwStats.c; https://github.com/deeptools/pyBigWig/issues/146).
PYBIGWIG_STATS_TYPES = ("mean", "min", "max", "std", "coverage")
# IGV asks for zoom data at chromosome length / (2^zoom × 700) bases per pixel at zoom 0, 1, 2,
# ...; the deepest zoom queried, though the queries stop sooner, at the first zoom where IGV
# finds no zoom level fine enough and summarises raw data itself.
IGV_MAX_ZOOM = 24


def tool(variable, default):
    """The command for a tool: its environment variable, or its name on PATH."""
    return os.environ.get(variable) or default


TOOLS = {
    "interop": tool("INTEROP_BIN", "bigwig-interop"),
    "bigWigInfo": tool("BIGWIGINFO", "bigWigInfo"),
    "bigWigToBedGraph": tool("BIGWIGTOBEDGRAPH", "bigWigToBedGraph"),
    "bigWigToWig": tool("BIGWIGTOWIG", "bigWigToWig"),
    "bigWigSummary": tool("BIGWIGSUMMARY", "bigWigSummary"),
    "bedGraphToBigWig": tool("BEDGRAPHTOBIGWIG", "bedGraphToBigWig"),
    "wigToBigWig": tool("WIGTOBIGWIG", "wigToBigWig"),
    "bigwiginfo": tool("BIGTOOLS_BIGWIGINFO", "bigwiginfo"),
    "bigwigtobedgraph": tool("BIGTOOLS_BIGWIGTOBEDGRAPH", "bigwigtobedgraph"),
    "bedgraphtobigwig": tool("BIGTOOLS_BEDGRAPHTOBIGWIG", "bedgraphtobigwig"),
    "js": tool("JS_RUNTIME", "node"),
    "bbi_dump": tool("BBI_DUMP", os.path.join(os.path.dirname(os.path.abspath(__file__)), "bbi_dump.mjs")),
    "java": tool("JAVA", "java"),
    "igv_classpath": tool("IGV_CLASSPATH", ""),
    "igvjs_dump": tool("IGVJS_DUMP",
                       os.path.join(os.path.dirname(os.path.abspath(__file__)), "igvjs_dump.mjs")),
    "igvjs_dir": tool("IGVJS_DIR", ""),
}


class Mismatch(Exception):
    """A check found a difference; the message says where."""


# ---------------------------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------------------------


@dataclass
class Result:
    tool: str
    check: str
    fixture: str
    ok: bool
    required: bool
    detail: str


@dataclass
class Report:
    results: list = field(default_factory=list)

    def run(self, tool_name, check, fixture, required, function, *args):
        """Runs one check, recording and printing its outcome. `function` returns a short
        description of what it compared, or raises Mismatch."""
        try:
            detail = function(*args)
            ok = True
        except Mismatch as error:
            detail, ok = str(error), False
        except Exception as error:  # A crash in a tool or in the harness fails the check.
            detail, ok = f"{type(error).__name__}: {error}", False
        result = Result(tool_name, check, fixture, ok, required, detail or "")
        self.results.append(result)
        status = "PASS" if ok else ("FAIL" if required else "NOTE")
        print(f"{status}  {tool_name:<9} {check:<24} {fixture:<22} {result.detail}", flush=True)
        return ok

    def summarise(self):
        """Prints per-tool, per-check totals and, for each fixture that some reader is not
        required to read, what every reader made of it; returns the exit status."""
        print("\n==== Summary (pass / required failures / recorded failures) ====")
        totals = {}
        for result in self.results:
            counts = totals.setdefault((result.tool, result.check), [0, 0, 0])
            if result.ok:
                counts[0] += 1
            elif result.required:
                counts[1] += 1
            else:
                counts[2] += 1
        for (tool_name, check), (passed, failed, recorded) in sorted(totals.items()):
            print(f"{tool_name:<9} {check:<24} {passed:>4} {failed:>4} {recorded:>4}")

        for fixture in sorted({r.fixture for r in self.results if not r.required}):
            print(f"\n==== The {fixture} fixture, by tool ====")
            for tool_name in sorted({r.tool for r in self.results if r.fixture == fixture}):
                mine = [r for r in self.results if r.fixture == fixture and r.tool == tool_name]
                failed = [r for r in mine if not r.ok]
                verdict = "reads it" if not failed else f"fails ({failed[0].check}: {failed[0].detail})"
                print(f"{tool_name:<9} {verdict}")

        failures = [r for r in self.results if r.required and not r.ok]
        print(f"\n{len(self.results)} checks, {len(failures)} required failures")
        for result in failures:
            print(f"  FAIL {result.tool} {result.check} {result.fixture}: {result.detail}")
        return 1 if failures else 0


# ---------------------------------------------------------------------------------------------
# Running tools and reading their output
# ---------------------------------------------------------------------------------------------


def run(command, check=True):
    """Runs `command`, returning the completed process; raises Mismatch on a non-zero exit when
    `check`."""
    process = subprocess.run(command, capture_output=True, text=True)
    if check and process.returncode != 0:
        stderr = process.stderr.strip().splitlines()
        raise Mismatch(f"{os.path.basename(command[0])} exited {process.returncode}: "
                       f"{stderr[-1] if stderr else ''}")
    return process


def dump(*args):
    """Runs `bigwig-interop dump ...`, returning its rows as dicts (no rows for empty output)."""
    output = run([TOOLS["interop"], "dump", *args]).stdout
    return list(csv.DictReader(io.StringIO(output), delimiter="\t"))


def optional_float(text):
    """A float column, or None where the dump left it empty."""
    return None if text == "" else float(text)


def item_tuples(rows):
    """Item rows as (chrom, start, end, value) tuples."""
    return [(r["chrom"], int(r["start"]), int(r["end"]), float(r["value"])) for r in rows]


# ---------------------------------------------------------------------------------------------
# Number comparisons
# ---------------------------------------------------------------------------------------------


def f32_ulp(value):
    """The spacing of f32 values at `value`."""
    if value == 0:
        return 2.0 ** -149
    exponent = math.frexp(abs(value))[1]
    return 2.0 ** (exponent - 24)


def same_f32(text, value):
    """Whether `text`, a shortest round-trip printing of an f32 (bigtools), denotes `value`, the
    exact double of an f32: the decimal is closer to `value` than to any other f32."""
    return abs(float(text) - value) <= f32_ulp(value) / 2


def close(a, b, rel):
    """Relative closeness, with NaN equal to NaN."""
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return abs(a - b) <= rel * max(abs(a), abs(b))


def printed_close(text, value, rel, half_unit):
    """Whether `text`, a value printed with limited digits, matches `value` to a relative
    `rel` plus `half_unit`, half a unit in the last printed digit."""
    if text.lower() in ("nan", "-nan"):
        return math.isnan(value)
    if math.isnan(value):
        return False
    return abs(float(text) - value) <= half_unit + rel * abs(value)


def g_matches(text, value, rel):
    """Whether kent's `%g` printing `text` matches `value`: the same printed text, or within
    `rel` plus half a unit in the sixth significant digit."""
    if text == "%g" % value:
        return True
    if value == 0 or math.isnan(value):
        return printed_close(text, value, rel, 0.0)
    half_unit = 0.5 * 10.0 ** (math.floor(math.log10(abs(value))) - 5)
    return printed_close(text, value, rel, half_unit)


def f_matches(text, value, rel):
    """Whether kent's `%f` printing `text` matches `value` (six decimal places)."""
    return printed_close(text, value, rel, 5e-7)


# ---------------------------------------------------------------------------------------------
# Zoom-level selection rules
# ---------------------------------------------------------------------------------------------


def kent_level(reductions, bin_width):
    """kent's bbiBestZoom: the largest reduction strictly below bin_width // 2, or None when
    that is at most 1 or no level qualifies."""
    desired = bin_width // 2
    if desired <= 1:
        return None
    below = [(r, i) for i, r in enumerate(reductions) if r < desired]
    return max(below)[1] if below else None


def libbigwig_level(reductions, bin_width):
    """libBigWig's determineZoomLevel: the largest reduction at most bin_width // 2."""
    desired = bin_width // 2
    below = [(r, i) for i, r in enumerate(reductions) if r <= desired]
    return max(below)[1] if below else None


def bbi_level(reductions, scale):
    """@gmod/bbi's getView: scanning from the last header, the first level whose reduction is at
    most 2 / scale."""
    for index in range(len(reductions) - 1, -1, -1):
        if reductions[index] <= 2 / scale:
            return index
    return None


def igv_level(reductions, bp_per_pixel, tolerance):
    """IGV's BBFile.zoomLevelForScale: the coarsest level whose reduction is below
    `bp_per_pixel`; failing that, the finest, if its reduction // `tolerance` is below it; else
    None, and IGV summarises raw data itself."""
    below = [(r, i) for i, r in enumerate(reductions) if r < bp_per_pixel]
    if below:
        return max(below)[1]
    if not reductions:
        return None
    finest = min(range(len(reductions)), key=reductions.__getitem__)
    return finest if reductions[finest] // tolerance < bp_per_pixel else None


def igvjs_level(reductions, bp_per_pixel):
    """igv.js's zoomLevelForScale: the coarsest level whose reduction is below `bp_per_pixel`."""
    below = [(r, i) for i, r in enumerate(reductions) if r < bp_per_pixel]
    return max(below)[1] if below else None


# ---------------------------------------------------------------------------------------------
# Our view of a fixture, and the queries asked of every tool
# ---------------------------------------------------------------------------------------------


@dataclass
class Query:
    """A summary query: `bins` bins over [start, end) of `chrom`. `zoom_level` is the header
    index of the level kent's rule picks, or None for raw data; `aligned` means the region and
    bins fall on that level's record grid."""

    chrom: str
    start: int
    end: int
    bins: int
    zoom_level: object
    aligned: bool


class Ours:
    """What our reader sees in one fixture, and the regions and queries derived from it."""

    def __init__(self, path):
        self.path = path
        self.info = dump("info", path)[0]
        self.chroms = [(r["name"], int(r["chrom_id"]), int(r["length"])) for r in dump("chroms", path)]
        self.lengths = {name: length for name, _, length in self.chroms}
        self.reductions = [int(r) for r in self.info["zoom_reductions"].split(",") if r]
        self.items = item_tuples(dump("items", path))
        self.items_by_chrom = {}
        for item in self.items:
            self.items_by_chrom.setdefault(item[0], []).append(item)
        self.big_endian = self.info["byte_order"] == "big"
        self.sampled = self._sample_chroms()
        self.empty_chrom = next((name for name, _, _ in self.chroms if name not in self.items_by_chrom), None)
        self.windows = {chrom: self._window(chrom) for chrom in self.sampled}
        self.queries = [q for chrom in self.sampled for q in self._queries(chrom)]
        # The regions every reader's raw data is read over: each sampled chromosome whole and its
        # window, and the chromosome without data.
        self.raw_regions = [(chrom, 0, self.lengths[chrom]) for chrom in self.sampled]
        self.raw_regions += [(chrom, *self.windows[chrom]) for chrom in self.sampled]
        if self.empty_chrom:
            self.raw_regions.append((self.empty_chrom, 0, self.lengths[self.empty_chrom]))

    def _sample_chroms(self):
        """Up to MAX_SAMPLED_CHROMS chromosomes with data: the first, middle and last by
        chromId."""
        with_data = [name for name, _, _ in self.chroms if name in self.items_by_chrom]
        if len(with_data) <= MAX_SAMPLED_CHROMS:
            return with_data
        return [with_data[0], with_data[len(with_data) // 2], with_data[-1]]

    def _window(self, chrom):
        """A region of at most WINDOW_LENGTH bases inside the data of `chrom`, starting a third
        of the way in at an offset that is on no grid."""
        items = self.items_by_chrom[chrom]
        first, last = items[0][1], items[-1][2]
        start = min(first + (last - first) // 3 + 7, last - 1)
        return start, min(start + WINDOW_LENGTH, self.lengths[chrom])

    def _queries(self, chrom):
        """One raw query, one grid-aligned query per zoom level and one unaligned zoom query
        on `chrom`."""
        length = self.lengths[chrom]
        window_start, window_end = self.windows[chrom]
        queries = []

        # Raw: bins at most twice the finest reduction wide, so kent's rule reads raw data.
        finest = min(self.reductions) if self.reductions else WINDOW_LENGTH
        raw_length = min(window_end - window_start, RAW_BINS * 2 * finest)
        if raw_length >= RAW_BINS:
            query = Query(chrom, window_start, window_start + raw_length, RAW_BINS, None, False)
            if kent_level(self.reductions, raw_length // RAW_BINS) is None:
                queries.append(query)

        # One aligned query per level, with bins a multiple of the reduction wide that kent and
        # libBigWig both resolve to this level, so their results are directly comparable.
        for level, reduction in enumerate(self.reductions):
            widths = [m * reduction for m in range(3, 9)
                      if kent_level(self.reductions, m * reduction) == level
                      and libbigwig_level(self.reductions, m * reduction) == level]
            if not widths:
                continue
            width = widths[0]
            start = window_start // reduction * reduction
            bins = min(MAX_ZOOM_BINS, (length - start) // width)
            if bins < 1:
                start = 0
                bins = min(MAX_ZOOM_BINS, length // width)
            if bins >= 1:
                queries.append(Query(chrom, start, start + bins * width, bins, level, True))

        # Unaligned: bins of 3.x times the finest reduction from an off-grid start, so records
        # straddle bin edges and kent pro-rates them.
        if self.reductions:
            width = 3 * self.reductions[0] + 17
            level = kent_level(self.reductions, width)
            bins = min(5, (length - window_start) // width)
            if level is not None and bins >= 1:
                queries.append(Query(chrom, window_start, window_start + bins * width, bins, level, False))
        return queries

    def expect_items(self, chrom, found):
        """Guards against vacuous passes: a region queried on a chromosome with data must hold
        something (every region is chosen inside the data), and one on the empty chromosome
        nothing."""
        if chrom == self.empty_chrom and found:
            raise Mismatch(f"our reader found {len(found)} items on {chrom}, which has no data")
        if chrom != self.empty_chrom and not found:
            raise Mismatch(f"no data in the region queried on {chrom}; nothing was compared")

    def summaries(self, query, source):
        """Our summaries for `query` from `source` ("raw" or "best-zoom")."""
        rows = dump("summaries", self.path, query.chrom, str(query.start), str(query.end),
                    str(query.bins), "--source", source)
        expected_source = "raw" if source == "raw" or query.zoom_level is None \
            else f"zoom:{self.reductions[query.zoom_level]}"
        for row in rows:
            if row["source"] != expected_source:
                raise Mismatch(f"our BestZoom read {row['source']}, kent's rule reads {expected_source}")
        return rows


def describe(query):
    """A short label for a summary query."""
    kind = "raw" if query.zoom_level is None else f"zoom{query.zoom_level}"
    grid = "" if query.zoom_level is None or query.aligned else "-unaligned"
    return f"{query.chrom}:{query.start}-{query.end}/{query.bins} {kind}{grid}"


# ---------------------------------------------------------------------------------------------
# kent
# ---------------------------------------------------------------------------------------------


def parse_info(text):
    """bigWigInfo (kent or bigtools) output with -chroms and -zooms: the key: value fields, the
    zoom reductions, and the chromosomes as (name, id, size)."""
    fields, zooms, chroms, section = {}, [], [], None
    for line in text.splitlines():
        if line.startswith("\t"):
            parts = line.split()
            if section == "zoomLevels":
                zooms.append(int(parts[0]))
            elif section == "chromCount":
                chroms.append((parts[0], int(parts[1]), int(parts[2])))
            continue
        key, _, value = line.partition(":")
        fields[key] = value.strip()
        section = key
    return fields, zooms, chroms


def check_info_common(ours, fields, zooms, chroms):
    """The checks bigWigInfo and bigtools' bigwiginfo share: byte order, zoom levels and
    chromosomes as written."""
    expected_swapped = "1" if ours.big_endian else "0"
    if fields.get("isSwapped") != expected_swapped:
        raise Mismatch(f"isSwapped {fields.get('isSwapped')}, expected {expected_swapped}")
    if int(fields["zoomLevels"]) != len(ours.reductions) or zooms != ours.reductions:
        raise Mismatch(f"zoom levels {zooms}, expected {ours.reductions}")
    if int(fields["chromCount"]) != len(ours.chroms) or sorted(chroms) != sorted(ours.chroms):
        raise Mismatch(f"{len(chroms)} chromosomes differ from our {len(ours.chroms)}")
    if int(fields["basesCovered"].replace(",", "")) != int(float(ours.info["bases_covered"])):
        raise Mismatch(f"basesCovered {fields['basesCovered']}, ours {ours.info['bases_covered']}")


def kent_info(ours, fixture):
    """kent bigWigInfo parses the file; compression, chromosomes, zoom levels and the total
    summary match our info()."""
    fields, zooms, chroms = parse_info(run([TOOLS["bigWigInfo"], "-chroms", "-zooms", ours.path]).stdout)
    check_info_common(ours, fields, zooms, chroms)
    # kent reports a file as compressed when uncompressBufSize > 0, which every compressed file
    # must declare, the empty one included.
    compressed = fixture not in UNCOMPRESSED_FIXTURES
    if fields["isCompressed"] != ("yes" if compressed else "no") \
            or (int(ours.info["uncompress_buf_size"]) > 0) != compressed:
        raise Mismatch(f"isCompressed {fields['isCompressed']}, uncompressBufSize "
                       f"{ours.info['uncompress_buf_size']}, written {'' if compressed else 'un'}compressed")
    for key in ("min", "max"):
        if fields[key] != "%f" % float(ours.info[key]):
            raise Mismatch(f"{key} {fields[key]}, ours {ours.info[key]}")
    mean, std = float(ours.info["mean"]), float(ours.info["std"])
    if not f_matches(fields["mean"], mean, RAW_REL_TOL):
        raise Mismatch(f"mean {fields['mean']}, ours {mean}")
    # With no bases covered kent's calcStdFromSums divides 0 by 0 and prints nan; our
    # Summary::std is documented as 0 when no more than one base is covered.
    if float(ours.info["bases_covered"]) == 0 and fields["std"] in ("nan", "-nan"):
        std = math.nan
    if not f_matches(fields["std"], std, RAW_REL_TOL):
        raise Mismatch(f"std {fields['std']}, ours {std}")
    return f"{len(chroms)} chroms, zooms {zooms}, {fields['basesCovered']} bases"


def compare_text_items(theirs, ours, value_matches, label):
    """Compares (chrom, start, end, value-text) rows with our (chrom, start, end, value) items,
    both sorted by chromosome and start, since tools list chromosomes in name order and we in
    chromId order."""
    theirs = sorted(theirs, key=lambda item: (item[0], item[1]))
    ours = sorted(ours, key=lambda item: (item[0], item[1]))
    if len(theirs) != len(ours):
        raise Mismatch(f"{len(theirs)} {label}, ours {len(ours)}")
    for their, our in zip(theirs, ours):
        if their[:3] != our[:3] or not value_matches(their[3], our[3]):
            raise Mismatch(f"{their} differs from our {our}")
    return f"{len(ours)} {label}"


def expect_data(fixture, items):
    """Fails when a fixture that has data produced nothing to compare."""
    if fixture not in EMPTY_FIXTURES and not items:
        raise Mismatch("nothing to compare")


def f32_of(value):
    """The f32 nearest `value`, a number or text such as a shortest round-trip printing of an f32
    (as `convert` writes values), as an exact double."""
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


def region_args(region, kent):
    """Command-line arguments restricting kent's tools (`kent`) or our dump to `region`, a
    (chrom, start, end) tuple or None for the whole file."""
    if region is None:
        return []
    chrom, start, end = region
    if kent:
        return [f"-chrom={chrom}", f"-start={start}", f"-end={end}"]
    return ["--chrom", chrom, "--start", str(start), "--end", str(end)]


def describe_region(region):
    """A short label for a region, or `whole file`."""
    return "whole file" if region is None else f"{region[0]}:{region[1]}-{region[2]}"


def kent_bedgraph(ours, fixture, region):
    """kent bigWigToBedGraph's output equals convert::bigwig_to_bedgraph with merge_equal over
    the whole file or a region (both clip to the region), with our shortest-form values printed
    as kent's %g would print them."""
    lines = run([TOOLS["bigWigToBedGraph"], *region_args(region, True), ours.path, "stdout"]).stdout
    theirs = [(c, int(s), int(e), v) for c, s, e, v in (line.split("\t") for line in lines.splitlines())]
    text = run([TOOLS["interop"], "dump", "bedgraph", ours.path, *region_args(region, False)]).stdout
    mine = [(c, int(s), int(e), f32_of(v)) for c, s, e, v in (line.split("\t") for line in text.splitlines())]
    expect_data(fixture, mine)
    detail = compare_text_items(theirs, mine, lambda text, value: text == "%g" % value, "runs")
    return f"{describe_region(region)}: {detail}"


@dataclass
class Wig:
    """Wig text parsed: the items with their block kind, the blocks, and the
    `#bedGraph section` comments with the extent of the lines that follow each."""

    items: list = field(default_factory=list)
    blocks: list = field(default_factory=list)
    comments: list = field(default_factory=list)


def parse_wig(text):
    """Wig text (kent's or ours) parsed. Items are (chrom, start, end, value-text, kind) with
    kind "bedGraph", "variableStep" or ("fixedStep", step); blocks are dicts with the kind,
    chrom, span, first start and item count; comments are [chrom, start, end, lines, extent]."""
    wig, block, comment = Wig(), None, None
    for line in text.splitlines():
        if line.startswith("#bedGraph section"):
            chrom, _, extent = line.split()[2].rpartition(":")
            start, end = extent.split("-")
            comment = [chrom, int(start), int(end), 0, None]
            wig.comments.append(comment)
            block = {"kind": "bedGraph", "count": 0}
            wig.blocks.append(block)
        elif line.startswith(("variableStep", "fixedStep")):
            keys = dict(word.split("=") for word in line.split()[1:])
            kind = "variableStep" if line.startswith("variableStep") else ("fixedStep", int(keys["step"]))
            block = {"kind": kind, "chrom": keys["chrom"], "span": int(keys["span"]), "count": 0,
                     "next": int(keys["start"]) - 1 if "start" in keys else None}
            wig.blocks.append(block)
            comment = None
        else:
            fields = line.split()
            if block["kind"] == "bedGraph":
                item = (fields[0], int(fields[1]), int(fields[2]), fields[3], "bedGraph")
                if comment is not None:
                    comment[3] += 1
                    first = comment[4][0] if comment[4] else item[1]
                    comment[4] = (first, item[2])
            elif block["kind"] == "variableStep":
                start = int(fields[0]) - 1
                item = (block["chrom"], start, start + block["span"], fields[1], "variableStep")
            else:
                start = block["next"]
                item = (block["chrom"], start, start + block["span"], fields[0], block["kind"])
                block["next"] += block["kind"][1]
            block.setdefault("first", item[1])
            block["last"] = item[1]
            block["count"] += 1
            wig.items.append(item)
    return wig


def check_our_wig_structure(wig):
    """convert::bigwig_to_wig's documented layout: no header without a value after it,
    consecutive sections of one fixedStep or variableStep run under one header, and each
    `#bedGraph section chrom:start-end` comment giving the extent of the lines after it."""
    previous = None
    for block in wig.blocks:
        if block["count"] == 0:
            raise Mismatch(f"empty {block['kind']} block")
        if previous is not None and previous["kind"] == block["kind"] and block["kind"] != "bedGraph" \
                and previous["chrom"] == block["chrom"] and previous["span"] == block["span"]:
            continues = block["kind"] == "variableStep" or \
                block["first"] == previous["last"] + block["kind"][1]
            if continues:
                raise Mismatch(f"a {block['kind']} run on {block['chrom']} at {block['first']} "
                               f"is split across headers")
        previous = block
    for chrom, start, end, lines, extent in wig.comments:
        if lines == 0 or extent != (start, end):
            raise Mismatch(f"#bedGraph section {chrom}:{start}-{end} heads lines spanning {extent}")


def clip_items(items, region):
    """Wig items clipped to `region`, dropping those outside it, as convert clips them."""
    if region is None:
        return items
    _, region_start, region_end = region
    clipped = []
    for chrom, start, end, value, kind in items:
        start, end = max(start, region_start), min(end, region_end)
        if start < end:
            clipped.append((chrom, start, end, value, kind))
    return clipped


def kent_wig(ours, fixture, region):
    """kent bigWigToWig and convert::bigwig_to_wig hold the same items, each in a block of the
    same kind (fixedStep with the same step), with our values printed as kent's %g would print
    them. The layouts differ: kent writes one header per stored section, including empty ones,
    and does not clip items to the region, while ours merges runs, writes no empty header and
    clips. Over the whole file both write one `#bedGraph section` comment per bedGraph section,
    with the same extents."""
    kent = parse_wig(run([TOOLS["bigWigToWig"], *region_args(region, True), ours.path, "stdout"]).stdout)
    mine = parse_wig(run([TOOLS["interop"], "dump", "wig", ours.path, *region_args(region, False)]).stdout)
    check_our_wig_structure(mine)
    theirs = sorted(clip_items(kent.items, region), key=lambda item: (item[0], item[1]))
    items = sorted(mine.items, key=lambda item: (item[0], item[1]))
    expect_data(fixture, items)
    if len(theirs) != len(items):
        raise Mismatch(f"{describe_region(region)}: kent {len(theirs)} items, ours {len(items)}")
    for their, my in zip(theirs, items):
        same = their[:3] == my[:3] and their[4] == my[4] and their[3] == "%g" % f32_of(my[3])
        if not same:
            raise Mismatch(f"{describe_region(region)}: kent {their} differs from our {my}")
    if region is None:
        kent_comments = sorted(tuple(c[:3]) for c in kent.comments)
        our_comments = sorted(tuple(c[:3]) for c in mine.comments)
        if kent_comments != our_comments:
            raise Mismatch(f"{len(kent_comments)} kent bedGraph comments differ from our {len(our_comments)}")
    return (f"{describe_region(region)}: {len(items)} items; {len(kent.blocks)} kent blocks, "
            f"{len(mine.blocks)} ours")


def variance(covered, total, total_squares):
    """The sample variance from sums, as kent's calcVarianceFromSums computes it: possibly
    slightly negative when rounding cancels the two terms."""
    value = total_squares - total * total / covered
    return value / (covered - 1) if covered > 1 else value


def variance_floor(covered, total_squares):
    """The absolute accuracy of a variance computed from sums: 1e-12 of the mean square.
    `sum_squares - sum²/n` cancels to the rounding noise of its terms when the values are
    nearly constant, so a standard deviation near zero is only known to about 1e-6 of the RMS
    value in any reader that uses sums (ours, kent's, libBigWig's zoom path)."""
    return 1e-12 * total_squares / covered


def std_matches(their_std, expected_std, expected_variance, floor, rel):
    """Whether a reported std matches: relatively close, or with variances within `floor` of
    each other. A NaN (kent prints n/a, pyBigWig None) is a variance rounded below zero, which
    matches when the expected variance is within `floor` of zero."""
    if their_std is None or math.isnan(their_std):
        return expected_variance <= floor
    return close(their_std, expected_std, rel) or abs(their_std ** 2 - expected_variance) <= floor


def kent_counts(ours, query, rows):
    """kent's pro-rated base count for each bin, summed in kent's order and arithmetic from our
    items (raw path, bbiIntervalSlice: size × (overlap / size) per item) or our zoom records
    (zoom path, bbiSummarySlice: validCount × overlap / extent per record)."""
    if query.zoom_level is None:
        # kent clips the items to the query region before pro-rating them (bwgQuery.c).
        items = dump("items", ours.path, *region_args((query.chrom, query.start, query.end), False))
        clipped = [(max(int(r["start"]), query.start), min(int(r["end"]), query.end)) for r in items]
        pieces = [(start, end, end - start) for start, end in clipped if start < end]
    else:
        records = dump("zoom", ours.path, str(query.zoom_level), query.chrom, str(query.start), str(query.end))
        pieces = [(int(r["start"]), int(r["end"]), float(r["covered_bases"])) for r in records]
    counts = []
    for row in rows:
        bin_start, bin_end, count = int(row["start"]), int(row["end"]), 0.0
        for start, end, weight in pieces:
            overlap = min(end, bin_end) - max(start, bin_start)
            if overlap > 0:
                count += weight * (overlap / (end - start))
        counts.append(count)
    return counts


def kent_expected(row, kent_count, statistic, query):
    """What kent's bigWigSummary should print for one of our summary bins, and for std the
    expected variance and its accuracy floor.

    kent rounds each bin's pro-rated base count up to a whole number and scales the sums to
    match (normalizeCount in lib/bbiRead.c). That changes zoom-based coverage and std where
    records straddle bin edges, and even raw coverage by one base where the count's rounding
    error lands just above a whole number, so the expectation repeats kent's arithmetic on our
    sums. kent's coverage divides by the mean bin width, not the bin's own width."""
    rounded = math.ceil(kent_count)
    scale = rounded / kent_count
    if statistic == "coverage":
        return rounded * query.bins / (query.end - query.start)
    if statistic == "std":
        total, total_squares = float(row["sum"]) * scale, float(row["sum_squares"]) * scale
        expected_variance = variance(rounded, total, total_squares)
        return (math.sqrt(max(expected_variance, 0.0)), expected_variance,
                variance_floor(rounded, total_squares))
    return float(row[statistic])


def kent_summary(ours, query, require_data=True):
    """kent bigWigSummary agrees with our summaries for every statistic, raw or from the zoom
    level kent's rule picks (our BestZoom), at kent's %g precision. With `require_data`, a
    query where no bin has data fails, since nothing was compared."""
    rows = ours.summaries(query, "best-zoom")
    counts = kent_counts(ours, query, rows)
    tolerance = RAW_REL_TOL if query.zoom_level is None else ZOOM_REL_TOL
    label = describe(query)
    for index, (row, count) in enumerate(zip(rows, counts)):
        covered = optional_float(row["covered_bases"])
        if covered is not None and not close(covered, count, 1e-9):
            raise Mismatch(f"{label} bin {index}: our {covered} bases, kent's arithmetic {count}")
    compared = 0
    for statistic in KENT_SUMMARY_TYPES:
        process = run([TOOLS["bigWigSummary"], f"-type={statistic}", ours.path, query.chrom,
                       str(query.start), str(query.end), str(query.bins)], check=False)
        if process.returncode != 0:
            if "no data in region" in process.stderr and all(r["covered_bases"] == "" for r in rows):
                continue
            raise Mismatch(f"{label}: bigWigSummary -type={statistic} exited {process.returncode}: "
                           f"{process.stderr.strip()}")
        printed = process.stdout.split()
        if len(printed) != query.bins:
            raise Mismatch(f"{label}: {len(printed)} {statistic} values for {query.bins} bins")
        for index, (text, row, count) in enumerate(zip(printed, rows, counts)):
            if row["covered_bases"] == "":
                if text != "n/a":
                    raise Mismatch(f"{label} bin {index} {statistic}: kent {text}, ours has no data")
                continue
            expected = kent_expected(row, count, statistic, query)
            if statistic == "std":
                expected_std, expected_variance, floor = expected
                # kent prints n/a for the NaN its unguarded sqrt gives a negative variance.
                their = math.nan if text == "n/a" else float(text)
                ok = std_matches(their, expected_std, expected_variance, floor, tolerance) or \
                    (text != "n/a" and g_matches(text, expected_std, tolerance))
            elif statistic in ("min", "max"):
                ok = text == "%g" % expected
            else:
                ok = g_matches(text, expected, tolerance)
            if not ok:
                raise Mismatch(f"{label} bin {index} {statistic}: kent {text}, ours {expected!r}")
            compared += 1
    if compared == 0 and require_data:
        raise Mismatch(f"{label}: no bin has data; nothing was compared")
    return f"{label}: {compared} values"


# ---------------------------------------------------------------------------------------------
# pyBigWig
# ---------------------------------------------------------------------------------------------


def pybigwig_checks(report, ours, fixture, required):
    """pyBigWig's chroms, header, intervals, values and stats, exact and approximate."""
    import pyBigWig

    handle = {}

    def open_file():
        handle["bw"] = pyBigWig.open(ours.path)
        if handle["bw"] is None:
            raise Mismatch("pyBigWig.open returned None")
        return "opened"

    if not report.run("pyBigWig", "open", fixture, required, open_file):
        return
    bw = handle["bw"]

    def chroms():
        theirs = bw.chroms()
        expected = {name: length for name, _, length in ours.chroms}
        if theirs != expected:
            raise Mismatch(f"{len(theirs)} chromosomes differ from our {len(expected)}")
        return f"{len(theirs)} chroms"

    def header():
        # pyBigWig truncates the header's floats to integers.
        theirs = bw.header()
        expected = {
            "nBasesCovered": int(float(ours.info["bases_covered"])),
            "minVal": int(float(ours.info["min"])),
            "maxVal": int(float(ours.info["max"])),
            "sumData": int(float(ours.info["sum"])),
            "sumSquared": int(float(ours.info["sum_squares"])),
            "nLevels": len(ours.reductions),
        }
        for key, value in expected.items():
            if theirs[key] != value:
                raise Mismatch(f"{key} {theirs[key]}, ours {value}")
        return "total summary"

    report.run("pyBigWig", "chroms", fixture, required, chroms)
    report.run("pyBigWig", "header", fixture, required, header)

    for chrom, start, end in ours.raw_regions:
        def intervals(chrom=chrom, start=start, end=end):
            theirs = [(chrom, s, e, v) for s, e, v in (bw.intervals(chrom, start, end) or ())]
            mine = item_tuples(dump("items", ours.path, "--chrom", chrom, "--start", str(start),
                                    "--end", str(end)))
            ours.expect_items(chrom, mine)
            if theirs != mine:
                raise Mismatch(f"{len(theirs)} intervals differ from our {len(mine)}")
            return f"{chrom}:{start}-{end}: {len(mine)} intervals"

        report.run("pyBigWig", "intervals", fixture, required, intervals)

    for chrom in ours.sampled:
        def values(chrom=chrom):
            start, end = ours.windows[chrom]
            theirs = bw.values(chrom, start, end)
            mine = [float(r["value"]) for r in dump("values", ours.path, chrom, str(start), str(end))]
            if len(theirs) != len(mine):
                raise Mismatch(f"{len(theirs)} values, ours {len(mine)}")
            if all(math.isnan(value) for value in mine):
                raise Mismatch(f"no data in {chrom}:{start}-{end}; nothing was compared")
            for offset, (their, my) in enumerate(zip(theirs, mine)):
                if not (their == my or (math.isnan(their) and math.isnan(my))):
                    raise Mismatch(f"{chrom}:{start + offset}: {their}, ours {my}")
            return f"{chrom}:{start}-{end}: {len(mine)} values"

        report.run("pyBigWig", "values", fixture, required, values)

    for query in ours.queries:
        report.run("pyBigWig", "stats-exact", fixture, required, pybigwig_stats, bw, ours, query, True)
        if query.aligned:
            report.run("pyBigWig", "stats-approx", fixture, required, pybigwig_stats, bw, ours, query, False)
    bw.close()


def pybigwig_stats(bw, ours, query, exact):
    """pyBigWig stats over a query's bins: exact against our raw summaries; approximate (on
    grid-aligned queries, where libBigWig and kent pick the same level) against our BestZoom."""
    rows = ours.summaries(query, "raw" if exact else "best-zoom")
    tolerance = RAW_REL_TOL if exact else ZOOM_REL_TOL
    compared = 0
    failed_std = False
    for statistic in PYBIGWIG_STATS_TYPES:
        try:
            theirs = bw.stats(query.chrom, query.start, query.end, type=statistic, nBins=query.bins, exact=exact)
        except RuntimeError:
            # libBigWig's zoom std takes the square root of a variance rounded below -1e-8
            # (bwStats.c). glibc's sqrt then sets errno, which libBigWig takes as a failure of
            # the whole call; elsewhere it gives NaN. Either is a variance near zero, so the call
            # may fail only when some bin's variance is within the floor of zero.
            if statistic == "std" and not exact and any(variance_near_zero(row) for row in rows):
                failed_std = True
                continue
            raise
        for index, (their, row) in enumerate(zip(theirs, rows)):
            if row["covered_bases"] == "":
                if their is not None:
                    raise Mismatch(f"bin {index} {statistic}: {their}, ours has no data")
                continue
            covered = float(row["covered_bases"])
            if statistic == "coverage":
                expected = covered / (int(row["end"]) - int(row["start"]))
            else:
                expected = float(row[statistic])
            if statistic == "std":
                # libBigWig's exact std is two-pass, so exactly 0 for constant values, while its
                # zoom std is NaN (None) for at most one covered base or a variance rounded
                # below -1e-8 (bwStats.c). Ours comes from sums, so near zero all agree only to
                # the variance floor.
                total, total_squares = float(row["sum"]), float(row["sum_squares"])
                ok = std_matches(their, expected, variance(covered, total, total_squares),
                                 variance_floor(covered, total_squares), tolerance)
            elif their is None:
                ok = False
            elif statistic in ("min", "max"):
                ok = their == expected
            else:
                ok = close(their, expected, tolerance)
            if not ok:
                raise Mismatch(f"{describe(query)} bin {index} {statistic}: {their!r}, ours {expected!r}")
            compared += 1
    if compared == 0:
        raise Mismatch(f"{describe(query)}: no bin has data; nothing was compared")
    if failed_std:
        return f"{describe(query)}: {compared} values; std failed on a variance near zero"
    return f"{describe(query)}: {compared} values"


def variance_near_zero(row):
    """Whether a summary row has data and a variance within the floor of zero."""
    if row["covered_bases"] == "":
        return False
    covered, total, total_squares = float(row["covered_bases"]), float(row["sum"]), float(row["sum_squares"])
    return variance(covered, total, total_squares) <= variance_floor(covered, total_squares)


# ---------------------------------------------------------------------------------------------
# bigtools
# ---------------------------------------------------------------------------------------------


def bigtools_info(ours):
    """bigtools bigwiginfo parses the file and agrees with our info() at full precision."""
    fields, zooms, chroms = parse_info(run([TOOLS["bigwiginfo"], "--chroms", "--zooms", ours.path]).stdout)
    check_info_common(ours, fields, zooms, chroms)
    for key in ("min", "max"):
        if float(fields[key]) != float(ours.info[key]):
            raise Mismatch(f"{key} {fields[key]}, ours {ours.info[key]}")
    for key in ("mean", "std"):
        theirs, mine = float(fields[key]), float(ours.info[key])
        # bigtools prints NaN for the std of an empty file, where ours is 0 (Summary::std).
        if key == "std" and float(ours.info["bases_covered"]) == 0 and math.isnan(theirs):
            continue
        if not close(theirs, mine, RAW_REL_TOL):
            raise Mismatch(f"{key} {theirs}, ours {mine}")
    return f"{len(chroms)} chroms, zooms {zooms}"


def bigtools_bedgraph(ours, fixture):
    """bigtools bigwigtobedgraph prints exactly our items (it does not merge them)."""
    lines = run([TOOLS["bigwigtobedgraph"], "-t", "1", ours.path, "-"]).stdout.splitlines()
    theirs = [(c, int(s), int(e), v) for c, s, e, v in (line.split("\t") for line in lines)]
    expect_data(fixture, ours.items)
    return compare_text_items(theirs, ours.items, same_f32, "items")


# ---------------------------------------------------------------------------------------------
# @gmod/bbi
# ---------------------------------------------------------------------------------------------


def bbi_checks(report, ours, fixture, required):
    """@gmod/bbi's header, raw getFeatures, and getFeatures at a scale picking each zoom level."""
    # A scale whose 2 / scale is half a base above a level's reduction makes bbi pick exactly
    # that level, since reductions are whole numbers.
    zoom_queries = [(chrom, level, *ours.windows[chrom], 2 / (reduction + 0.5))
                    for chrom in ours.sampled for level, reduction in enumerate(ours.reductions)]
    queries = [{"chrom": c, "start": s, "end": e} for c, s, e in ours.raw_regions]
    queries += [{"chrom": c, "start": s, "end": e, "scale": scale} for c, _, s, e, scale in zoom_queries]

    output = {}

    def header():
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as request:
            json.dump({"file": os.path.abspath(ours.path), "queries": queries}, request)
        try:
            process = run([TOOLS["js"], TOOLS["bbi_dump"], request.name])
        finally:
            os.unlink(request.name)
        output.update(json.loads(process.stdout))
        if "error" in output:
            raise Mismatch(output["error"])
        theirs = output["header"]
        if [tuple(r) for r in theirs["refsByNumber"]] != ours.chroms:
            raise Mismatch(f"{len(theirs['refsByNumber'])} chromosomes differ from ours")
        if theirs["zoomLevels"] != ours.reductions:
            raise Mismatch(f"zoom levels {theirs['zoomLevels']}, ours {ours.reductions}")
        summary = theirs["totalSummary"]
        pairs = [("basesCovered", "bases_covered"), ("scoreMin", "min"), ("scoreMax", "max"),
                 ("scoreSum", "sum"), ("scoreSumSquares", "sum_squares")]
        for their_key, our_key in pairs:
            if float(summary[their_key]) != float(ours.info[our_key]):
                raise Mismatch(f"{their_key} {summary[their_key]}, ours {ours.info[our_key]}")
        return f"{len(theirs['refsByNumber'])} chroms, zooms {theirs['zoomLevels']}"

    if not report.run("bbi", "header", fixture, required, header):
        return
    results = output["results"]

    for (chrom, start, end), result in zip(ours.raw_regions, results):
        def raw(chrom=chrom, start=start, end=end, result=result):
            if "error" in result:
                raise Mismatch(result["error"])
            theirs = sorted((chrom, s, e, v) for s, e, v in result["features"])
            mine = sorted(item_tuples(dump("items", ours.path, "--chrom", chrom, "--start", str(start),
                                           "--end", str(end))))
            ours.expect_items(chrom, mine)
            if theirs != mine:
                raise Mismatch(f"{len(theirs)} features differ from our {len(mine)} items")
            return f"{chrom}:{start}-{end}: {len(mine)} features"

        report.run("bbi", "getFeatures-raw", fixture, required, raw)

    for (chrom, level, start, end, scale), result in zip(zoom_queries, results[len(ours.raw_regions):]):
        def zoom(chrom=chrom, level=level, start=start, end=end, scale=scale, result=result):
            if bbi_level(ours.reductions, scale) != level:
                raise Mismatch(f"scale {scale} picks level {bbi_level(ours.reductions, scale)}")
            if "error" in result:
                raise Mismatch(result["error"])
            theirs = sorted(tuple(f) for f in result["features"])
            records = dump("zoom", ours.path, str(level), chrom, str(start), str(end))
            mine = sorted((int(r["start"]), int(r["end"]), float(r["sum"]) / float(r["covered_bases"]),
                           float(r["min"]), float(r["max"])) for r in records)
            ours.expect_items(chrom, mine)
            if theirs != mine:
                raise Mismatch(f"{len(theirs)} zoom features differ from our {len(mine)} records")
            return f"{chrom}:{start}-{end} level {level}: {len(mine)} records"

        report.run("bbi", "getFeatures-zoom", fixture, required, zoom)


def dropped(counts):
    """A note of how many records a check dropped as known differences, by reason, or ''."""
    notes = [f"{count} {reason}" for reason, count in counts.items() if count]
    return f" (dropped {', '.join(notes)})" if notes else ""


# ---------------------------------------------------------------------------------------------
# IGV desktop
# ---------------------------------------------------------------------------------------------


@dataclass
class IgvAnswer:
    """What IgvDump printed for one request: its rows, the `zoomsel` fields (bases per pixel,
    the reduction IGV picks or `none`, and for the whole-genome view the home chromosome), the
    chromosomes the whole-genome view spans, and its error, if any."""

    rows: list = field(default_factory=list)
    chosen: list = None
    wgchroms: list = field(default_factory=list)
    error: str = None


def parse_igv(text):
    """IgvDump's output: the header fields; the chromosomes as (name, id, the id IGV's lookup of
    the name finds); the error that stopped IGV opening the file, if any; and the IgvAnswer to
    each request, by index."""
    header, chroms, error, answers = None, [], None, {}
    for line in text.splitlines():
        kind, *fields = line.split("\t")
        if kind == "header":
            header = fields
        elif kind == "chrom":
            chroms.append((fields[1], int(fields[0]), fields[2]))
        elif kind == "error" and fields[0] == "open":
            error = fields[1]
        else:
            answer = answers.setdefault(int(fields[0]), IgvAnswer())
            if kind == "raw":
                answer.rows.append((int(fields[1]), int(fields[2]), float(fields[3])))
            elif kind == "zoom":
                answer.rows.append((fields[1], int(fields[2]), int(fields[3]), *map(float, fields[4:])))
            elif kind == "genome":
                answer.rows.append((int(fields[1]), int(fields[2]), *map(float, fields[3:])))
            elif kind == "zoomsel":
                answer.chosen = fields[1:]
            elif kind == "wgchrom":
                answer.wgchroms.append(fields[1])
            elif kind == "done":
                pass  # The answer is complete, and recorded even if it has no rows.
            elif kind == "error":
                answer.error = fields[1]
            else:
                raise ValueError(f"unexpected IgvDump line: {line}")
    return header, chroms, error, answers


def igv_checks(report, ours, fixture, required):
    """IGV desktop's reading of our file through its own classes (IgvDump.java): the header and
    chromosome lookups; getRawData; every zoom level's records; what IGV draws at each of its
    zooms, from the zoom level it picks itself; and the whole-genome view."""
    requests = [("raw", region) for region in ours.raw_regions]
    requests += [("level", (chrom, *ours.windows[chrom], reduction))
                 for chrom in ours.sampled for reduction in ours.reductions]
    # IGV shows a whole chromosome at zoom 0, and narrower regions deeper in.
    for chrom in ours.sampled:
        for zoom in range(IGV_MAX_ZOOM + 1):
            start, end = ours.windows[chrom] if zoom else (0, ours.lengths[chrom])
            requests.append(("zoom", (chrom, start, end, zoom)))
            if igv_level(ours.reductions, igv_bp_per_pixel(ours, chrom, zoom), 2) is None:
                break
    requests.append(("genome", ()))
    output = {}

    def header():
        with tempfile.TemporaryDirectory() as work:
            sizes, requests_path = os.path.join(work, "chrom.sizes"), os.path.join(work, "requests.tsv")
            with open(sizes, "w") as out:
                out.writelines(f"{name}\t{length}\n" for name, _, length in ours.chroms)
            with open(requests_path, "w") as out:
                out.writelines("\t".join(map(str, (kind, *args))) + "\n" for kind, args in requests)
            # A home directory of its own, so that nothing IGV might write lands in the user's.
            process = run([TOOLS["java"], "-Djava.awt.headless=true", f"-Duser.home={work}", "-cp",
                           TOOLS["igv_classpath"], "IgvDump", ours.path, sizes, requests_path])
        output["header"], output["chroms"], error, output["answers"] = parse_igv(process.stdout)
        if error:
            raise Mismatch(error)
        version, reductions, buffer_size, *summary = output["header"]
        if int(version) != int(ours.info["version"]):
            raise Mismatch(f"version {version}, ours {ours.info['version']}")
        if int(buffer_size) != int(ours.info["uncompress_buf_size"]):
            raise Mismatch(f"uncompressBufSize {buffer_size}, ours {ours.info['uncompress_buf_size']}")
        # IGV orders the zoom levels coarsest first.
        theirs = [int(r) for r in reductions.split(",") if r]
        if theirs != sorted(ours.reductions, reverse=True):
            raise Mismatch(f"zoom levels {theirs}, ours {ours.reductions}")
        for value, key in zip(summary, ("bases_covered", "min", "max", "sum", "sum_squares")):
            if float(value) != float(ours.info[key]):
                raise Mismatch(f"total summary {key} {value}, ours {ours.info[key]}")
        return f"version {version}, zooms {theirs}"

    if not report.run("igv", "header", fixture, required, header):
        return
    chrom_ids = {name: chrom_id for name, chrom_id, _ in ours.chroms}

    def chroms():
        theirs = [(name, chrom_id) for name, chrom_id, _ in output["chroms"]]
        if theirs != [(name, chrom_id) for name, chrom_id, _ in ours.chroms]:
            raise Mismatch(f"{len(theirs)} chromosomes differ from our {len(ours.chroms)}")
        wrong = [(name, found) for name, chrom_id, found in output["chroms"] if found != str(chrom_id)]
        if wrong:
            raise Mismatch(f"looking up {len(wrong)} names finds the wrong id, first {wrong[0]}")
        return f"{len(theirs)} chroms, each found by name"

    report.run("igv", "chroms", fixture, required, chroms)

    for index, (kind, args) in enumerate(requests):
        answer = output["answers"].get(index, IgvAnswer(error="no output"))
        if kind == "raw":
            report.run("igv", "raw", fixture, required, igv_raw, ours, *args, answer)
        elif kind == "level":
            report.run("igv", "zoom-level", fixture, required, igv_zoom_level, ours, chrom_ids, *args, answer)
        elif kind == "zoom":
            report.run("igv", "zoom", fixture, required, igv_zoom, ours, chrom_ids, *args, answer)
        else:
            report.run("igv", "whole-genome", fixture, required, igv_genome, ours, fixture, answer)


def igv_bp_per_pixel(ours, chrom, zoom):
    """The bases per pixel IGV reads zoom data for at IGV zoom `zoom` on `chrom`
    (BBDataSource.getPrecomputedSummaryScores)."""
    return ours.lengths[chrom] / (2 ** zoom * 700)


def igv_raw(ours, chrom, start, end, answer):
    """getRawData's items in a region are ours exactly, but for items ending at its start, which
    IGV returns: BBFile.decodeWigData skips only items ending before the start."""
    if answer.error:
        raise Mismatch(answer.error)
    theirs = sorted((chrom, s, e, v) for s, e, v in answer.rows if e != start)
    mine = sorted(item_tuples(dump("items", ours.path, "--chrom", chrom, "--start", str(start),
                                   "--end", str(end))))
    ours.expect_items(chrom, mine)
    if theirs != mine:
        raise Mismatch(f"{chrom}:{start}-{end}: {len(theirs)} items differ from our {len(mine)}")
    touching = len(answer.rows) - len(theirs)
    return f"{chrom}:{start}-{end}: {len(mine)} items{dropped({'ending at the start': touching})}"


def igv_zoom_records(answer, chrom, chrom_id, start, end):
    """IGV's zoom records for a region on `chrom`, as sorted (chrom, start, end, min, mean, max),
    without the records BBFile.decodeZoomData returns outside the region, and how many there were:
    - it skips only records ending before the start, so it returns those ending at it;
    - it filters by region only `if (chrIdx > 0)`, so for chromId 0 it returns every record of the
      blocks it reads, other chromosomes' included, since zoom blocks span chromosomes.
    Any other record outside the region is kept, and fails the comparison."""
    if answer.error:
        raise Mismatch(answer.error)
    kept, touching, unfiltered = [], 0, 0
    for contig, s, e, *values in answer.rows:
        if contig == chrom and e == start:
            touching += 1
        elif chrom_id == 0 and (contig != chrom or e < start or s >= end):
            unfiltered += 1
        else:
            kept.append((contig, s, e, *values))
    return sorted(kept), {"ending at the start": touching, "outside the region on chromId 0": unfiltered}


def igv_expected_records(ours, level, chrom, start, end):
    """Our records of zoom level `level` in a region as IGV shows them, as sorted (chrom, start,
    end, min, mean, max). IGV's mean is the f32 quotient sumData / validCount, the count converted
    to f32 first (BBFile.decodeZoomData); dividing in f64 and rounding to f32 gives exactly that."""
    return sorted((chrom, int(r["start"]), int(r["end"]), float(r["min"]),
                   f32_of(float(r["sum"]) / f32_of(r["covered_bases"])), float(r["max"]))
                  for r in dump("zoom", ours.path, str(level), chrom, str(start), str(end)))


def igv_zoom_level(ours, chrom_ids, chrom, start, end, reduction, answer):
    """One zoom level's records in a region, read as IGV reads the level it picks."""
    theirs, counts = igv_zoom_records(answer, chrom, chrom_ids[chrom], start, end)
    mine = igv_expected_records(ours, ours.reductions.index(reduction), chrom, start, end)
    ours.expect_items(chrom, mine)
    if theirs != mine:
        raise Mismatch(f"{chrom}:{start}-{end} reduction {reduction}: {len(theirs)} records differ "
                       f"from our {len(mine)}")
    return f"{chrom}:{start}-{end} reduction {reduction}: {len(mine)} records{dropped(counts)}"


def igv_zoom(ours, chrom_ids, chrom, start, end, zoom, answer):
    """What IGV draws at IGV zoom `zoom`: the records of the zoom level its rule picks for that
    resolution, or nothing where no level is fine enough."""
    if answer.error:
        raise Mismatch(answer.error)
    bp_per_pixel = igv_bp_per_pixel(ours, chrom, zoom)
    level = igv_level(ours.reductions, bp_per_pixel, 2)
    expected = "none" if level is None else str(ours.reductions[level])
    if answer.chosen[1] != expected:
        raise Mismatch(f"{chrom} zoom {zoom} ({bp_per_pixel:.1f} bp/pixel): IGV picks {answer.chosen[1]}, "
                       f"expected {expected}")
    if level is None:
        if answer.rows:
            raise Mismatch(f"{chrom} zoom {zoom}: no zoom level, yet {len(answer.rows)} records")
        return f"{chrom} zoom {zoom}: {bp_per_pixel:.1f} bp/pixel, where IGV reads no zoom level"
    theirs, counts = igv_zoom_records(answer, chrom, chrom_ids[chrom], start, end)
    mine = igv_expected_records(ours, level, chrom, start, end)
    ours.expect_items(chrom, mine)
    if theirs != mine:
        raise Mismatch(f"{chrom}:{start}-{end} zoom {zoom}: {len(theirs)} records differ from our "
                       f"{len(mine)}")
    return f"{chrom}:{start}-{end} zoom {zoom}: reduction {expected}, {len(mine)} records{dropped(counts)}"


def igv_genome(ours, fixture, answer):
    """The whole-genome view: the records of the zoom level IGV's rule picks for its resolution,
    on every chromosome it spans, at genome coordinates in kbp (getWholeGenomeScores). IGV has no
    such view of a genome with only one long chromosome."""
    if answer.error:
        raise Mismatch(answer.error)
    bp_per_pixel, chosen, home = answer.chosen
    if home != "All":
        if answer.rows:
            raise Mismatch(f"home chromosome {home}, yet {len(answer.rows)} whole-genome records")
        return f"home chromosome {home}: no whole-genome view"
    level = igv_level(ours.reductions, float(bp_per_pixel), 1000)
    expected = "none" if level is None else str(ours.reductions[level])
    if chosen != expected:
        raise Mismatch(f"{float(bp_per_pixel):.0f} bp/pixel: IGV picks {chosen}, expected {expected}")
    mine, offset = [], 0
    for name in answer.wgchroms:
        if level is not None and name in ours.items_by_chrom:
            for _, s, e, *values in igv_expected_records(ours, level, name, 0, ours.lengths[name]):
                mine.append(((offset + s) // 1000, (offset + e) // 1000, *values))
        offset += ours.lengths[name]
    expect_data(fixture, mine)
    if sorted(answer.rows) != sorted(mine):
        raise Mismatch(f"{len(answer.rows)} records differ from our {len(mine)}")
    return (f"{len(answer.wgchroms)} chromosomes, {float(bp_per_pixel):.0f} bp/pixel: reduction {expected}, "
            f"{len(mine)} records")


# ---------------------------------------------------------------------------------------------
# igv.js
# ---------------------------------------------------------------------------------------------


def igvjs_checks(report, ours, fixture, required):
    """igv.js's reading of our file through its own reader (igvjs_dump.mjs): the header, the
    chromosome lookups, raw features, and features at a resolution picking each zoom level, over
    each window and each whole sampled chromosome, where a query spans several blocks."""
    # bpPerPixel half a base above a level's reduction makes igv.js pick exactly that level, since
    # reductions are whole numbers.
    zoom_queries = [(chrom, level, start, end, reduction + 0.5)
                    for chrom in ours.sampled
                    for start, end in (ours.windows[chrom], (0, ours.lengths[chrom]))
                    for level, reduction in enumerate(ours.reductions)]
    queries = [{"chrom": c, "start": s, "end": e} for c, s, e in ours.raw_regions]
    queries += [{"chrom": c, "start": s, "end": e, "bpPerPixel": b} for c, _, s, e, b in zoom_queries]
    names = [name for name, _, _ in ours.chroms]
    output = {}

    def header():
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as request:
            json.dump({"file": os.path.abspath(ours.path), "chroms": names, "queries": queries}, request)
        try:
            process = run([TOOLS["js"], TOOLS["igvjs_dump"], TOOLS["igvjs_dir"], request.name])
        finally:
            os.unlink(request.name)
        output.update(json.loads(process.stdout))
        if "error" in output:
            raise Mismatch(output["error"])
        theirs = output["header"]
        # igv.js leaves the type unset for a magic number it does not recognise.
        if theirs.get("type") != "bigwig":
            raise Mismatch(f"type {theirs.get('type')}: igv.js does not recognise the magic number")
        if theirs["littleEndian"] == ours.big_endian:
            raise Mismatch(f"littleEndian {theirs['littleEndian']}, ours is {ours.info['byte_order']}-endian")
        if theirs["version"] != int(ours.info["version"]):
            raise Mismatch(f"version {theirs['version']}, ours {ours.info['version']}")
        if theirs["uncompressBuffSize"] != int(ours.info["uncompress_buf_size"]):
            raise Mismatch(f"uncompressBuffSize {theirs['uncompressBuffSize']}, "
                           f"ours {ours.info['uncompress_buf_size']}")
        if theirs["chromCount"] != len(ours.chroms):
            raise Mismatch(f"{theirs['chromCount']} chromosomes, ours {len(ours.chroms)}")
        # igv.js orders the zoom levels coarsest first.
        if theirs["zoomLevels"] != sorted(ours.reductions, reverse=True):
            raise Mismatch(f"zoom levels {theirs['zoomLevels']}, ours {ours.reductions}")
        summary = theirs["totalSummary"]
        pairs = [("basesCovered", "bases_covered"), ("minVal", "min"), ("maxVal", "max"),
                 ("sumData", "sum"), ("sumSquares", "sum_squares")]
        for their_key, our_key in pairs:
            if float(summary[their_key]) != float(ours.info[our_key]):
                raise Mismatch(f"{their_key} {summary[their_key]}, ours {ours.info[our_key]}")
        return f"{theirs['chromCount']} chroms, zooms {theirs['zoomLevels']}"

    if not report.run("igvjs", "header", fixture, required, header):
        return
    found_ids = dict(zip(names, output["chromIds"]))
    our_ids = {name: chrom_id for name, chrom_id, _ in ours.chroms}

    def chrom_ids():
        # The dumper looks every name up before any feature query, since after one a name the
        # tree search misses can still resolve (see igvjs_dump.mjs).
        if len(output["chromIds"]) != len(names) + 1:
            raise Mismatch(f"{len(output['chromIds'])} lookups for {len(names)} names and an absent one")
        wrong = [(name, our_ids[name], found) for name, found in found_ids.items() if found != our_ids[name]]
        if wrong:
            raise Mismatch(f"{len(wrong)} of {len(names)} names resolve wrongly, first {wrong[:3]}")
        if output["chromIds"][-1] is not None:
            raise Mismatch(f"a name the file does not hold resolves to {output['chromIds'][-1]}")
        return f"all {len(names)} names resolve to our ids, and an absent name to none"

    report.run("igvjs", "chrom-ids", fixture, required, chrom_ids)
    results = output["results"]

    for (chrom, start, end), result in zip(ours.raw_regions, results):
        def raw(chrom=chrom, start=start, end=end, result=result):
            # igv.js returns no features for a name it cannot find, as for a region without data,
            # so the chromosome without data must be found for its empty result to mean anything.
            if found_ids[chrom] != our_ids[chrom]:
                raise Mismatch(f"igv.js resolves {chrom} to {found_ids[chrom]}, not our {our_ids[chrom]}")
            features, touching = igvjs_features(result, chrom, start)
            theirs = sorted((chrom, *f) for f in features)
            mine = sorted(item_tuples(dump("items", ours.path, "--chrom", chrom, "--start", str(start),
                                           "--end", str(end))))
            ours.expect_items(chrom, mine)
            if theirs != mine:
                raise Mismatch(f"{chrom}:{start}-{end}: {len(theirs)} features differ from our "
                               f"{len(mine)} items")
            return f"{chrom}:{start}-{end}: {len(mine)} features{dropped({'ending at the start': touching})}"

        report.run("igvjs", "raw", fixture, required, raw)

    for (chrom, level, start, end, bp_per_pixel), result in zip(zoom_queries, results[len(ours.raw_regions):]):
        def zoom(chrom=chrom, level=level, start=start, end=end, bp_per_pixel=bp_per_pixel, result=result):
            picked = igvjs_level(ours.reductions, bp_per_pixel)
            if picked != level:
                raise Mismatch(f"bpPerPixel {bp_per_pixel} picks level {picked}")
            features, touching = igvjs_features(result, chrom, start)
            theirs = sorted(features)
            records = dump("zoom", ours.path, str(level), chrom, str(start), str(end))
            mine = sorted((int(r["start"]), int(r["end"]), float(r["sum"]) / float(r["covered_bases"]),
                           float(r["min"]), float(r["max"])) for r in records)
            ours.expect_items(chrom, mine)
            if theirs != mine:
                raise Mismatch(f"{chrom}:{start}-{end} level {level}: {len(theirs)} zoom features differ "
                               f"from our {len(mine)} records")
            note = dropped({"ending at the start": touching})
            return f"{chrom}:{start}-{end} level {level}: {len(mine)} records{note}"

        report.run("igvjs", "zoom", fixture, required, zoom)


def igvjs_features(result, chrom, start):
    """igv.js's features for a query on `chrom` that starts at `start`, as tuples without the
    chromosome, which must be `chrom`, and without those ending at `start`; and how many of
    those there were. igv.js returns them because its decoders skip only features ending before
    the start."""
    if "error" in result:
        raise Mismatch(result["error"])
    others = {f[0] for f in result["features"]} - {chrom}
    if others:
        raise Mismatch(f"features on {sorted(others)[:3]} from a query on {chrom}")
    features = [tuple(f[1:]) for f in result["features"]]
    kept = [f for f in features if f[1] != start]
    return kept, len(features) - len(kept)


# ---------------------------------------------------------------------------------------------
# Their files, read by us
# ---------------------------------------------------------------------------------------------


def read_chrom_sizes(path):
    """chrom.sizes as a list of (name, length), in file order."""
    with open(path) as sizes:
        return [(name, int(length)) for name, length in (line.split() for line in sizes)]


def parse_wig_blocks(path):
    """Our generated wig (or bedGraph) text as blocks: ("fixed", chrom, start, step, span,
    values), ("variable", chrom, span, starts, values) or ("bedgraph", chrom, starts, ends,
    values), with 0-based starts. bedGraph lines are grouped while the chromosome is
    unchanged."""
    blocks = []
    with open(path) as text:
        for line in text:
            fields = line.split()
            if fields[0] in ("fixedStep", "variableStep"):
                keys = dict(word.split("=") for word in fields[1:])
                if fields[0] == "fixedStep":
                    blocks.append(("fixed", keys["chrom"], int(keys["start"]) - 1, int(keys["step"]),
                                   int(keys["span"]), []))
                else:
                    blocks.append(("variable", keys["chrom"], int(keys["span"]), [], []))
            elif len(fields) == 4:
                if not blocks or blocks[-1][0] != "bedgraph" or blocks[-1][1] != fields[0]:
                    blocks.append(("bedgraph", fields[0], [], [], []))
                blocks[-1][2].append(int(fields[1]))
                blocks[-1][3].append(int(fields[2]))
                blocks[-1][4].append(float(fields[3]))
            elif blocks[-1][0] == "fixed":
                blocks[-1][5].append(float(fields[0]))
            else:
                blocks[-1][3].append(int(fields[0]) - 1)
                blocks[-1][4].append(float(fields[1]))
    return blocks


def write_with_pybigwig(text_path, sizes, out_path):
    """Writes `text_path` with pyBigWig addEntries, one call per block, in the form matching
    the block: fixedStep, variableStep or bedGraph."""
    import pyBigWig

    bw = pyBigWig.open(out_path, "w")
    bw.addHeader(sizes)
    for block in parse_wig_blocks(text_path):
        if block[0] == "fixed":
            _, chrom, start, step, span, values = block
            bw.addEntries(chrom, start, values=values, span=span, step=step)
        elif block[0] == "variable":
            _, chrom, span, starts, values = block
            bw.addEntries(chrom, starts, values=values, span=span)
        else:
            _, chrom, starts, ends, values = block
            bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
    bw.close()


def their_items(bigwig, expected_path, sizes):
    """Our reader reproduces exactly the items the text input holds, and the chromosomes it
    lists have their chrom.sizes lengths."""
    with open(expected_path) as expected_file:
        expected = item_tuples(csv.DictReader(expected_file, delimiter="\t"))
    lengths = dict(sizes)
    for row in dump("chroms", bigwig):
        if lengths.get(row["name"]) != int(row["length"]):
            raise Mismatch(f"chromosome {row['name']} has length {row['length']}")
    mine = item_tuples(dump("items", bigwig))
    expected.sort(key=lambda item: (item[0], item[1]))
    mine.sort(key=lambda item: (item[0], item[1]))
    if len(mine) != len(expected):
        raise Mismatch(f"{len(mine)} items, input has {len(expected)}")
    for my, their in zip(mine, expected):
        if my != their:
            raise Mismatch(f"our {my} differs from input {their}")
    return f"{len(mine)} items"


def their_zooms(bigwig):
    """Every zoom level reads, and over each whole chromosome its records' min and max equal
    the raw data's, which holds for every writer (kent's wigToBigWig over-counts coverage,
    libBigWig zeroes some sums and bigtools emits phantom records over gaps, but none of those
    changes min or max)."""
    reductions = [int(r) for r in dump("info", bigwig)[0]["zoom_reductions"].split(",") if r]
    by_chrom = {}
    for item in item_tuples(dump("items", bigwig)):
        by_chrom.setdefault(item[0], []).append(item[3])
    lengths = {r["name"]: int(r["length"]) for r in dump("chroms", bigwig)}
    for level in range(len(reductions)):
        for chrom, values in by_chrom.items():
            records = dump("zoom", bigwig, str(level), chrom, "0", str(lengths[chrom]))
            if not records:
                raise Mismatch(f"level {level} has no records on {chrom}")
            low, high = min(float(r["min"]) for r in records), max(float(r["max"]) for r in records)
            if (low, high) != (min(values), max(values)):
                raise Mismatch(f"level {level} {chrom}: min/max {low}/{high}, raw {min(values)}/{max(values)}")
    return f"zooms {reductions}" if reductions else "no zoom levels"


def their_files(report, text_dir, work_dir):
    """Converts the text inputs with every other writer and checks what our reader sees."""
    sizes_path = os.path.join(text_dir, "chrom.sizes")
    sizes = read_chrom_sizes(sizes_path)
    os.makedirs(work_dir, exist_ok=True)

    def text(name):
        return os.path.join(text_dir, name)

    def out(name):
        return os.path.join(work_dir, name)

    writers = [
        ("kent", "bedGraphToBigWig", "input", out("kent_bedgraph.bw"),
         lambda o: run([TOOLS["bedGraphToBigWig"], text("input.bedGraph"), sizes_path, o])),
    ]
    for name in ("fixed", "variable", "mixed"):
        writers.append(("kent", f"wigToBigWig-{name}", name, out(f"kent_{name}.bw"),
                        lambda o, name=name: run([TOOLS["wigToBigWig"], text(f"{name}.wig"), sizes_path, o])))
    # The text lists chromosomes in BAM order; bigtools' default requires byte order, so tell it
    # only starts are sorted.
    writers.append(("bigtools", "bedgraphtobigwig", "input", out("bigtools_bedgraph.bw"),
                    lambda o: run([TOOLS["bedgraphtobigwig"], "--sorted", "start", text("input.bedGraph"),
                                   sizes_path, o])))
    for name, source in (("bedgraph", "input.bedGraph"), ("variable", "variable.wig"),
                         ("fixed", "fixed.wig"), ("mixed", "mixed.wig")):
        expected = "input" if name == "bedgraph" else name
        writers.append(("pyBigWig", f"addEntries-{name}", expected, out(f"pybigwig_{name}.bw"),
                        lambda o, source=source: write_with_pybigwig(text(source), sizes, o)))

    for tool_name, writer, expected, bigwig, write in writers:
        fixture = os.path.basename(bigwig)

        def convert(write=write, bigwig=bigwig):
            write(bigwig)
            return f"wrote {os.path.getsize(bigwig)} bytes"

        if report.run(tool_name, writer, fixture, True, convert):
            report.run("ours", "read-items", fixture, True, their_items, bigwig,
                       text(f"{expected}.items.tsv"), sizes)
            report.run("ours", "read-zooms", fixture, True, their_zooms, bigwig)


# ---------------------------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------------------------


def our_files(report, fixture_dir, selected):
    """Every check of our fixtures by the other readers, or of those named in `selected`."""
    with open(os.path.join(fixture_dir, "manifest.tsv")) as manifest:
        fixtures = list(csv.DictReader(manifest, delimiter="\t"))
    for entry in fixtures:
        if selected and entry["fixture"] not in selected:
            continue
        fixture, path = entry["fixture"], os.path.join(fixture_dir, entry["file"])
        ours = Ours(path)
        # kent and our reader are required everywhere. The other readers are not required to read
        # the big-endian file: pyBigWig and @gmod/bbi read only little-endian files, and IGV and
        # igv.js misread big-endian ones (IGV's BBFile.loadBinaryBuffer ignores the byte order it
        # is given; igv.js compares the magic number with a mistyped constant, igvteam/igv.js#786).
        required = fixture not in BIG_ENDIAN_FIXTURES

        report.run("kent", "bigWigInfo", fixture, True, kent_info, ours, fixture)
        regions = [None] + [(chrom, *ours.windows[chrom]) for chrom in ours.sampled]
        for region in regions:
            report.run("kent", "bigWigToBedGraph", fixture, True, kent_bedgraph, ours, fixture, region)
            report.run("kent", "bigWigToWig", fixture, True, kent_wig, ours, fixture, region)
        for query in ours.queries:
            report.run("kent", "bigWigSummary", fixture, True, kent_summary, ours, query)
        if fixture in EMPTY_FIXTURES:
            # kent must fail with "no data in region" where we report no data.
            chrom, _, length = ours.chroms[0]
            level = kent_level(ours.reductions, length // RAW_BINS)
            report.run("kent", "bigWigSummary", fixture, True, kent_summary, ours,
                       Query(chrom, 0, length, RAW_BINS, level, False), False)

        pybigwig_checks(report, ours, fixture, required)
        report.run("bigtools", "bigwiginfo", fixture, required, bigtools_info, ours)
        report.run("bigtools", "bigwigtobedgraph", fixture, required, bigtools_bedgraph, ours, fixture)
        bbi_checks(report, ours, fixture, required)
        igv_checks(report, ours, fixture, required)
        # Nor is igv.js required to read uncompressed files: it decodes every block of a query
        # that reads several from the start of the buffer holding them all.
        igvjs_checks(report, ours, fixture, required and fixture not in UNCOMPRESSED_FIXTURES)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    fixture_dir, selected = sys.argv[1], set(sys.argv[2:])
    report = Report()
    our_files(report, fixture_dir, selected)
    if not selected or "theirs" in selected:
        their_files(report, os.path.join(fixture_dir, "text"), os.path.join(fixture_dir, "theirs"))
    sys.exit(report.summarise())


if __name__ == "__main__":
    main()
