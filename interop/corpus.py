#!/usr/bin/env python3
"""One-time validation against public bigWig files, run by corpus.sh.

Usage: corpus.py <manifest.tsv> <download-dir> <filled-manifest.tsv> <results.md>
                 [--annotations <tsv>]...
       corpus.py --fetch-only <manifest.tsv> <download-dir> <filled-manifest.tsv>

The manifest (interop/corpus.tsv) is read by column name. It must have the columns name, url,
sha256, size, source_issue and expectation; writer, when present, names the claimed writer in the
results, and every other column is carried into the filled manifest unchanged. A row whose
expectation starts with `skip:` is recorded as skipped and not fetched.

Every other row is fetched to <download-dir>/<name>/<the URL's file name>, unless it is already
there. A file under CORPUS_LOCAL_DIRS (colon-separated) whose name ends with the URL's file name
and whose sha256 is the manifest's is hard-linked instead (copied across file systems). Downloads
are capped at CORPUS_MAX_FILE_BYTES per file (default 1 GiB) and CORPUS_BUDGET_BYTES in total
(default 4 GB); `none` lifts either cap, and a row whose size is over a cap is skipped. Whatever is
fetched must match the manifest's sha256, where one is given.

With --fetch-only, that is all: the filled manifest is written, and the exit status is 1 only if a
checksum did not match. Otherwise zip, gzip and tar downloads are unpacked, and for each bigWig:

- `scan`: it opens, every item is read, the total summary is recomputed, and the items are
  rewritten with our writer and read back identically;
- on random regions holding data: kent bigWigToBedGraph and bigWigToWig against our convert
  output; our `values` against kent's bedGraph; our raw summaries against kent's bedGraph; and
  kent bigWigSummary against our summaries on a raw-path and a zoom-path query, with compare.py's
  tolerances and encoded kent quirks;
- nothing panics (exit status 101 from our binary is reported as a panic).

The filled manifest (sha256 and size of what was fetched) and a Markdown table of results are
written. Annotation TSVs (columns url plus writer or writer_if_known) override the manifest's
claimed writers.
"""

import collections
import csv
import gzip
import hashlib
import os
import re
import shutil
import struct
import subprocess
import sys
import tarfile
import zipfile

import compare
from compare import Mismatch, Query, kent_bedgraph, kent_level, kent_summary, kent_wig, run

# Manifest columns corpus.py needs; any others are carried through.
REQUIRED_COLUMNS = ("name", "url", "sha256", "size", "source_issue", "expectation")
# A row's name is a directory under <download-dir>.
NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._+-]*")
# Default largest file fetched, and default total download budget.
MAX_FILE_BYTES = 1 << 30
DEFAULT_BUDGET_BYTES = 4 * 10**9
# bigWig magic in each byte order, and other formats worth naming.
MAGICS = {b"\x26\xfc\x8f\x88": "bigWig", b"\x88\x8f\xfc\x26": "bigWig (big-endian)",
          b"\xeb\xf2\x89\x87": "bigBed", b"\x87\x89\xf2\xeb": "bigBed (big-endian)"}
# Random regions with data checked per file, and their sizes.
REGIONS_PER_FILE = 4
VALUES_WINDOW = 20_000
SUMMARY_WINDOW = 200_000
SUMMARY_BINS = 10
MAX_RAW_BINS = 1_000
# Seconds allowed for one scan (which rewrites the whole file) and for any other command.
SCAN_TIMEOUT = 3_600


class CorpusFile:
    """One bigWig under test, with what compare.py's kent checks need of it."""

    summaries = compare.Ours.summaries

    def __init__(self, path):
        self.path = path
        self.info = compare.dump("info", path)[0]
        self.chroms = [(r["name"], int(r["length"])) for r in compare.dump("chroms", path)]
        self.reductions = [int(r) for r in self.info["zoom_reductions"].split(",") if r]


class Limits:
    """The download caps, in bytes, from CORPUS_MAX_FILE_BYTES and CORPUS_BUDGET_BYTES. None is
    no cap; the budget shrinks as files are downloaded."""

    def __init__(self):
        self.max_file = byte_cap("CORPUS_MAX_FILE_BYTES", MAX_FILE_BYTES)
        self.budget = byte_cap("CORPUS_BUDGET_BYTES", DEFAULT_BUDGET_BYTES)

    def refusal(self, size):
        """Why a file of `size` bytes (None if unknown) may not be downloaded, or None."""
        if size is not None and self.max_file is not None and size > self.max_file:
            return f"skip: {size} bytes is over the per-file cap of {self.max_file} (CORPUS_MAX_FILE_BYTES)"
        if self.budget is not None and (self.budget <= 0 or (size or 0) > self.budget):
            return f"skip: over what is left of the total budget, {max(self.budget, 0)} bytes (CORPUS_BUDGET_BYTES)"
        return None

    def curl_cap(self):
        """The most bytes one download may take, or None for no cap."""
        caps = [cap for cap in (self.max_file, self.budget) if cap is not None]
        return min(caps) if caps else None

    def spend(self, size):
        if self.budget is not None:
            self.budget -= size


def byte_cap(variable, default):
    """A cap in bytes from an environment variable: a positive count, `none` for no cap (None), or
    `default` when unset."""
    value = os.environ.get(variable, "").strip()
    if not value:
        return default
    if value.lower() == "none":
        return None
    if not value.isdigit() or int(value) == 0:
        sys.exit(f"corpus.py: {variable} must be a positive number of bytes or 'none', not {value!r}")
    return int(value)


def read_manifest(path):
    """The manifest's columns and rows; exits naming the problem if it lacks a required column, a
    row has a stray field, or a name is missing, unsafe as a directory name, or repeated."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        rows = list(reader)
        columns = reader.fieldnames or []
    missing = [c for c in REQUIRED_COLUMNS if c not in columns]
    if missing:
        sys.exit(f"corpus.py: {path} lacks the column(s) {', '.join(missing)}")
    names = set()
    for line, row in enumerate(rows, start=2):
        if None in row:
            sys.exit(f"corpus.py: {path} line {line} has more fields than the header")
        for column in columns:
            row[column] = row[column] or ""
        if not NAME_PATTERN.fullmatch(row["name"]) or row["name"] in names:
            sys.exit(f"corpus.py: {path} line {line}: name {row['name']!r} is empty, repeated or "
                     "not of the form [A-Za-z0-9][A-Za-z0-9._+-]*")
        names.add(row["name"])
    return columns, rows


def write_manifest(path, columns, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n",
                                quoting=csv.QUOTE_NONE)
        writer.writeheader()
        writer.writerows(rows)


def layout(path):
    """The writer-revealing facts of a file's fixed header, parsed directly."""
    with open(path, "rb") as f:
        header = f.read(64)
    order = "<" if header[:4] == b"\x26\xfc\x8f\x88" else ">"
    (_, version, zoom_levels, chrom_tree, full_data, full_index, _, _, _, total_summary,
     uncompress_buf_size, _) = struct.unpack(order + "IHHQQQHHQQIQ", header)
    facts = []
    if chrom_tree > full_data:
        facts.append("chrom tree after data")
    if total_summary == 0:
        facts.append("no total summary")
    elif total_summary == 64 + 24 * 10 and zoom_levels < 10:
        facts.append("10 reserved zoom slots")
    elif total_summary == 64 + 24 * zoom_levels:
        facts.append("summary after the zoom headers")
    facts.append("uncompressed" if uncompress_buf_size == 0 else f"bufSize {uncompress_buf_size}")
    return version, "big" if order == ">" else "little", "; ".join(facts)


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def local_files():
    """Every file under the CORPUS_LOCAL_DIRS directories."""
    return [os.path.join(root, f)
            for directory in filter(None, os.environ.get("CORPUS_LOCAL_DIRS", "").split(":"))
            for root, _, files in os.walk(directory) for f in files]


def find_local(candidates, file_name, sha256):
    """The first of `candidates` whose name ends with `file_name` and whose digest is `sha256`, or
    None. Matching the end of the name finds prefixed copies, such as `07-test.bw` or
    `pyBigWig_test.bw` for `test.bw`; the digest makes any match exact."""
    if not sha256:
        return None
    for path in candidates:
        if os.path.basename(path).endswith(file_name) and sha256_of(path) == sha256:
            return path
    return None


def fetch(row, download_dir, limits, candidates):
    """Fetches one manifest row to <download_dir>/<name>/<the URL's file name>, reusing a local copy
    from `candidates` where one matches; returns (path, status), path being None if nothing was
    fetched."""
    url = row["url"]
    file_name = os.path.basename(url.split("?")[0]) or "download"
    directory = os.path.join(download_dir, row["name"])
    path = os.path.join(directory, file_name)
    if os.path.exists(path):
        return path, "present"
    local = find_local(candidates, file_name, row["sha256"])
    if local:
        os.makedirs(directory, exist_ok=True)
        try:
            os.link(local, path)
            return path, f"linked from {local}"
        except OSError:
            shutil.copyfile(local, path)
            return path, f"copied from {local}"
    refusal = limits.refusal(int(row["size"]) if row["size"] else None)
    if refusal:
        return None, refusal
    os.makedirs(directory, exist_ok=True)
    partial = path + ".part"
    cap = limits.curl_cap()
    process = subprocess.run(
        ["curl", "-sSL", "--fail", "--retry", "2", "--max-time", "3600", "-A", "Mozilla/5.0",
         *(["--max-filesize", str(cap)] if cap is not None else []),
         "-w", "%{http_code}", "-o", partial, url],
        capture_output=True, text=True)
    if process.returncode != 0:
        if os.path.exists(partial):
            os.remove(partial)
        if not os.listdir(directory):
            os.rmdir(directory)
        http_code = process.stdout.strip() or "-"
        # curl's exit status 63 is a response larger than --max-filesize, which may be an error page.
        if process.returncode == 63 and http_code.startswith("2"):
            return None, f"skip: larger than the download cap of {cap} bytes"
        return None, f"unreachable: HTTP {http_code}, curl exit {process.returncode}"
    os.rename(partial, path)
    limits.spend(os.path.getsize(path))
    return path, "downloaded"


def unpack(path):
    """The files to examine from a download: archive members extracted next to it, a gunzipped
    copy, or the file itself."""
    with open(path, "rb") as f:
        head = f.read(4)
    target = path + ".d"
    if head[:2] == b"\x1f\x8b" and tarfile.is_tarfile(path):
        os.makedirs(target, exist_ok=True)
        with tarfile.open(path) as archive:
            archive.extractall(target, filter="data")
        return [os.path.join(root, f) for root, _, files in os.walk(target) for f in sorted(files)]
    if head[:2] == b"\x1f\x8b":
        out = path[:-3] if path.endswith(".gz") else path + ".unzipped"
        if not os.path.exists(out):
            with gzip.open(path, "rb") as source, open(out, "wb") as sink:
                shutil.copyfileobj(source, sink)
        return [out]
    if head == b"PK\x03\x04":
        os.makedirs(target, exist_ok=True)
        with zipfile.ZipFile(path) as archive:
            archive.extractall(target)
        return [os.path.join(root, f) for root, _, files in os.walk(target)
                for f in sorted(files) if "__MACOSX" not in root]
    return [path]


def kind_of(path):
    """The format a file's magic number names, or None."""
    with open(path, "rb") as f:
        return MAGICS.get(f.read(4))


def our_command(*args, timeout=600):
    """Runs our binary; raises Mismatch naming a panic or an error."""
    process = subprocess.run([compare.TOOLS["interop"], *args], capture_output=True, text=True,
                             timeout=timeout)
    if process.returncode == 101 or "panicked at" in process.stderr:
        raise Mismatch(f"PANIC: {process.stderr.strip().splitlines()[0]}")
    if process.returncode != 0:
        raise Mismatch(f"our {args[0]} failed: {process.stderr.strip()}")
    return process.stdout


def pick_regions(bigwig, seed):
    """Up to REGIONS_PER_FILE windows of VALUES_WINDOW bases, as (chrom, start, end), each
    centred on an item sampled uniformly from the whole file, so every window holds data."""
    lengths = dict(bigwig.chroms)
    out = our_command("sample", bigwig.path, str(REGIONS_PER_FILE), str(seed))
    regions = []
    for row in csv.DictReader(out.splitlines(), delimiter="\t"):
        chrom, item_start = row["chrom"], int(row["start"])
        start = max(0, min(item_start - VALUES_WINDOW // 2, lengths[chrom] - VALUES_WINDOW))
        regions.append((chrom, start, min(lengths[chrom], start + VALUES_WINDOW)))
    return regions


def values_vs_kent(bigwig, region):
    """Our per-base values equal kent bigWigToBedGraph's runs over the region, expanded, at %g."""
    chrom, start, end = region
    text = run([compare.TOOLS["bigWigToBedGraph"], *compare.region_args(region, True), bigwig.path,
                "stdout"]).stdout
    expected = [None] * (end - start)
    for line in text.splitlines():
        _, run_start, run_end, value = line.split("\t")
        for position in range(int(run_start), int(run_end)):
            expected[position - start] = value
    values = [float(r["value"]) for r in compare.dump("values", bigwig.path, chrom, str(start), str(end))]
    for offset, (text_value, value) in enumerate(zip(expected, values)):
        if text_value is None:
            if value == value:
                raise Mismatch(f"{chrom}:{start + offset}: ours {value}, kent has no data")
        elif text_value != "%g" % value:
            raise Mismatch(f"{chrom}:{start + offset}: ours {value!r}, kent {text_value}")
    return f"{sum(v is not None for v in expected)} bases"


def raw_summaries_vs_kent_bedgraph(bigwig, query):
    """Our raw summaries equal statistics from kent bigWigToBedGraph's runs over each bin:
    bases covered exactly, min and max at %g, and the mean to within %g's rounding of the
    values."""
    rows = bigwig.summaries(query, "raw")
    text = run([compare.TOOLS["bigWigToBedGraph"], *compare.region_args((query.chrom, query.start, query.end), True),
                bigwig.path, "stdout"]).stdout
    runs = [(int(s), int(e), v) for _, s, e, v in (line.split("\t") for line in text.splitlines())]
    for index, row in enumerate(rows):
        bin_start, bin_end = int(row["start"]), int(row["end"])
        covered, total, magnitude, low, high = 0, 0.0, 0.0, None, None
        for run_start, run_end, value in runs:
            overlap = min(run_end, bin_end) - max(run_start, bin_start)
            if overlap > 0:
                number = float(value)
                covered += overlap
                total += number * overlap
                magnitude += abs(number) * overlap
                low = value if low is None or number < float(low) else low
                high = value if high is None or number > float(high) else high
        if row["covered_bases"] == "":
            if covered:
                raise Mismatch(f"bin {index}: kent covers {covered} bases, ours none")
            continue
        if covered != float(row["covered_bases"]):
            raise Mismatch(f"bin {index}: kent covers {covered} bases, ours {row['covered_bases']}")
        if low != "%g" % float(row["min"]) or high != "%g" % float(row["max"]):
            raise Mismatch(f"bin {index}: kent min/max {low}/{high}, ours {row['min']}/{row['max']}")
        # %g keeps six significant digits, so each value is off by at most 5e-6 of itself.
        if abs(total / covered - float(row["mean"])) > 5e-6 * magnitude / covered + 1e-9:
            raise Mismatch(f"bin {index}: kent mean {total / covered}, ours {row['mean']}")
    return f"{compare.describe(query)}: {len(rows)} bins"


def check_bigwig(path, label, digest):
    """Every check of one bigWig; returns (facts, checks), checks being (name, ok, detail).
    The file's `digest` seeds the choice of regions, so reruns check the same ones."""
    checks = []

    def record(name, function, *args):
        try:
            checks.append((name, True, function(*args)))
        except Mismatch as error:
            checks.append((name, False, str(error)))
        except subprocess.TimeoutExpired as error:
            checks.append((name, False, f"timed out after {error.timeout} s"))
        except Exception as error:
            checks.append((name, False, f"{type(error).__name__}: {error}"))

    version, order, facts = layout(path)
    scan = {}

    def scan_file():
        out = our_command("scan", path, "--round-trip", path + ".roundtrip.bw", timeout=SCAN_TIMEOUT)
        rows = list(csv.DictReader(out.splitlines(), delimiter="\t"))
        scan.update(rows[0])
        read = (f"{scan['items']} items on {scan['chroms_with_data']} of {scan['chrom_count']} chroms "
                f"in {float(scan['seconds']):.1f} s")
        if scan["round_trip"] == "identical":
            return f"{read}; round trip identical"
        # Our writer rejects overlapping items and non-finite values, which lax writers store, so
        # such files cannot be rewritten; that is not a failure.
        overlapping = int(scan["overlapping_items"])
        if scan["round_trip"].startswith("rejected") and (
                (overlapping and "overlapping" in scan["round_trip"])
                or "non-finite" in scan["round_trip"].lower()):
            return (f"{read}; round trip not possible, as designed: {overlapping} overlapping items; "
                    f"{scan['round_trip']}")
        raise Mismatch(f"round trip {scan['round_trip']}")

    record("scan", scan_file)
    summary = {"version": version, "order": order, "layout": facts, "scan": scan}
    if not scan:
        return summary, checks
    bigwig = CorpusFile(path)
    regions = pick_regions(bigwig, int(digest[:15], 16)) if int(scan["items"]) else []
    for region in regions:
        record("bigWigToBedGraph", kent_bedgraph, bigwig, label, region)
        record("bigWigToWig", kent_wig, bigwig, label, region)
        record("values", values_vs_kent, bigwig, region)
        chrom, start, end = region
        length = dict(bigwig.chroms)[chrom]
        # Bins at most twice the finest reduction wide make kent read raw data.
        if bigwig.reductions:
            raw_bins = min(MAX_RAW_BINS, end - start, max(1, -(-(end - start) // (2 * min(bigwig.reductions)))))
        else:
            raw_bins = min(100, end - start)
        raw_query = Query(chrom, start, end, raw_bins, kent_level(bigwig.reductions, (end - start) // raw_bins), False)
        summary_start = max(0, start - (SUMMARY_WINDOW - VALUES_WINDOW) // 2)
        summary_end = min(length, summary_start + SUMMARY_WINDOW)
        width = (summary_end - summary_start) // SUMMARY_BINS
        if width >= 1:
            zoom_query = Query(chrom, summary_start, summary_end, SUMMARY_BINS,
                               kent_level(bigwig.reductions, width), False)
            record("bigWigSummary-wide", kent_summary, bigwig, zoom_query)
            record("summaries-raw", raw_summaries_vs_kent_bedgraph, bigwig, zoom_query)
        if raw_query.bins <= end - start:
            record("bigWigSummary-narrow", kent_summary, bigwig, raw_query)
    if int(scan["items"]) and not regions:
        checks.append(("regions", False, "no item sampled"))
    return summary, checks


def read_annotations(paths):
    """Claimed writers by URL from the annotation TSVs."""
    writers = {}
    for path in paths:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                writer = row.get("writer") or row.get("writer_if_known") or ""
                if row.get("url") and writer:
                    writers[row["url"]] = writer
    return writers


def markdown_cell(text):
    return str(text).replace("|", "\\|").replace("\n", " ")


def write_results(results_path, manifest_path, results, passed, not_checked, writers):
    """Writes the Markdown results: a row per bigWig checked, and a row per entry or archive
    member not checked, with the reason."""
    with open(results_path, "w") as out:
        out.write("# Public-file corpus validation\n\n")
        out.write(f"Generated by `interop/corpus.sh` from `{os.path.basename(manifest_path)}`. "
                  f"{len(results)} bigWigs checked: {passed} passed, {len(results) - passed} failed. "
                  f"{len(not_checked)} manifest entries or archive members not checked.\n\n")
        out.write("## Checked files\n\n")
        out.write("| File | Source | Claimed writer | Version | Byte order | Layout | Chroms (with data) "
                  "| Zoom levels | Items | Result | Notes |\n")
        out.write("|---|---|---|---|---|---|---|---|---|---|---|\n")
        for row, label, facts, checks in results:
            scan = facts["scan"]
            failed = [c for c in checks if not c[1]]
            result = f"pass ({len(checks)} checks)" if not failed else f"FAIL ({len(failed)} of {len(checks)})"
            notes = "; ".join(f"{name}: {detail}" for name, _, detail in failed) or row["expectation"]
            chroms = f"{scan.get('chrom_count', '?')} ({scan.get('chroms_with_data', '?')})"
            zooms = len([r for r in scan.get("zoom_reductions", "").split(",") if r]) if scan else "?"
            cells = [label, source(row), writers.get(row["url"], ""), facts["version"], facts["order"],
                     facts["layout"], chroms, zooms, scan.get("items", "?"), result, notes]
            out.write("| " + " | ".join(markdown_cell(c) for c in cells) + " |\n")
        out.write("\n## Not checked\n\n| URL | Source | Reason |\n|---|---|---|\n")
        for row, reason in not_checked:
            out.write(f"| {markdown_cell(row['url'])} | {markdown_cell(source(row))} | {markdown_cell(reason)} |\n")


def source(row):
    """Where a manifest entry came from: its issue, or `curated` for a file chosen directly."""
    return row["source_issue"] or "curated"


def main():
    args, annotations, fetch_only = [], [], False
    argv = iter(sys.argv[1:])
    for arg in argv:
        if arg == "--annotations":
            annotations.append(next(argv))
        elif arg == "--fetch-only":
            fetch_only = True
        else:
            args.append(arg)
    if len(args) != (3 if fetch_only else 4):
        sys.exit(__doc__)
    manifest_path, download_dir, filled_path = args[:3]
    columns, manifest = read_manifest(manifest_path)
    os.makedirs(download_dir, exist_ok=True)

    limits, candidates = Limits(), local_files()
    fetched, not_fetched, outcomes = [], [], collections.Counter()
    for row in manifest:
        if row["expectation"].startswith("skip:"):
            path, status = None, row["expectation"]
        else:
            path, status = fetch(row, download_dir, limits, candidates)
        outcome = status.split()[0].rstrip(":")
        outcome = "skipped" if outcome == "skip" else outcome
        if path is not None:
            digest, size = sha256_of(path), os.path.getsize(path)
            if row["sha256"] and row["sha256"] != digest:
                path, outcome = None, "mismatched"
                status = f"checksum mismatch: got {digest}, expected {row['sha256']}"
        outcomes[outcome] += 1
        if path is None:
            not_fetched.append((row, status))
            print(f"{'BAD ' if outcome == 'mismatched' else 'SKIP'}  {row['name']}: {status}", flush=True)
            continue
        row["sha256"], row["size"] = digest, str(size)
        fetched.append((row, path, digest))
        print(f"GOT   {row['name']}: {status}, {size} bytes", flush=True)
    write_manifest(filled_path, columns, manifest)
    print(f"\n{len(fetched)} of {len(manifest)} entries fetched; "
          + ", ".join(f"{count} {outcome}" for outcome, count in outcomes.items()), flush=True)
    if fetch_only:
        sys.exit(1 if outcomes["mismatched"] else 0)

    writers = {row["url"]: row.get("writer", "") for row in manifest} | read_annotations(annotations)
    results, not_checked = [], list(not_fetched)
    for row, path, digest in fetched:
        members = [m for m in unpack(path) if kind_of(m)]
        if not members:
            not_checked.append((row, "no bigWig in the download"))
            print(f"NONE  {row['name']}: no bigWig in the download", flush=True)
            continue
        for member in members:
            label = os.path.relpath(member, download_dir)
            if not kind_of(member).startswith("bigWig"):
                not_checked.append((row, f"{label} is a {kind_of(member)}"))
                continue
            facts, checks = check_bigwig(member, label, digest)
            failed = [c for c in checks if not c[1]]
            outcome = "pass" if not failed else "FAIL"
            print(f"{outcome}  {label}: {len(checks)} checks" + (f"; {failed[0][0]}: {failed[0][2]}" if failed else ""), flush=True)
            for name, ok, detail in checks:
                print(f"      {'ok  ' if ok else 'FAIL'} {name}: {detail}", flush=True)
            results.append((row, label, facts, checks))

    results_path = args[3]
    passed = sum(all(c[1] for c in checks) for _, _, _, checks in results)
    write_results(results_path, manifest_path, results, passed, not_checked, writers)
    print(f"\n{len(results)} bigWigs checked: {passed} passed, {len(results) - passed} failed; "
          f"{len(not_checked)} not checked; results in {results_path}")
    sys.exit(0 if passed == len(results) else 1)


if __name__ == "__main__":
    main()
