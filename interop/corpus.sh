#!/usr/bin/env bash
# One-time validation of our reader and writer against public bigWig files.
# Run by hand, never in CI: downloads are too brittle to gate on.
#
# Usage: interop/corpus.sh <manifest.tsv> <download-dir> <filled-manifest.tsv> <results.md> \
#            [--annotations <tsv>]...
#        interop/corpus.sh --fetch-only <manifest.tsv> <download-dir> <filled-manifest.tsv>
#
# The manifest is interop/corpus.tsv (its columns are described in interop/README.md). Downloads
# are verified against sha256 where given; the filled manifest records the digest and size of
# everything fetched, and the results table what every check found. --fetch-only stops after the
# downloads and checksums, and needs only curl and Python: no kent tools and no Rust build.
#
# Tools are found on PATH unless overridden, as in run.sh: KENT_DIR (or BIGWIGINFO,
# BIGWIGTOBEDGRAPH, BIGWIGTOWIG, BIGWIGSUMMARY), PYTHON (default python3) and INTEROP_BIN.
# CORPUS_LOCAL_DIRS (colon-separated) names directories whose files, matched by name and sha256,
# are hard-linked instead of downloaded. CORPUS_MAX_FILE_BYTES caps each download (1 GiB) and
# CORPUS_BUDGET_BYTES the total (4 GB); `none` lifts either cap.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The path of a tool: `$dir/$name` when `dir` is set, otherwise `name` for a PATH lookup.
tool_in() {
    local dir="$1" name="$2"
    if [[ -n "${dir}" ]]; then echo "${dir}/${name}"; else echo "${name}"; fi
}

export BIGWIGINFO="${BIGWIGINFO:-$(tool_in "${KENT_DIR:-}" bigWigInfo)}"
export BIGWIGTOBEDGRAPH="${BIGWIGTOBEDGRAPH:-$(tool_in "${KENT_DIR:-}" bigWigToBedGraph)}"
export BIGWIGTOWIG="${BIGWIGTOWIG:-$(tool_in "${KENT_DIR:-}" bigWigToWig)}"
export BIGWIGSUMMARY="${BIGWIGSUMMARY:-$(tool_in "${KENT_DIR:-}" bigWigSummary)}"
python="${PYTHON:-python3}"

fetch_only=0
for arg in "$@"; do
    if [[ "${arg}" == "--fetch-only" ]]; then
        fetch_only=1
    fi
done

required=("${python}" curl)
if (( ! fetch_only )); then
    required+=("${BIGWIGTOBEDGRAPH}" "${BIGWIGTOWIG}" "${BIGWIGSUMMARY}")
fi
missing=0
for command in "${required[@]}"; do
    if ! command -v "${command}" > /dev/null; then
        echo "corpus.sh: ${command} not found" >&2
        missing=1
    fi
done
if (( missing )); then
    exit 2
fi

if (( ! fetch_only )) && [[ -z "${INTEROP_BIN:-}" ]]; then
    cargo build --release --locked --quiet --manifest-path "${here}/Cargo.toml"
    export INTEROP_BIN="${here}/target/release/bigwig-interop"
fi

"${python}" "${here}/corpus.py" "$@"
