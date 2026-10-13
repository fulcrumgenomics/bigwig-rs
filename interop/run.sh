#!/usr/bin/env bash
# Runs the cross-implementation interop matrix: builds the interop binary, generates
# the fixtures, then has compare.py check them against kent, pyBigWig, bigtools, @gmod/bbi, IGV
# and igv.js in both directions. Exits non-zero if any required check fails. Arguments, if any,
# name the fixtures to check (see compare.py); with none, everything runs.
#
# Every tool is found on PATH unless overridden:
#   KENT_DIR       directory holding the kent utilities (bigWigInfo, wigToBigWig, ...)
#   BIGTOOLS_DIR   directory holding the bigtools binaries (bigwiginfo, ...)
#   PYTHON         a Python with pyBigWig (default python3)
#   JS_RUNTIME     node, or bun (default node); interop/node_modules must hold @gmod/bbi
#   JAVA           a Java 21 java, with javac beside it (default java)
#   IGV_LIB        directory holding IGV's jars (default $CONDA_PREFIX/lib/igv, where bioconda's
#                  igv package puts them)
#   IGVJS_DIR      an igv.js source tree laid out as igvjs_dump.mjs describes; without it, the
#                  pinned release is fetched from GitHub once, into interop/target/igvjs-<version>
#   SEED           fixture seed (default 1)
#   OUT_DIR        where fixtures and converted files go (default interop/target/interop-run)
#   INTEROP_BIN    a prebuilt bigwig-interop binary; without it, run.sh builds one
# Individual tools can also be set directly: BIGWIGINFO, BIGWIGTOBEDGRAPH, BIGWIGTOWIG,
# BIGWIGSUMMARY, BEDGRAPHTOBIGWIG, WIGTOBIGWIG, BIGTOOLS_BIGWIGINFO, BIGTOOLS_BIGWIGTOBEDGRAPH,
# BIGTOOLS_BEDGRAPHTOBIGWIG.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# igv.js v3.8.9 and igv-utils v1.8.4, the igv-utils its package.json names, by commit. The hash is
# of the files kept (see content_sha256), not of the tarballs, whose bytes GitHub does not promise
# to keep stable.
IGVJS_VERSION=3.8.9
IGVJS_COMMIT=e6bc7403caed0767c88a5d197fa40092d0441e5a
IGV_UTILS_COMMIT=d8c3eb7efc5bc435b1a294146f943ba6f2d22bda
IGVJS_CONTENT_SHA256=7fb2a0fe3e6dabcb51363d71af0fef08877aa41fe39bd5570a39cf294284ba84

# The path of a tool: `$dir/$name` when `dir` is set, otherwise `name` for a PATH lookup.
tool_in() {
    local dir="$1" name="$2"
    if [[ -n "${dir}" ]]; then echo "${dir}/${name}"; else echo "${name}"; fi
}

# The sha256 of a listing of every file's sha256 and path under a directory, in byte order of
# path.
content_sha256() {
    (cd "$1" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 shasum -a 256) \
        | shasum -a 256 | cut -d' ' -f1
}

# Fetches the pinned igv.js source and the igv-utils source it imports into `dest`, keeping only
# what the bigWig reader needs, with igv-utils' src/ nested where bwReader.js imports it from.
# Nothing is installed and no package scripts run. The tree is assembled and checked beside
# `dest` and moved into place only once its hash matches, so a failed fetch leaves nothing to
# reuse.
fetch_igvjs() {
    local dest="$1" actual
    # Global, so the trap can still name it if the script exits here.
    fetch_dir="$(mktemp -d "${dest}.fetch.XXXXXX")"
    trap 'rm -rf "${fetch_dir}"' EXIT
    curl -fsSL --retry 3 "https://codeload.github.com/igvteam/igv.js/tar.gz/${IGVJS_COMMIT}" \
        | tar -xz -C "${fetch_dir}"
    curl -fsSL --retry 3 "https://codeload.github.com/igvteam/igv-utils/tar.gz/${IGV_UTILS_COMMIT}" \
        | tar -xz -C "${fetch_dir}"
    mkdir -p "${fetch_dir}/tree/node_modules/igv-utils"
    cp -R "${fetch_dir}/igv.js-${IGVJS_COMMIT}/js" "${fetch_dir}/igv.js-${IGVJS_COMMIT}/package.json" \
        "${fetch_dir}/tree/"
    cp -R "${fetch_dir}/igv-utils-${IGV_UTILS_COMMIT}/src" "${fetch_dir}/tree/node_modules/igv-utils/"
    actual="$(content_sha256 "${fetch_dir}/tree")"
    if [[ "${actual}" != "${IGVJS_CONTENT_SHA256}" ]]; then
        echo "run.sh: the igv.js files fetched hash to ${actual}, not ${IGVJS_CONTENT_SHA256}" >&2
        exit 2
    fi
    mv "${fetch_dir}/tree" "${dest}"
    rm -rf "${fetch_dir}"
    trap - EXIT
    echo "Fetched igv.js ${IGVJS_VERSION} into ${dest}"
}

export BIGWIGINFO="${BIGWIGINFO:-$(tool_in "${KENT_DIR:-}" bigWigInfo)}"
export BIGWIGTOBEDGRAPH="${BIGWIGTOBEDGRAPH:-$(tool_in "${KENT_DIR:-}" bigWigToBedGraph)}"
export BIGWIGTOWIG="${BIGWIGTOWIG:-$(tool_in "${KENT_DIR:-}" bigWigToWig)}"
export BIGWIGSUMMARY="${BIGWIGSUMMARY:-$(tool_in "${KENT_DIR:-}" bigWigSummary)}"
export BEDGRAPHTOBIGWIG="${BEDGRAPHTOBIGWIG:-$(tool_in "${KENT_DIR:-}" bedGraphToBigWig)}"
export WIGTOBIGWIG="${WIGTOBIGWIG:-$(tool_in "${KENT_DIR:-}" wigToBigWig)}"
export BIGTOOLS_BIGWIGINFO="${BIGTOOLS_BIGWIGINFO:-$(tool_in "${BIGTOOLS_DIR:-}" bigwiginfo)}"
export BIGTOOLS_BIGWIGTOBEDGRAPH="${BIGTOOLS_BIGWIGTOBEDGRAPH:-$(tool_in "${BIGTOOLS_DIR:-}" bigwigtobedgraph)}"
export BIGTOOLS_BEDGRAPHTOBIGWIG="${BIGTOOLS_BEDGRAPHTOBIGWIG:-$(tool_in "${BIGTOOLS_DIR:-}" bedgraphtobigwig)}"
export JS_RUNTIME="${JS_RUNTIME:-node}"
export BBI_DUMP="${here}/bbi_dump.mjs"
export IGVJS_DUMP="${here}/igvjs_dump.mjs"
export JAVA="${JAVA:-java}"
igv_lib="${IGV_LIB:-${CONDA_PREFIX:+${CONDA_PREFIX}/lib/igv}}"
igvjs_dir="${IGVJS_DIR:-${here}/target/igvjs-${IGVJS_VERSION}}"
python="${PYTHON:-python3}"
seed="${SEED:-1}"
out_dir="${OUT_DIR:-${here}/target/interop-run}"

# Fail early, naming every missing tool, rather than as a cascade of failed checks.
missing=0
required=("${BIGWIGINFO}" "${BIGWIGTOBEDGRAPH}" "${BIGWIGTOWIG}" "${BIGWIGSUMMARY}"
    "${BEDGRAPHTOBIGWIG}" "${WIGTOBIGWIG}" "${BIGTOOLS_BIGWIGINFO}" "${BIGTOOLS_BIGWIGTOBEDGRAPH}"
    "${BIGTOOLS_BEDGRAPHTOBIGWIG}" "${JS_RUNTIME}" "${python}" "${JAVA}")
# javac sits beside java in a JDK's bin/, and so do the links to them in conda's bin/ or /usr/bin.
if command -v "${JAVA}" > /dev/null; then
    javac="$(dirname "$(command -v "${JAVA}")")/javac"
    required+=("${javac}")
fi
if [[ -z "${IGVJS_DIR:-}" && ! -d "${igvjs_dir}" ]]; then
    required+=(curl tar shasum)
fi
for command in "${required[@]}"; do
    if ! command -v "${command}" > /dev/null; then
        echo "run.sh: ${command} not found" >&2
        missing=1
    fi
done
if ! "${python}" -c "import pyBigWig" 2> /dev/null; then
    echo "run.sh: ${python} cannot import pyBigWig" >&2
    missing=1
fi
if [[ ! -d "${here}/node_modules/@gmod/bbi" ]]; then
    echo "run.sh: @gmod/bbi is not installed; run 'npm ci' in ${here}" >&2
    missing=1
fi
if [[ -z "${igv_lib}" ]]; then
    echo "run.sh: IGV_LIB is not set, and there is no conda environment to find IGV in" >&2
    missing=1
elif [[ ! -f "${igv_lib}/igv.jar" ]]; then
    echo "run.sh: IGV is not installed: no igv.jar in ${igv_lib}" >&2
    missing=1
fi
if [[ -n "${IGVJS_DIR:-}" && ! -f "${IGVJS_DIR}/js/bigwig/bwSource.js" ]]; then
    echo "run.sh: IGVJS_DIR ${IGVJS_DIR} is not an igv.js source tree" >&2
    missing=1
fi
if (( missing )); then
    exit 2
fi

if [[ ! -d "${igvjs_dir}" ]]; then
    mkdir -p "$(dirname "${igvjs_dir}")"
    fetch_igvjs "${igvjs_dir}"
fi
export IGVJS_DIR="${igvjs_dir}"

if [[ -z "${INTEROP_BIN:-}" ]]; then
    cargo build --release --locked --quiet --manifest-path "${here}/Cargo.toml"
    export INTEROP_BIN="${here}/target/release/bigwig-interop"
fi

rm -rf "${out_dir}"
"${INTEROP_BIN}" generate "${seed}" "${out_dir}"
echo "Fixtures for seed ${seed} in ${out_dir}"

# Compiled once here, rather than run as a source file, which would compile it again for every
# fixture.
"${javac}" -d "${out_dir}/igv-dump" -cp "${igv_lib}/*" "${here}/IgvDump.java"
export IGV_CLASSPATH="${out_dir}/igv-dump:${igv_lib}/*"

"${python}" "${here}/compare.py" "${out_dir}" "$@"
