// Reads one bigWig with igv.js's own reader (js/bigwig/bwSource.js, the source behind a wig
// track) and prints what it sees as JSON.
//
// Usage: node igvjs_dump.mjs <igv.js-dir> <queries.json>
//
// <igv.js-dir> is an igv.js source tree with igv-utils' src/ at node_modules/igv-utils, where
// bwReader.js imports it from by relative path; run.sh fetches one. The npm package can't be used:
// it ships only a bundle, which needs a browser's `document`.
//
// queries.json is {"file": "<bigwig>", "chroms": [name, ...], "queries": [{"chrom", "start",
// "end", "bpPerPixel"?}, ...]}. The output is {"header", "chromIds", "results"}, or {"error"} if
// the header can't be read:
// - chromIds: the id igv.js's chromosome lookup finds for each of "chroms", in order, then for a
//   name no file holds; null where it finds none.
// - results: one {"features": [...]} or {"error"} per query, in order. A query without bpPerPixel
//   reads raw data, as [chr, start, end, value]. One with bpPerPixel reads the zoom level igv.js
//   picks for that resolution, as [chr, start, end, mean, min, max] from its mean, min and max
//   window functions.
//
// The file is handed to igv.js as a File, as a browser's local-file picker does, so it is read
// with File.slice() rather than over HTTP; everything above that is the code a URL goes through.
// Plain ESM, so it runs under node (20 or later, for the global File) and bun alike.
import { readFileSync } from 'node:fs'
import { basename, resolve } from 'node:path'
import { pathToFileURL } from 'node:url'

const [igvDir, queriesPath] = process.argv.slice(2)
if (!igvDir || !queriesPath) {
  console.error('usage: igvjs_dump.mjs <igv.js-dir> <queries.json>')
  process.exit(2)
}
const { file, chroms, queries } = JSON.parse(readFileSync(queriesPath, 'utf8'))
const { default: BWSource } = await import(pathToFileURL(resolve(igvDir, 'js/bigwig/bwSource.js')).href)

const ABSENT_CHROM = '__igvjs_dump_absent__'
const WINDOW_FUNCTIONS = ['mean', 'min', 'max']

const out = {}
try {
  const source = new BWSource({ url: new File([readFileSync(file)], basename(file)), format: 'bigwig' })
  const reader = source.reader
  const header = await source.getHeader()
  out.header = {
    type: reader.type,
    littleEndian: reader.littleEndian,
    version: header.bwVersion,
    uncompressBuffSize: header.uncompressBuffSize,
    chromCount: reader.chromTree.getItemCount(),
    // Coarsest first, as igv.js stores them.
    zoomLevels: reader.zoomLevelHeaders.map(z => z.reductionLevel),
    totalSummary: reader.totalSummary,
  }

  // Before any feature query: decoding features looks names up by id, which caches the names of
  // every leaf it walks, and a name cached that way resolves even where the tree search fails.
  out.chromIds = []
  for (const name of [...chroms, ABSENT_CHROM]) {
    out.chromIds.push((await reader.getIdForChr(name)) ?? null)
  }

  out.results = []
  for (const { chrom, start, end, bpPerPixel } of queries) {
    try {
      if (bpPerPixel === undefined) {
        const features = await source.getFeatures({ chr: chrom, start, end, windowFunction: 'none' })
        out.results.push({ features: features.map(f => [f.chr, f.start, f.end, f.value]) })
      } else {
        const perFunction = []
        for (const windowFunction of WINDOW_FUNCTIONS) {
          perFunction.push(await source.getFeatures({ chr: chrom, start, end, bpPerPixel, windowFunction }))
        }
        out.results.push({ features: zipWindowFunctions(...perFunction) })
      }
    } catch (e) {
      out.results.push({ error: String(e) })
    }
  }
} catch (e) {
  out.error = String(e)
}
process.stdout.write(JSON.stringify(out))

// One feature per zoom record from each window function, in the same order, zipped by position
// rather than by interval, which would hide a duplicated record.
function zipWindowFunctions(means, mins, maxes) {
  const sameIntervals = (a, b) => a.length === b.length && a.every((f, i) => f.start === b[i].start && f.end === b[i].end)
  if (!sameIntervals(means, mins) || !sameIntervals(means, maxes)) {
    throw Error('the mean, min and max window functions returned different records')
  }
  return means.map((f, i) => [f.chr, f.start, f.end, f.value, mins[i].value, maxes[i].value])
}
