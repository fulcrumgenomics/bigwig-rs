// Reads one bigWig with @gmod/bbi (the JBrowse 2 reader) and prints what it sees as JSON.
//
// Usage: node bbi_dump.mjs <queries.json>
//
// queries.json is {"file": "<bigwig>", "queries": [{"chrom", "start", "end", "scale"?}, ...]}.
// The output is {"header": {...}} or {"error": "..."} if the header can't be read, and otherwise
// "results": one {"features": [...]} or {"error": "..."} per query, in order. A query with a
// scale reads the zoom level bbi picks for it; one without reads raw data.
//
// Plain ESM with no runtime-specific APIs, so it runs under node and bun alike.
import { readFileSync } from 'node:fs'

import { BigWig } from '@gmod/bbi'
import { LocalFile } from 'generic-filehandle2'

const [queriesPath] = process.argv.slice(2)
if (!queriesPath) {
  console.error('usage: bbi_dump.mjs <queries.json>')
  process.exit(2)
}
const { file, queries } = JSON.parse(readFileSync(queriesPath, 'utf8'))

const out = {}
try {
  const bw = new BigWig({ filehandle: new LocalFile(file) })
  const header = await bw.getHeader()
  out.header = {
    zoomLevels: header.zoomLevels.map(z => z.reductionLevel),
    totalSummary: header.totalSummary,
    uncompressBufSize: header.uncompressBufSize,
    refsByNumber: header.refsByNumber.map(r => [r.name, r.id, r.length]),
  }
  out.results = []
  for (const q of queries) {
    try {
      const opts = q.scale === undefined ? {} : { scale: q.scale }
      const features = await bw.getFeatures(q.chrom, q.start, q.end, opts)
      out.results.push({
        features: features.map(f =>
          f.summary
            ? [f.start, f.end, f.score, f.minScore, f.maxScore]
            : [f.start, f.end, f.score],
        ),
      })
    } catch (e) {
      out.results.push({ error: String(e) })
    }
  }
} catch (e) {
  out.error = String(e)
}
// bigints (from 64-bit header fields) become numbers; every value here fits in a double.
process.stdout.write(JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? Number(v) : v)))
