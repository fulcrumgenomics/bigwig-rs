// Reads one bigWig with IGV desktop's own classes, headless, and prints what IGV sees as TSV.
//
// Usage, with IGV 2.19.x's jars (a release's lib/ directory, or $CONDA_PREFIX/lib/igv):
//
//   javac -d <classes> -cp "$IGV_LIB/*" IgvDump.java
//   java -Djava.awt.headless=true -cp "<classes>:$IGV_LIB/*" IgvDump <bigwig> <chrom.sizes> <requests.tsv>
//
// The bigWig is opened as IGV opens a bigWig track: a BBFile and a BBDataSource, against a Genome
// loaded from <chrom.sizes>. Each line of <requests.tsv> is one request, answered in order:
//
//   raw    <chrom> <start> <end>              BBDataSource.getRawData: what IGV draws zoomed in
//   level  <chrom> <start> <end> <reduction>  the records of the zoom level with that reduction,
//                                             decoded as BBDataSource decodes a level it chose
//   zoom   <chrom> <start> <end> <igvZoom>    BBDataSource.getPrecomputedSummaryScores at IGV
//                                             zoom <igvZoom>, so IGV chooses the zoom level
//   genome                                    the whole-genome view, IGV's first view of a genome
//                                             with several chromosomes
//
// Output lines are tab-separated, each starting with its kind. <query> is the request's 0-based
// line number, and a zoom record's values are what IGV shows under its min, mean and max window
// functions:
//
//   header  <version> <reductions, coarsest first> <uncompressBufSize> <basesCovered> <min> <max> <sum> <sumSquares>
//   chrom   <id> <name> <the id IGV's lookup of that name finds>
//   raw     <query> <start> <end> <value>
//   zoom    <query> <chrom> <start> <end> <min> <mean> <max>    (level and zoom requests)
//   genome  <query> <start kbp> <end kbp> <min> <mean> <max>
//   wgchrom <query> <name>                       (the chromosomes the whole-genome view spans)
//   zoomsel <query> <bp per pixel> <reduction IGV picks, or none> [<home chromosome>, genome only]
//   done    <query> <records>
//   error   <query, or "open"> <message>
//
// Values are f32 widened to f64 and printed with Double.toString, so they read back exactly.

import org.broad.igv.Globals;
import org.broad.igv.data.DataTile;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.load.ChromsizesLoader;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ucsc.bb.BBDataSource;
import org.broad.igv.ucsc.bb.BBFile;
import org.broad.igv.ucsc.bb.BBHeader;
import org.broad.igv.ucsc.bb.BBTotalSummary;
import org.broad.igv.ucsc.bb.BBZoomHeader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

public class IgvDump {

    static final WindowFunction[] WINDOW_FUNCTIONS = {WindowFunction.min, WindowFunction.mean, WindowFunction.max};

    private final BBFile file;
    private final Genome genome;
    private final Source source;
    private final PrintWriter out;

    public static void main(String[] args) throws IOException {
        if (args.length != 3) {
            System.err.println("usage: IgvDump <bigwig> <chrom.sizes> <requests.tsv>");
            System.exit(2);
        }
        Globals.setHeadless(true);
        PrintWriter out = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out, StandardCharsets.UTF_8)));
        try {
            IgvDump dump;
            try {
                dump = new IgvDump(args[0], args[1], out);
                dump.header();
            } catch (Throwable e) {
                print(out, "error", "open", describe(e));
                return;
            }
            List<String> requests = Files.readAllLines(Path.of(args[2]), StandardCharsets.UTF_8);
            for (int query = 0; query < requests.size(); query++) {
                try {
                    dump.answer(query, requests.get(query).split("\t", -1));
                } catch (Throwable e) {
                    print(out, "error", query, describe(e));
                }
            }
        } finally {
            out.flush();
        }
    }

    IgvDump(String bigwig, String chromSizes, PrintWriter out) throws IOException {
        this.genome = new ChromsizesLoader(chromSizes).loadGenome();
        this.file = new BBFile(bigwig, genome);
        this.source = new Source(file, genome);
        this.out = out;
    }

    /** Prints the header fields, the zoom levels and total summary, and every chromosome with
     *  the id IGV's name lookup finds for it. */
    void header() throws ReflectiveOperationException {
        BBHeader header = file.getHeader();
        StringJoiner reductions = new StringJoiner(",");
        for (BBZoomHeader zoom : Internals.zoomHeaders(file)) {
            reductions.add(Long.toString(zoom.reductionLevel));
        }
        BBTotalSummary total = file.getTotalSummary();
        print(out, "header", header.version, reductions, header.uncompressBuffSize, total.basesCovered, total.minVal,
                total.maxVal, total.sumData, total.sumSquares);
        String[] names = file.getChromosomeNames();
        for (int id = 0; id < names.length; id++) {
            print(out, "chrom", id, names[id], Internals.ID_FOR_CHR.invoke(file, names[id]));
        }
    }

    void answer(int query, String[] request) throws ReflectiveOperationException {
        int records = switch (request[0]) {
            case "raw" -> raw(query, request[1], Integer.parseInt(request[2]), Integer.parseInt(request[3]));
            case "level" -> level(query, request[1], Integer.parseInt(request[2]), Integer.parseInt(request[3]),
                    Long.parseLong(request[4]));
            case "zoom" -> zoom(query, request[1], Integer.parseInt(request[2]), Integer.parseInt(request[3]),
                    Integer.parseInt(request[4]));
            case "genome" -> genome(query);
            default -> throw new IllegalArgumentException("unknown request " + request[0]);
        };
        print(out, "done", query, records);
    }

    int raw(int query, String chrom, int start, int end) {
        DataTile tile = source.raw(chrom, start, end);
        int n = tile == null ? 0 : tile.getStartLocations().length;
        for (int i = 0; i < n; i++) {
            print(out, "raw", query, tile.getStartLocations()[i], tile.getEndLocations()[i], (double) tile.getValues()[i]);
        }
        return n;
    }

    /** The zoom level's records, read as getPrecomputedSummaryScores reads the level it picks, so
     *  that every level is read whether or not IGV would pick it. */
    @SuppressWarnings("unchecked")
    int level(int query, String chrom, int start, int end, long reduction) throws ReflectiveOperationException {
        BBZoomHeader header = null;
        for (BBZoomHeader zoom : Internals.zoomHeaders(file)) {
            if (zoom.reductionLevel == reduction) header = zoom;
        }
        if (header == null) throw new IllegalArgumentException("no zoom level with reduction " + reduction);
        Integer chromId = (Integer) Internals.ID_FOR_CHR.invoke(file, chrom);
        List<byte[]> chunks = (List<byte[]>) Internals.LEAF_CHUNKS.invoke(file, chrom, start, chrom, end, header.indexOffset);
        return printZoom("zoom", query, function -> {
            List<LocusScore> scores = new ArrayList<>();
            for (byte[] chunk : chunks) {
                Internals.DECODE_ZOOM_DATA.invoke(file, chunk, chromId, start, end, function, scores);
            }
            return scores;
        });
    }

    /** What IGV draws at IGV zoom `igvZoom`, where it asks for chromosome length / (2^zoom × 700)
     *  bases per pixel, and the zoom level it picks for that (getPrecomputedSummaryScores's own
     *  arithmetic). */
    int zoom(int query, String chrom, int start, int end, int igvZoom) throws ReflectiveOperationException {
        double bpPerPixel = genome.getChromosome(chrom).getLength() / (Math.pow(2, igvZoom) * 700);
        BBZoomHeader chosen = (BBZoomHeader) Internals.ZOOM_LEVEL_FOR_SCALE.invoke(file, bpPerPixel, 2);
        print(out, "zoomsel", query, bpPerPixel, chosen == null ? "none" : chosen.reductionLevel);
        return printZoom("zoom", query, function -> {
            source.setWindowFunction(function);
            return source.zoom(chrom, start, end, igvZoom);
        });
    }

    /** The whole-genome view's records, in genome coordinates (kbp), and the zoom level it reads
     *  (getWholeGenomeScores's own arithmetic). */
    int genome(int query) throws ReflectiveOperationException {
        for (String name : genome.getLongChromosomeNames()) {
            print(out, "wgchrom", query, name);
        }
        double bpPerPixel = genome.getWGLength() / 1000;
        BBZoomHeader chosen = (BBZoomHeader) Internals.ZOOM_LEVEL_FOR_SCALE.invoke(file, bpPerPixel, 1000);
        print(out, "zoomsel", query, bpPerPixel, chosen == null ? "none" : chosen.reductionLevel,
                genome.getHomeChromosome());
        return printZoom("genome", query, function -> {
            source.setWindowFunction(function);
            return source.zoom(Globals.CHR_ALL, 0, Integer.MAX_VALUE, 0);
        });
    }

    /** Reads zoom records under the min, mean and max window functions and prints each record once
     *  with its three values, as `zoom <query> <chrom> <start> <end> <min> <mean> <max>`, or for
     *  the whole-genome view, whose records have no chromosome, `genome <query> <start> <end> ...`.
     *  Every read must return the same records in the same order. */
    int printZoom(String kind, int query, ZoomRead read) throws ReflectiveOperationException {
        List<List<LocusScore>> reads = new ArrayList<>();
        for (WindowFunction function : WINDOW_FUNCTIONS) {
            List<LocusScore> scores = read.apply(function);
            reads.add(scores == null ? List.of() : scores);
        }
        List<LocusScore> records = reads.get(0);
        for (List<LocusScore> other : reads) {
            boolean same = other.size() == records.size();
            for (int i = 0; same && i < records.size(); i++) {
                same = other.get(i).getStart() == records.get(i).getStart() && other.get(i).getEnd() == records.get(i).getEnd();
            }
            if (!same) throw new IllegalStateException("the window functions read different records");
        }
        for (int i = 0; i < records.size(); i++) {
            LocusScore record = records.get(i);
            List<Object> fields = new ArrayList<>(List.of(kind, query));
            if (kind.equals("zoom")) fields.add(record.getContig());
            fields.addAll(List.of(record.getStart(), record.getEnd()));
            for (List<LocusScore> scores : reads) fields.add((double) scores.get(i).getScore());
            print(out, fields.toArray());
        }
        return records.size();
    }

    /** A read of zoom records under one window function. */
    interface ZoomRead {
        List<LocusScore> apply(WindowFunction function) throws ReflectiveOperationException;
    }

    /** The root cause of a failure, where it was thrown, on one line. */
    static String describe(Throwable e) {
        while (e.getCause() != null) e = e.getCause();
        String where = e.getStackTrace().length > 0 ? " at " + e.getStackTrace()[0] : "";
        return (e.getClass().getSimpleName() + ": " + e.getMessage() + where).replace('\t', ' ').replace('\n', ' ');
    }

    static void print(PrintWriter out, Object... fields) {
        StringJoiner line = new StringJoiner("\t");
        for (Object field : fields) line.add(String.valueOf(field));
        out.println(line);
    }

    /** Exposes the protected methods IGV's track rendering calls on a bigWig data source. */
    static class Source extends BBDataSource {
        Source(BBFile file, Genome genome) throws IOException {
            super(file, genome);
        }

        DataTile raw(String chrom, int start, int end) {
            return getRawData(chrom, start, end);
        }

        List<LocusScore> zoom(String chrom, int start, int end, int igvZoom) {
            return getPrecomputedSummaryScores(chrom, start, end, igvZoom);
        }
    }

    /**
     * The private and package-private BBFile members of IGV 2.19.x this harness relies on, reached
     * reflectively; that works because every class is on the class path, in the unnamed module.
     * Together with BBDataSource's protected getRawData and getPrecomputedSummaryScores (see
     * Source), these are what an IGV upgrade must port. IGV 3.0 renames the package org.broad.igv
     * to org.igv and changes several of these signatures.
     */
    static class Internals {
        /** The zoom levels, coarsest first. */
        static final Field ZOOM_HEADERS = field("zoomHeaders");
        /** The chromosome B+ tree lookup of a name. */
        static final Method ID_FOR_CHR = method("getIdForChr", String.class);
        /** The zoom level IGV reads for a resolution, given a tolerance for the finest level. */
        static final Method ZOOM_LEVEL_FOR_SCALE = method("zoomLevelForScale", double.class, int.class);
        /** The R-tree walk: the data blocks overlapping a region of an index. */
        static final Method LEAF_CHUNKS = method("getLeafChunks", String.class, int.class, String.class, int.class,
                long.class);
        /** Decodes one zoom data block's records for a region under a window function. */
        static final Method DECODE_ZOOM_DATA = method("decodeZoomData", byte[].class, int.class, int.class, int.class,
                WindowFunction.class, List.class);

        static BBZoomHeader[] zoomHeaders(BBFile file) throws IllegalAccessException {
            return (BBZoomHeader[]) ZOOM_HEADERS.get(file);
        }

        static Field field(String name) {
            try {
                Field field = BBFile.class.getDeclaredField(name);
                field.setAccessible(true);
                return field;
            } catch (NoSuchFieldException e) {
                throw new IllegalStateException("BBFile has no field " + name, e);
            }
        }

        static Method method(String name, Class<?>... parameters) {
            try {
                Method method = BBFile.class.getDeclaredMethod(name, parameters);
                method.setAccessible(true);
                return method;
            } catch (NoSuchMethodException e) {
                throw new IllegalStateException("BBFile has no method " + name, e);
            }
        }
    }
}
