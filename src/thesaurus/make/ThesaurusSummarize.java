/*
 * Copyright 2013 Tomasz Konopka.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package thesaurus.make;

import thesaurus.util.KilobaseIntArray;
import thesaurus.util.ThesaurusIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.file.RleWriter;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePosition;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.regions.GenomeBitSet;

/**
 * This tool reads files output by ThesaurusWrite. The function will output some
 * bed files showing regions that are described in the thesaurus with penalty 0,
 * 1, 2, etc.
 *
 *
 * @author tkonopka
 */
public class ThesaurusSummarize extends ThesaurusMapTool {

    private String output = "stdout";
    private String what = "align";
    private final ArrayList<File> input = new ArrayList<File>();
    private GenomeInfo ginfo;

    private void printWriteHelp() {
        System.out.println("GeneticThesaurus summarize: extract simple bed from thesaurus files");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar summarize ");
        System.out.println();

        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--thesaurus <File,File,...>", "thesaurus files, comma separated");
        ThesaurusIO.printHelpItem("--what <String>", "either 'align' or 'origin' or 'coverage'");

        System.out.println();
    }

    /**
     * Parse command line parameters.
     *
     * @param args
     * @return
     *
     * true if everything parsed correctly. false if any information is
     * missing/unreadable.
     *
     */
    private boolean parseWriteParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // change input genome
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // input and output 
        prs.accepts("thesaurus").withRequiredArg().ofType(String.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("what").withRequiredArg().ofType(String.class);

        // now use OptionSet to parse the command line
        OptionSet options;

        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        if (options.has("genome")) {
            genome = (File) options.valueOf("genome");
            if (!genome.exists() || !genome.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist");
                return false;
            }
        }

        // get prefix for output tracks
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        } else {
            System.out.println("Missing required argument --output");
            return false;
        }

        // get comma-separated input files
        if (options.has("thesaurus")) {
            String tempinput = (String) options.valueOf("thesaurus");
            String[] tempsplit = tempinput.split(",");
            for (int i = 0; i < tempsplit.length; i++) {
                File ff = new File(tempsplit[i]);
                if (!ff.exists() || !ff.canRead()) {
                    System.out.println("Cannot read thesaurus file, or file does not exist: " + tempsplit[i]);
                    return false;
                }
                input.add(ff);
            }
        } else {
            System.out.println("Missing required argument --thesaurus");
            return false;
        }

        //
        if (options.has("what")) {
            what = (String) options.valueOf("what");
            if (what.equalsIgnoreCase("align")) {
                what = "align";
            } else if (what.equalsIgnoreCase("origin")) {
                what = "origin";
            } else if (what.equalsIgnoreCase("coverage")) {
                what = "coverage";
            } else {
                System.out.println("unrecognized value for --what: " + what);
                return false;
            }
        }

        return true;
    }

    public ThesaurusSummarize(String[] args) {
        if (args.length == 0) {
            printWriteHelp();
            ginfo = null;
            return;
        }

        penalty = 0;
        super.loadDefaults();
        // parse parameters
        super.setOk(parseWriteParameters(args));
        if (!super.isOk()) {
            ginfo = null;
            return;
        }

        try {
            // get chromosome order and prevent further work if it does not work
            ginfo = new GenomeInfo(genome);
        } catch (IOException ex) {
            System.out.println("Failed reading genome information: " + ex.getMessage());
        }

    }

    private void runCoverage() {
        HashMap<String, KilobaseIntArray> cov = new HashMap<String, KilobaseIntArray>();
        for (int i = 0; i < ginfo.getNumChromosomes(); i++) {
            cov.put(ginfo.getChrName(i), new KilobaseIntArray(ginfo.getChrLength(i)));
        }
        for (int i = 0; i < input.size(); i++) {
            try {
                processCoverage(input.get(i), ginfo, cov);
            } catch (Exception ex) {
                System.out.println("Error processing coverage for file: " + input.get(i).getAbsolutePath());
            }
        }
        try {
            outputGenomeCov(cov, output);
        } catch (Exception ex) {
            System.out.println("Error saving coverage: " + ex.getMessage());
        }
    }

    private void runBed() {

        // make genome bitset
        GenomeBitSet genomebitset = null;
        // align bed or origin bed can be done in a simpler approach.
        genomebitset = new GenomeBitSet(ginfo);

        // process one file at a time. 
        try {
            if (what.equals("align")) {
                for (int i = 0; i < input.size(); i++) {
                    processAlign(input.get(i), genomebitset);
                }
            } else if (what.equals("origin")) {
                for (int i = 0; i < input.size(); i++) {
                    processOrigin(input.get(i), genomebitset);
                }
            }
        } catch (Exception ex) {
            System.out.println("something went wrong: " + ex.getMessage());
        }


        // make output streams
        OutputStream outstream;
        try {
            outstream = OutputStreamMaker.makeOutputStream(output);
        } catch (Exception ex) {
            System.out.println("Something went wrong when opening outstream: " + ex.getMessage());
            Logger.getLogger(ThesaurusSummarize.class.getName()).log(Level.SEVERE, null, ex);
            return;
        } 

        // output the bitset as a bed file
        outputGenomeBitSet(genomebitset, outstream);

        if (outstream != System.out) {
            try {
                outstream.close();
            } catch (IOException ex) {
                System.out.println("Something went wrong when closing outstream: " + ex.getMessage());
            }
        }
    }

    @Override
    public void runTool() {
        // coverage calculations
        if (what.equals("coverage")) {
            runCoverage();
        } else {
            runBed();
        }
    }

    private void processAlign(File infile, GenomeBitSet genomebitset) throws IOException {
        // open the reader
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(infile);

        String ss;
        while ((ss = br.readLine()) != null) {
            if (!ss.startsWith("#") && !ss.startsWith("Align.chr")) {
                String[] tokens = ss.split("\t", 4);
                int nowstart = Integer.parseInt(tokens[1]);
                int nowend = Integer.parseInt(tokens[2]);
                genomebitset.set(tokens[0], nowstart - 1, nowend, true);
            }
        }

        br.close();
    }

    private void processOrigin(File infile, GenomeBitSet genomebitset) throws IOException {
        // open the reader
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(infile);

        String ss;
        while ((ss = br.readLine()) != null) {
            if (!ss.startsWith("#") && !ss.startsWith("Align.chr")) {
                String[] tokens = ss.split("\t", 7);
                int nowstart = Integer.parseInt(tokens[4]);
                int nowend = Integer.parseInt(tokens[5]);
                genomebitset.set(tokens[3], nowstart - 1, nowend, true);
            }
        }

        br.close();
    }

    /**
     * Outputs all bitsets in bed format. Order of chromsomes in output is taken
     * from the arraylist.
     *
     * @param genomebitset
     * @param chrorder
     * @param outstream
     */
    private void outputGenomeBitSet(GenomeBitSet genomebitset, OutputStream outstream) {

        for (int i = 0; i < ginfo.getNumChromosomes(); i++) {
            StringBuilder sb = new StringBuilder();
            String nowchr = ginfo.getChrName(i);
            int chrlen = ginfo.getChrLength(i);
            BitSet nowbitset = genomebitset.get(nowchr, 0, chrlen);

            // look for set/unset regions and output them in bed format
            int nowpos = 0;
            while (nowbitset.nextSetBit(nowpos) >= 0) {
                int nextset = nowbitset.nextSetBit(nowpos);
                int nextunset = nowbitset.nextClearBit(nextset);
                if (nextunset < 0) {
                    nextunset = nowbitset.length();
                }
                sb.append(nowchr).append("\t").append(nextset).append("\t").append(nextunset).append("\n");
                nowpos = nextunset;
            }

            // output the chromsome bed to the stream
            try {
                if (sb.length() > 0) {
                    outstream.write(sb.toString().getBytes());
                }
            } catch (IOException ex) {
                Logger.getLogger(ThesaurusSummarize.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    /**
     * Output coverage of thesaurus into one file per chromosome
     *
     * @param cov
     * @param outdir
     * @throws FileNotFoundException
     * @throws IOException
     */
    private void outputGenomeCov(HashMap<String, KilobaseIntArray> cov, String outdir) throws FileNotFoundException, IOException {
        File outdirfile = new File(outdir);
        outdirfile.mkdirs();
        for (java.util.Map.Entry<String, KilobaseIntArray> entry : cov.entrySet()) {
            String nowchr = entry.getKey();
            // convert from kilobase array into normal array
            KilobaseIntArray kbarray = entry.getValue();
            int nowlen = kbarray.size();
            int[] array = new int[nowlen];
            for (int i = 0; i < nowlen; i++) {
                array[i] = kbarray.get(i);
            }
            // output the normal array as rle object
            OutputStream nowstream = OutputStreamMaker.makeOutputStream(new File(outdirfile, nowchr + ".txt.gz"));
            RleWriter.write(nowstream, array, true);
            nowstream.close();
        }
    }

    private void processCoverage(File infile, GenomeInfo ginfo, HashMap<String, KilobaseIntArray> cov) throws IOException {

        ThesaurusRegionsMap tmap = new ThesaurusRegionsMap(infile, ginfo, 0);
        GenomePositionComparator gcomp = new GenomePositionComparator();

        // look at one chromosome at a time
        int numchr = ginfo.getNumChromosomes();
        for (int i = 0; i < numchr; i++) {
            String nowchr = ginfo.getChrName(i);
            System.out.println(nowchr);
            int nowlen = ginfo.getChrLength(i);
            KilobaseIntArray chrcov = cov.get(nowchr);

            // compute coverage for each position
            for (int j = 0; j < nowlen; j++) {
                if (j % 10000000 == 0) {
                    System.out.println("  " + (j / 10000000));
                }
                // get thesaurus entries
                ThesaurusEntry[] tentries = tmap.lookup(nowchr, j + 1, 0);

                int herecov = 0;
                // get GenomicPositions of origin sites

                if (tentries != null) {
                    GenomePosition[] allpos = new GenomePosition[tentries.length];
                    for (int k = 0; k < tentries.length; k++) {
                        ThesaurusEntry te = tentries[k];
                        if (te.alignStrand == '+') {
                            allpos[k] = new GenomePosition(te.originChrIndex, (j - te.alignStart) + te.originStart);
                        } else {
                            allpos[k] = new GenomePosition(te.originChrIndex, te.originEnd - (j - te.alignStart));
                        }
                    }
                    // look for unique origin chr positions
                    Arrays.sort(allpos, gcomp);
                    for (int k = 0; k < allpos.length; k++) {
                        if (k == 0) {
                            herecov++;
                        } else {
                            if (gcomp.compare(allpos[k], allpos[k - 1]) != 0) {
                                herecov++;
                            }
                        }
                    }
                }

                if (herecov > 0) {
                    chrcov.set(j, herecov);
                }
            }
        }

    }

    private void processCoverageFast(File infile, GenomeInfo ginfo, HashMap<String, KilobaseIntArray> cov) throws IOException {

        BufferedReader inbr = BufferedReaderMaker.makeBufferedReader(infile);
        //GenomePositionComparator gcomp = new GenomePositionComparator();

        String s;
        while ((s = inbr.readLine()) != null) {
            if (!s.startsWith("#") && !s.startsWith("Align.chr")) {
                String[] tokens = s.split("\t", 4);
                String nowchr = tokens[0];
                int nowstart = Integer.parseInt(tokens[1]);
                int nowend = Integer.parseInt(tokens[2]);

                KilobaseIntArray nowcov = cov.get(nowchr);
                for (int i = nowstart - 1; i < nowend; i++) {
                    nowcov.increment(i);
                }
            }
        }
    }
}
