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

import thesaurus.util.ThesaurusLog;
import thesaurus.util.ThesaurusIO;
import thesaurus.util.ThesaurusSAMRecord;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.variants.VCFEntry;
import jsequtils.variants.VCFEntrySet;
import net.sf.samtools.SAMRecord;

/**
 * Helper utility. Reads a vcf file or otherwise a table with genomic locations.
 * Output a table with chr, position, and counts for A, T, C, G, N
 *
 *
 * @author tkonopka
 */
public class ThesaurusDetails extends ThesaurusMapTool {

    private String output = "stdout";
    private String validate = "STRICT";
    private File vcffile = null;
    private File bamfile = null;
    private int minmapqual = -1;
    private boolean verbose = true;
    private final ThesaurusLog mylog;

    /**
     * Data structure that holds a locus and counts for all bases
     */
    class VCFEntryDetails {

        VCFEntry entry;
        int countA, countT, countC, countG, countN;

        public VCFEntryDetails(VCFEntry entry) {
            this.entry = entry;
            countA = 0;
            countT = 0;
            countC = 0;
            countG = 0;
            countN = 0;
        }

        public int getPosition() {
            return entry.getPosition();
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(entry.getChr()).append("\t").append(entry.getPosition())
                    .append("\t").append(countA)
                    .append("\t").append(countT)
                    .append("\t").append(countC)
                    .append("\t").append(countG)
                    .append("\t").append(countN).append("\n");
            return sb.toString();
        }

        public String getHeader() {
            return "chr\tposition\tA\tT\tC\tG\tN\n";
        }
    }

    private void printDetailsHelp() {
        System.out.println("GeneticThesaurus details: get counts for ATCG bases at genomic loci");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar details ");
        System.out.println();
        ThesaurusIO.printHelpItem("--bam <File>", "alignment files matching --vcf");
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--minmapqual <int>", "minimum mapping quality");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--vcf <File>", "variant call file matching --bam");

        System.out.println();
    }

    /**
     *
     * @param args
     * @return
     */
    private boolean parseDetailsParameters(String[] args) {

        OptionParser prs = new OptionParser();

        // change input genome
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // input and output         
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("minmapqual").withRequiredArg().ofType(Integer.class);
        prs.accepts("vcf").withRequiredArg().ofType(File.class);
        prs.accepts("bam").withRequiredArg().ofType(File.class);
        prs.accepts("validate").withRequiredArg().ofType(String.class);

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

        // get prefix for output tracks
        if (options.has("validate")) {
            validate = (String) options.valueOf("validate");
            if (!validate.equals("STRICT") && !validate.equals("LENIENT") && !validate.equals("SILENT")) {
                System.out.println("Option --validate must be set to STRICT, LENIENT, or SILENT");
                return false;
            }
        }

        // get input vcf file
        if (options.has("vcf")) {
            vcffile = (File) options.valueOf("vcf");
            if (!vcffile.exists() || !vcffile.canRead()) {
                System.out.println("Cannot read vcf file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument --vcf");
            return false;
        }

        // get minimum mapping quality
        if (options.has("minmapqual")) {
            minmapqual = (Integer) options.valueOf("minmapqual");
            if (minmapqual < 0) {
                System.out.println("Minimum mapping quality must be >= 0");
                return false;
            }
        }

        if (options.has("bam")) {
            bamfile = (File) options.valueOf("bam");
            if (!bamfile.exists() || !bamfile.canRead()) {
                System.out.println("Cannot read bam file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument --bam");
            return false;
        }

        return true;
    }

    public ThesaurusDetails(String[] args) {
        if (args.length == 0) {
            printDetailsHelp();
            mylog = null;
            return;
        }
        super.loadDefaults();
        super.setOk(parseDetailsParameters(args));
        mylog = new ThesaurusLog(System.out);
        mylog.setVerbose(verbose);
    }

    @Override
    void runTool() {
        // set up input and output streams
        OutputStream outdetails; // output 

        try {
            outdetails = OutputStreamMaker.makeOutputStream(output + ".details.txt.gz");
        } catch (Exception ex) {
            System.out.println("Something went wrong during stream setup: " + ex.getMessage());
            return;
        }

        // read all the variants into memory        
        mylog.log("Reading variants: " + vcffile.getAbsolutePath());
        GenomeInfo ginfo;
        try {
            ginfo = new GenomeInfo(genome);
        } catch (Exception ex) {
            System.out.println("Could not look up genome info: " + ex.getMessage());
            return;
        }
        VCFEntrySet allvariants = new VCFEntrySet(vcffile, ginfo, false);

        // perform the filtering in another function
        try {
            mylog.log("Looking at the details");
            getDetails(allvariants, bamfile, outdetails);
        } catch (Exception ex) {
            System.out.println("Something went wrong during detail extraction: " + ex.getMessage());
        }

        // finish up, close all the streams
        try {
            outdetails.close();
        } catch (Exception ex) {
            System.out.println("Something went wrong during stream closing: " + ex.getMessage());
        }

        mylog.log("done");
    }

    private void getDetails(VCFEntrySet variants, File bamfile, OutputStream outdetails) throws IOException {

        BamRegionsMap bamregions = new BamRegionsMap(bamfile, new GenomeInfo(genome), "SILENT", minmapqual);

        outdetails.write((new VCFEntryDetails(null)).getHeader().getBytes());

        int numvariants = variants.size();
        System.out.println("Number of variants " + numvariants);
        for (int i = 0; i < numvariants; i++) {
            VCFEntry entry = variants.getVariant(i);

            SAMRecord[] varinbam = bamregions.lookup(entry.getChr(), entry.getPosition());
            VCFEntryDetails entrydetails = new VCFEntryDetails(entry);

            // collect information about the locus from each bam record
            countDetails(entrydetails, varinbam);

            outdetails.write(entrydetails.toString().getBytes());

            if ((i % 200000) == 0) {
                mylog.log(entry.getChr() + ":" + entry.getPosition());
            }
        }

    }

    private void countDetails(VCFEntryDetails entrydetails, SAMRecord[] varinbam) {

        if (varinbam == null) {
            return;
        }

        int vsize = varinbam.length;
        for (int i = 0; i < vsize; i++) {
            ThesaurusSAMRecord tr = new ThesaurusSAMRecord(varinbam[i]);
            byte nowbase = tr.getBaseAtGenomicPosition(entrydetails.getPosition());
            switch (nowbase) {
                case 'A':
                    entrydetails.countA++;
                    break;
                case 'T':
                    entrydetails.countT++;
                    break;
                case 'C':
                    entrydetails.countC++;
                    break;
                case 'G':
                    entrydetails.countG++;
                    break;
                case 'N':
                    entrydetails.countN++;
                    break;
                default:
                    break;
            }
        }
    }
}
