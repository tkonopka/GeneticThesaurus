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
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.variants.VCFEntry;
import jsequtils.variants.VCFEntrySet;

/**
 * Tool to compare VCF files, one of which is annotated via a thesaurus object
 *
 *
 * @author tkonopka
 */
public class ThesaurusCompare extends ThesaurusMapTool {

    // prefix for output files
    private String output = "stdout";
    // "actual" variants file
    private File refvcffile = null;
    // predicted variants with thesaurus
    private File vcffile = null;
    private File synonymsfile = null;
    private boolean verbose = true;
    private final ThesaurusLog mylog;

    private void printCompareHelp() {
        System.out.println("GeneticThesaurus compare: compare two VCF files, one of which is annotated with thesaurus");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar compare ");
        System.out.println();
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--ref <File>", "vcf file with true variants");
        ThesaurusIO.printHelpItem("--vcf <File>", "variant call file");
        ThesaurusIO.printHelpItem("--vtf <File>", "synonyms for loci in vcf file (optional)");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");

        System.out.println();
    }

    /**
     *
     * @param args
     * @return
     */
    private boolean parseCompareParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // change input genome
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // input and output                 
        prs.accepts("ref").withRequiredArg().ofType(File.class);
        prs.accepts("vcf").withRequiredArg().ofType(File.class);
        prs.accepts("vtf").withRequiredArg().ofType(File.class);
        // output
        prs.accepts("output").withRequiredArg().ofType(String.class);

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

        // get input vcf file
        if (options.has("vcf")) {
            vcffile = (File) options.valueOf("vcf");
            if (!vcffile.exists() || !vcffile.canRead()) {
                System.out.println("Cannot read vcf file, or file does not exist: " + vcffile.getAbsolutePath());
                return false;
            }
        } else {
            System.out.println("Missing required argument --vcf");
            return false;
        }

        // get input vcf file
        if (options.has("vtf")) {
            synonymsfile = (File) options.valueOf("vtf");
            if (!synonymsfile.exists() || !synonymsfile.canRead()) {
                System.out.println("Cannot read synonyms file, or file does not exist");
                return false;
            }
        }

        // get reference vcf file
        if (options.has("ref")) {
            refvcffile = (File) options.valueOf("ref");
            if (!refvcffile.exists() || !refvcffile.canRead()) {
                System.out.println("Cannot read  reference vcf file, or file does not exist: " + refvcffile.getAbsolutePath());
                return false;
            }
        } else {
            System.out.println("Missing required argument --ref");
            return false;
        }

        return true;
    }

    public ThesaurusCompare(String[] args) {
        if (args.length == 0) {
            printCompareHelp();
            mylog = null;
            return;
        }
        super.loadDefaults();
        super.setOk(parseCompareParameters(args));
        mylog = new ThesaurusLog(System.out);
        mylog.setVerbose(verbose);
    }

    @Override
    void runTool() {
        // load the reference variant set
        GenomePositionComparator vcfcomp = null;
        GenomeInfo ginfo = null;
        try {
            ginfo = new GenomeInfo(genome);
            vcfcomp = new GenomePositionComparator();
        } catch (Exception ex) {
            System.out.println("Failed to create genomic position comparator: " + ex.getMessage());
            return;
        }
        VCFEntrySet refset = new VCFEntrySet(refvcffile, ginfo, false);

        try {
            compareRefVcf(refset, vcffile, synonymsfile, output, ginfo);
        } catch (Exception ex) {
            System.out.println("Something went wrong in compareRefVcf: " + ex.getMessage());
        }
    }

    private String vcfEntryString(VCFEntry entry) {
        return entry.getChr() + ":" + entry.getPosition();
    }

    private String logExplanation() {
        StringBuilder sb = new StringBuilder();
        sb.append("#Log for Thesaurus compare\n");
        sb.append("#chr and position identify location of variant annotated with thesaurus filter code\n");
        sb.append("#NumSynonyms shows number of alternative mapping sites.\n");
        sb.append("#   Value 1 indicates variant can be at the chr/position locus, or one other location.\n");
        sb.append("#TrueSynonym shows the index of the true synonym.\n");
        sb.append("#   Value 0 indicates the declared position is the true site.\n");
        sb.append("#   Values 1 or greater indicate the chr/position locus is false, but the synonym with indicated index is true.\n");
        sb.append("#   Value -1 indicates the whole set if false.\n");
        sb.append("chr\tposition\tNumSynonyms\tTrueSynonym\n");
        return sb.toString();
    }

    /**
     *
     * @param filterfield
     *
     * the string in the filter column in a vcf file, i.e. a series of filter
     * names separated by commas.
     *
     * @param filtername
     *
     * a certain substring
     *
     * @return
     *
     * true if the filterfield contain filtername as a filter
     *
     */
    private boolean hasFilter(String filterfield, String filtername) {
        String[] tokens = filterfield.split(";");
        for (int i = 0; i < tokens.length; i++) {
            if (filtername.equals(tokens[i])) {
                return true;
            }
        }
        return false;
    }

    public void compareRefVcf(VCFEntrySet refset, File vcffile, File synonymsfile, String output, GenomeInfo ginfo) throws IOException {

        // make a hashmap with all the refset positions, this will keep track of the false negatives
        // at the beginning of the calculations, all variants are false negatives
        HashMap<String, Boolean> FNmap = new HashMap<String, Boolean>(refset.size() + (refset.size() / 5));
        int rssize = refset.size();
        for (int i = 0; i < rssize; i++) {
            FNmap.put(vcfEntryString(refset.getVariant(i)), true);
        }

        int countTP = 0, countFP = 0, countFN = 0, countTPthes = 0;

        String vcffileheader = "";
        try {
            vcffileheader = getVcfHeader(vcffile);
        } catch (Exception ex) {
            System.out.println("Error getting vcf file header: " + ex.getMessage());
            return;
        }

        // readers for the vcf files and associated synonyms file
        BufferedReader vcfreader = BufferedReaderMaker.makeBufferedReader(vcffile);
        String vcfline = skipHeader(vcfreader);

        BufferedReader synonymsreader = null;
        String synline = null;
        if (synonymsfile != null) {
            synonymsreader = BufferedReaderMaker.makeBufferedReader(synonymsfile);
            synline = skipHeader(synonymsreader);
        }

        // these streams will show variants of the TP,FP,etc types
        OutputStream outTP = OutputStreamMaker.makeOutputStream(output + ".compare.TP.vcf.gz");
        OutputStream outFP = OutputStreamMaker.makeOutputStream(output + ".compare.FP.vcf.gz");
        // the TPthesaurus will have variants that are "true" only by nature of their synonyms
        OutputStream outTPthes = OutputStreamMaker.makeOutputStream(output + ".compare.TPthesaurus.vcf.gz");
        OutputStream outlog = OutputStreamMaker.makeOutputStream(output + ".compare.log.gz");

        // write out some header for the output files
        outTP.write(vcffileheader.getBytes());
        outFP.write(vcffileheader.getBytes());
        outTPthes.write(vcffileheader.getBytes());
        outlog.write(logExplanation().getBytes());

        while (vcfline != null) {
            VCFEntry nowentry = new VCFEntry(vcfline, ginfo);
            if (!nowentry.isIndel()) {

                // check if it is a TP
                boolean nowTP = false;
                if (refset.containsPosition(nowentry)) {
                    outTP.write(nowentry.toString().getBytes());
                    FNmap.put(vcfEntryString(nowentry), false);
                    nowTP = true;
                    countTP++;
                }

                // check for thesaurus status anyway (to keep synonyms on track)                
                if (hasFilter(nowentry.getFilter(), "thesaurus") && synonymsfile != null) {
                    // check that the synonyms are synced with the vcf file
                    String entrycode = vcfEntryString(nowentry);
                    String[] synonyms = synline.split("\t");
                    if (!synonyms[0].equals(entrycode)) {
                        System.out.println("Desynced vcf and synonyms files: " + synonyms[0] + " " + entrycode);
                        return;
                    }

                    // try if the synonyms are in the reference set
                    int synlen = synonyms.length;
                    int syntrue = -1;
                    boolean nowTPthes = false;
                    for (int i = 0; i < synlen; i++) {
                        if (refset.containsPosition(synonyms[i])) {
                            FNmap.put(synonyms[i], false);
                            syntrue = i;
                            nowTPthes = true;
                            i = synlen;
                        }
                    }

                    if (!nowTP) {
                        if (nowTPthes) {
                            outTPthes.write(nowentry.toString().getBytes());
                            countTPthes++;
                        } else {
                            outFP.write(nowentry.toString().getBytes());
                            countFP++;
                        }
                    }

                    String temp = nowentry.getChr() + "\t" + nowentry.getPosition() + "\t" + (synlen - 1) + "\t" + syntrue + "\n";
                    outlog.write(temp.getBytes());

                    synline = synonymsreader.readLine();

                } else {
                    if (!nowTP) {
                        outFP.write(nowentry.toString().getBytes());
                        countFP++;
                    }
                }
            }
            // read in the next item from the vcf file
            vcfline = vcfreader.readLine();
        }

        outTP.close();
        outFP.close();
        outlog.close();
        outTPthes.close();
        vcfreader.close();
        if (synonymsreader != null) {
            synonymsreader.close();
        }

        // now write out the false negatives. These are the items that are still marked
        // as true in the FNmap
        OutputStream outFN = OutputStreamMaker.makeOutputStream(output + ".compare.FN.vcf.gz");
        outFN.write(vcffileheader.getBytes());
        for (int i = 0; i < rssize; i++) {
            VCFEntry nowentry = refset.getVariant(i);
            boolean isFN = FNmap.get(vcfEntryString(nowentry));
            if (isFN) {
                outFN.write(nowentry.toString().getBytes());
                countFN++;
            }
        }
        outFN.close();

        // print out a summary of the TP, FP, etc
        System.out.println("TP:\t" + countTP);
        System.out.println("FP:\t" + countFP);
        System.out.println("FN:\t" + countFN);
        System.out.println("TPthes:\t" + countTPthes);

    }

    /**
     * looks into a buffered and skips the header rows. Returns the first
     * non-header line.
     *
     * @param br
     * @return
     */
    private String skipHeader(BufferedReader br) throws IOException {
        String s = br.readLine();
        while (s != null) {
            if (!s.startsWith("#")) {
                return s;
            }
            s = br.readLine();
        }
        return null;
    }

    private String getVcfHeader(File f) throws IOException {

        StringBuilder sb = new StringBuilder();
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(f);

        String s = br.readLine();
        while (s != null && s.startsWith("#")) {
            sb.append(s).append("\n");
            s = br.readLine();
        }
        return sb.toString();
    }
}
