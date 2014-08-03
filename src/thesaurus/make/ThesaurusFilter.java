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

import thesaurus.GeneticThesaurus;
import thesaurus.util.ThesaurusLog;
import thesaurus.util.ThesaurusIO;
import thesaurus.util.SNVPosition;
import thesaurus.util.SNVPositionDetails;
import thesaurus.util.ThesaurusSAMRecord;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.genome.GenomePositionInterface;
import jsequtils.sequence.FastaReader;
import jsequtils.variants.VCFEntry;
import jsequtils.variants.VCFEntrySet;
import net.sf.samtools.SAMRecord;

/**
 * One of main programs under GeneticThesaurus. Supposed to concurrently scan a
 * thesaurus file, a vcf file, and a bam file with the purpose to annotate
 * variants.
 *
 * @author tkonopka
 */
public class ThesaurusFilter extends ThesaurusMapTool {

    // constants describing how to deal with information from bam
    private static final int METHOD_FULL = 0;
    private static final int METHOD_BAF = 1;
    // variables collected from the command line
    private String output = "stdout";
    private String validate = "STRICT";
    private File thesaurusfile = null;
    private File vcffile = null;
    private File bamfile = null;
    // method - full of baf
    // set to full to use bam file when computing synonyms
    // set to baf to use bam file only to compute the baf
    private int method = METHOD_FULL;
    private int minmapqual = -1;
    private boolean verbose = true;
    private final ThesaurusLog mylog;
    // thresholds used during filtering
    private int tolerance = 1;
    private int maxtolerance = 2;
    private int softclip = 5;
    private int hitcount = 3;
    private double hitproportion = 0.5;
    private int many = 100;
    private int toomany = 400;
    private int insertsize = 200;

    private void printFilterHelp() {
        System.out.println("GeneticThesaurus filter: use a thesaurus to filter a VCF file");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar filter ");
        System.out.println();
        System.out.println("Core options:");
        ThesaurusIO.printHelpItem("--bam <File>", "alignment files matching --vcf");
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        //ThesaurusIO.printHelpItem("--method <String>", "choose method to use information from bam file (values baf, full; default full)");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--thesaurus <File>", "thesaurus file");
        ThesaurusIO.printHelpItem("--vcf <File>", "variant call file matching --bam");
        System.out.println("\nFiltering details:");
        ThesaurusIO.printHelpItem("--insertsize <int>", "insert size for paired-end reads [default " + insertsize + "]");
        ThesaurusIO.printHelpItem("--clip <int>", "number of bases to clip from read ends [default " + softclip + "]");
        ThesaurusIO.printHelpItem("--minmapqual <int>", "minimum mappinq quality [default " + minmapqual + "]");
        ThesaurusIO.printHelpItem("--readlen <int>", "read length [default " + readlen + "]");
        ThesaurusIO.printHelpItem("--tolerance <int>", "number of errors in read to tolerate [default " + tolerance + "]");
        ThesaurusIO.printHelpItem("--maxtolerance <int>", "maximum errors in read to tolerate [default " + maxtolerance + "]");
        ThesaurusIO.printHelpItem("--hitcount <int>", "number of reads supporting thesaurus link [default " + hitcount + "]");
        ThesaurusIO.printHelpItem("--hitproportion <int>", "proportion of reads supporting variant that also thesaurus link [default " + hitproportion + "]");
        ThesaurusIO.printHelpItem("--many <int>", "decrease tolerance in areas with too many hits [default " + many + "]");
        ThesaurusIO.printHelpItem("--toomany <int>", "abandon variants with too many hits [default " + toomany + "]");
        System.out.println();
    }

    /**
     *
     * @param args
     * @return
     */
    private boolean parseFilterParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // change input genome
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // input and output 
        prs.accepts("thesaurus").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("method").withRequiredArg().ofType(String.class);
        prs.accepts("minmapqual").withRequiredArg().ofType(Integer.class);
        prs.accepts("vcf").withRequiredArg().ofType(File.class);
        prs.accepts("bam").withRequiredArg().ofType(File.class);

        prs.accepts("tolerance").withRequiredArg().ofType(Integer.class);
        prs.accepts("insertsize").withRequiredArg().ofType(Integer.class);
        prs.accepts("readlen").withRequiredArg().ofType(Integer.class);
        prs.accepts("maxtolerance").withRequiredArg().ofType(Integer.class);
        prs.accepts("many").withRequiredArg().ofType(Integer.class);
        prs.accepts("toomany").withRequiredArg().ofType(Integer.class);
        prs.accepts("hitcount").withRequiredArg().ofType(Integer.class);
        prs.accepts("hitproportion").withRequiredArg().ofType(Double.class);
        prs.accepts("validate").withRequiredArg().ofType(String.class);
        prs.accepts("clip").withRequiredArg().ofType(Integer.class);

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

        if (options.has("minmapqual")) {
            minmapqual = (Integer) options.valueOf("minmapqual");
            if (minmapqual < 0) {
                System.out.println("Parameter minmapqual must be >= 0");
                return false;
            }
        }

        if (options.has("softclip")) {
            softclip = (Integer) options.valueOf("clip");
            if (softclip < 0) {
                System.out.println("Parameter clip must be >= 0");
                return false;
            }
        }
        if (options.has("readlen")) {
            readlen = (Integer) options.valueOf("readlen");
            if (readlen < 0) {
                System.out.println("Parameter readlen must be >= 0");
                return false;
            }
        }
        if (options.has("insertsize")) {
            insertsize = (Integer) options.valueOf("insertsize");
            if (insertsize < readlen) {
                System.out.println("Parameter insertsize must be >= readlen");
                return false;
            }
        }


        if (options.has("tolerance")) {
            tolerance = (Integer) options.valueOf("tolerance");
            if (tolerance < 0) {
                System.out.println("Parameter tolerance must be >= 0");
                return false;
            }
        }
        if (options.has("maxtolerance")) {
            maxtolerance = (Integer) options.valueOf("maxtolerance");
            if (maxtolerance < tolerance) {
                System.out.println("Parameter maxtolerance must be >= value of parameter tolerance");
                return false;
            }
        }
        if (options.has("many")) {
            many = (Integer) options.valueOf("many");
            if (many < 2) {
                System.out.println("Parameter many must be > 1");
                return false;
            }
        }
        if (options.has("toomany")) {
            toomany = (Integer) options.valueOf("toomany");
            if (toomany < 2) {
                System.out.println("Parameter toomany must be > 1");
                return false;
            }
        }

        // thresholds deciding which thesaurus links to accept
        if (options.has("hitcount")) {
            hitcount = (Integer) options.valueOf("hitcount");
            if (hitcount < 1) {
                System.out.println("Parameter hitcount must be >= 1");
                return false;
            }
        }
        if (options.has("hitproportion")) {
            hitproportion = (Double) options.valueOf("hitproportion");
            if (hitproportion < 0) {
                System.out.println("Parameter hitproportion must be >= 0");
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

        if (options.has("bam")) {
            bamfile = (File) options.valueOf("bam");
            if (!bamfile.exists() || !bamfile.canRead()) {
                System.out.println("Cannot read bam file, or file does not exist");
                return false;
            }
        }

        // get input vcf file
        if (options.has("thesaurus")) {
            thesaurusfile = (File) options.valueOf("thesaurus");
            if (!thesaurusfile.exists() || !thesaurusfile.canRead()) {
                System.out.println("Cannot read thesaurus file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument --thesaurus");
            return false;
        }

        if (options.has("method")) {
            String methodString = (String) options.valueOf("method");
            if (methodString.equalsIgnoreCase("full")) {
                method = METHOD_FULL;
            } else if (methodString.equalsIgnoreCase("baf")) {
                method = METHOD_BAF;
            } else {
                System.out.println("Unrecognized value for parameter method - " + methodString);
                System.out.println("Allowed values are baf and full");
                return false;
            }
        }

        return true;
    }

    public ThesaurusFilter(String[] args) {
        if (args.length == 0) {
            printFilterHelp();
            mylog = null;
            return;
        }
        super.loadDefaults();
        super.setOk(parseFilterParameters(args));
        mylog = new ThesaurusLog(System.out);
        mylog.setVerbose(verbose);
    }

    @Override
    void runTool() {

        // set up input and output streams
        OutputStream outvcf; // output variant call file file (vcf)
        OutputStream outvtf; // output variant thesaurus file (vtf)          

        // ------------------ PASS 1
        try {
            outvcf = OutputStreamMaker.makeOutputStream(output + ".vcf.gz");
            outvtf = OutputStreamMaker.makeOutputStream(output + ".vtf.gz");
        } catch (Exception ex) {
            System.out.println("Something went wrong during stream setup: " + ex.getMessage());
            return;
        }

        // read all the variants into memory        
        mylog.log("Reading variants");
        GenomeInfo ginfo = null;
        GenomePositionComparator vcomp = new GenomePositionComparator();
        try {
            ginfo = new GenomeInfo(genome);
        } catch (Exception ex) {
            System.out.println("Failed to load genomic position comparator: " + ex.getMessage());
            return;
        }
        VCFEntrySet allvariants = new VCFEntrySet(vcffile, ginfo, true);
        SNVPositionNetwork allnetwork = null;

        // perform the filtering in another function
        try {
            mylog.log("Filtering variants using thesaurus");
            if (bamfile == null) {
                mylog.log("Filtering without bam file - not specified on command line");
            }
            allnetwork = filterVcfWithBamFile(allvariants, ginfo, outvcf, outvtf);
        } catch (Exception ex) {
            System.out.println("Something went wrong during filtering: " + ex.getMessage());
        }

        // finish up, close all the streams
        try {
            outvcf.close();
            outvtf.close();
        } catch (Exception ex) {
            System.out.println("Something went wrong during stream closing: " + ex.getMessage());
        }

        // ----------------- PASS 2
        // need to read the variants and synonyms        
        OutputStream outBAF;

        try {
            outBAF = OutputStreamMaker.makeOutputStream(output + ".baf.tsv.gz");
        } catch (Exception ex) {
            System.out.println("Something went wrong during BAF stream setup: " + ex.getMessage());
            return;
        }

        // Scan the bam file again and collect allele frequency information about the variants.
        try {
            mylog.log("Computing BAF for all variants, with and without thesaurus annotation");
            getBAFs(allvariants, allnetwork, ginfo, bamfile, outBAF);
        } catch (Exception ex) {
            System.out.println("Something went wrong during BAF computation: " + ex.getMessage());
        }

        // finish up, close all the streams
        try {
            outBAF.close();
        } catch (Exception ex) {
            System.out.println("Something went wrong during BAF stream closing: " + ex.getMessage());
        }

        mylog.log("done");
    }

    /**
     * Function collects coverage information on variants/loci of interest and
     * computes the BAFs for those positions
     *
     * @param variants
     * @param synonyms
     * @param vcomp
     * @param bamfile
     * @param outBAF
     * @throws IOException
     */
    private void getBAFs(VCFEntrySet variants,
            SNVPositionNetwork synonyms, GenomeInfo ginfo,
            File bamfile, OutputStream outBAF) throws IOException {

        // make a list of all variats in the synonyms map
        ArrayList<SNVPositionDetails> allloci = makeListOfAllLoci(variants, synonyms, ginfo);

        // in the constructor for bamregions, I had the minimum mapping quality set at "minmapqual"
        // but this only worked well for minmapqual<=1. So now I set this to zero regardless of 
        // the mapping quality used in variant calling.
        mylog.log("Collecting coverage information from bam file");
        BamRegionsMap bamregions = new BamRegionsMap(bamfile, new GenomeInfo(genome), "SILENT", 0);
        int oldchr = -1;
        int numloci = allloci.size();
        for (int i = 0; i < numloci; i++) {
            SNVPositionDetails entry = allloci.get(i);
            int nowchr = entry.getChrIndex();
            // output a log line at every chromosome
            if (nowchr != oldchr) {
                oldchr = nowchr;
            }
            // look up the position in the bam
            SAMRecord[] varinbam = bamregions.lookup(ginfo.getChrName(nowchr), entry.getPosition());
            // collect information about the locus from each bam record
            countDetails(entry, varinbam);
        }
        bamregions.close();

        GenomePositionComparator vcomp = new GenomePositionComparator();

        mylog.log("Writing baf table");
        prepareBAFStream(outBAF);
        int vs = variants.size();
        for (int i = 0; i < vs; i++) {
            VCFEntry entry = variants.getVariant(i);
            SNVPosition entrypos = new SNVPosition(entry, ginfo);

            // collect information about this location only
            SNVPositionDetails nowGPD = getGPD(allloci, entry, vcomp);

            String ref = entry.getRef();
            String alt = entry.getAlt();

            StringBuilder sb = new StringBuilder();
            sb.append(entry.getChr()).append("\t").append(entry.getPosition()).append("\t").append(ref).append("\t").append(alt);

            // collect information about this location only                
            sb.append("\t").append(((double) nowGPD.getCountAlt()) / (double) nowGPD.getCountATCG());

            // add in all the synonym regions
            ArrayList<SNVPosition> nowsynonyms = synonyms.getNeighborNodes(entrypos);
            int numsynonyms;

            if (nowsynonyms != null) {
                numsynonyms = nowsynonyms.size();
                for (int k = 0; k < numsynonyms; k++) {
                    nowGPD.incrementCounts(getGPD(allloci, nowsynonyms.get(k), vcomp));
                }
            } else {
                numsynonyms = 0;
            }

            // output the BAF including the synonyms
            sb.append("\t").append(((numsynonyms + 1) * (double) nowGPD.getCountAlt()) / (double) nowGPD.getCountATCG());

            // output the total coverage
            sb.append("\t").append(nowGPD.getCountATCG());

            // output the number of synonyms
            if (nowsynonyms == null) {
                sb.append("\t0\n");
            } else {
                sb.append("\t").append(nowsynonyms.size()).append("\n");
            }

            outBAF.write(sb.toString().getBytes());
        }

    }

    /**
     * find an entry in an array using binary search. This is just a wrapper for
     * Collections.binarySearch
     *
     * @param gpdlist
     *
     * @param gpdfind
     *
     * @param vcomp
     *
     * @return
     *
     */
    private SNVPositionDetails getGPD(ArrayList<SNVPositionDetails> gpdlist, GenomePositionInterface gpdfind, GenomePositionComparator vcomp) {
        int a = Collections.binarySearch(gpdlist, gpdfind, vcomp);
        if (a < 0) {
            return null;
        } else {
            return gpdlist.get(a);
        }
    }

    /**
     * Convert between a structure of positions in a network form to an array.
     *
     * @param variants
     * @param network
     * @param ginfo
     * @return
     */
    private ArrayList<SNVPositionDetails> makeListOfAllLoci(VCFEntrySet variants,
            SNVPositionNetwork network, GenomeInfo ginfo) {

        // make array of loci  from the network. 
        // Important: getAllNodes() will return the positions in sorted order.
        ArrayList<SNVPosition> allloci = network.getAllNodes();
        int allsize = allloci.size();

        // make array with details
        ArrayList<SNVPositionDetails> alldetails = new ArrayList<SNVPositionDetails>(allloci.size());
        for (int i = 0; i < allsize; i++) {
            alldetails.add(new SNVPositionDetails(allloci.get(i)));
        }

        return alldetails;
    }

    /**
     * runs through the thesaurus regions and records alignment regions into a
     * bitset.
     *
     *
     * @param chrbitset
     * @param thesregions
     */
    private void updateBitSetWithThesEntries(BitSet chrbitset, ThesaurusEntry[] thesregions) {
        for (int i = 0; i < thesregions.length; i++) {
            ThesaurusEntry te = thesregions[i];
            chrbitset.set(te.alignStart - 1, te.alignEnd);
        }
    }

    /**
     * Perform filtering using bam file as a help for inferring alternate
     * variant loci
     *
     * @param calledvariants
     *
     * set of variants to be filtered
     *
     * @param vcomp
     *
     * comparator for the current genome
     *
     * @param outvcf
     *
     * output vcf will contain the same variants as the called variants, except
     * that some will be annotated with a thesaurus field and an extra label in
     * the sample info field.
     *
     * @param outvtf
     *
     * output vtf will be an auxilliary file holding location of alternate loci
     *
     * @return
     *
     * a hashmap with chromosome as key. The values will contain objects
     * describing all the called variants and all the alternate loci.
     *
     * @throws IOException
     */
    private SNVPositionNetwork filterVcfWithBamFile(VCFEntrySet calledvariants,
            GenomeInfo ginfo,
            OutputStream outvcf, OutputStream outvtf) throws IOException {

        // copy the header from the input to the output vcf and add filter line
        prepareVcfStreams(calledvariants, outvcf);
        prepareVtfStream(outvtf);

        // open readers for all the input files        
        final ThesaurusRegionsMap thesregions = new ThesaurusRegionsMap(thesaurusfile, ginfo, insertsize);
        final BamRegionsMap bamregions = new BamRegionsMap(bamfile, ginfo, "SILENT", minmapqual);
        final FastaReader genomereader = new FastaReader(genome);

        //
        SNVPositionNetwork synnetwork = new SNVPositionNetwork(calledvariants, ginfo);
        String oldchr = "NA";

        // a chromosome bitset
        BitSet chrbitset = null;

        int numvariants = calledvariants.size();
        for (int i = 0; i < numvariants; i++) {
            VCFEntry entry = calledvariants.getVariant(i);

            if (entry.isIndel()) {
                // skip processing indels, but still change the filter and format fields
                adjustVcfEntryFilter(entry, true, "thesaurus");
                adjustVcfEntryFormatGenotype(entry, 0);
            } else {
                // variant is a substitutions - this is the interesting case here

                // write a log line for every new chromosome
                if (!entry.getChr().equals(oldchr)) {
                    oldchr = entry.getChr();
                    mylog.log(oldchr + " (synonyms network: " + synnetwork.size() + ")");

                    // load sequences from the genome until find sequence for the current chromosome
                    while (genomereader.hasNext()) {
                        genomereader.readNext();
                        if (genomereader.getChromosomeName().equals(oldchr)) {
                            break;
                        }
                    }
                    if (!genomereader.getChromosomeName().equals(oldchr)) {
                        mylog.log("Chromosome not in genome fasta (check order of chromosomes?)");
                        break; // this breaks the for loop over variants
                    }
                    // make a new bitset for the chromosome
                    chrbitset = new BitSet(genomereader.getChromosomeLength());
                }


                // get entries in a wider window to update the bitset
                ThesaurusEntry[] thesLinesOnLocus = thesregions.lookup(entry.getChr(), entry.getPosition(), insertsize);
                if (thesLinesOnLocus != null) {
                    updateBitSetWithThesEntries(chrbitset, thesLinesOnLocus);
                }

                // get entries in a narrower window around the variant - this will be used in analysis
                thesLinesOnLocus = thesregions.lookup(entry.getChr(), entry.getPosition(), readlen);
                // get also the bam file information
                SAMRecord[] readsOnLocus = bamregions.lookup(entry.getChr(), entry.getPosition());

                // in output vcf file, signal there are synonymous entries via the filter field
                if (thesLinesOnLocus == null) {
                    adjustVcfEntryFilter(entry, true, "thesaurus");
                    adjustVcfEntryFormatGenotype(entry, 0);
                } else {

                    // collect some information about the synonyms (number and loci)
                    ThesaurusSynonyms tsyns = new ThesaurusSynonyms(thesLinesOnLocus, entry, ginfo, calledvariants);
                    ThesaurusSAMRecord[] tbamrecords = null;
                    // when method is full, make an array of bam 
                    if (method == METHOD_FULL && readsOnLocus != null) {
                        // convert the SAMRecords into ThesaurusSAMRecords                    
                        tbamrecords = makeClippedRecords(readsOnLocus, genomereader, softclip, entry.getPosition());
                    }
                    // calculate the alternate loci
                    ArrayList<SNVPosition> synonyms;
                    try {
                        synonyms = tsyns.findVariants(tbamrecords, hitcount, hitproportion, tolerance, maxtolerance, toomany, chrbitset);
                    } catch (Exception ex) {
                        System.out.println("Error while finding synonyms for " + entry.getChr() + ":" + entry.getPosition() + " " + ex.getMessage());
                        break;
                    }

                    // synonyms can be either null (difficult) or not null, in which case could process
                    if (synonyms == null) {
                        adjustVcfEntryFilter(entry, false, "thesaurushard");
                        adjustVcfEntryFormatGenotype(entry, 0);
                    } else {
                        int numsynonyms = synonyms.size();
                        // if there are too many, give up
                        if (numsynonyms > toomany || containsLinkToUnaligned(synonyms)) {
                            adjustVcfEntryFilter(entry, numsynonyms == 0, "thesaurusmany");
                            adjustVcfEntryFormatGenotype(entry, numsynonyms);
                        } else {
                            adjustVcfEntryFilter(entry, numsynonyms == 0, "thesaurus");
                            adjustVcfEntryFormatGenotype(entry, numsynonyms);
                            if (numsynonyms > 0) {
                                StringBuilder sb = new StringBuilder();
                                sb.append(entry.getChr()).append(":").append(entry.getPosition());
                                sb.append(ThesaurusSynonyms.chainSynonyms(synonyms, ginfo)).append("\n");
                                outvtf.write(sb.toString().getBytes());
                                synnetwork.addLinks(new SNVPosition(entry, ginfo), synonyms);
                            }
                        }
                    }
                }
            }

            outvcf.write(entry.toString().getBytes());
        }

        thesregions.close();
        bamregions.close();
        genomereader.close();
        System.gc();

        return synnetwork;
    }

    /**
     *
     * screens candidate sites for SNVpositions that do not point to a real
     * position on the genome (i.e. to an "unaligned" position)
     *
     * @param synonyms
     *
     * @return
     *
     * true if one of the elements in the array is unaligned (negative
     * chromosome index)
     *
     */
    private boolean containsLinkToUnaligned(ArrayList<SNVPosition> synonyms) {
        int ss = synonyms.size();
        for (int i = 0; i < ss; i++) {
            if (synonyms.get(i).getChrIndex() < 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * Serves two purposes. 1 - to wrap the SAMRecord objects into
     * ThesaurusSAMRecords 2 - to clip the records if they have mismatches near
     * the ends. (But if the clipped records do not support the variant, the
     * clipping is not carried out)
     *
     *
     * @param varinbam
     * @param genomereader
     * @param cliplength
     * @return
     */
    private ThesaurusSAMRecord[] makeClippedRecords(SAMRecord[] varinbam, FastaReader genomereader,
            int cliplength, int variantposition) {

        ThesaurusSAMRecord[] tbamrecords = new ThesaurusSAMRecord[varinbam.length];
        for (int k = 0; k < varinbam.length; k++) {
            SAMRecord record = varinbam[k];

            byte[] refsequence = genomereader.getSequenceBase1(record.getAlignmentStart(), record.getAlignmentEnd());

            ThesaurusSAMRecord trecord = new ThesaurusSAMRecord(record, refsequence, cliplength);
            // if the clipping cancel the variant position, revert the clipping                        
            if (trecord.getAlignmentStart() > variantposition || trecord.getAlignmentEnd() < variantposition) {
                trecord = new ThesaurusSAMRecord(record, refsequence, 0);
            }

            tbamrecords[k] = trecord;
        }
        return tbamrecords;
    }

    private void countDetails(SNVPositionDetails details, SAMRecord[] varinbam) {

        if (varinbam == null) {
            return;
        }

        int vsize = varinbam.length;
        for (int i = 0; i < vsize; i++) {
            ThesaurusSAMRecord tr = new ThesaurusSAMRecord(varinbam[i]);
            byte nowbase = tr.getBaseAtGenomicPosition(details.getPosition());
            switch (nowbase) {
                case 'A':
                    details.incrementA();
                    break;
                case 'T':
                    details.incrementT();
                    break;
                case 'C':
                    details.incrementC();
                    break;
                case 'G':
                    details.incrementG();
                    break;
                case 'N':
                    details.incrementN();
                    break;
                default:
                    break;
            }
        }
    }

    /**
     * Function returns nothing, but change the filter field in the entry object
     *
     * @param entry
     * @param pass
     *
     * true if the entry is ok, false if the entry is annotated via the
     * thesaurus
     */
    private void adjustVcfEntryFilter(VCFEntry entry, boolean pass, String filter) {
        String oldfilter = entry.getFilter();
        if (pass) {
            if (oldfilter.equals(".")) {
                entry.setFilter("PASS");
            }
        } else {
            if (oldfilter.equals(".") || oldfilter.equals("PASS")) {
                entry.setFilter(filter);
            } else {
                entry.setFilter(oldfilter + ";" + filter);
            }
        }
    }

    /**
     * annotates a vcf entry by the number of synonyms found by the thesaurus.
     *
     * If the TS field is not present in the format, this field is appended and
     * the number of synonyms is added.
     *
     * If the TS field is present, the value in the genotype column is changed
     * to reflect the value numsynonyms
     *
     * @param entry
     *
     * vcf entry to annotate
     *
     *
     * @param numsynonyms
     *
     * number of synonyms found for this entry
     *
     */
    private void adjustVcfEntryFormatGenotype(VCFEntry entry, int numsynonyms) {

        String format = entry.getFormat();
        String genotype = entry.getGenotype();

        // check if the TS is already defined
        String[] formattokens = format.split("\t");
        int TSindex = -1;
        for (int i = 0; i < formattokens.length; i++) {
            if (formattokens[i].equals("TS")) {
                TSindex = i;
            }
        }

        if (TSindex >= 0) {
            String[] genotypetokens = genotype.split("\t");
            genotypetokens[TSindex] = "" + numsynonyms;

            StringBuilder sb = new StringBuilder();
            sb.append(genotypetokens[0]);
            for (int i = 0; i < genotypetokens.length; i++) {
                sb.append(":").append(genotypetokens[i]);
            }
            entry.setGenotype(sb.toString());
        } else {
            genotype = genotype + ":" + numsynonyms;
            format = format + ":TS";
            entry.setFormat(format);
            entry.setGenotype(genotype);
        }

    }

    /**
     * When this function finishes and the vcf file is well formed, the
     * vcfreader should start producing actual data table lines when readLines()
     * is called.
     *
     *
     * @param vcfreader
     * @param outvcf
     * @throws IOException
     */
    private void prepareVcfStreams(VCFEntrySet allvariants, OutputStream outvcf) throws IOException {

        // read and copy vcf file header from the original file
        StringBuilder oldheader = new StringBuilder();
        oldheader.append(allvariants.getHeader());

        // add the filter message line
        String formatmessage = "##FORMAT=<ID=TS,Number=1,Type=Integer,Description=\"Number of Thesaurus synonyms\">\n";
        oldheader.append(formatmessage);
        String filtermessage = "##FILTER=<ID=thesaurus,Description=\"Annotated with thesaurus " + thesaurusfile.getAbsoluteFile() + "\">\n";
        String filtermessage2 = "##FILTER=<ID=thesaurushard,Description=\"Difficult to annotated with thesaurus " + thesaurusfile.getAbsoluteFile() + "\">\n";
        String filtermessage3 = "##FILTER=<ID=thesaurusmany,Description=\"More than " + toomany + " alternative sites found with thesaurus\">\n";
        oldheader.append(filtermessage);
        oldheader.append(filtermessage2);
        oldheader.append(filtermessage3);

        // add the #CHROM line 
        oldheader.append(allvariants.getColDefLine());
        outvcf.write(oldheader.toString().getBytes());
    }

    /**
     * When this function finishes and the vcf file is well formed, the
     * vcfreader should start producing actual data table lines when readLines()
     * is called.
     *
     *
     * @param vcfreader
     * @param outvcf
     * @throws IOException
     */
    private void prepareVtfStream(OutputStream outvtf) throws IOException {

        // read and copy vcf file header from the original file
        StringBuilder header = new StringBuilder();
        header.append("##Variant Thesaurus File\n");
        header.append("##Processed by GeneticThesaurus ").append(GeneticThesaurus.getVersion()).append("\n");
        header.append("##Matching filtered variant call file ").append(output).append(".vcf\n");
        header.append("##Format: position in vcf file, followed by thesaurus synonyms, separated by tabs\n");
        outvtf.write(header.toString().getBytes());
    }

    /**
     * write a header suitable for a list of BAF
     *
     * @param allvariants
     * @param outvcf
     */
    private void prepareBAFStream(OutputStream outvcf) throws IOException {
        // read and copy vcf file header from the original file
        StringBuilder oldheader = new StringBuilder();
        oldheader.append("## BAF for variants\n");
        oldheader.append("chr\tposition\tref\talt\tnaive.BAF\tthesaurus.BAF\tthesaurus.cov\tthesaurus.synonyms\n");
        outvcf.write(oldheader.toString().getBytes());
    }
}
