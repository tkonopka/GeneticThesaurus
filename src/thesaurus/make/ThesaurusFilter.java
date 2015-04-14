/*
 * Copyright 2013-2015 Tomasz Konopka.
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
import thesaurus.util.SNVPositionComparator;
import thesaurus.util.SNVPositionDetailsList;
import thesaurus.util.ThesaurusSAMRecord;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
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

    // variables collected from the command line
    private String output = "stdout";
    private String validate = "SILENT";
    private File thesaurusfile = null;
    private File vcffile = null;
    private File bamfile = null;
    private int minmapqual = -1;
    private boolean verbose = true;
    private final ThesaurusLog mylog;
    // for comparisons
    private ArrayList<File> bamfiles = new ArrayList<>(2);
    private ArrayList<String> labels = new ArrayList<>(2);
    // thresholds used during filtering
    private int tolerance = 1;
    private int maxtolerance = 2;
    private int softclip = 5;
    private int hitcount = 3;
    private double hitproportion = 0.5;
    private int many = 100;
    private int toomany = 400;
    private int insertsize = 200;
    // for outputing 
    private final DecimalFormat BAFformat = new DecimalFormat("0.0000");
    //private final int specialpos = 18637725;
    //private final int specialpos2 = 21357558;    

    private void printFilterHelp() {
        System.out.println("GeneticThesaurus filter: use a thesaurus to filter a VCF file");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar filter ");
        System.out.println();
        System.out.println("Core options:");
        ThesaurusIO.printHelpItem("--bam <File>", "alignment files matching --vcf");
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--thesaurus <File>", "thesaurus file");
        ThesaurusIO.printHelpItem("--vcf <File>", "variant call file matching --bam");
        System.out.println("\nComparisons:");
        ThesaurusIO.printHelpItem("--comparebam <File>", "additional alignemnt files for comparisons of BAFs; can specify multiple times");
        ThesaurusIO.printHelpItem("--comparelabel <String>", "labels matching --comparebam; specify one for each --comparebam in order");
        ThesaurusIO.printHelpItem("--label <String>", "label for alignment --bam");
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
        ThesaurusIO.printHelpItem("--validate <String>", "validation level for SAM records [default " + validate + "]");
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
        prs.accepts("minmapqual").withRequiredArg().ofType(Integer.class);
        prs.accepts("vcf").withRequiredArg().ofType(File.class);
        prs.accepts("bam").withRequiredArg().ofType(File.class);

        // for comparisons
        prs.accepts("comparebam").withRequiredArg().ofType(File.class);
        prs.accepts("comparelabel").withRequiredArg().ofType(String.class);
        prs.accepts("label").withRequiredArg().ofType(String.class);

        // thresholds for filtering
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
        boolean ok = true;

        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        try {
            if (options.has("genome")) {
                genome = (File) options.valueOf("genome");
                if (!genome.exists() || !genome.canRead()) {
                    System.out.println("Cannot read genome file, or file does not exist");
                    ok = false;
                }
            }

            if (options.has("minmapqual")) {
                minmapqual = (Integer) options.valueOf("minmapqual");
                if (minmapqual < 0) {
                    System.out.println("Parameter minmapqual must be >= 0");
                    ok = false;
                }
            }

            if (options.has("softclip")) {
                softclip = (Integer) options.valueOf("clip");
                if (softclip < 0) {
                    System.out.println("Parameter clip must be >= 0");
                    ok = false;
                }
            }
            if (options.has("readlen")) {
                readlen = (Integer) options.valueOf("readlen");
                if (readlen < 0) {
                    System.out.println("Parameter readlen must be >= 0");
                    ok = false;
                }
            }
            if (options.has("insertsize")) {
                insertsize = (Integer) options.valueOf("insertsize");
                if (insertsize < readlen) {
                    System.out.println("Parameter insertsize must be >= readlen");
                    ok = false;
                }
            }

            if (options.has("tolerance")) {
                tolerance = (Integer) options.valueOf("tolerance");
                if (tolerance < 0) {
                    System.out.println("Parameter tolerance must be >= 0");
                    ok = false;
                }
            }
            if (options.has("maxtolerance")) {
                maxtolerance = (Integer) options.valueOf("maxtolerance");
                if (maxtolerance < tolerance) {
                    System.out.println("Parameter maxtolerance must be >= value of parameter tolerance");
                    ok = false;
                }
            }
            if (options.has("many")) {
                many = (Integer) options.valueOf("many");
                if (many < 2) {
                    System.out.println("Parameter many must be > 1");
                    ok = false;
                }
            }
            if (options.has("toomany")) {
                toomany = (Integer) options.valueOf("toomany");
                if (toomany < 2) {
                    System.out.println("Parameter toomany must be > 1");
                    ok = false;
                }
            }

            // thresholds deciding which thesaurus links to accept
            if (options.has("hitcount")) {
                hitcount = (Integer) options.valueOf("hitcount");
                if (hitcount < 1) {
                    System.out.println("Parameter hitcount must be >= 1");
                    ok = false;
                }
            }
            if (options.has("hitproportion")) {
                hitproportion = (Double) options.valueOf("hitproportion");
                if (hitproportion < 0) {
                    System.out.println("Parameter hitproportion must be >= 0");
                    ok = false;
                }
            }

            // get prefix for output tracks
            if (options.has("output")) {
                output = (String) options.valueOf("output");
            } else {
                System.out.println("Missing required argument --output");
                ok = false;
            }

            // get prefix for output tracks
            if (options.has("validate")) {
                validate = (String) options.valueOf("validate");
                if (!validate.equals("STRICT") && !validate.equals("LENIENT") && !validate.equals("SILENT")) {
                    System.out.println("Option --validate must be set to STRICT, LENIENT, or SILENT");
                    ok = false;
                }
            }

            // get input vcf file
            if (options.has("vcf")) {
                vcffile = (File) options.valueOf("vcf");
                if (!vcffile.exists() || !vcffile.canRead()) {
                    System.out.println("Cannot read vcf file, or file does not exist: "+vcffile.getAbsolutePath());
                    ok = false;
                }
            } else {
                System.out.println("Missing required argument --vcf");
                ok = false;
            }

            if (options.has("bam")) {
                bamfile = (File) options.valueOf("bam");
                if (!bamfile.exists() || !bamfile.canRead()) {
                    System.out.println("Cannot read bam file, or file does not exist");
                    ok = false;
                }
                bamfiles.add(bamfile);
            }

            // get input vcf file
            if (options.has("thesaurus")) {
                thesaurusfile = (File) options.valueOf("thesaurus");
                if (!thesaurusfile.exists() || !thesaurusfile.canRead()) {
                    System.out.println("Cannot read thesaurus file, or file does not exist");
                    ok = false;
                }
            } else {
                System.out.println("Missing required argument --thesaurus");
                ok = false;
            }

            // Settings for sample comparisons
            // Get label for the primary sample               
            if (options.has("label")) {
                labels.add((String) options.valueOf("label"));
            } else {
                labels.add("sample");
            }

            // get labels for comparison samplse
            if (options.has("comparelabel")) {
                Object[] temp = options.valuesOf("comparelabel").toArray();
                for (int i = 0; i < temp.length; i++) {
                    labels.add((String) temp[i]);
                }
            }

            if (options.has("comparebam")) {
                Object[] temp = options.valuesOf("comparebam").toArray();
                // copy each file from temp into an array, print warnings 
                for (int i = 0; i < temp.length; i++) {
                    File nowbam = (File) temp[i];
                    bamfiles.add(nowbam);
                    if (!nowbam.exists() || !nowbam.canRead()) {
                        System.out.println("Cannot read bam file, or file does not exist\n" + nowbam.getAbsolutePath());
                        ok = false;
                    }
                }
            }

            // final checks on comparisons - number of comparison bams and labels must match
            if (labels.size() != bamfiles.size()) {
                System.out.println("Number of values for comparebam and comparelabel must match");
                System.out.println("Got " + labels.size() + " labels and " + bamfiles.size() + " bams");
                ok = false;
            }

        } catch (Exception ex) {
            System.out.println("Error parsing command line arguments\n"
                    + "This error often appears when a parameter is unexpected set more than once");
            ok = false;
        }

        return ok;
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

        // more setup - need a comparator object for chromosomes and positions
        GenomeInfo ginfo = null;        
        try {
            ginfo = new GenomeInfo(genome);
        } catch (Exception ex) {
            System.out.println("Failed to load genomic position comparator: " + ex.getMessage());
            System.out.println("Is the genome set on the command line or in the defaults?");
            return;
        }

        // read all the variants into memory                
        mylog.log("Reading variants");
        VCFEntrySet allvariants = new VCFEntrySet(vcffile, ginfo, true);
        allvariants.separateMultiSNVs();
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
            getBAFs(allvariants, allnetwork, ginfo, bamfiles, labels, outBAF);
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
            ArrayList<File> bamfiles, ArrayList<String> labels, OutputStream outBAF) throws IOException {

        // find out how many bam files to compare to
        int numbams = bamfiles.size();

        // make a list of all variats in the synonyms map
        ArrayList<SNVPosition> allpositions = synonyms.getAllNodes();

        // collect coverage and BAF details from each bam file                
        SNVPositionDetailsList[] comparedetails = new SNVPositionDetailsList[numbams];
        for (int i = 0; i < numbams; i++) {
            mylog.log("Collecting coverage information from bam file: " + labels.get(i));
            comparedetails[i] = collectBAFDetailsFromBam(allpositions, bamfiles.get(i), ginfo);
        }

        SNVPositionComparator vcomp = new SNVPositionComparator();

        mylog.log("Writing baf table");
        prepareBAFStream(outBAF);
        int vs = variants.size();
        for (int i = 0; i < vs; i++) {
            // get the variant  
            VCFEntry entry = variants.getVariant(i);            

            // collect some basic information about what the variant is and thesaurus filters
            String ref = entry.getRef();
            String alt = entry.getAlt();
            SNVPosition entrypos = new SNVPosition(entry, ginfo);
            ArrayList<SNVPosition> nowsynonyms = synonyms.getNeighborNodes(entrypos);
           
            // start a BAF file entry with information about the varint
            StringBuilder sb = new StringBuilder();
            sb.append(entry.getChr()).append("\t").append(entry.getPosition()).append("\t").append(ref).append("\t").append(alt);
            // output the number of synonyms
            if (nowsynonyms == null) {
                sb.append("\t0");
            } else {
                sb.append("\t").append(nowsynonyms.size());
            }

            // add in information about this variant from each bam file
            for (int j = 0; j < numbams; j++) {
                sb.append(getBAFentryBlock(entrypos, comparedetails[j], nowsynonyms, vcomp, ginfo));
            }
            sb.append("\n");
            
            outBAF.write(sb.toString().getBytes());
        }
    }

    /**
     *
     * @param snvposition
     *
     * a genomic position with a called variant
     *
     * @param posdetailslist
     *
     * a sorted list of positions with details
     *
     * @param nowsynonyms
     *
     * synonyms for the snvposition site
     *
     * @param vcomp
     *
     * a comparator (used to search for positions in the posdetailslist object)
     *
     * @return
     *
     * A string that can be output in a BAF file (block of values separated by
     * tabs). Format should match definition in the BAF file header. Always has
     * a leading tab.
     */
    private String getBAFentryBlock(SNVPosition snvposition,
            SNVPositionDetailsList posdetailslist, ArrayList<SNVPosition> nowsynonyms,
            SNVPositionComparator vcomp, GenomeInfo ginfo) {

        // get "naive" estimates using only the posdetails object
        SNVPositionDetails posdetails = posdetailslist.find(snvposition, vcomp);
        long naivecov = posdetails.getCountATCG();
        double naiveBAF = ((double) posdetails.getCountAlt()) / ((double) naivecov);
        if (naivecov == 0) {
            naiveBAF = 0.0;
        }
        
        // get "thesaurus" estimates 
        SNVPositionDetails thesdetails = new SNVPositionDetails(posdetails);        
        int numsynonyms = 0;
        if (nowsynonyms != null) {
            numsynonyms = nowsynonyms.size();
            for (int i = 0; i < nowsynonyms.size(); i++) {                
                thesdetails.incrementCounts(posdetailslist.find(nowsynonyms.get(i), vcomp));                
            }
        }

        long thescov = thesdetails.getCountATCG();
        double thesBAF = ((numsynonyms + 1) * (double) thesdetails.getCountAlt()) / ((double) thescov);
        if (thescov == 0) {
            thesBAF = 0.0;
        }

        // create the block with naive and thesaurus estimates, separated by tabs
        StringBuilder sb = new StringBuilder();
        sb.append("\t").append(BAFformat.format(naiveBAF)).append("\t").append(BAFformat.format(thesBAF));
        sb.append("\t").append(naivecov).append("\t").append(thescov);
        return sb.toString();
    }

    /**
     * Collects information about SNV positions in a given alignment file.
     *
     *
     * @param snvpositions
     *
     * a set of genomic position with base changes
     *
     * @param bamf
     *
     * file to scan to fill in coverage information
     *
     * @param ginfo
     * @return
     *
     * a list of same length as snvposition in the input. This new list includes
     * counters for each base in the alphabet. The function scans the alignment
     * to fill in the counters.
     *
     * @throws IOException
     */
    private SNVPositionDetailsList collectBAFDetailsFromBam(ArrayList<SNVPosition> snvpositions,
            File bamf, GenomeInfo ginfo) throws IOException {

        // from a list of positions, create a list that includes base/coverage counters
        SNVPositionDetailsList alldetails = makeDetailsList(snvpositions);

        // in the constructor for bamregions, I had the minimum mapping quality set at "minmapqual"
        // but this only worked well for minmapqual<=1. So now I set this to zero regardless of 
        // the mapping quality used in variant calling.        
        BamRegionsMap bamregions = new BamRegionsMap(bamf, ginfo, "SILENT", 0);

        int numloci = alldetails.size();
        for (int i = 0; i < numloci; i++) {
            SNVPositionDetails entry = alldetails.get(i);
            int nowchr = entry.getChrIndex();
            // look up the position in the bam
            SAMRecord[] varinbam = bamregions.lookup(ginfo.getChrName(nowchr), entry.getPosition());
            // collect information about the locus from each bam record
            countDetailsFromRecords(entry, varinbam);
        }
        bamregions.close();

        return (alldetails);
    }

    /**
     * Convert between a list of positions to an array holding also counters for
     * each base.
     *
     * @param variants
     *
     *
     * @param allloci
     *
     * all SNVPositions in sorted order
     *
     * @param ginfo
     *
     * @return
     */
    private SNVPositionDetailsList makeDetailsList(ArrayList<SNVPosition> allloci) {
        // make array with details
        int allsize = allloci.size();
        SNVPositionDetailsList alldetails = new SNVPositionDetailsList(allsize);
        for (int i = 0; i < allsize; i++) {
            alldetails.add(new SNVPositionDetails(allloci.get(i)));
        }                
        alldetails.sortList();        
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
            //System.out.print(entry.toString(ginfo));

            if (entry.isIndel()) {
                // skip processing indels, but still change the filter and format fields
                adjustVcfEntryFilter(entry, true, "thesaurus");
                adjustVcfEntryFormatGenotype(entry, 0);
            } else {
                // variant is a substitution - this is the interesting case here

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
                    if (readsOnLocus != null) {
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

    /**
     * Increment count values stored in details using the sequences in SAM
     * records. Function does not return anything, but the counters in the
     * details object are edited/updated.
     *
     * @param details
     * @param varinbam
     */
    private void countDetailsFromRecords(SNVPositionDetails details, SAMRecord[] varinbam) {

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
        // prepare a header for a file with BAF estimates
        StringBuilder bafheader = new StringBuilder();
        bafheader.append("## BAF for variants\n");
        bafheader.append("chr\tposition\tref\talt\tthesaurus.synonyms\t");

        for (int i = 0; i < labels.size(); i++) {
            String nowlabel = labels.get(i);
            bafheader.append("\t" + nowlabel + ".naive.BAF"
                    + "\t" + nowlabel + ".thesaurus.BAF"
                    + "\t" + nowlabel + ".naive.cov"
                    + "\t" + nowlabel + ".thesarus.cov");
        }

        bafheader.append("\n");
        outvcf.write(bafheader.toString().getBytes());
    }
}
