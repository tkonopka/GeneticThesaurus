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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.sequence.SequenceComplementer;
import jsequtils.sequence.SequenceMap;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * Write a thesaurus of genomic regions. Uses as input alignment file(s) with
 * many non-primary records. Outputs a table showing regions of origin and
 * closely related regions.
 *
 * Works in two passes.
 *
 * Pass one output several temporary files, one for each alignment chromosome.
 * This section can work in parallel, each thread reading from a distinct input
 * file.
 *
 * Pass two scans the chromosome files and compiles the merged thesaurus table.
 * This part in particular requires a lot of memory to hold/sort/merge all the
 * small segments.
 *
 *
 * @author tkonopka
 */
public class ThesaurusWrite extends ThesaurusMapTool {

    //private File bamfile = null;
    private ArrayList<File> bamfiles = new ArrayList<File>();
    // prefix for output files
    private String output = "stdout";
    // some output files can be stored on a temporary disk
    private String tempoutput = null;
    private final int minmapqual;
    // maxpenalty and errorrate will determine how many mismatches are reported in the 
    // final thesaurus table
    private int maxpenalty = 2;
    private final double errorrate;
    // settings for multi-threading, this is essentially for running stream compression in parallel
    private int numthreads = 4;
    // objects holding summaries of the genome and enabling sorting/comparison of thesaurus entries
    private final GenomeInfo ginfo;
    private final ThesaurusEntryAlignComparator teac;
    private final ThesaurusEntryMergingComparator temc;
    //advanced features - allows to skip pass 1 or pass 2 (useful when one fails)
    private boolean skip1 = false, skip2 = false;
    ThesaurusLog theslog = new ThesaurusLog();
    // bucketspacing is used when sorting chromosome files
    // records are distributed into buckets of these many bases
    private final int bucketspacing = 1000000;

    private void printWriteHelp() {
        System.out.println("GeneticThesaurus write: write a thesaurus of genetic variation");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar write ");
        System.out.println();
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta");
        ThesaurusIO.printHelpItem("--bam <File>", "alignment (can specify more than once)");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--tempdisk <String>", "prefix for temporary files (use a fast disk)");
        ThesaurusIO.printHelpItem("--maxpenalty <int>", "maximum penalty value in output thesaurus table [default 2]");
        ThesaurusIO.printHelpItem("--readlen <int>", "readlength of aligned reads");
        ThesaurusIO.printHelpItem("--skip1", "skip pass one, go directly to pass 2");
        ThesaurusIO.printHelpItem("--skip2", "skip pass two");
        ThesaurusIO.printHelpItem("--threads <int>", "number of threads used in pass 1");
        System.out.println();
    }

    /**
     * The constructor for this tool parses command line parameters and
     * initializes some variables. The computing only begins when the Runnable
     * is executed.
     *
     * @param args
     */
    public ThesaurusWrite(String[] args) {
        if (args.length == 0) {
            printWriteHelp();
            ginfo = null;
            teac = null;
            temc = null;
            errorrate = 0.0;
            minmapqual = 0;
            return;
        }

        super.loadDefaults();
        // parse parameters
        super.setOk(parseWriteParameters(args));
        if (!super.isOk()) {
            ginfo = null;
            teac = null;
            temc = null;
            errorrate = 0.0;
            minmapqual = 0;
            return;
        }

        GenomeInfo tempgenomeinfo;
        try {
            tempgenomeinfo = new GenomeInfo(genome);
        } catch (Exception ex) {
            System.out.println("Error reading genome information: " + ex.getMessage());
            super.setOk(false);
            teac = null;
            temc = null;
            tempgenomeinfo = null;
            ginfo = null;
            errorrate = 0.0;
            minmapqual = 0;
            return;
        }

        ginfo = tempgenomeinfo;
        teac = new ThesaurusEntryAlignComparator();
        temc = new ThesaurusEntryMergingComparator();

        // minmapqual is related to the read length and maximal number of mismatches expected
        minmapqual = readlen - maxpenalty;
        errorrate = (double) maxpenalty / readlen;
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
        prs.accepts("bam").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("tempdisk").withRequiredArg().ofType(String.class);

        // tuning for how the thesaurus should be written
        prs.accepts("maxpenalty").withRequiredArg().ofType(Integer.class);
        prs.accepts("readlen").withRequiredArg().ofType(Integer.class);
        prs.accepts("threads").withRequiredArg().ofType(Integer.class);

        prs.accepts("skip1");
        prs.accepts("skip2");

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

        if (options.has("bam")) {
            bamfiles.addAll((List<File>) options.valuesOf("bam"));
            // check that files exists
            boolean failbamtest = false;
            for (int i = 0; i < bamfiles.size(); i++) {
                File nowfile = bamfiles.get(i);
                if (!nowfile.exists() | !nowfile.canRead()) {
                    //System.out.println("Cannot read input bam file: " + nowfile.getAbsolutePath());
                    //failbamtest = true;
                }
            }
            if (failbamtest) {
                return false;
            }
        }

        if (options.has("output")) {
            output = (String) options.valueOf("output");
        } else {
            System.out.println("Missing required argument output");
            return false;
        }

        if (options.has("tempdisk")) {
            tempoutput = (String) options.valueOf("tempdisk");
        } else {
            tempoutput = output;
        }

        if (options.has("maxpenalty")) {
            maxpenalty = (Integer) options.valueOf("maxpenalty");
            if (maxpenalty < 0) {
                System.out.println("maxpenalty must be >= 0");
                return false;
            }
        }

        if (options.has("readlen")) {
            readlen = (Integer) options.valueOf("readlen");
            if (readlen < 0) {
                System.out.println("readlen must be >= 0");
                return false;
            }
        }

        if (options.has("threads")) {
            numthreads = (Integer) options.valueOf("threads");
            if (numthreads < 1) {
                System.out.println("threads must be >= 1");
                return false;
            }
        }

        skip1 = options.has("skip1");
        skip2 = options.has("skip2");

        return true;
    }

    @Override
    public void runTool() {

        // load all the genome
        SequenceMap genomeSeqMap;
        theslog.log(true, "Loading reference genome map");
        try {
            genomeSeqMap = new SequenceMap(genome, true);
        } catch (IOException ex) {
            System.out.println("Error reading genome: " + ex.getMessage());
            return;
        }

        // Pass one of thesaurus. Write a large file with all information
        try {
            if (!skip1) {
                // make temporary files, one per chromosome
                HashMap<String, OutputStream> tempfiles = makePassOneTempFiles(output);
                // read all the bamfiles one by one 
                theslog.log(true, "Pass one: " + bamfiles.size() + " alignment files");
                writeThesaurusPassOne(genomeSeqMap, tempfiles);
                // close all the temporary files
                for (Map.Entry<String, OutputStream> ee : tempfiles.entrySet()) {
                    ee.getValue().close();
                }
            } else {
                theslog.log(true, "Skipping pass one");
            }
        } catch (IOException ex) {
            System.out.println("Something went wrong in pass one: " + ex.getMessage());
        }

        // Pass two of the thesaurus. Read thesaurus tables and try to compress them
        // This involves a sort step for each chromosome, followed by merging steps
        try {
            if (!skip2) {
                theslog.log(true, "Pass two");
                writeThesaurusPassTwo(output, genomeSeqMap);
            } else {
                theslog.log(true, "Skipping pass two");
            }
        } catch (IOException ex) {
            System.out.println("Something went wrong in pass two: " + ex.getMessage());
        }

    }

    /**
     * PassTwo is a sorting pass. Sort each chromosome by first copying entries
     * into buckets and then processing each bucket separately.
     *
     * @param output
     * @throws IOException
     */
    private void writeThesaurusPassTwo(String output, SequenceMap seqmap) throws IOException {

        // create output stream for the main output file
        OutputStream outstream;
        try {
            outstream = OutputStreamMaker.makeOutputStream(output + ".tsv.gz");
        } catch (Exception ex) {
            System.out.println("Error creating output stream in pass two: " + ex.getMessage());
            return;
        }
        // write the header again to the outstream
        StringBuilder header = new StringBuilder();
        header.append("##maxpenalty=").append(maxpenalty).append("\n");
        header.append("##readlen=").append(readlen).append("\n");
        header.append("##genome=").append(genome.getAbsolutePath());
        writeThesaurusHeader(outstream, minmapqual, header.toString());

        // process each chromosome, one at a time in series
        int numchr = ginfo.getNumChromosomes();
        for (int i = 0; i < numchr; i++) {
            String nowchr = ginfo.getChrName(i);
            theslog.log(true, nowchr);
            sortAndProcessChromosome(output, nowchr, outstream, seqmap);
        }

        outstream.close();
    }

    /**
     * sort a single chromosome.
     *
     * Note: the "raw" file will be deleted at the end.
     *
     *
     * @param output
     *
     * @param chr
     *
     * @throws IOException
     *
     */
    private void sortAndProcessChromosome(String output, String chr, OutputStream outstream, SequenceMap seqmap) throws IOException {


        // find out the number of buckets for this chromosome
        int chrlen = ginfo.getChrLength(chr);
        int numbuckets = 1 + (chrlen / bucketspacing);

        // keep track of how many items were put into each bucket
        int[] bucketsize = new int[numbuckets];

        // --- Part 1, read from one file and distribute into buckets        

        // for debugging

        // make array of output streams
        OutputStream[] outbuckets = new OutputStream[numbuckets];
        for (int i = 0; i < numbuckets; i++) {
            // compress some and keep others not compressed?            
            outbuckets[i] = OutputStreamMaker.makeOutputStream(tempoutput + "." + chr + "." + i + ".tsv");
            writeThesaurusHeader(outbuckets[i], minmapqual, "##GeneticThesaurus " + chr + " bucket " + i);
        }

        // open the chromosome file and distribute into buckets, perform the writing 
        // using multiple threads to speed up output
        ExecutorService service = Executors.newFixedThreadPool(numthreads);
        BlockingQueue<StringBucket> queue = new LinkedBlockingQueue<StringBucket>(32 * numthreads);
        for (int i = 0; i < numthreads; i++) {
            service.submit(new DistributeRecordsPassTwoRunnable(outbuckets, queue, bucketspacing));
        }

        File rawfile = new File(output + "." + chr + ".raw.tsv.gz");
        BufferedReader chrreader = BufferedReaderMaker.makeBufferedReader(rawfile);
        String ss = null;
        try {
            while ((ss = chrreader.readLine()) != null) {
                if (!ss.startsWith("#") && !ss.startsWith("Align.chr")) {
                    StringBucket entry = new StringBucket(ss);
                    queue.put(entry);
                    int nowbucket = entry.val / bucketspacing;
                    bucketsize[nowbucket]++;
                }
            }
            // add poison pills into the queue (ss is a null string at this point)            
            for (int i = 0; i < numthreads + 1; i++) {
                queue.put(new StringBucket(ss));
            }
        } catch (Exception ex) {
            System.out.println("Error reading/distributiong entries: " + ex.getMessage());
        }

        try {
            service.shutdown();
            service.awaitTermination(1000, java.util.concurrent.TimeUnit.DAYS);
            queue = null;
        } catch (Exception ex) {
            System.out.println("Error shutting down service: " + ex.getMessage());
        }

        // close the chromosome and buckets streams from part 1
        chrreader.close();
        for (int i = 0; i < numbuckets; i++) {
            outbuckets[i].close();
        }
        System.gc();

        // --- Part 2, read from the buckets, sort, merge, and output
        for (int i = 0; i < numbuckets; i++) {
            File bucketfile;
            bucketfile = new File(tempoutput + "." + chr + "." + i + ".tsv");
            int minalignstart = 1 + (i * bucketspacing);
            summarizeFromBucket(bucketfile, bucketsize[i], minalignstart, outstream, seqmap);
            bucketfile.delete();
            System.gc();
        }
    }

    /**
     * helper function that processes entries from a bucket of known size
     *
     * @param bucketfile
     *
     * File holding a set of thesaurus entries
     *
     * @param bucketsize
     *
     * number of the thesaurus entries to read from the file (allows to use
     * fixed size object array rather than collection)
     *
     * @param minalignstart
     *
     * when making output, do not allow to extend thesaurus entries to aligns
     * start position lower than this value.
     *
     * @param outstream
     *
     * output
     *
     * @param seqmap
     *
     * helper object containing genome sequence. This is used while extending
     * thesaurus entries while keeping track of mismatches.
     *
     *
     * @throws IOException
     */
    private void summarizeFromBucket(File bucketfile, int bucketsize, int minalignstart,
            OutputStream outstream, SequenceMap seqmap) throws IOException {
        // set up array to hold entries
        ThesaurusEntry[] bucketentries = new ThesaurusEntry[bucketsize];

        // read from the bucket into array, then close and delete the bucket        
        BufferedReader bucketreader = BufferedReaderMaker.makeBufferedReader(bucketfile);
        int counter = 0;
        String ss;
        while ((ss = bucketreader.readLine()) != null) {
            if (!ss.startsWith("#") && !ss.startsWith("Align.chr")) {
                bucketentries[counter] = new ThesaurusEntry(ss, ginfo);
                counter++;
            }
        }
        bucketreader.close();

        // simplify/merge/collapse the entries within the bucket            
        ArrayList<ThesaurusEntry> mergedentries = mergeEntries(bucketentries, seqmap, minalignstart);

        // outut the entries to the final out stream
        int nummerged = mergedentries.size();
        for (int j = 0; j < nummerged; j++) {
            outstream.write(mergedentries.get(j).toString().getBytes());
        }
    }

    boolean containsMine(ThesaurusEntry te, String chr, int pos) {
        if (!te.getAlignChr().equals(chr)) {
            return false;
        }
        //return (te.alignStart <= pos && te.alignEnd >= pos);
        return false;
    }

    /**
     * simplifies a collection of thesaurus entries. This is done by first
     * special sorting the entries, then combining neighboring ones together.
     *
     * @param tempentries
     *
     * @param seqmap
     *
     * @param minalignstart
     *
     * minimum number of the align start position. This is used during the
     * extension step, which prevents extending an entry to beyond the bucket
     * boundary.
     *
     * @return
     *
     * A list of merged and expanded entries, sorted by alignment, ready to be
     * output.
     *
     */
    private ArrayList<ThesaurusEntry> mergeEntries(ThesaurusEntry[] tempentries, SequenceMap seqmap, int minalignstart) {

        // make an array to hold the merged entries
        int numentries = tempentries.length;
        ArrayList<ThesaurusEntry> ans = new ArrayList<ThesaurusEntry>(numentries);

        if (tempentries.length < 1) {
            return ans;
        }

        // sort the elements of this bucket so that they can be easily merged                        
        Arrays.sort(tempentries, temc);

        ThesaurusEntry lastentry = tempentries[0];
        lastentry.extendLeft(seqmap, errorrate, 2, minalignstart);
        lastentry.extendRight(seqmap, errorrate, 2);
        ans.add(lastentry);

        for (int i = 1; i < numentries; i++) {
            ThesaurusEntry entry = tempentries[i];
            if (!lastentry.mergeWith(entry)) {
                // try to make last entry as long as possible                                                
                lastentry.extendLeft(seqmap, errorrate, 2, minalignstart);
                lastentry.extendRight(seqmap, errorrate, 2);

                // try to make current entry as long as possible, and symmetrical
                ThesaurusEntry entrycopyR = new ThesaurusEntry(entry);
                entrycopyR.extendRight(seqmap, errorrate, 2);
                entry.extendLeft(seqmap, errorrate, 2, minalignstart);
                entry.mergeWith(entrycopyR);

                // last attempt to merge with previous entry
                if (!lastentry.mergeWith(entry)) {
                    ans.add(entry);
                    lastentry = entry;
                }
            }
        }

        //System.out.print(" -> " + ans.size());
        // make sure the output is properly sorted by alignment position        
        int anssize = ans.size();
        tempentries = new ThesaurusEntry[anssize];
        for (int i = 0; i < anssize; i++) {
            tempentries[i] = ans.get(i);
        }
        Arrays.sort(tempentries, teac);

        // do another round of merging, this time on teac sorted
        ans = new ArrayList<ThesaurusEntry>(anssize);
        lastentry = tempentries[0];
        ans.add(lastentry);
        for (int i = 1; i < anssize; i++) {
            ThesaurusEntry entry = tempentries[i];
            if (!lastentry.mergeWith(entry)) {
                // add this entry to the list                
                ans.add(entry);
                lastentry = entry;
            }
        }

        // do another round just to make sure
        Collections.sort(ans, temc);
        ArrayList<ThesaurusEntry> ans2 = new ArrayList<ThesaurusEntry>(ans.size());
        lastentry = ans.get(0);
        ans2.add(lastentry);
        anssize = ans.size();
        for (int i = 1; i < anssize; i++) {
            ThesaurusEntry entry = ans.get(i);
            if (!lastentry.mergeWith(entry)) {
                // add this entry to the list
                ans2.add(entry);
                lastentry = entry;
            }
        }

        // sort the collection now, 
        // It should already be almost sorted, but the merging can disturb
        // the order a little bit, so resort to be sure        
        Collections.sort(ans2, teac);

        return ans2;
    }

    /**
     * makes a set of files/buckets for alignments, one per chromosome defined
     * in the genome info
     *
     * @param outbase
     * @param seqmap
     * @return
     */
    private HashMap<String, OutputStream> makePassOneTempFiles(String outbase) throws FileNotFoundException, IOException {
        HashMap<String, OutputStream> ans = new HashMap<String, OutputStream>();
        int numchr = ginfo.getNumChromosomes();
        for (int i = 0; i < numchr; i++) {
            String nowchr = ginfo.getChrName(i);
            OutputStream outstream = OutputStreamMaker.makeOutputStream(outbase + "." + nowchr + ".raw.tsv.gz");
            writeThesaurusHeader(outstream, minmapqual, "##GeneticThesaurus " + nowchr);
            ans.put(nowchr, outstream);
        }
        return ans;
    }

    /**
     * One of the core functions of this class. It compares reads in the
     * alignment with the genome sequence and outputs lines of the thesaurus.
     * Output is in raw format, i.e. one line per read.
     *
     * @param genomereader
     * @param inputSam
     * @param output
     * @throws IOException
     */
    private void writeThesaurusPassOne(SequenceMap seqmap, HashMap<String, OutputStream> outstreams) throws IOException {

        // process all the bam files in separate runnables/threads
        ExecutorService service = Executors.newFixedThreadPool(numthreads);
        for (File bamfile : bamfiles) {
            service.submit(new ProcessRecordsPassOneRunnable(bamfile, seqmap, outstreams, minmapqual, theslog, ginfo));
        }
        // wait for all the computation to finish and exit
        try {
            service.shutdown();
            service.awaitTermination(1000, TimeUnit.DAYS);
        } catch (Exception ex) {
            System.out.println("Exception shutting down service: " + ex.getMessage());
        }

    }

    static void writeThesaurusHeader(OutputStream outstream, int minmapqual, String custom) throws IOException {
        StringBuilder header = new StringBuilder();
        header.append("##GeneticThesaurus\n");
        header.append("##minmapqual=").append(minmapqual).append("\n");
        header.append(custom).append("\n");
        header.append(ThesaurusEntry.getHeader());
        outstream.write(header.toString().getBytes());
    }
}

/**
 * Class holding a string and an integer. The integer is evaluated from the
 * second column of the string (tab separated) This allows holding the alignment
 * position of a thesaurus entry without full parsing of the entry.
 *
 * @author tkonopka
 */
class StringBucket {

    final String text;
    final int val;
    static int tabint = (int) (byte) '\t';

    public StringBucket(String s) {
        if (s == null) {
            this.text = null;
            this.val = 0;
        } else {
            this.text = s;
            // parse out the align start position, 
            // knowing that it is in the second tab-separated column
            int tab1 = 1 + s.indexOf(tabint, 0);
            int tab2 = s.indexOf(tabint, tab1);
            this.val = Integer.parseInt(s.substring(tab1, tab2));
        }
    }
}

/**
 * Runnable class that reads an alignment file and outputs thesaurus records to
 * output stream by chromosome.
 *
 *
 */
class ProcessRecordsPassOneRunnable implements Runnable {

    private final File bamfile;
    private final SequenceMap seqmap;
    private final HashMap<String, OutputStream> outstreams;
    private final int minmapqual;
    private final ThesaurusLog theslog;
    private final GenomeInfo ginfo;

    public ProcessRecordsPassOneRunnable(File bamfile, SequenceMap seqmap,
            HashMap<String, OutputStream> outstreams, int minmapqual, ThesaurusLog theslog, GenomeInfo ginfo) {
        this.bamfile = bamfile;
        this.seqmap = seqmap;
        this.outstreams = outstreams;
        this.minmapqual = minmapqual;
        this.theslog = theslog;
        this.ginfo = ginfo;
    }

    @Override
    public void run() {

        theslog.log(true, "Pass one: " + bamfile.getAbsolutePath());

        // set up input/output objects using SAM library
        SAMFileReader inputSam = new SAMFileReader(bamfile);
        inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        for (final SAMRecord record : inputSam) {
            int mapqual = record.getMappingQuality();
            // only consider reads with a minimum mapping quality
            // (this here is not a "real" mapping quality, but a proxy for number of mismatches)
            if (mapqual >= minmapqual && !record.getReadUnmappedFlag()) {
                // check that read is aligned                
                processOneRecord(record);
            }
        }

        // and close the input file stream
        inputSam.close();
    }

    /**
     * converts a SAM record into a thesaurus entry
     *
     * @param record
     * @param seqmap
     * @param outsream
     */
    private void processOneRecord(SAMRecord record) {
        // create a new entry from the record, but give up if correct alignment
        ThesaurusEntry entry = new ThesaurusEntry(record, ginfo);
        if (entry.isTrivial()) {
            return;
        }

        if (isAlignmentClean(record)) {
            completeThesaurusEntry(entry, record);
            outputThesaurusEntry(entry, outstreams.get(entry.getAlignChr()));
        }

    }

    /**
     *
     * @return
     *
     * true if the cigar contains only one element that describes a sequence
     * match, e.g. 100M
     *
     *
     */
    private boolean isAlignmentClean(SAMRecord record) {

        Cigar cigar = record.getCigar();
        int numcigel = cigar.numCigarElements();
        if (numcigel > 1) {
            return false;
        }
        CigarOperator co = cigar.getCigarElement(0).getOperator();
        if (co != CigarOperator.M) {
            return false;
        }

        return true;
    }

    /**
     * updates a ThesaurusEntry with information about anchor points.
     *
     * This is done by comparing sequences stored in the SAM record and the
     * genome fasta file.
     *
     * @param entry
     * @param record
     * @param genomereader
     */
    private void completeThesaurusEntry(ThesaurusEntry entry, SAMRecord record) {

        // get sequences of the read and the genome alignment
        byte[] readseq = record.getReadBases();
        byte[] genomeseq = seqmap.getSequenceBase1(entry.getAlignChr(), entry.alignStart, entry.alignEnd);

        // if clipping has occured, need to redefine readseq
        int numclip = readseq.length - genomeseq.length;
        if (numclip != 0) {
            Cigar cig = record.getCigar();
            if (cig.getCigarElement(0).getOperator() == CigarOperator.S) {
                // clip from the start
                readseq = Arrays.copyOfRange(readseq, numclip, readseq.length);
            } else {
                // clip from the end
                readseq = Arrays.copyOfRange(readseq, 0, readseq.length - numclip);
            }
        }

        // compare the sequences to count substitutions 
        int substitutions = 0;
        int nowreadlen = readseq.length;
        for (int i = 0; i < nowreadlen; i++) {
            if (readseq[i] != genomeseq[i]) {
                substitutions++;
            }
        }

        // reloop and output
        if (substitutions == 0) {
            entry.penalty = 0;
            return;
        }

        // loop and collect information about the type of substitution
        for (int i = 0; i < nowreadlen; i++) {
            if (readseq[i] != genomeseq[i]) {
                ThesaurusAnchor anchor = new ThesaurusAnchor();

                if (record.getReadNegativeStrandFlag()) {
                    anchor.originPosition = entry.originEnd - i;
                    byte[] temp = new byte[1];
                    temp[0] = readseq[i];
                    anchor.originRef = (char) (SequenceComplementer.complement(temp)[0]);
                    temp[0] = genomeseq[i];
                    anchor.originAlt = (char) (SequenceComplementer.complement(temp)[0]);

                    anchor.alignPosition = entry.alignStart + i;
                    anchor.alignRef = (char) genomeseq[i];
                    anchor.alignAlt = (char) readseq[i];

                } else {
                    anchor.originPosition = (entry.originStart + i);
                    anchor.originRef = (char) readseq[i];
                    anchor.originAlt = (char) genomeseq[i];

                    anchor.alignPosition = (entry.alignStart + i);
                    anchor.alignRef = (char) genomeseq[i];
                    anchor.alignAlt = (char) readseq[i];
                }
                entry.addAnchor(anchor);
            }
        }
    }

    /**
     * Outputs one line of the thesaurus
     *
     * The write operation is synchronized on the output stream
     *
     * @param outstream
     * @param entry
     *
     */
    private void outputThesaurusEntry(ThesaurusEntry entry, OutputStream outstream) {
        try {
            byte[] towrite = entry.toString().getBytes();
            synchronized (outstream) {
                outstream.write(towrite);
            }
        } catch (IOException ex) {
            System.out.println("Error writing thesaurus entry: " + ex.getMessage());
        }
    }
}

/**
 * Runnable class that looks for thesaurus entries and writes them into output
 * stream by chromosome.
 *
 * The write operation is synchronized on the output stream.
 */
class DistributeRecordsPassTwoRunnable implements Runnable {

    private final OutputStream[] outstreams;
    private final BlockingQueue<StringBucket> queue;
    private final int bucketspacing;

    public DistributeRecordsPassTwoRunnable(OutputStream[] outstreams,
            BlockingQueue<StringBucket> queue, int bucketspacing) {
        this.outstreams = outstreams;
        this.queue = queue;
        this.bucketspacing = bucketspacing;
    }

    @Override
    public void run() {
        try {
            StringBucket entry;
            while (true) {
                entry = queue.take();
                if (entry.text != null) {
                    int bucket = entry.val / bucketspacing;
                    byte[] towrite = (entry.text + "\n").getBytes();
                    synchronized (outstreams[bucket]) {
                        outstreams[bucket].write(towrite);
                    }
                } else {
                    break;
                }
            }
        } catch (Exception ex) {
            System.out.println("Exception in readoutput: " + ex.getMessage());
        }
    }
}