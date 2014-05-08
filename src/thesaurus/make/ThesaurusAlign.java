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

import thesaurus.make.ThesaurusEarIndex;
import thesaurus.util.ThesaurusLog;
import thesaurus.util.ThesaurusIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.sequence.FastaReader;
import jsequtils.sequence.SequenceComplementer;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Tool that will ``align'' fasta reads onto a genome and report all possible
 * alignment in a bam file.
 *
 * The ``alignment'' is done with a custom method, i.e. internally in this tool
 * rather than using a third-party aligner like blat.
 *
 * @author tkonopka
 */
public class ThesaurusAlign extends ThesaurusMapTool {

    // earlength is a distance on the left/right of each read that MUST match the genome for an alignment to be triggered
    private int earlength = 5;
    // maxmismatches is the number of mismatches allowed in an alignment
    private int maxmismatches = 5;
    // input/output 
    private File input = null;
    private String output = "";
    // number of processor threads
    private int numthreads = 4;
    // for logging and progress tracking
    private int reportevery = 100000;
    ThesaurusLog theslog = new ThesaurusLog();

    private void printAlignHelp() {
        System.out.println("GeneticThesaurus align: produce alignments");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar align ");
        System.out.println();
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--input <String>", "input fasta file");
        ThesaurusIO.printHelpItem("--mismatches <int>", "maximum number of mismatches allowed in alignments [default " + maxmismatches + "]");
        ThesaurusIO.printHelpItem("--earlength <int>", "ear length [default " + earlength + "]");
        ThesaurusIO.printHelpItem("--reportevery <int>", "print out progress every so many reads [default " + reportevery + "]");
        ThesaurusIO.printHelpItem("--threads <int>", "number of processor threads");
        System.out.println();
    }

    private boolean parseAlignParameters(String[] args) {

        OptionParser prs = new OptionParser();

        // genome 
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // settings for input/output
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("input").withRequiredArg().ofType(File.class);

        // tuning blat runtime and 
        //prs.accepts("readlen").withRequiredArg().ofType(Integer.class);
        prs.accepts("mismatches").withRequiredArg().ofType(Integer.class);
        prs.accepts("earlength").withRequiredArg().ofType(Integer.class);
        prs.accepts("threads").withRequiredArg().ofType(Integer.class);
        prs.accepts("reportevery").withRequiredArg().ofType(Integer.class);

        // now use OptionSet to parse the command line
        OptionSet options;

        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // extract value for each parameter
        if (options.has("genome")) {
            genome = (File) options.valueOf("genome");
            if (!genome.exists() || !genome.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist");
                return false;
            }
        }

        if (options.has("mismatches")) {
            maxmismatches = Math.abs((int) (Integer) options.valueOf("mismatches"));
        }

        if (options.has("threads")) {
            numthreads = Math.max(1, Math.abs((int) (Integer) options.valueOf("threads")));
        }

        if (options.has("earlength")) {
            earlength = Math.max(1, Math.abs((int) (Integer) options.valueOf("earlength")));
        }

        if (options.has("reportevery")) {
            reportevery = Math.abs((int) (Integer) options.valueOf("reportevery"));
        }

        // check for input/output
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        } else {
            System.out.println("Missing required argument output");
            return false;
        }

        if (options.has("input")) {
            input = (File) options.valueOf("input");
            if (!input.exists() || !input.canRead()) {
                System.out.println("Cannot read input file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument input");
            return false;
        }

        return true;
    }

    /**
     * Set up the alignment program. This just loads and parses the command line
     * parameters.
     *
     * @param args
     */
    public ThesaurusAlign(String[] args) {
        if (args.length == 0) {
            printAlignHelp();
            return;
        }
        loadDefaults();
        if (!parseAlignParameters(args)) {
            return;
        }
        super.setOk(true);
    }

    @Override
    void runTool() {

        FastaReader fr = null;
        try {
            fr = new FastaReader(genome);
        } catch (Exception ex) {
            System.out.println("Exception while creating genome fasta reader: " + ex.getMessage());
            return;
        }

        SAMFileHeader samheader = null;
        try {
            samheader = makeGenomeSAMFileHeader(genome, "thesalign");
            samheader.addComment("Alignment by GenThesaurus - up to " + maxmismatches + " mismatches");
        } catch (Exception ex) {
            System.out.println("Error making SAM file header");
            return;
        }

        SAMFileWriter alignSAM;
        alignSAM = new SAMFileWriterFactory().makeSAMOrBAMWriter(
                samheader, true, new File(output + ".bam"));

        theslog.log(true, "Starting alignment");

        while (fr.hasNext()) {
            try {
                fr.readNext(true);
            } catch (IOException ex) {
                System.out.println("Exception while reading the genome: " + ex.getMessage());
            }
            processOneChrom(fr.getChromosomeName(),
                    fr.getSequenceBase0(0, fr.getChromosomeLength()),
                    alignSAM, samheader);
        }

        fr.close();
        alignSAM.close();
        theslog.log(true, "done");
    }

    /**
     * create a header object complete with a HD and SQ lines
     *
     * @param seqinfo
     * @return
     */
    private SAMFileHeader makeGenomeSAMFileHeader(File genome, String readgroup) throws IOException {

        // get sequence information from the genome fai file
        File genomefai = new File(genome.getAbsoluteFile() + ".fai");
        BufferedReader faireader = BufferedReaderMaker.makeBufferedReader(genomefai);
        String s;

        SAMFileHeader header = new SAMFileHeader();
        header.setTextHeader("@HD	VN:1.0 SO:unsorted");
        while ((s = faireader.readLine()) != null) {
            String[] ssplit = s.split("\t");
            header.addSequence(new SAMSequenceRecord(ssplit[0], Integer.parseInt(ssplit[1])));
        }

        // add a read group line to the header with the string readgroup specifying the 
        // the actual RG name, the sample, and library
        SAMReadGroupRecord rgr = new SAMReadGroupRecord(readgroup);
        rgr.setPlatform("ILLUMINA");
        rgr.setLibrary(readgroup);
        rgr.setSample(readgroup);
        header.addReadGroup(new SAMReadGroupRecord(readgroup, rgr));

        return header;
    }

    /**
     *
     * This function manages a producer/processor/reporter chain of threads that
     * map reads from an input file onto a chromosome.
     *
     *
     * @param chromname
     *
     * name of current chromosome
     *
     * @param chromseq
     *
     * sequence of choromosome
     *
     * @param alignSAM
     *
     * output SAM file writer
     *
     * @param samheader
     *
     * header for SAM file (necessary for creating new SAM records)
     *
     */
    private void processOneChrom(String chromname, byte[] chromseq,
            SAMFileWriter alignSAM, SAMFileHeader samheader) {
        // build an index for this chromsome
        int chromlen = chromseq.length;

        theslog.log(true, chromname);
        //theslog.log(true, "Building index for chromosome " + chromname);
        ThesaurusEarIndex index = new ThesaurusEarIndex(chromname, chromseq, earlength);
        //theslog.log(true, "Starting alignment");

        // queues for the various stages of the pipeline
        BlockingQueue<OneRead> q1 = new LinkedBlockingQueue<OneRead>(8 * numthreads);
        BlockingQueue<OneReadWithAlignment> q2 = new LinkedBlockingQueue<OneReadWithAlignment>(8 * numthreads);

        // create service for producing, processing, and reporting alignments
        ExecutorService service = Executors.newFixedThreadPool(2 + numthreads);
        // add processes for reading fasta and outputing alignment, then fill the remaining
        // threads with threads that use the index to align
        service.execute(new OneReadProducer(input, q1, reportevery));
        for (int i = 0; i < numthreads; i++) {
            service.execute(new OneReadAligner(index, q1, q2, maxmismatches));
        }
        service.execute(new OneReadWithAlignmentReporter(alignSAM, samheader, q2, numthreads));

        // let all the threads finish and 
        try {
            service.shutdown();
            service.awaitTermination(1000, java.util.concurrent.TimeUnit.DAYS);
        } catch (Exception ex) {
            System.out.println("Exception shutting down service: " + ex.getMessage());
        }

    }
}

/**
 * Runnable for parallel processing, which gets reads from an input file and
 * pushes them onto a queue for other threads to process
 *
 * @author tkonopka
 */
class OneReadProducer implements Runnable {

    private FastaReadReader frr = null;
    private final BlockingQueue<OneRead> queue;
    private final int reportevery;

    /**
     *
     * @param input
     *
     * input file to scane
     *
     * @param queue
     *
     * queue on which reads will be put
     *
     * @param numpoisonpills
     *
     * number of poison pills to put onto the queue (to make sure all consumers
     * get at least one poison pill)
     *
     */
    public OneReadProducer(File input, BlockingQueue<OneRead> queue, int reportevery) {
        this.queue = queue;
        this.reportevery = reportevery;
        try {
            frr = new FastaReadReader(input);
        } catch (Exception ex) {
            System.out.println("Exception while creating fasta read reader");
            frr = null;
        }
    }

    @Override
    public void run() {

        int counter = 0;

        try {
            // avoid work if fasta did not initialize properly
            if (frr == null) {
                queue.put(new OneRead());
                return;
            }

            // get all reads from file and push them onto the queue            
            while (true) {
                OneRead read = frr.getRead();
                if (read == null) {
                    break;
                }
                queue.put(read);

                counter++;
                if (counter % reportevery == 0) {
                    System.out.println(read.getReadname());
                }
            }

            // finish up, signal to the queue that no more reads are coming, and close the stream             
            queue.put(new OneRead());

        } catch (Exception ex) {
            System.out.println("Exception in producing reads: " + ex.getMessage());
        }

        frr.close();
    }
}

/**
 * Runnable in a thread, which looks on a queue for reads with alignments and
 * output them into a SAM file.
 *
 * @author tkonopka
 */
class OneReadWithAlignmentReporter implements Runnable {

    private final BlockingQueue<OneReadWithAlignment> queue;
    private final SAMFileHeader samheader;
    private final SAMFileWriter alignSAM;
    private final int poisonpills;

    /**
     *
     * @param alignSAM
     * @param samheader
     * @param queue
     * @param poisonpills
     *
     * runnable will wait for more than one poisonpills before it really stops
     *
     */
    public OneReadWithAlignmentReporter(SAMFileWriter alignSAM, SAMFileHeader samheader,
            BlockingQueue<OneReadWithAlignment> queue, int poisonpills) {
        this.alignSAM = alignSAM;
        this.samheader = samheader;
        this.queue = queue;
        this.poisonpills = poisonpills;
    }

    @Override
    public void run() {

        // count the number of poison pills received in queue
        int pills = 0;

        try {
            OneReadWithAlignment one;
            while (true) {
                one = queue.take();
                if (one.isOk()) {
                    SAMRecord rec = one.getSAMRecord(samheader);
                    alignSAM.addAlignment(rec);
                } else {
                    pills++;
                    if (pills >= poisonpills) {
                        break;
                    }
                }
            }
        } catch (Exception ex) {
            System.out.println("Exception in readoutput: " + ex.getMessage());
        }

    }
}

/**
 * Runnable that will look at individual reads and align them using an index of
 * ears and sequence of the chromosome
 *
 * @author tkonopka
 */
class OneReadAligner implements Runnable {

    private final BlockingQueue<OneRead> inqueue;
    private final BlockingQueue<OneReadWithAlignment> outqueue;
    private final ThesaurusEarIndex index;
    private final int maxmismatches;

    public OneReadAligner(ThesaurusEarIndex index, BlockingQueue<OneRead> inqueue,
            BlockingQueue<OneReadWithAlignment> outqueue, int maxmismatches) {
        this.inqueue = inqueue;
        this.outqueue = outqueue;
        this.index = index;
        this.maxmismatches = maxmismatches;
    }

    @Override
    public void run() {
        try {
            OneRead one;
            while (true) {
                one = inqueue.take();
                if (one.isOk()) {
                    alignAndOutputRead(one);
                } else {
                    // put it back! so that maybe other aligners can see it and also break                    
                    inqueue.put(one);
                    // Also send a poison pill to the report thread listing on outqueue
                    outqueue.put(new OneReadWithAlignment());
                    break;
                }
            }
        } catch (Exception ex) {
            System.out.println("Exception in readoutput: " + ex.getMessage());
        }
    }

    /**
     * uses the index to produce candidate alignment positions. Then creates
     * objects holding the read sequence and alignment setup, and pushes them
     * onto the output queue.
     *
     * @param read
     */
    public void alignAndOutputRead(OneRead read) {

        // get some book-keeping variables
        int readlen = read.readsequence.length;
        int earlength = index.getEarlength();

        // look at the start/end ear sequences
        int startcode = index.getCode(read.readsequence, 0, earlength);
        int endcode = index.getCode(read.readsequence, readlen - earlength, readlen);

        if (startcode < 0 || endcode < 0) {
            return;
        }

        // look for alignments of the read on the positive strand
        ArrayList<Integer> alignplus = index.alignSequence(read.readsequence, maxmismatches);
        ArrayList<Integer> alignminus = index.alignSequence(SequenceComplementer.complement(read.readsequence), maxmismatches);

        // output the alignments onto the queue
        try {
            for (int i = 0; i < alignplus.size(); i += 2) {
                int astart = alignplus.get(i);
                int mapqual = alignplus.get(i + 1);
                OneReadWithAlignment onewitha = new OneReadWithAlignment(read, index.getSeqname(), astart, mapqual, false);
                outqueue.put(onewitha);
            }
            for (int i = 0; i < alignminus.size(); i += 2) {
                int astart = alignminus.get(i);
                int mapqual = alignminus.get(i + 1);
                OneReadWithAlignment onewitha = new OneReadWithAlignment(read, index.getSeqname(), astart, mapqual, true);
                outqueue.put(onewitha);
            }
        } catch (Exception ex) {
            System.out.println("Exception while putting alignments onto queue: " + ex.getMessage());
        }

    }
}

/**
 * Class that can get reads from an input file. Assumes reads are fasta. All
 * read sequence fits on one line.
 *
 * @author tkonopka
 */
class FastaReadReader {

    private final BufferedReader br;

    public FastaReadReader(File infile) throws IOException {
        br = BufferedReaderMaker.makeBufferedReader(infile);
    }

    public OneRead getRead() {
        try {
            String readname = br.readLine();
            byte[] readseq = br.readLine().getBytes();
            return new OneRead(readname, readseq);
        } catch (Exception ex) {
            return null;
        }
    }

    public void close() {
        try {
            br.close();
        } catch (IOException ex) {
            System.out.println("Exception while close FastaReadReader");
        }
    }
}

/**
 * Simple container for a read name and sequence
 *
 * @author tkonopka
 */
class OneRead {

    final String readname;
    final byte[] readsequence;

    public OneRead() {
        readname = null;
        readsequence = null;
    }

    public OneRead(String readname, byte[] readsequence) {
        this.readname = readname;
        this.readsequence = readsequence;
    }

    public String getReadname() {
        return readname.substring(1);
    }

    public boolean isOk() {
        return (readname != null);
    }
}

/**
 * a class that holds a read name, sequence, and attributes for alignment
 *
 * @author tkonopka
 */
class OneReadWithAlignment extends OneRead {

    final private String chrom;
    final private int alignpos;
    final private int mappingquality;
    final private boolean negativestrand;

    public OneReadWithAlignment() {
        super();
        chrom = null;
        alignpos = -1;
        mappingquality = -1;
        negativestrand = true;
    }

    public OneReadWithAlignment(OneRead read, String chrom, int alignpos,
            int mappingquality, boolean negativestrand) {
        super(read.readname, read.readsequence);
        this.chrom = chrom;
        this.alignpos = alignpos;
        this.mappingquality = mappingquality;
        this.negativestrand = negativestrand;
    }

    private String makeBaseQuality(int size) {
        byte[] ans = new byte[size];
        for (int i = 0; i < size; i++) {
            ans[i] = 's';
        }
        return new String(ans);
    }

    /**
     * checks if the alignment is primary or not using read name (Warning - only
     * properly formatted thesaurus reads will match this)
     *
     * @return
     */
    private boolean isPrimary() {
        String[] nametokens = getReadname().split(";|:|-");
        if (nametokens[1].equals(chrom)) {
            if ((int) Integer.parseInt(nametokens[2]) == alignpos + 1) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    public SAMRecord getSAMRecord(SAMFileHeader hh) {
        SAMRecord sr = new SAMRecord(hh);

        // set flags for a neutral read
        sr.setFlags(0);
        sr.setNotPrimaryAlignmentFlag(!isPrimary());
        sr.setReadName(getReadname());
        sr.setReadBases(readsequence);
        sr.setAlignmentStart(alignpos + 1);
        sr.setCigarString(readsequence.length + "M");
        sr.setReferenceName(chrom);
        sr.setMappingQuality(mappingquality);
        if (negativestrand) {
            sr.setReadNegativeStrandFlag(true);
            sr.setReadBases(SequenceComplementer.complement(readsequence));
        } else {
            sr.setReadNegativeStrandFlag(false);
            sr.setReadBases(readsequence);
        }
        sr.setBaseQualityString(makeBaseQuality(readsequence.length));
        sr.setAttribute("RG", "thesalign");

        return sr;
    }
}