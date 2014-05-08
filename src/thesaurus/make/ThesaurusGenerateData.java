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

import thesaurus.util.ThesaurusIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.sequence.FastaReader;
import jsequtils.sequence.SequenceComplementer;

/**
 * Systematically generate fasta reads from a genome and place them into files.
 *
 * This program works in a consumer/producer design and runs using two threads.
 *
 *
 * @author tkonopka
 */
public class ThesaurusGenerateData extends ThesaurusMapTool {

    // get some defaults values from preset preferences        
    // setup for generating reads
    private int maxperfile = 3000000;
    private String output = "stdout";
    private File bedfile = null;
    HashMap<String, ArrayList<TwoInt>> bedcontents = null;
    // for paired-end reads
    private int insertSize = 0;

    private void printGenerateDataHelp() {
        System.out.println("GeneticThesaurus generate: generate perfect reads from a genome");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar generate ");
        System.out.println();
        ThesaurusIO.printHelpItem("--bed <File>", "bed file specifying regions of genome to consider");
        ThesaurusIO.printHelpItem("--divisor <int>", "generate reads at regular intervals");
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--insertsize <int>", "insert size for paired reads [default " + insertSize + " for single-end reads]");
        ThesaurusIO.printHelpItem("--maxperfile <int>", "maximum number of reads to place in one output file [default " + maxperfile + "]");
        ThesaurusIO.printHelpItem("--offset <int>", "generate reads at regular intervals (coordinate/divisor)+offset [default " + offset + "]");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--readlen <int>", "read length");
        System.out.println();
    }

    public ThesaurusGenerateData(String[] args) {
        if (args == null || args.length == 0) {
            printGenerateDataHelp();
            return;
        }
        super.loadDefaults();
        super.setOk(parseParameters(args));
    }

    private boolean parseParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // genome - input genome fasta
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // settings for generating reads
        prs.accepts("divisor").withRequiredArg().ofType(Integer.class);
        prs.accepts("offset").withRequiredArg().ofType(Integer.class);
        prs.accepts("readlen").withRequiredArg().ofType(Integer.class);
        prs.accepts("maxperfile").withRequiredArg().ofType(Integer.class);
        prs.accepts("bed").withRequiredArg().ofType(File.class);
        prs.accepts("insertsize").withRequiredArg().ofType(Integer.class);

        // settings for producing output
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

        if (options.has("bed")) {
            bedfile = (File) options.valueOf("bed");
            if (!bedfile.exists() || !bedfile.canRead()) {
                System.out.println("Cannot read bed file, or file does not exist");
                return false;
            }
        }

        if (options.has("divisor")) {
            divisor = Math.abs((int) (Integer) options.valueOf("divisor"));
        }

        if (options.has("offset")) {
            offset = Math.abs((int) (Integer) options.valueOf("offset"));
        }

        if (options.has("readlen")) {
            readlen = Math.abs((int) (Integer) options.valueOf("readlen"));
        }

        if (options.has("maxperfile")) {
            maxperfile = Math.abs((int) (Integer) options.valueOf("maxperfile"));
        }

        if (options.has("insertsize")) {
            insertSize = (int) (Integer) options.valueOf("insertsize");
        }

        if (options.has("output")) {
            output = (String) options.valueOf("output");
        } else {
            System.out.println("Missing required argument output");
            return false;
        }

        return true;
    }

    /**
     * Function reads a bed file and stores its content in a hashmap.
     *
     * @throws IOException
     */
    private void readBedFile() throws IOException {

        // create a non-null hashmap
        bedcontents = new HashMap<String, ArrayList<TwoInt>>();

        // copy contents of the file into the hashmap
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(bedfile);
        String s;
        while ((s = br.readLine()) != null) {
            String[] ssplit = s.split("\t");
            TwoInt twoint = new TwoInt(Integer.parseInt(ssplit[1]), Integer.parseInt(ssplit[2]));
            String nowchrom = ssplit[0];
            if (!bedcontents.containsKey(nowchrom)) {
                bedcontents.put(nowchrom, new ArrayList<TwoInt>());
            }

            ArrayList<TwoInt> aa = bedcontents.get(nowchrom);
            aa.add(twoint);
        }
        br.close();
    }

    @Override
    public void runTool() {

        // perhaps read in all the bed information
        if (bedfile != null) {
            try {
                readBedFile();
            } catch (IOException ex) {
                Logger.getLogger(ThesaurusGenerateData.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        // create a consumer/producer pair
        BlockingQueue<ReadLabelSequence> qq = new LinkedBlockingQueue<ReadLabelSequence>(128);
        ExecutorService service = Executors.newFixedThreadPool(2);
        service.execute(new ThesaurusGeneratorProducer(genome, bedcontents, readlen, divisor, offset, qq, insertSize));
        service.execute(new ThesaurusGeneratorConsumer(output, maxperfile, qq, insertSize > 0));
        try {
            service.shutdown();
            service.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
        } catch (InterruptedException ex) {
            System.out.println("Error while waiting for termination");
        }
    }
}

/**
 * Container class holding a read label and a sequence
 *
 * @author tkonopka
 */
class ReadLabelSequence {

    String label;
    String sequence1;
    String sequence2;

    public ReadLabelSequence(String label, String sequence1, String sequence2) {
        this.label = label;
        this.sequence1 = sequence1;
        this.sequence2 = sequence2;
    }

    public ReadLabelSequence() {
        this.label = null;
        this.sequence1 = null;
        this.sequence2 = null;
    }

    public boolean isOk() {
        return (label != null);
    }
}

/**
 * Container class for two integers
 */
class TwoInt {

    int start;
    int end;

    public TwoInt(int start, int end) {
        this.start = start;
        this.end = end;
    }
}

/**
 * Creates reads and puts them onto a queue
 *
 *
 * @author tkonopka
 */
class ThesaurusGeneratorProducer implements Runnable {

    final BlockingQueue<ReadLabelSequence> queue;
    private final File genome;
    private final HashMap<String, ArrayList<TwoInt>> bedcontents;
    final int readlen;
    final int divisor;
    final int offset;
    final int insertSize;

    public ThesaurusGeneratorProducer(File genome, HashMap<String, ArrayList<TwoInt>> bedcontents,
            int readlen, int divisor, int offset,
            BlockingQueue<ReadLabelSequence> queue, int insertSize) {
        this.bedcontents = bedcontents;
        this.queue = queue;
        this.genome = genome;
        this.readlen = readlen;
        this.divisor = divisor;
        this.offset = offset;
        this.insertSize = insertSize;
    }

    /**
     * counts the number of Ns in the sequence
     *
     * @param nowsequence
     * @return
     */
    private int countNs(byte[] nowsequence) {
        int nn = 0;
        int nowlen = nowsequence.length;
        for (int i = 0; i < nowlen; i++) {
            if (nowsequence[i] == 'N') {
                nn++;
            }
        }
        return nn;
    }

    private void generateForChrom(FastaReader freader) throws InterruptedException {

        // record the details of this chromosome
        String chromname = freader.getChromosomeName();
        int chromlength = freader.getChromosomeLength();

        // make a bitset that will indicate where to generate reads for
        BitSet generateBitset = new BitSet(chromlength);
        if (bedcontents == null) {
            // if user did not specify any regions at all, bedcontets will be null.
            // This interpreted as: generate for the full chromosome
            generateBitset.set(0, chromlength, true);
        } else {
            // if bedcontents does have some information, check if it has something 
            // for this chromosome
            ArrayList<TwoInt> aa = bedcontents.get(chromname);
            if (aa == null) {
                return;
            }
            // record the bed regions into the bitset
            for (int i = 0; i < aa.size(); i++) {
                TwoInt nowtwoint = aa.get(i);
                generateBitset.set(nowtwoint.start, nowtwoint.end, true);
            }
        }

        String nameprefix = ">T" + readlen + ";" + chromname + ":";

        // below generate reads for single or paired-end data
        // in both cases, first check if both mates fall within the bitset requested by user,
        // then extract read sequence from the genome, and finally put the read onto the queue for output
        if (insertSize > 0) {
            // here paired-end data
            int imax = chromlength - insertSize;
            for (int i = offset; i < imax; i += divisor) {
                if (generateBitset.get(i) && generateBitset.get(i + readlen - 1)
                        && generateBitset.get(i + insertSize - readlen) && generateBitset.get(i + insertSize - 1)) {
                    byte[] nowsequence1 = freader.getSequenceBase0(i, i + readlen);
                    byte[] nowsequence2 = SequenceComplementer.complement(freader.getSequenceBase0(i + insertSize - readlen, i + insertSize));
                    boolean goodpair = (countNs(nowsequence1) == 0) && (countNs(nowsequence2) == 0);
                    if (goodpair) {
                        String nowname = nameprefix + (i + 1) + "-" + (i + readlen) + ";"
                                + chromname + ":" + (i + insertSize - readlen + 1) + "-" + (i + insertSize);
                        ReadLabelSequence rls = new ReadLabelSequence(
                                nowname, new String(nowsequence1), new String(nowsequence2));
                        queue.put(rls);
                    }

                }
            }
        } else {
            // here single-end data
            int imax = chromlength - readlen;
            for (int i = offset; i < imax; i += divisor) {
                // check if need to generate at this location
                // check that both the first and last position are in the recorded bitset
                if (generateBitset.get(i) && generateBitset.get(i + readlen - 1)) {
                    // prepare a read/length pair and dump it onto the queue
                    byte[] nowsequence = freader.getSequenceBase0(i, i + readlen);
                    boolean goodpair = countNs(nowsequence) == 0;
                    if (goodpair) {
                        String nowname = nameprefix + (i + 1) + "-" + (i + readlen);
                        ReadLabelSequence rls = new ReadLabelSequence(
                                nowname, new String(nowsequence), null);
                        queue.put(rls);
                    }

                }
            }
        }
    }

    @Override
    public void run() {
        try {
            FastaReader fr = new FastaReader(BufferedReaderMaker.makeBufferedReader(genome));
            while (fr.hasNext()) {
                fr.readNext();
                generateForChrom(fr);
            }
            fr.close();

            // put a poison pill onto the queue to allow the consumer to stop            
            queue.put(new ReadLabelSequence());

        } catch (Exception ex) {
            Logger.getLogger(ThesaurusGeneratorProducer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

/**
 * Gets reads off a queue and puts them into a file. Nuance is that files get
 * only up to a certain number of reads. Once the limit is reached, excesss
 * items are put into a different file.
 *
 * @author tkonopka
 */
class ThesaurusGeneratorConsumer implements Runnable {

    private final BlockingQueue<ReadLabelSequence> queue;
    private final String outprefix;
    private final String outsuffix;
    private final int maxperfile;
    private final boolean paired;
    private OutputStream os1 = null;
    private OutputStream os2 = null;

    public ThesaurusGeneratorConsumer(String outprefix, int maxperfile,
            BlockingQueue<ReadLabelSequence> queue, boolean paired) {
        this.queue = queue;
        this.outprefix = outprefix;
        this.paired = paired;
        outsuffix = ".fa.gz";
        this.maxperfile = Math.max(1, maxperfile);
    }

    @Override
    public void run() {

        int readcount = 0;
        int filecount = 0;
        ReadLabelSequence nowdata;

        // read data from the queue and write it into output files
        try {
            if (paired) {
                while ((nowdata = queue.take()).isOk()) {
                    byte[] nowbyte1 = (nowdata.label + "\n" + nowdata.sequence1 + "\n").getBytes();
                    byte[] nowbyte2 = (nowdata.label + "\n" + nowdata.sequence2 + "\n").getBytes();
                    if (readcount >= maxperfile) {
                        readcount = 0;
                        filecount++;
                        os1.close();
                        os1 = null;
                        os2.close();
                        os2 = null;
                    }
                    if (os1 == null || os2 == null) {
                        os1 = OutputStreamMaker.makeOutputStream(outprefix + "." + filecount + ".1" + outsuffix);
                        os2 = OutputStreamMaker.makeOutputStream(outprefix + "." + filecount + ".2" + outsuffix);
                    }
                    os1.write(nowbyte1);
                    os2.write(nowbyte2);
                    readcount++;
                }
            } else {
                while ((nowdata = queue.take()).isOk()) {
                    byte[] nowbyte = (nowdata.label + "\n" + nowdata.sequence1 + "\n").getBytes();
                    if (readcount >= maxperfile) {
                        readcount = 0;
                        filecount++;
                        os1.close();
                        os1 = null;
                    }
                    if (os1 == null) {
                        os1 = OutputStreamMaker.makeOutputStream(outprefix + "." + filecount + outsuffix);
                    }
                    os1.write(nowbyte);
                    readcount++;
                }
            }

        } catch (Exception ex) {
            Logger.getLogger(ThesaurusGeneratorConsumer.class.getName()).log(Level.SEVERE, null, ex);
        }

        // make sure the output stream is closed if not done already
        try {
            if (os1 != null) {
                os1.close();
            }
            if (os2 != null) {
                os2.close();
            }
        } catch (IOException ex) {
            Logger.getLogger(ThesaurusGeneratorConsumer.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
