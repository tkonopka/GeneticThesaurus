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

//import Thesaurus.ThesaurusConfigure;
import thesaurus.util.ThesaurusIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.sequence.SequenceComplementer;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Tool that will align fasta reads onto a genome with blat, and then convert
 * blat output to bam.
 *
 *
 * @author tkonopka
 *
 * @deprecated use ThesaurusAlign instead
 *
 */
@Deprecated
public class ThesaurusBlat extends ThesaurusMapTool {

    // variables for running the blat process
    protected Process process = null;
    protected ProcessBuilder pb = null;
    // setting to tune input/output etc.        
    private File input = null;
    private String output = "";
    // for internal use
    private final PslComparator pec = new PslComparator();

    private void printBlatHelp() {
        System.out.println("GeneticThesaurus blat: produce alignments");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar blat ");
        System.out.println();
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--input <String>", "input fasta file");
        ThesaurusIO.printHelpItem("--readlen <String>", "input read length");
        ThesaurusIO.printHelpItem("--blatpath <String>", "path blat");
        ThesaurusIO.printHelpItem("--blatoptions <File>", "path to file containing string with blat options");
        ThesaurusIO.printHelpItem("--keeppsl <boolean>", "set to true to keep psl result from blat, false to delete it after computation");
        System.out.println();
    }

    private boolean parseBlatParameters(String[] args) {

        OptionParser prs = new OptionParser();

        // genome 
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // settings for input/output
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("input").withRequiredArg().ofType(File.class);

        // tuning blat runtime and 
        prs.accepts("readlen").withRequiredArg().ofType(Integer.class);
        //prs.accepts("penalty").withRequiredArg().ofType(Integer.class);
        prs.accepts("blatpath").withRequiredArg().ofType(String.class);
        prs.accepts("blatoptions").withRequiredArg().ofType(File.class);
        prs.accepts("keeppsl").withRequiredArg().ofType(Boolean.class);

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

        if (options.has("readlen")) {
            readlen = Math.abs((int) (Integer) options.valueOf("readlen"));
        }

        //if (options.has("penalty")) {
        //    penalty = Math.abs((int) (Integer) options.valueOf("penalty"));
        // }

        // check for input/output
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        } else {
            System.out.println("Missing required argument output");
            return false;
        }

        if (options.has("keeppsl")) {
            keeppsl = (Boolean) options.valueOf("keeppsl");
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

        // check for blat tuning
        if (options.has("blatpath")) {
            blatpath = (String) options.valueOf("blatpath");
        }
        if (options.has("blatoptions")) {
            BufferedReader optionsBR = null;
            try {
                File optionsfile = (File) options.valueOf("blatoptions");
                StringBuilder optionsSB = new StringBuilder();
                optionsBR = BufferedReaderMaker.makeBufferedReader(optionsfile);
                String s;
                while ((s = optionsBR.readLine()) != null) {
                    optionsSB.append(s);
                }
                blatoptions = optionsSB.toString();
            } catch (IOException ex) {
                System.out.println("Error reading blatoptions: " + ex.getMessage());
            } finally {
                try {
                    optionsBR.close();
                } catch (IOException ex) {
                    System.out.println("Error reading blatoptions: " + ex.getMessage());
                }
            }
        }

        return true;
    }

    public ThesaurusBlat(String[] newargs) {

        if (newargs.length == 0) {
            printBlatHelp();
            return;
        }
        loadDefaults();
        if (!parseBlatParameters(newargs)) {
            return;
        }

        // if reached here, set up the blat process via the builder   
        String cmd = blatpath + " " + genome + " " + input.getAbsolutePath() + " " + blatoptions + " " + output + ".psl";
        // make sure the array that is passed on to the process builder does not contain any blank items
        ArrayList<String> cmdsplit = new ArrayList<String>(Arrays.asList(cmd.split(" ")));
        ArrayList<String> cmdsplit2 = new ArrayList<String>();
        for (int i = 0; i < cmdsplit.size(); i++) {
            if (!cmdsplit.get(i).isEmpty()) {
                cmdsplit2.add(cmdsplit.get(i));
            }
        }
        pb = new ProcessBuilder(cmdsplit2);
        super.setOk(true);
    }

    @Override
    public void runTool() {

        // check if the setup was successful, i.e. if a process builder was created
        if (pb == null) {
            return;
        }

        // make sure the output file, if it exists, is deleted first
        File outputfile = new File(output + ".psl");
        if (outputfile.exists()) {
            outputfile.delete();
        }

        Process p = null;
        int blatresult = 0;
        try {
            p = pb.start();
            blatresult = p.waitFor();
        } catch (IOException ex) {
            Logger.getLogger(ThesaurusBlat.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(ThesaurusBlat.class.getName()).log(Level.SEVERE, null, ex);
        }

        if (blatresult != 0) {
            System.out.println("blat did not finish normally");
            return;
        }

        // if reached here and file exists, process the blat output
        outputfile = new File(output + ".psl");
        if (outputfile.exists()) {
            try {
                processBlatOutput(input, outputfile);
            } catch (IOException ex1) {
                Logger.getLogger(ThesaurusBlat.class.getName()).log(Level.SEVERE, null, ex1);
            }
        }

        // get rid of the blat psl output if not needed
        if (!keeppsl) {
            outputfile.delete();
        }

    }

    /**
     * Reads a fasta file and blat psl file in tandem, then converts the psl
     * output to sam format and outputs the reads into bam files indicating the
     * mapping qualities.
     *
     *
     *
     *
     * @param inputreads
     * @param blatoutput
     * @throws IOException
     */
    private void processBlatOutput(File inputreads, File blatoutput) throws IOException {

        // create readers for the inputs
        BufferedReader inputreader = BufferedReaderMaker.makeBufferedReader(inputreads);
        BufferedReader blatreader = BufferedReaderMaker.makeBufferedReader(blatoutput);

        // create streams for the outputs. The good alignments and unmapped readnames
        SAMFileWriter goodSam, badSam, indelSam;
        SAMFileHeader goodHeader = makeGenomeSAMFileHeader(genome, "good");
        SAMFileHeader badHeader = makeGenomeSAMFileHeader(genome, "bad");
        SAMFileHeader indelHeader = makeGenomeSAMFileHeader(genome, "indel");
        goodHeader.addComment("Alignments by blat - unique best matches only");
        badHeader.addComment("Alignments by blat - multiple and sub-optimal matches");
        indelHeader.addComment("Alignments by blat - multiple and sub-optimal matches with indels");

        goodSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
                goodHeader, true, new File(output + ".good.bam"));
        badSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
                badHeader, true, new File(output + ".bad.bam"));
        indelSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
                indelHeader, true, new File(output + ".bad.indel.bam"));

        OutputStream unmapped = OutputStreamMaker.makeOutputStream(output + ".unmapped");

        // read from blat output one line at a times
        String blatline;
        String inputread = inputreader.readLine().substring(1);
        String inputsequence = inputreader.readLine();
        ArrayList<PslEntry> alternatives = new ArrayList<PslEntry>();

        while ((blatline = blatreader.readLine()) != null) {
            PslEntry blatentry = new PslEntry(blatline);

            // check if this blatentry corresponds with the line from the input            
            if (inputread != null) {
                if (!inputread.equals(blatentry.Qname)) {
                    // process the blat output collected so far
                    Collections.sort(alternatives, pec);
                    decideAndOutput(inputsequence, alternatives, goodSam, badSam, indelSam, goodHeader);
                    // make sure the alternatives are cleared for next read.
                    alternatives.clear();

                    // move onto next input read                    
                    inputread = inputreader.readLine().substring(1);
                    inputsequence = inputreader.readLine();
                    // check for equality with current blat query
                    while (inputread != null && !inputread.equals(blatentry.Qname)) {
                        unmapped.write((inputread + "\n").getBytes());
                        inputread = inputreader.readLine().substring(1);
                        inputsequence = inputreader.readLine();
                    }

                }

                if (inputread == null) {
                    System.out.println("Somethign went wrong");
                }
            }

            alternatives.add(blatentry);
        }

        // make sure that any remaining items in the alternatives array are taken care of
        if (alternatives.size() > 0) {
            Collections.sort(alternatives, pec);
            decideAndOutput(inputsequence, alternatives, goodSam, badSam, indelSam, goodHeader);
        }

        // close all the streams
        unmapped.close();
        inputreader.close();
        blatreader.close();
        goodSam.close();
        badSam.close();
        indelSam.close();
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
     * Converts a blat psl entry into a SAM record.
     *
     * @param entry
     * @param bases
     * @param hh
     * @param ns
     *
     * number of blat entries with the same score
     *
     * @param sd
     *
     * rank-dependent interpretation. For first rank, distance to next best
     * match. For other ranks distance from first match.
     *
     * @param bs
     *
     * number of matches in best score
     *
     * @param matchrankA
     *
     * rank of this blat entry in a sorted list of entries (low rank means
     * better alignment, modulo equal matching entries)
     *
     * @param matchrankB
     *
     * alternative rank (optimistic/pessimistic rank)
     *
     * @param readgroup
     *
     * @return
     */
    private SAMRecord createSAMRecord(PslEntry entry, String bases, SAMFileHeader hh,
            int ns, int sd, int bs, int matchrankA, int matchrankB, String readgroup) {

        //System.out.println("in create pair");

        SAMRecord sr = new SAMRecord(hh);

        // set flags for a neutral read
        sr.setFlags(0);
        sr.setAlignmentStart(entry.Tstart + 1);
        sr.setMateAlignmentStart(entry.Tstart + 1);
        sr.setCigarString(entry.makeCigar());
        sr.setReferenceName(entry.Tname);
        sr.setReadName(entry.Qname);
        sr.setMappingQuality(Math.min(254, entry.match));
        if (entry.strand == '-') {
            sr.setReadNegativeStrandFlag(true);
            sr.setReadBases(SequenceComplementer.complement(bases.getBytes()));
        } else {
            sr.setReadNegativeStrandFlag(false);
            sr.setReadBases(bases.getBytes());
        }
        sr.setBaseQualityString(makeBaseQuality(bases.length()));
        sr.setAttribute("ns", ns);
        sr.setAttribute("sd", sd);
        sr.setAttribute("bs", bs);
        sr.setAttribute("ra", matchrankA);
        sr.setAttribute("rb", matchrankB);
        sr.setAttribute("RG", readgroup);

        return sr;
    }

    /**
     * get a fake base-quality string of a certain length
     *
     * @param size
     * @return
     */
    private String makeBaseQuality(int size) {
        byte[] ans = new byte[size];
        for (int i = 0; i < size; i++) {
            ans[i] = 's';
        }
        return new String(ans);
    }

    private void decideAndOutput(String inputsequence, ArrayList<PslEntry> blatentries,
            SAMFileWriter goodSam, SAMFileWriter badSam, SAMFileWriter indelSam,
            SAMFileHeader header) {

        // deal with easiest cases first, where there is nothing or one item in the array
        if (blatentries.isEmpty()) {
            return;
        }
        PslEntry firstentry = blatentries.get(0);
        int firstmatch = firstentry.match;
        if (blatentries.size() == 1) {
            goodSam.addAlignment(createSAMRecord(blatentries.get(0), inputsequence, header,
                    1, 0, firstmatch, 1, 1, "good"));
            return;
        }

        // if reached here, there are multiple possible alignments
        int bes = blatentries.size();

        // get an array of all the match values
        int[] matches = new int[bes];
        for (int i = 0; i < bes; i++) {
            matches[i] = blatentries.get(i).match;
        }

        // get ranks for those matches
        int[] ranksA = getRanks(matches, true);
        int[] ranksB = getRanks(matches, false);

        // then get multiplicities associated with rank values
        int[] multiplicitiesA = new int[bes];
        for (int i = 0; i < bes; i++) {
            multiplicitiesA[ranksA[i] - 1]++;
        }

        // if the best match is separated from the sub-optimal ones, report it as good one
        // Other wise all the entries will go into the bad file.        
        // In practice, check if the an entry with rank 1 exists
        // Reads with indels always go to the special indel bam file

        if (ranksA[0] == 1) {
            int secondmatch = blatentries.get(1).match;
            SAMRecord record = createSAMRecord(firstentry, inputsequence, header,
                    multiplicitiesA[ranksA[0] - 1], 0, firstmatch, ranksA[0], ranksB[0], "bad");
            if (hasLongIndel(record, 0)) {
                indelSam.addAlignment(record);
            } else {
                goodSam.addAlignment(createSAMRecord(firstentry, inputsequence, header,
                        1, firstmatch - secondmatch, firstmatch, ranksA[0], ranksB[0], "good"));
            }

        } else {
            SAMRecord record = createSAMRecord(firstentry, inputsequence, header,
                    multiplicitiesA[ranksA[0] - 1], 0, firstmatch, ranksA[0], ranksB[0], "bad");
            if (hasLongIndel(record, 0)) {
                indelSam.addAlignment(record);
            } else {
                badSam.addAlignment(record);
            }
        }

        // output the remainder of the records
        for (int i = 1; i < bes; i++) {
            PslEntry nowentry = blatentries.get(i);
            SAMRecord record = createSAMRecord(nowentry, inputsequence, header,
                    multiplicitiesA[ranksA[i] - 1], nowentry.match - firstmatch, firstmatch, ranksA[i], ranksB[i], "bad");

            // determine whether to output as a "normal" bad item or as an indel-containing bad item
            if (hasLongIndel(record, 0)) {
                indelSam.addAlignment(record);
            } else {
                badSam.addAlignment(record);
            }
        }
    }

    /**
     * same as hasLongIndel(record, maxlength), with maxlength set to 15.
     *
     * @param record
     * @return
     */
    private boolean hasLongIndel(SAMRecord record) {
        return (hasLongIndel(record, 15));
    }

    /**
     *
     * @param record
     * @param maxlength
     * @return
     *
     * true if the cigar of the record has an indel longer than maxlength.
     *
     */
    private boolean hasLongIndel(SAMRecord record, int maxlength) {
        Cigar cc = record.getCigar();
        for (int k = 0; k < cc.numCigarElements(); k++) {
            CigarElement cele = cc.getCigarElement(k);
            if ((cele.getOperator() == CigarOperator.D || cele.getOperator() == CigarOperator.I)
                    && cc.getCigarElement(k).getLength() > maxlength) {
                return true;
            }
        }
        return false;
    }

    /**
     * get Ranks for a sorted array.
     *
     *
     * @param values
     *
     * Assumed sorted array of values.
     *
     * @param conservative
     *
     * this determines how ranks of duplicate values are recorded. If set to
     * true, function will report pessimistic/conservative ranks, otherwise will
     * report optimistic ranks.
     *
     * e.g. Consider values 2,4,4,10.
     *
     * Set to true to obtain: 1,3,3,4. Set to false to obtain: 1,2,2,4.
     *
     * @return
     *
     * another array with ranks, i.e. numbers in range [1, values.length]
     *
     */
    private int[] getRanks(int[] values, boolean conservative) {
        int[] rank = new int[values.length];

        // get an initial set of ranks, which will work assuming all values are distinct and already ordered
        for (int i = 0; i < values.length; i++) {
            rank[i] = i + 1;
        }


        if (conservative) {
            // now go from the end and decrease ranks on duplicate values
            for (int i = values.length - 2; i >= 0; i--) {
                if (values[i] == values[i + 1]) {
                    rank[i] = rank[i + 1];
                }
            }
        } else {
            // go from the beginning and lower ranks on duplicate values
            for (int i = 1; i < values.length; i++) {
                if (values[i] == values[i - 1]) {
                    rank[i] = rank[i - 1];
                }
            }
        }

        return rank;
    }
}
