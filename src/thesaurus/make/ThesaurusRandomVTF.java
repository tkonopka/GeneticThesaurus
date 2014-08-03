/*
 * Copyright 2014 Tomasz Konopka.
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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePosition;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.sequence.SequenceMap;
import jsequtils.variants.VCFEntry;
import thesaurus.GeneticThesaurus;
import thesaurus.util.ThesaurusIO;

/**
 * Whereas ThesaurusFilter will produce an annotated vcf and vtf files that
 * should contain useful links between genomic positions, this tool will produce
 * vtf files with random links.
 *
 *
 * @author tkonopka
 */
public class ThesaurusRandomVTF extends ThesaurusMapTool {

    private File bedfile = null;
    private File vcffile = null;
    private File genomefile = null;
    private int seed = 0;
    private String output = "stdout";
    private GenomePositionComparator gcomp = new GenomePositionComparator();
        
    private void printRandomVTFHelp() {
        System.out.println("GeneticThesaurus randomvtf: generate random vtf files");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar randomvtf ");
        System.out.println();
        ThesaurusIO.printHelpItem("--bed <File>", "bed file with possible alternative positions");
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--vcf <File>", "vcf file output by Thesaurus filter");
        ThesaurusIO.printHelpItem("--output <File>", "output with random vtf links");
        ThesaurusIO.printHelpItem("--seed <int>", "seed for random number generation [default 0]");
        System.out.println();
    }

    /**
     *
     * @param args
     * @return
     */
    private boolean parseRandomVTFParameters(String[] args) {

        OptionParser prs = new OptionParser();

        prs.accepts("bed").withRequiredArg().ofType(File.class);
        prs.accepts("vcf").withRequiredArg().ofType(File.class);
        prs.accepts("genome").withRequiredArg().ofType(File.class);
        prs.accepts("seed").withRequiredArg().ofType(Integer.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);

        // now use OptionSet to parse the command line
        OptionSet options;

        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // get genome
        if (options.has("seed")) {
            seed = (Integer) options.valueOf("seed");
            seed = Math.abs(seed);
        }

        // get bed regions of potential link sites
        if (options.has("bed")) {
            bedfile = (File) options.valueOf("bed");
            if (!bedfile.exists() || !bedfile.canRead()) {
                System.out.println("Cannot read bedfile, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument --bed");
            return false;
        }

        // get definition of genome
        if (options.has("genome")) {
            genomefile = (File) options.valueOf("genome");
            if (!genomefile.exists() || !genomefile.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument --genome");
            return false;
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

        // get prefix for output tracks
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }

        return true;
    }

    public ThesaurusRandomVTF(String[] args) {
        if (args.length == 0) {
            printRandomVTFHelp();
            return;
        }
        super.loadDefaults();        
        super.setOk(parseRandomVTFParameters(args));
    }

    @Override
    void runTool() {

        // load the bed information
        GenomeInfo ginfo;
        SimpleBedSet bedset;
        SequenceMap seqmap;
        try {
            ginfo = new GenomeInfo(genomefile);
            bedset = new SimpleBedSet(bedfile, ginfo);
            //bedset.print();
            seqmap = new SequenceMap(genomefile, true);
        } catch (Exception ex) {
            System.out.println("Error loading genome and bed set: " + ex.getLocalizedMessage());
            return;
        }

        // create random number generator
        Random RNG = new Random(seed);

        // set up input and output streams
        OutputStream outstream;
        BufferedReader br;
        try {
            outstream = OutputStreamMaker.makeOutputStream(output);                            
            outstream.write(getOutComment().getBytes());
                        
            br = BufferedReaderMaker.makeBufferedReader(vcffile);
            
            String s;
            while ((s = br.readLine()) != null) {
                if (!s.startsWith("#")) {
                    makeRandomLinksForVariant(new VCFEntry(s, ginfo), bedset, seqmap, RNG, outstream);
                }
            }

            if (outstream != System.out) {
                outstream.close();
            }
            br.close();
        } catch (Exception ex) {
            System.out.println("Something went wrong during random link generation: " + ex.getMessage());
            return;
        }

    }

    /**
     * Generates a comment that is output into the vtf file.
     * 
     * @return 
     */
    private String getOutComment () {
        StringBuilder sb = new StringBuilder();
        sb.append("##Variant Thesaurus File\n");
        sb.append("##Random links generated by GeneticThesaurus v").append(GeneticThesaurus.getVersion()).append("\n");
        sb.append("##Matching filtered variant call file ").append(vcffile.getAbsolutePath()).append("\n");
        sb.append("##Links drawn from bed intervals file ").append(bedfile.getAbsolutePath()).append("\n");
        sb.append("##Seed ").append(seed).append("\n");
        return sb.toString();                
    }
    
    private void makeRandomLinksForVariant(VCFEntry entry, SimpleBedSet bset,
            SequenceMap seqmap, Random RNG, OutputStream out) throws IOException {

        // first check if the entry contains a filter field "thesaurus"
        String[] efilter = entry.getFilter().split(";");
        boolean makerandom = false;
        for (int i = 0; i < efilter.length; i++) {
            if (efilter[i].equals("thesaurus")) {
                makerandom = true;
            }
        }

        // avoid further work if variant is not supposed to have any links
        if (!makerandom) {
            return;
        }

        // at this stage, we know the variant should be thesaurus annotated
        // figure out how many links to make
        int numlinks = 0;
        String[] eformat = entry.getFormat().split(":");
        String[] egeno = entry.getGenotype().split(":");
        for (int i = 0; i < eformat.length; i++) {
            if (eformat[i].equals("TS")) {
                numlinks = Integer.parseInt(egeno[i]);
            }
        }

        if (numlinks == 0) {
            System.out.println("Warning: site has thesaurus filter status but no/zero TS annotation:\n" + entry.toString());
        }

        ArrayList<GenomePosition> links = new ArrayList<GenomePosition>();        
        for (int i = 0; i < numlinks; i++) {
            links.add(getRandomPosition(bset, seqmap, RNG));
        }
        Collections.sort(links, gcomp);

        StringBuilder sb = new StringBuilder();
        sb.append(entry.getChr()).append(":").append(entry.getPosition());
        for (int i=0; i<numlinks; i++) {
            GenomePosition nowpos = links.get(i);
            sb.append("\t").append(nowpos.getChr(bset.ginfo)).append(":").append(nowpos.getPosition());
        }
        sb.append("\n");
        out.write(sb.toString().getBytes());
        
    }

    /**
     * Picks a random position from the SimpleBedSet. This function also checks
     * whether the position is
     *
     * @param bset
     * @param seqmap
     * @param RNG
     * @return
     */
    private GenomePosition getRandomPosition(SimpleBedSet bset, SequenceMap seqmap, Random RNG) {
        while (true) {
            // get a candidate position
            GenomePosition gpos;
            try {
                gpos = bset.getRandomPosition(RNG);
            } catch (Exception ex) {
                System.out.println("AAA\t" + ex.getMessage());
                return null;
            }

            try {
                // check that this position has an ATCG character
                byte gposchar = seqmap.getSequenceBase1(gpos.getChr(bset.ginfo), gpos.getPosition());
                if (gposchar != 'N') {
                    //System.out.println("Got a hit at " + gpos.toString() + "\t" + gpos.toString(bset.ginfo));
                    return gpos;
                } else {
                    //System.out.println("Got an N at " + gpos.toString() + "\t" + gpos.toString(bset.ginfo));
                }
            } catch (Exception ex) {
                System.out.println("BBB\t" + ex.getMessage());
                System.out.println(gpos.toString());
                return null;
            }
        }
    }
}

/**
 * object stores a bed file verbatim (does not collapse repeat intervals)
 *
 * @author tkonopka
 */
class SimpleBedSet {

    int[] chrs;
    int[] start;
    int[] end;
    int[] offset;
    int totsize;
    GenomeInfo ginfo;

    public SimpleBedSet(File bedfile, GenomeInfo ginfo) {

        this.ginfo = ginfo;

        // read through the bedfile and count lines
        int numlines = 0;
        try {
            BufferedReader br = BufferedReaderMaker.makeBufferedReader(bedfile);
            String s;
            while ((s = br.readLine()) != null) {
                if (!s.startsWith("#")) {
                    numlines++;
                }
            }
            br.close();
        } catch (Exception ex) {
            System.out.println("Error counting : " + ex.getLocalizedMessage());
        }

        chrs = new int[numlines];
        start = new int[numlines];
        end = new int[numlines];
        offset = new int[numlines];

        int nowoffset = 0;
        int nowline = 0;

        // read through the bedfile and set the bitsets        
        try {
            BufferedReader br = BufferedReaderMaker.makeBufferedReader(bedfile);
            String s;
            while ((s = br.readLine()) != null) {
                if (!s.startsWith("#")) {
                    String tokens[] = s.split("\t", 5);
                    int nowstart = Integer.parseInt(tokens[1]);
                    int nowend = Integer.parseInt(tokens[2]);
                    int intervallen = Math.abs(Integer.parseInt(tokens[2]) - Integer.parseInt(tokens[1]));
                    // do not record intervals that have zero length
                    if (intervallen > 0) {
                        chrs[nowline] = ginfo.getChrIndex(tokens[0]);
                        start[nowline] = Math.min(nowstart, nowend);
                        end[nowline] = Math.max(nowstart, nowend);
                        nowoffset += Math.abs(end[nowline] - start[nowline]);
                        offset[nowline] = nowoffset;
                        totsize = nowoffset;
                        nowline++;
                    }
                }
            }
            br.close();
        } catch (Exception ex) {
            System.out.println("Error loading bed regions: " + ex.getLocalizedMessage());
        }
    }

    void print() {
        for (int i = 0; i < chrs.length; i++) {
            System.out.println(i + "\t" + ginfo.getChrName(chrs[i]) + "\t" + start[i] + "\t" + end[i] + "\t" + offset[i]);
        }
    }

    /**
     * Selects a random position within this bedset
     *
     * @param RNG
     * @return
     */
    GenomePosition getRandomPosition(Random RNG) {
        // pick a random position
        int nowoffset = RNG.nextInt(totsize);

        // find the row in the bedset that corresponds to the position
        int findpos = Arrays.binarySearch(offset, nowoffset);

        if (findpos < 0) {
            findpos = -findpos - 1;
        }        

        // the +1 is to convert from bed format into a GenomePosition (1-based)
        // start[findpos] is the start position of the interval
        // nowoffset-offset[findpos] is the offset within the interval                
        if (findpos == 0) {
            return (new GenomePosition(chrs[findpos], 1 + start[findpos] + nowoffset));
        } else {
            return (new GenomePosition(chrs[findpos], 1 + start[findpos] + nowoffset - offset[findpos - 1]));
        }
    }
}