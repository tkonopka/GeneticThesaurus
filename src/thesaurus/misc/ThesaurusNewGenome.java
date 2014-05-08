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
package thesaurus.misc;

import thesaurus.GeneticThesaurus;
import thesaurus.util.ThesaurusIO;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.Preferences;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.sequence.FastaReader;

/**
 * This tool constructs a reference Genome fasta file by introducing point
 * substitution mutations into an existing reference genome.
 *
 * Mutations are randomly placed with a set density. Mutations occur equally
 * likely between AT,AC,AG, etc. combinations.
 *
 *
 * @author tkonopka
 */
public class ThesaurusNewGenome implements Runnable {

    private File genome;
    private String output;
    private double density = 0.001;
    private int seed = 0;
    boolean isok = false;

    private void printNewGenomeHelp() {
        System.out.println("GeneticThesaurus newgenome: create a new genome fasta file with point substitutions");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar newgenomes ");
        System.out.println();
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--density <double>", "density of random substitutions. Use this to create a genome with randomly placed SNVs.");
        ThesaurusIO.printHelpItem("--seed <integer>", "seed for random number generation");

        System.out.println();
    }

    /**
     * Extract parameter values from the command line
     *
     * @param args
     * @return
     */
    private boolean parseNewGenomeParameters(String[] args) {
        OptionParser prs = new OptionParser();


        // input and output 
        //prs.accepts("vcf").withRequiredArg().ofType(File.class);
        prs.accepts("genome").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("density").withRequiredArg().ofType(Double.class);
        prs.accepts("seed").withRequiredArg().ofType(Integer.class);

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
            // check that file exists
            if (!genome.exists() | !genome.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist: " + genome.getAbsolutePath());
                return false;
            }
        }

        if (options.has("density")) {
            density = Math.abs((Double) options.valueOf("density"));
        } else {
            System.out.println("Missing required argument density");
            return false;
        }

        if (options.has("seed")) {
            seed = Math.abs((Integer) options.valueOf("seed"));
        }

        if (options.has("output")) {
            output = (String) options.valueOf("output");
        } else {
            System.out.println("Missing required argument output");
            return false;
        }

        return true;
    }

    public ThesaurusNewGenome(String[] args) {
        if (args == null || args.length == 0) {
            printNewGenomeHelp();
            return;
        }

        // get a default genome first
        Preferences prefs = Preferences.userNodeForPackage(GeneticThesaurus.class);
        genome = new File(prefs.get("genome", GeneticThesaurus.DEFAULT_GENOME));

        isok = parseNewGenomeParameters(args);
    }

    @Override
    public void run() {
        if (!isok) {
            return;
        }

        // check/create the seed for random number generation
        if (seed == 0) {
            seed = Math.abs((new Random()).nextInt());
        }
        Random rng = new Random(seed);

        // build the new genome one chromosome at a time
        try {
            FastaReader freader = new FastaReader(BufferedReaderMaker.makeBufferedReader(genome));
            OutputStream vcfstream = OutputStreamMaker.makeOutputStream(output + ".vcf");
            OutputStream newgenomestream = OutputStreamMaker.makeOutputStream(output + ".fa");

            // output some lines to the vcfstream
            writeVcfHeader(vcfstream, density, seed);

            while (freader.hasNext()) {
                freader.readNext();
                makeChromosomeWithRandomSNVs(freader, vcfstream, newgenomestream, density, rng);
            }
            freader.close();
            vcfstream.close();
            newgenomestream.close();
        } catch (IOException ex) {
            System.out.println("exception: " + ex.getMessage());
            Logger.getLogger(ThesaurusNewGenome.class.getName()).log(Level.SEVERE, null, ex);
        }


    }

    /**
     * Writes some header lines into a vcf file. The header indicates the
     * density at which the variants were introduced as well as the seed for
     * random number generation.
     *
     * @param vcfstream
     * @param density
     * @param seed
     * @throws IOException
     */
    private void writeVcfHeader(OutputStream vcfstream, double density, int seed) throws IOException {
        StringBuilder sb = new StringBuilder();

        sb.append("##fileformat=VCFv4.1\n");
        sb.append("##random.SNVs.density=").append(density).append("\n");
        sb.append("##random.SNVs.seed=").append(seed).append("\n");
        // write out the canonical line signaling start of the variant list
        sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        vcfstream.write(sb.toString().getBytes());
    }

    /**
     * @param rng
     * @return
     *
     * a random letter from alphabet ATCG
     *
     */
    private byte getRandomBase(Random rng) {
        int xx = rng.nextInt(4);
        switch (xx) {
            case 0:
                return 'A';
            case 1:
                return 'T';
            case 2:
                return 'C';
            case 3:
                return 'G';
            default:
                return 'N';
        }
    }

    /**
     * Get a radom letter from alphabet ATCG that is not equal to x
     *
     * @param rng
     * @param x
     * @return
     */
    private byte getRandomBaseNotX(Random rng, byte x) {
        byte y = getRandomBase(rng);
        while (y == x) {
            y = getRandomBase(rng);
        }
        return y;
    }

    /**
     * One of the main functions of this tool. It looks at the currently loaded
     * chromosome and introduces the point substitutions.
     *
     * After introducing the mutations, the chromosome is written to a stream.
     * The introduced mutations are also logged into another output stream.
     *
     *
     * @param freader
     * @param vcfstream
     * @param newgenomestream
     * @param density
     * @param rng
     * @throws IOException
     */
    private void makeChromosomeWithRandomSNVs(FastaReader freader, OutputStream vcfstream,
            OutputStream newgenomestream, double density, Random rng) throws IOException {

        String chromname = freader.getChromosomeName();
        int chromlen = freader.getChromosomeLength();
        int numSNVs = (int) Math.round(density * (double) chromlen);

        // copy the sequence of the chromsome
        byte[] fullchrom = freader.getSequenceBase0(0, chromlen);

        // randomly mutate sequences and keep track of the locations
        TreeMap<Integer, String> changepositions = new TreeMap<Integer, String>();
        int numchanged = 0;
        while (numchanged != numSNVs) {
            int pos = rng.nextInt(chromlen);
            if (!changepositions.containsKey(pos)) {
                numchanged++;
                byte oldchar = fullchrom[pos];
                if (oldchar != 'N') {
                    byte newchar = getRandomBaseNotX(rng, oldchar);
                    fullchrom[pos] = newchar;
                    changepositions.put(pos, "" + (char) oldchar + "" + (char) newchar);
                }
            }
        }


        // write lines of the vcf file
        StringBuilder sb = new StringBuilder();
        String threedots = "\t.\t.\t.\n";
        for (Map.Entry pairs : changepositions.entrySet()) {
            int pos = (Integer) pairs.getKey();
            String refalt = (String) pairs.getValue();
            sb.append(chromname).append("\t").append(pos + 1).append("\t.\t").
                    append(refalt.substring(0, 1)).append("\t").
                    append(refalt.substring(1, 2)).append(threedots);
        }
        vcfstream.write(sb.toString().getBytes());


        // write out the new genome
        sb = new StringBuilder();
        sb.append(">").append(freader.getChromosomeName()).append("\n");
        int nowpos = 0;
        while (nowpos < chromlen) {
            if (nowpos + 50 < chromlen) {
                sb.append(new String(Arrays.copyOfRange(fullchrom, nowpos, nowpos + 50)));
            } else {
                sb.append(new String(Arrays.copyOfRange(fullchrom, nowpos, chromlen)));
            }
            sb.append("\n");
            nowpos += 50;
        }
        newgenomestream.write(sb.toString().getBytes());
    }
}
