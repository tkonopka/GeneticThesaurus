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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.regions.GenomeBitSet;

/**
 * Tool that extract lines of thesaurus for desired "align" regions. For
 * example, to create an exome-thesaurus from a whole genome thesaurus.
 *
 *
 * @author tkonopka
 */
public class ThesaurusSubset implements Runnable {

    private File inputfile = null;
    private File bedfile = null;
    private File genomefile = null;
    private String output = "stdout";
    private final ArrayList<String> bedregions = new ArrayList<String>();
    private boolean isok = false;

    private void printSubsetHelp() {
        System.out.println("GeneticThesaurus subset: extract lines from thesaurus matching regions of interest");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar subset ");
        System.out.println();
        System.out.println("  --thesaurus <File>       - thesaurus file");
        System.out.println("  --bed <File>             - bed file with regions of interest");
        System.out.println("  --region <String>        - interval of interest (eg. chr1:1000-2000)");
        System.out.println("  --genome <File>          - file with genome");
        System.out.println("  --output <File>          - where to save output [default stdout]");
        System.out.println(" ");
    }

    private boolean parseSubsetParameters(String[] args) {

        OptionParser prs = new OptionParser();

        prs.accepts("thesaurus").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("bed").withRequiredArg().ofType(File.class);
        prs.accepts("region").withRequiredArg().ofType(String.class);
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // extract input/output settings
        if (options.has("thesaurus")) {
            inputfile = (File) options.valueOf("thesaurus");
            if (!inputfile.exists() || !inputfile.canRead()) {
                System.out.println("Cannot read thesaurus file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required parameter: thesaurus");
            return false;
        }

        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }

        if (options.has("bed")) {
            bedfile = (File) options.valueOf("bed");
            if (!bedfile.exists() || !bedfile.canRead()) {
                System.out.println("Cannot read bed file, or file does not exist");
                return false;
            }
        }

        if (options.has("region")) {
            bedregions.addAll((List<String>) options.valuesOf("region"));
        }
        if (options.has("genome")) {
            genomefile = (File) options.valueOf("genome");
            if (!genomefile.exists() || !genomefile.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required parameter: genome");
            return false;
        }

        // check that at least one region or bed was specified
        if (bedregions.isEmpty() && bedfile == null) {
            System.out.println("Must specify at least one region or bed file");
            return false;
        }

        return true;
    }

    public ThesaurusSubset(String[] args) {
        if (args == null || args.length == 0) {
            printSubsetHelp();
            return;
        }
        isok = parseSubsetParameters(args);
    }

    /**
     *
     */
    @Override
    public void run() {
        if (!isok) {
            return;
        }

        try {
            subsetThesaurus(inputfile, output, genomefile, bedfile, bedregions);
        } catch (IOException ex) {
            System.out.println("Something went wrong when subsetting thesaurus");
        }

    }

    private void subsetThesaurus(File thesaurusfile, String output, File genomefile,
            File bedfile, ArrayList<String> bedregions) throws IOException {

        // create a genome bitset        
        GenomeBitSet bgs = new GenomeBitSet(new GenomeInfo(genomefile));

        // get regions of interest from file or from the arraylist
        try {
            if (bedfile != null) {
                BufferedReader br = BufferedReaderMaker.makeBufferedReader(bedfile);
                String s;
                while ((s = br.readLine()) != null) {
                    String[] tokens = s.split("\t");
                    bgs.set(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
                }
                br.close();
            }

        } catch (IOException ex) {
            System.out.println("Could not read bed: " + ex.getMessage());
            return;
        }

        try {
            for (int i = 0; i < bedregions.size(); i++) {
                String[] tokens = bedregions.get(i).split(":|-");
                bgs.set(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
            }
        } catch (Exception ex) {
            System.out.println("Could not parse the regions of interest: " + ex.getMessage());
            return;
        }

        BufferedReader tr = BufferedReaderMaker.makeBufferedReader(thesaurusfile);
        OutputStream os = OutputStreamMaker.makeOutputStream(output);

        // get the header from the thesaurus
        String s = tr.readLine();
        StringBuilder sb = new StringBuilder();
        while (s != null && s.startsWith("#")) {
            sb.append(s).append("\n");
            s = tr.readLine();
        }
        // ad this point s should contain the header line in the thesaurus

        // add one line to the header
        if (bedfile != null) {
            sb.append("##Subset bed file: ").append(bedfile.getAbsolutePath()).append("\n");
        }
        for (int i = 0; i < bedregions.size(); i++) {
            sb.append("##Subset region: ").append(bedregions.get(i)).append("\n");
        }
        sb.append(s).append("\n");
        os.write(sb.toString().getBytes());

        // read lines from the thesaurus one by one and check them against the thesaurus                
        while ((s = tr.readLine()) != null) {
            String[] tokens = s.split("\t", 4);
            BitSet tempset = bgs.get(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
            if (tempset.cardinality() > 0) {
                os.write((s + "\n").getBytes());
            }
        }

        tr.close();
        os.close();
    }
}
