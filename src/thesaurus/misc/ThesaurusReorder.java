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
package thesaurus.misc;

import thesaurus.util.ThesaurusLog;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;

/**
 * Small tool that reorders chromosomes in a thesaurus table. This is useful for
 * users who want to download a thesaurus file and the apply it on their data.
 * Such users will have to reorder the thesaurus table to match their own
 * conventions.
 *
 *
 * @author tkonopka
 */
public class ThesaurusReorder implements Runnable {

    private File inputfile = null;
    private File genomefile = null;
    private String output = "stdout";
    private boolean isok = false;
    private boolean verbose = false;
    private final ThesaurusLog mylog = new ThesaurusLog();
    private final int PROGRESS_INTERVAL = 10000000;

    private void printReorderHelp() {
        System.out.println("GeneticThesaurus reorder: reorder chromosomes in a thesaurus table");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar reorder ");
        System.out.println();
        System.out.println("  --thesaurus <File>       - thesaurus file");
        System.out.println("  --genome <File>          - file with genome");
        System.out.println("  --output <File>          - where to save output [default stdout]");
        System.out.println("  --verbose                - print some progress information");
        System.out.println(" ");
    }

    private boolean parseReorderParameters(String[] args) {

        OptionParser prs = new OptionParser();

        prs.accepts("thesaurus").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("genome").withRequiredArg().ofType(File.class);
        prs.accepts("verbose");

        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        verbose = options.has("verbose");

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

        return true;
    }

    public ThesaurusReorder(String[] args) {
        if (args == null || args.length == 0) {
            printReorderHelp();
            return;
        }
        isok = parseReorderParameters(args);
    }

    @Override
    public void run() {
        if (!isok) {
            return;
        }

        // get information about the genome
        GenomeInfo ginfo = null;
        try {
            ginfo = new GenomeInfo(genomefile);
        } catch (IOException ex) {
            System.out.println("Error reading genome information: " + ex.getMessage());
        }

        // create file names for chromosomes
        String[] chrfiles = new String[ginfo.getNumChromosomes()];
        for (int i = 0; i < ginfo.getNumChromosomes(); i++) {
            chrfiles[i] = output + ".temp." + ginfo.getChrName(i) + ".txt.gz";
        }

        try {
            // get the existing header and append a new line annotating new genome
            mylog.log(verbose, "Reorder reading header");
            String header = getNewThesaurusHeader(inputfile, genomefile);
            // Pass 1: read from thesaurus, redistribute to temporary files
            mylog.log(verbose, "Reorder pass 1/2");
            ReorderPassOne(inputfile, ginfo, chrfiles);
            // Pass 2: read from temporary files and collect into one new thesaurus table
            mylog.log(verbose, "Reorder pass 2/2");
            ReorderPassTwo(output, ginfo, chrfiles, header);
            mylog.log(verbose, "Done");
        } catch (Exception ex) {
            System.out.println("Error in reorder: " + ex.getMessage());
        }
    }

    /**
     * Reads first lines from a thesaurus tha
     *
     * @param infile
     * @return
     * @throws IOException
     */
    private String getNewThesaurusHeader(File infile, File genomefile) throws IOException {

        StringBuilder sb = new StringBuilder();
        // read the lines from the input that start with ##
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(infile);
        String s;
        while ((s = br.readLine()) != null) {
            if (s.startsWith("#")) {
                sb.append(s).append("\n");
            } else {
                break;
            }
        }
        br.close();

        // add a new comment with new genome
        sb.append("##reorder.genome=").append(genomefile.getAbsolutePath()).append("\n");
        // append last line read from file, should be "Align.chr ..."
        sb.append(s).append("\n");

        return sb.toString();
    }

    /**
     * read from a thesaurus table and distribute data into multiple files
     *
     * @param infile
     * @param ginfo
     * @param chrfiles
     * @param headerfile
     */
    private void ReorderPassOne(File infile, GenomeInfo ginfo, String[] chrfiles) throws FileNotFoundException, IOException {

        // make an array of output files
        OutputStream[] outchr = new OutputStream[chrfiles.length];
        for (int i = 0; i < chrfiles.length; i++) {
            outchr[i] = OutputStreamMaker.makeOutputStream(chrfiles[i]);
        }

        // start reading the input thesaurus
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(infile);
        String s;
        long counter = 0;
        while ((s = br.readLine()) != null) {
            if (!s.startsWith("#") && !s.startsWith("Align")) {
                String[] tokens = s.split("\t", 5);
                int alignindex = ginfo.getChrIndex(tokens[0]);
                int originindex = ginfo.getChrIndex(tokens[3]);
                if (alignindex >= 0 && originindex >= 0) {
                    outchr[alignindex].write((s + "\n").getBytes());
                }
            }
            counter++;
            if (counter % PROGRESS_INTERVAL == 0) {
                mylog.log(verbose, "progress: " + counter);
            }
        }
        br.close();

        // close the output streams
        for (int i = 0; i < chrfiles.length; i++) {
            outchr[i].close();
        }
    }

    /**
     * read from multiple files and collect into a single item
     *
     * @param outfile
     * @param ginfo
     * @param chrfiles
     * @param headerfile
     *
     */
    private void ReorderPassTwo(String outfile, GenomeInfo ginfo, String[] chrfiles, String header) throws FileNotFoundException, IOException {

        // make a stream for the output
        OutputStream os = OutputStreamMaker.makeOutputStream(outfile);

        // write the header
        os.write(header.getBytes());

        // read from all the chromosome files and write to single file
        for (int i = 0; i < chrfiles.length; i++) {
            mylog.log(verbose, ginfo.getChrName(i));
            BufferedReader chrreader = BufferedReaderMaker.makeBufferedReader(chrfiles[i]);
            String s;
            while ((s = chrreader.readLine()) != null) {
                os.write((s + "\n").getBytes());
            }
            chrreader.close();
            // delete the temporary file
            new File(chrfiles[i]).delete();
        }

        // clean up
        os.close();
    }
}
