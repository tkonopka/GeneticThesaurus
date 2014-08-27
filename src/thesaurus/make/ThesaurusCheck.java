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

import thesaurus.util.ThesaurusIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;

/**
 *
 * Utility for developer stage. Checks if a thesaurus file is ordered.
 *
 * @author tkonopka
 */
public class ThesaurusCheck extends ThesaurusMapTool {

    File thesaurusfile;
    String output = "stdout";

    enum checkWhatEnum {

        ORDER, TRIVIAL
    }
    checkWhatEnum checkwhat = null;

    private void printCheckOrderHelp() {
        System.out.println("GeneticThesaurus checkorder: check that a thesaurus file is properly ordered");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar checkorder ");
        System.out.println();
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--thesaurus <File>", "thesaurus file, merged and sorted");
        ThesaurusIO.printHelpItem("--check <String>", "either order or bad [default order]");
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

        // input and output 
        prs.accepts("check").withRequiredArg().ofType(String.class);


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

        // get prefix for output tracks
        if (options.has("output")) {
            output = (String) options.valueOf("output");
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

        if (options.has("check")) {
            String cw = (String) options.valueOf("check");
            if (cw.equalsIgnoreCase("order")) {
                checkwhat = checkWhatEnum.ORDER;
            } else if (cw.equalsIgnoreCase("trivial")) {
                checkwhat = checkWhatEnum.TRIVIAL;
            } else {
                System.out.println("Unrecognized value for argument --check: " + cw);
                return false;
            }
        }

        return true;
    }

    public ThesaurusCheck(String[] args) {
        if (args.length == 0) {
            printCheckOrderHelp();
            return;
        }
        super.loadDefaults();
        super.setOk(parseFilterParameters(args));
    }

    @Override
    void runTool() {
        // set up input and output streams
        OutputStream outstream;
        BufferedReader br;
        try {
            outstream = OutputStreamMaker.makeOutputStream(output);
            br = BufferedReaderMaker.makeBufferedReader(thesaurusfile);
        } catch (Exception ex) {
            System.out.println("Something went wrong during stream setup: " + ex.getMessage());
            return;
        }

        ThesaurusEntryAlignComparator teac;
        GenomeInfo ginfo = null;
        try {
            ginfo = new GenomeInfo(genome);
            teac = new ThesaurusEntryAlignComparator();
        } catch (Exception ex) {
            System.out.println("Failed to load genomic position comparator: " + ex.getMessage());
            return;
        }

        // read the thesaurus file        

        try {
            switch (checkwhat) {
                case ORDER:
                    runCheckOrder(br, teac, ginfo);
                case TRIVIAL:
                    runCheckBad(br, ginfo);
                default:
            }

        } catch (Exception ex) {
            System.out.println("Error reading thesaurus file: " + ex.getMessage());
        }

    }

    /**
     * runs through thesaurus and checks if the entries are ordered correctly,
     * i.e. chromosome by chromosome and in increasing value of the align start
     * position.
     *
     * @param br
     * @param outstream
     * @param teac
     * @throws IOException
     */
    private void runCheckOrder(BufferedReader br, ThesaurusEntryAlignComparator teac, GenomeInfo ginfo) throws IOException {
        String s;
        ThesaurusEntry last = null, current = null;
        long counter = 0;
        boolean lastbad = false;
        while ((s = br.readLine()) != null) {
            if (!s.startsWith("#") && !s.startsWith("Align")) {
                counter++;
                current = new ThesaurusEntry(s, ginfo);
                if (lastbad) {
                    System.out.println(current.toString());
                }
                if (last != null) {
                    if (teac.compare(last, current) > 0) {
                        System.out.println("Misordered at " + counter + ": ");
                        System.out.print(last.toString());
                        System.out.print(current.toString());
                        lastbad = true;
                    } else {
                        lastbad = false;
                    }
                }
                last = current;
            }
        }
    }

    /**
     * runs through thesaurus and checks if none of the entries are trivial,
     * i.e. map onto themselves
     *
     *
     * @param br
     * @param outstream
     * @throws IOException
     */
    private void runCheckBad(BufferedReader br, GenomeInfo ginfo) throws IOException {
        String s;
        ThesaurusEntry last = null, current = null;
        long counter = 0;
        boolean lastbad = false;
        while ((s = br.readLine()) != null) {
            if (!s.startsWith("#") && !s.startsWith("Align")) {
                counter++;
                current = new ThesaurusEntry(s, ginfo);
                if (last != null) {
                    if (current.isTrivial()) {
                        System.out.println("Trivial at " + counter + ": ");
                        System.out.print(current.toString());
                        lastbad = true;
                    } else {
                        lastbad = false;
                    }
                }
                last = current;
            }
        }
    }
}
