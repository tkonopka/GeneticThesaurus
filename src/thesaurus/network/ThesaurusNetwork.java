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
package thesaurus.network;

import thesaurus.util.ThesaurusLog;
import thesaurus.util.ThesaurusIO;
import java.io.File;
import java.io.OutputStream;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.FileExtensionGetter;
import jsequtils.file.OutputStreamMaker;

/**
 * Perform a network analysis based on a vcf file and its associated vtf file
 *
 * @author tkonopka
 */
public class ThesaurusNetwork implements Runnable {

    private final ThesaurusLog mylog;
    private boolean isok = false;
    private String output = "stdout";
    private File vcffile = null;
    private File vtffile = null;
    private File genome = null;

    private void printNetworkHelp() {
        System.out.println("GeneticThesaurus network: network analysis based on vcf and vtf calls");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar network ");
        System.out.println();
        System.out.println("Core options:");
        ThesaurusIO.printHelpItem("--output <String>", "prefix for output files");
        ThesaurusIO.printHelpItem("--genome <File>", "genome fasta file");
        ThesaurusIO.printHelpItem("--vcf <File>", "variant call file");
        ThesaurusIO.printHelpItem("--vtf <File>", "variant thesaurus file [if not specified, will be loaded automatically based on vcf file name]");
        System.out.println();
    }

    /**
     *
     * @param args
     * @return
     */
    private boolean parseNetworkParameters(String[] args) {
        OptionParser prs = new OptionParser();

        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("genome").withRequiredArg().ofType(File.class);
        prs.accepts("vcf").withRequiredArg().ofType(File.class);
        prs.accepts("vtf").withRequiredArg().ofType(File.class);

        // now use OptionSet to parse the command line
        OptionSet options;

        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }


        // get prefix for output tracks
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }

        if (options.has("genome")) {
            genome = (File) options.valueOf("genome");
            if (!genome.exists() || !genome.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist");
                return false;
            }
        }

        // get input vcf file
        if (options.has("vcf")) {
            vcffile = (File) options.valueOf("vcf");
            if (!vcffile.exists() || !vcffile.canRead()) {
                System.out.println("Cannot read vcf file, or file does not exist: " + vcffile.getAbsolutePath());
                return false;
            }
        } else {
            System.out.println("Missing required argument --vcf");
            return false;
        }

        // get matching vtf file
        if (options.has("vtf")) {
            vtffile = (File) options.valueOf("vtf");
        } else {
            // try to get vtf file by parsing vcf file name
            String[] vcftokens = FileExtensionGetter.getExtensionSplit(vcffile);
            vtffile = new File(vcftokens[0] + ".vtf.gz");
        }
        if (!vtffile.exists() || !vtffile.canRead()) {
            System.out.println("Cannot read vtf file, or file does not exist: " + vtffile.getAbsolutePath());
            return false;
        }

        return true;
    }

    /**
     *
     * @param args
     *
     */
    public ThesaurusNetwork(String[] args) {
        if (args.length == 0) {
            printNetworkHelp();
            mylog = null;
            return;
        }
        isok = parseNetworkParameters(args);
        mylog = new ThesaurusLog(System.out);
        mylog.setVerbose(true);
    }

    @Override
    public void run() {
        if (!isok) {
            return;
        }

        // setup for the network
        GenomePositionNetwork gpnet = new GenomePositionNetwork(vcffile, vtffile, genome);

        // compute summary statistics and output 
        try {
            OutputStream outstream = OutputStreamMaker.makeOutputStream(output);
            gpnet.outputClusterIds(outstream);
            if (outstream != System.out) {
                outstream.close();
            }
            //System.out.println(gpnet.size()+"\t"+gpnet.countDisconnectedComponents());
        } catch (Exception ex) {
            System.out.println("Error output results: " + ex.getMessage());
        }


    }
}
