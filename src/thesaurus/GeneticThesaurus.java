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
package thesaurus;

import thesaurus.make.ThesaurusAlign;
import thesaurus.misc.ThesaurusCountPatterns;
import thesaurus.misc.ThesaurusNewGenome;
import thesaurus.make.ThesaurusSummarize;
import thesaurus.make.ThesaurusGenerateData;
import thesaurus.make.ThesaurusCheck;
import thesaurus.make.ThesaurusCompare;
import thesaurus.make.ThesaurusDetails;
import thesaurus.make.ThesaurusFilter;
import thesaurus.make.ThesaurusWrite;
import thesaurus.make.ThesaurusMerge2Chrom;
import thesaurus.misc.ThesaurusReorder;
import thesaurus.network.ThesaurusNetwork;
import thesaurus.misc.ThesaurusSubset;

/**
 * Gateway to tools associated with generating a thesaurus of genetic variants
 *
 * @author tkonopka
 */
public class GeneticThesaurus {

    // set the version of this Genetic Thesaurus
    final static String version = "0.1.0";
    // set some default values
    public static final int DEFAULT_READLEN = 100;
    public static final int DEFAULT_DIVISOR = 10;
    public static final int DEFAULT_OFFSET = 0;
    public static final int DEFAULT_PENALTY = 10;
    public static final String DEFAULT_GENOME = "NA";
    //public static final String DEFAULT_BLATPATH = "blat";
    //public static final String DEFAULT_BLATOPTIONS = "-t=dna -tileSize=12 -minMatch=3 -minScore=90 -noHead -repMatch=1024";
    //public static final Boolean DEFAULT_KEEPPSL = false;

    public static String getVersion() {
        return version;
    }

    public static void printHelp() {
        System.out.println("GeneticThesaurus: construct and apply a thesaurus of genetic variation");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar TYPE options");
        System.out.println();
        System.out.println("Making a thesaurus:");
        System.out.println("  align                    - create basis for thesaurus by aligning reads onto genome");
        System.out.println("  generate                 - generate reads");
        System.out.println("  reorder                  - reorder chromosomes in a theausus");
        System.out.println("  write                    - use bam files to write a thesaurus");

        System.out.println("\nUsing a thesaurus:");
        System.out.println("  compare                  - compare vcf files annotated with thesaurus");
        System.out.println("  filter                   - annotate vcf file using thesaurus");
        System.out.println("  network                  - network analysis of thesaurus-annotated calls");

        System.out.println("\nDevelopment:");
        System.out.println("  check                    - perform checks of thesuarus entries");
        System.out.println("  configure                - set default values");
        System.out.println("  countpatterns            - count items matching pairs of patterns");
        System.out.println("  details                  - get details on loci");
        System.out.println("  merge2chrom              - merge many bam files and automatically split by chromosome");
        System.out.println("  subset                   - create a smaller thesaurus file");
        System.out.println("  summarize                - create bed or coverage files from a thesaurus file");
        System.out.println("  version                  - version");
        System.out.println();
        System.out.println("\nAuthor: Tomasz Konopka (tkonopka@cemm.oeaw.ac.be)");
        System.out.println();
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args == null || args.length == 0) {
            printHelp();
            return;
        }

        // copy all but one of the argument into a different arguments array
        String[] newargs = new String[args.length - 1];
        for (int i = 0; i < newargs.length; i++) {
            newargs[i] = args[i + 1];
        }

        // now call one of the tools in this toolkit
        if (args[0].equalsIgnoreCase("configure")) {
            new ThesaurusConfigure(newargs).run();
        } else if (args[0].equalsIgnoreCase("countpatterns")) {
            new ThesaurusCountPatterns(newargs).run();

            // The next few are related to mapping thesaurus
        } else if (args[0].equalsIgnoreCase("align")) {
            new ThesaurusAlign(newargs).run();
        } else if (args[0].equalsIgnoreCase("generate")) {
            new ThesaurusGenerateData(newargs).run();
        } else if (args[0].equalsIgnoreCase("write")) {
            new ThesaurusWrite(newargs).run();
        } else if (args[0].equalsIgnoreCase("merge2chrom")) {
            new ThesaurusMerge2Chrom(newargs).run();
        } else if (args[0].equalsIgnoreCase("summarize")) {
            new ThesaurusSummarize(newargs).run();
        } else if (args[0].equalsIgnoreCase("filter")) {
            new ThesaurusFilter(newargs).run();
        } else if (args[0].equalsIgnoreCase("network")) {
            new ThesaurusNetwork(newargs).run();
        } else if (args[0].equalsIgnoreCase("details")) {
            new ThesaurusDetails(newargs).run();
        } else if (args[0].equalsIgnoreCase("reorder")) {
            new ThesaurusReorder(newargs).run();
        } else if (args[0].equalsIgnoreCase("subset")) {
            new ThesaurusSubset(newargs).run();
        } else if (args[0].equalsIgnoreCase("compare")) {
            new ThesaurusCompare(newargs).run();
        } else if (args[0].equalsIgnoreCase("check")) {
            new ThesaurusCheck(newargs).run();

            // creates a new genome by adding variants into a reference genome
        } else if (args[0].equalsIgnoreCase("newgenome")) {
            new ThesaurusNewGenome(newargs).run();

        } else if (args[0].equalsIgnoreCase("version")) {
            System.out.println("GeneticThesaurus v" + version);
        } else if (args[0].equalsIgnoreCase("help")) {
            printHelp();
        } else {
            System.out.println("Unrecognized command: " + args[0]);
        }
    }
}
