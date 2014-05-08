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
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.regex.Pattern;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;

/**
 * This tool scans an input file and looks for pairs of patterns. It counts the
 * number of lines where the patterns are matched and outputs a summary table
 * with those counts.
 *
 * The output can be later used to perform, eg. a Fisher exact test on the
 * co-occurance of the two patterns.
 *
 * @author tkonopka
 */
public class ThesaurusCountPatterns implements Runnable {

    boolean isok = false;
    private File inputfile = null;
    private String output = "stdout";
    private String patternstringA = ".";
    private String patternstringB = ".";
    private String sep = "\t";
    private String ignore = "#";
    private int columnA = -1;
    private int columnB = -1;
    Pattern pA, pB;

    private void printCountFilteredVariantsHelp() {
        System.out.println("GeneticThesaurus countfiltered: count variants with certain filter status");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar countpatterns ");
        System.out.println();
        System.out.println("  --input <File>           - input file");
        System.out.println("  --pattern <String>       - regex pattern (use \\\\t for tab)");
        System.out.println("                             (specify once or twice)");
        System.out.println("  --column <int>           - column index for pattern search (first column is index 1)");
        System.out.println("                             (specify once or twice, set to -1 to ignore columns)");
        System.out.println("  --sep <String>           - character separating columns (use \\\\t for tab) [default \\t]");
        System.out.println("  --ignore <String>        - character marking comment line [default #]");
        System.out.println("  --output <File>          - where to save output [default stdout]");
        System.out.println(" ");
    }

    public ThesaurusCountPatterns(String[] newargs) {
        if (newargs == null || newargs.length == 0) {
            printCountFilteredVariantsHelp();
            return;
        }

        isok = parseCountFilteredParameters(newargs);
    }

    private boolean parseCountFilteredParameters(String[] args) {

        OptionParser prs = new OptionParser();

        // options for running blat        
        prs.accepts("pattern").withRequiredArg().ofType(String.class);
        prs.accepts("column").withRequiredArg().ofType(Integer.class);
        prs.accepts("sep").withRequiredArg().ofType(String.class);
        prs.accepts("ignore").withRequiredArg().ofType(String.class);
        prs.accepts("input").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);


        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // extract input/output settings
        if (options.has("input")) {
            inputfile = (File) options.valueOf("input");
            if (!inputfile.exists() || !inputfile.canRead()) {
                System.out.println("Cannot read input file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required parameter: input");
            return false;
        }
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }

        // extract patterns and matching column ids
        if (options.has("pattern")) {
            ArrayList<Object> pstrings = new ArrayList<Object>(options.valuesOf("pattern"));
            patternstringA = (String) pstrings.get(0);
            if (pstrings.size() > 1) {
                patternstringB = (String) pstrings.get(1);
            }
        }
        if (options.has("column")) {
            ArrayList<Object> cols = new ArrayList<Object>(options.valuesOf("column"));
            columnA = -1 + ((Integer) cols.get(0)).intValue();
            if (cols.size() > 1) {
                columnB = -1 + ((Integer) cols.get(1)).intValue();
            }
        }

        // try to create Patterns from the pattern strings
        try {
            pA = Pattern.compile(patternstringA);
        } catch (Exception ex) {
            System.out.println("Error while constructing pattern matcher with regex: " + patternstringA);
            System.out.println(ex.getMessage());
            return false;
        }
        try {
            pB = Pattern.compile(patternstringB);
        } catch (Exception ex) {
            System.out.println("Error while constructing pattern matcher with regex: " + patternstringB);
            System.out.println(ex.getMessage());
            return false;
        }

        // extract parameters tweaking input file parsing
        if (options.has("sep")) {
            sep = (String) options.valueOf("sep");
        }

        if (options.has("ignore")) {
            ignore = (String) options.valueOf("ignore");
        }

        return true;
    }

    private boolean matchesPattern(String row, int column, Pattern p, String sep) {
        // perhaps do the matching on the entire row
        if (column < 0) {
            return (p.matcher(row).find());
        }

        // if reached here, need to split the row into columns
        String[] tokens = row.split(sep);
        if (tokens.length <= column) {
            return false;
        } else {
            return (p.matcher(tokens[column]).find());
        }
    }

    @Override
    public void run() {
        if (!isok) {
            return;
        }

        OutputStream outstream = null;
        try {
            // scan the vcf file looking for hits
            BufferedReader vcfreader = BufferedReaderMaker.makeBufferedReader(inputfile);
            // keep track of variants that have the filter and that match the pattern
            int A0B0 = 0, A0B1 = 0, A1B0 = 0, A1B1 = 0;
            String ss;
            while ((ss = vcfreader.readLine()) != null) {
                if (ss.startsWith(ignore)) {
                    // skip this line
                } else {

                    if (matchesPattern(ss, columnA, pA, sep)) {
                        if (matchesPattern(ss, columnB, pB, sep)) {
                            A1B1++;
                        } else {
                            A1B0++;
                        }
                    } else {
                        if (matchesPattern(ss, columnB, pB, sep)) {
                            A0B1++;
                        } else {
                            A0B0++;
                        }
                    }
                }
            }
            vcfreader.close();

            // prepare the output
            StringBuilder sb = new StringBuilder();
            sb.append("##input=").append(inputfile.getAbsolutePath()).append("\n");
            sb.append("##patternA=").append(patternstringA).append("\n");
            sb.append("##columnA=").append(1 + columnA).append("\n");
            sb.append("##patternB=").append(patternstringB).append("\n");
            sb.append("##columnB=").append(1 + columnB).append("\n");
            sb.append("patternA\tpatternB\tcount\n");
            sb.append("0\t0\t").append(A0B0).append("\n");
            sb.append("0\t1\t").append(A0B1).append("\n");
            sb.append("1\t0\t").append(A1B0).append("\n");
            sb.append("1\t1\t").append(A1B1).append("\n");

            outstream = OutputStreamMaker.makeOutputStream(output);
            outstream.write(sb.toString().getBytes());
            outstream.close();

        } catch (Exception ex) {
            System.out.println("Something went wront: " + ex.getMessage());
        } finally {
            try {
                outstream.close();
            } catch (Exception ex) {
                System.out.println("Something went wrong: " + ex.getMessage());
            }
        }

    }
}
