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
import java.util.HashMap;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import thesaurus.util.ThesaurusIO;

/**
 * Create an annotated vtf file, i.e. a vtf file with additional dbSNP ids.
 *
 * @author tkonopka
 */
public class ThesaurusAnnotateVtf extends ThesaurusMapTool {

    private String output = "stdout";
    private File dbfile = null;
    private File vtffile = null;

    private void printAnnotateVtfHelp() {
        System.out.println("GeneticThesaurus annotatevtf: adds ids (eg. dbSNP) into a vtf file");
        System.out.println();
        System.out.println("Usage: java -jar GeneticThesaurus.jar annotatevtf ");
        System.out.println();
        System.out.println("Core options:");
        ThesaurusIO.printHelpItem("--vtf <File>", "input vtf file");
        ThesaurusIO.printHelpItem("--database <File>", "database file, e.g. dbSNP");
        ThesaurusIO.printHelpItem("--output <String>", "output file");
        System.out.println();
    }

    /**
     *
     * @param args
     * @return
     */
    private boolean parseAnnotateVtfParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // input and output 
        prs.accepts("vtf").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);
        prs.accepts("database").withRequiredArg().ofType(File.class);

        // now use OptionSet to parse the command line
        OptionSet options;

        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        if (options.has("database")) {
            dbfile = (File) options.valueOf("database");
            if (!dbfile.exists() || !dbfile.canRead()) {
                System.out.println("Cannot read database file, or file does not exist");
                return false;
            }
        }

        // get prefix for output tracks
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        } else {
            System.out.println("Missing required argument --output");
            return false;
        }

        // get input vcf file
        if (options.has("vtf")) {
            vtffile = (File) options.valueOf("vtf");
            if (!vtffile.exists() || !vtffile.canRead()) {
                System.out.println("Cannot read vtf file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument --vtf");
            return false;
        }

        return true;
    }

    public ThesaurusAnnotateVtf(String[] args) {
        if (args.length == 0) {
            printAnnotateVtfHelp();
            return;
        }
        super.loadDefaults();
        super.setOk(parseAnnotateVtfParameters(args));
    }

    @Override
    void runTool() {

        // make a hashmap that will hold "chr" position and dbsnp
        IDMap posids = new IDMap();

        // read the input vtf file
        try {
            getPosFromVtf(vtffile, posids);
        } catch (Exception ex) {
            System.out.println("Error reading/parsing vtf file: " + ex.getMessage());
            return;
        }


        // scan the database file and add ids into the map
        boolean readok = false;
        try {
            readok = learnIdsFromDatabase(dbfile, posids);
        } catch (Exception ex) {
            System.out.println("Error reading database: " + ex.getMessage());
            return;
        }

        if (!readok) {
            return;
        }

        // output a new vtf file
        try {
            writeAnnotatedVtf(vtffile, output, posids);
        } catch (Exception ex) {
            System.out.println("Error outputing annotated vtf file: " + ex.getMessage());
            return;
        }

    }

    /**
     * Scan a vcf dbfile (eg. dbSNP) and record the ids associated with a locus
     *
     * @param dbfile
     * @param posmap
     * @throws IOException
     */
    private boolean learnIdsFromDatabase(File dbfile, IDMap posmap) throws IOException {
        
        // open the database for reading
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(dbfile);
        
        // process each line in the database one by one
        String s;        
        while ((s = br.readLine()) != null) {           
            if (!s.startsWith("#")) {                
                String[] tokens = s.split("\t", 4);
                try {
                    posmap.set(tokens[0], Integer.parseInt(tokens[1]), tokens[2]);
                } catch (Exception ex) {
                    System.out.println("error processing databse entry:\n" + s);
                    return false;
                }
            }
        }

        br.close();

        // signal that learning finished ok
        return true;
    }

    /**
     * Read a vtf file and record all loci into a hashmap
     *
     * @param vtffile
     * @param posmap
     * @throws IOException
     */
    private void getPosFromVtf(File vtffile, IDMap posmap) throws IOException {
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(vtffile);
        String s;

        // read one line at a time, enter all positions encountered into a map object
        while ((s = br.readLine()) != null) {
            if (!s.startsWith("#")) {
                String[] tokens = s.split("\t");
                for (int i = 0; i < tokens.length; i++) {
                    posmap.put(tokens[i]);
                }
            }
        }

        br.close();
    }

    /**
     * Write an annotated vtf file. This works by reading a vtf file and
     * outputing a related version with all the same positions, but with the ids
     * stored withing the IDmap
     *
     * @param vtffile
     * @param output
     * @param posmap
     * @throws IOException
     */
    private void writeAnnotatedVtf(File vtffile, String output, IDMap posmap) throws IOException {

        OutputStream os = OutputStreamMaker.makeOutputStream(output);
        BufferedReader br = BufferedReaderMaker.makeBufferedReader(vtffile);
        String s;

        while ((s = br.readLine()) != null) {
            if (s.startsWith("#")) {
                // copy comments from old vtf into new vtf
                os.write((s + "\n").getBytes());
            } else {
                // look up ids for each element and output them in the same order 
                String[] tokens = s.split("\t");
                StringBuilder sb = new StringBuilder();
                sb.append(tokens[0]).append(":").append(posmap.get(tokens[0]));
                for (int i = 1; i < tokens.length; i++) {
                    sb.append("\t").append(tokens[i]).append(":").append(posmap.get(tokens[i]));
                }
                sb.append("\n");
                os.write(sb.toString().getBytes());
            }
        }

        os.close();
        br.close();
    }

    /**
     * A class holds genomic positions and a string. In this program the string
     * is used to hold a label/id from dbSNP.
     *
     */
    class IDMap {

        HashMap<String, HashMap<Integer, String>> idmap = new HashMap<String, HashMap<Integer, String>>();

        // adds a new element into the hashmap
        public void put(String poslabel) {
            String[] tokens = poslabel.split(":");
            String chr = tokens[0];
            int pos = Integer.parseInt(tokens[1]);

            if (!idmap.containsKey(chr)) {
                idmap.put(chr, new HashMap<Integer, String>());
            }

            HashMap<Integer, String> nowchrmap = idmap.get(chr);
            if (!nowchrmap.containsKey(pos)) {
                nowchrmap.put(pos, ".");
            }
        }

        // adds a new element into the hashmap
        public void set(String chr, int position, String id) {
            if (idmap.containsKey(chr)) {
                HashMap<Integer, String> nowchrmap = idmap.get(chr);
                if (nowchrmap.containsKey(position)) {                    
                    nowchrmap.put(position, id);
                }
            }
        }        

        public String get(String chr, int position) {
            if (idmap.containsKey(chr)) {
                HashMap<Integer, String> nowchrmap = idmap.get(chr);
                if (nowchrmap.containsKey(position)) {
                    return nowchrmap.get(position);
                }
            }
            return null;
        }

        public String get(String poslabel) {
            String[] tokens = poslabel.split(":");
            return get(tokens[0], Integer.parseInt(tokens[1]));
        }
    }
}
