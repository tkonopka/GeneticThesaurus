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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.Preferences;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;

/**
 * Small utility that allows the user to store some settings into the
 * Preferences object associated with the GenThesaurus class.
 *
 *
 * @author tkonopka
 */
class ThesaurusConfigure implements Runnable {

    private Preferences prefs = Preferences.userNodeForPackage(GeneticThesaurus.class);
    private boolean reset = false;
    private boolean show = false;

    private void printConfigureHelp() {
        System.out.println("GenThesaurus configure:");
        System.out.println();
        System.out.println("To view or reset the defaults:");
        System.out.println();
        System.out.println("  --reset                  - reset all values to defaults");
        System.out.println("  --show                   - print all currently set values");
        System.out.println();
        System.out.println("To change the configuration, run with one or more of the following:");
        System.out.println();
        System.out.println("  --genome <File>          - genome fasta");
        //System.out.println("  --blatpath <File>        - path to blat executable");
        //System.out.println("  --blatoptions <File>     - text file holding blat options");
        System.out.println("                             (is read immediately, not at runtime)");
        System.out.println("  --readlen <int>          - read length");
        System.out.println("  --divisor <int>          - generate reads every so many bases");
        System.out.println("  --offset <int>           - when using divisor > 1");
        System.out.println();
    }

    private boolean parseConfigureParameters(String[] args) {

        OptionParser prs = new OptionParser();

        prs.accepts("reset");
        prs.accepts("show");

        // options for running blat
        //prs.accepts("blatpath").withRequiredArg().ofType(String.class);
        //prs.accepts("blatoptions").withRequiredArg().ofType(File.class);
        //prs.accepts("keeppsl").withRequiredArg().ofType(Boolean.class);

        // genome fasta file
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // read generation settings
        prs.accepts("readlen").withRequiredArg().ofType(Integer.class);
        prs.accepts("divisor").withRequiredArg().ofType(Integer.class);
        prs.accepts("offset").withRequiredArg().ofType(Integer.class);

        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // check if user wants to just see/reset the default values
        if (options.has("reset") || options.has("show")) {
            if (options.has("reset")) {
                reset = true;
            }
            if (options.has("show")) {
                show = true;
            }
            return true;
        }


        // if not just seeing or resetting, change the values
        if (options.has("genome")) {
            File genome = (File) options.valueOf("genome");
            if (!genome.exists() || !genome.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist");
                return false;
            } else {
                prefs.put("genome", genome.getAbsolutePath());
            }
        }

        if (options.has("blatpath")) {
            prefs.put("blatpath", (String) options.valueOf("blatpath"));
        }

        // read options for blat from a file
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
                prefs.put("blatoptions", optionsSB.toString());
            } catch (IOException ex) {
                Logger.getLogger(ThesaurusConfigure.class.getName()).log(Level.SEVERE, null, ex);
            } finally {
                try {
                    optionsBR.close();
                } catch (IOException ex) {
                    Logger.getLogger(ThesaurusConfigure.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        if (options.has("readlen")) {
            prefs.put("readlen", "" + (Integer) options.valueOf("readlen"));
        }
        if (options.has("divisor")) {
            prefs.put("divisor", "" + (Integer) options.valueOf("divisor"));
        }
        if (options.has("offset")) {
            prefs.put("offset", "" + (Integer) options.valueOf("offset"));
        }

        if (options.has("keeppsl")) {
            prefs.putBoolean("keeppsl", (Boolean) options.valueOf("keeppsl"));
        }

        return true;
    }

    ThesaurusConfigure(String[] newargs) {
        if (newargs.length == 0) {
            printConfigureHelp();
            return;
        }
        parseConfigureParameters(newargs);
    }

    @Override
    public void run() {
        if (reset) {
            resetValues();
        }
        if (show) {
            showValues();
        }
    }

    private void resetValues() {
        prefs.put("genome", GeneticThesaurus.DEFAULT_GENOME);
        //prefs.put("blatpath", GenThesaurus.DEFAULT_BLATPATH);
        //prefs.put("blatoptions", GenThesaurus.DEFAULT_BLATOPTIONS);
        //prefs.putBoolean("keeppsl", GenThesaurus.DEFAULT_KEEPPSL);
        prefs.putInt("readlen", GeneticThesaurus.DEFAULT_READLEN);
        prefs.putInt("divisor", GeneticThesaurus.DEFAULT_DIVISOR);
        prefs.putInt("offset", GeneticThesaurus.DEFAULT_OFFSET);
        prefs.putInt("penalty", GeneticThesaurus.DEFAULT_PENALTY);
    }

    private void showValues() {
        System.out.println("Genetic Thesaurus v" + GeneticThesaurus.getVersion());
        System.out.println("genome:  \t" + prefs.get("genome", GeneticThesaurus.DEFAULT_GENOME));
        //System.out.println("blatpath:\t" + prefs.get("blatpath", GenThesaurus.DEFAULT_BLATPATH));
        //System.out.println("blatoptions:\t" + prefs.get("blatoptions", GenThesaurus.DEFAULT_BLATOPTIONS));
        //System.out.println("keeppsl:\t" + prefs.getBoolean("keeppsl", GenThesaurus.DEFAULT_KEEPPSL));
        System.out.println("readlen:  \t" + prefs.getInt("readlen", GeneticThesaurus.DEFAULT_READLEN));
        System.out.println("divisor: \t" + prefs.getInt("divisor", GeneticThesaurus.DEFAULT_DIVISOR));
        System.out.println("offset:  \t" + prefs.getInt("offset", GeneticThesaurus.DEFAULT_OFFSET));
        System.out.println("penalty:  \t" + prefs.getInt("offset", GeneticThesaurus.DEFAULT_PENALTY));
    }
}
