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
package thesaurus.make;

import thesaurus.GeneticThesaurus;
import java.io.File;
import java.util.prefs.Preferences;

/**
 *
 * This is just a means to keep track of settings that are saved in users'
 * preferences.
 *
 * @author tkonopka
 */
public abstract class ThesaurusMapTool implements Runnable {

    File genome;
    int divisor;
    int offset;
    int readlen;
    int penalty;
    boolean keeppsl;
    String blatpath;
    String blatoptions;
    // used in the runnable
    private boolean ok = false;

    public boolean isOk() {
        return ok;
    }

    public void setOk(boolean ok) {
        this.ok = ok;
    }

    /**
     * load the default values for the blat options
     */
    void loadDefaults() {
        Preferences prefs = Preferences.userNodeForPackage(GeneticThesaurus.class);
        genome = new File(prefs.get("genome", GeneticThesaurus.DEFAULT_GENOME));
        readlen = prefs.getInt("readlen", GeneticThesaurus.DEFAULT_READLEN);
        divisor = prefs.getInt("divisor", GeneticThesaurus.DEFAULT_DIVISOR);
        penalty = prefs.getInt("penalty", GeneticThesaurus.DEFAULT_PENALTY);
        offset = prefs.getInt("offset", GeneticThesaurus.DEFAULT_OFFSET);
        readlen = prefs.getInt("readlen", GeneticThesaurus.DEFAULT_READLEN);        
    }

    @Override
    public void run() {
        if (!ok) {
            return;
        }
        runTool();
    }

    abstract void runTool();
}
