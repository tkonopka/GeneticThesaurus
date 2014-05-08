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
package thesaurus.util;

import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Class basically copied from Bamformatics project
 *
 *
 * Basically, a wrapper for println. The logged string is written to a
 * previously specified logstream. The logged output also has some special
 * formating including the date.
 *
 *
 * @author tomasz
 */
public class ThesaurusLog {

    private final static SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private final PrintStream logstream;
    private boolean verbose = false;

    public ThesaurusLog(PrintStream logstream) {
        this.logstream = logstream;
    }

    public ThesaurusLog() {
        this.logstream = System.out;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    /**
     * output comment, but only if verbose level for this log has been set to
     * true
     *
     * @param s
     *
     * comment string to display in the log
     */
    public void log(String s) {
        if (verbose) {
            logstream.println("[Thesaurus][" + sdf.format(new Date()) + "] " + s);
            logstream.flush();
        }
    }

    /**
     * output comment
     *
     * @param verbose
     *
     * verbose level for this comment. will temporarily override verbose level
     * for the class.
     *
     * @param s
     *
     * comment string to display in the log
     *
     */
    public void log(boolean verbose, String s) {
        if (verbose) {
            logstream.println("[Thesaurus][" + sdf.format(new Date()) + "] " + s);
            logstream.flush();
        }
    }
}
