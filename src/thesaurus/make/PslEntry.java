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

import java.util.ArrayList;
import java.util.Comparator;

/**
 * A class storing a parsed psl entry from BLAT.
 *
 * This code is adapted from code from the Bamformatics project
 *
 *
 * @author tkonopka
 * @deprecated no longer needed because ThesaurusBlat is deprecated
 *
 */
@Deprecated
class PslEntry {

    int match = 0;
    int mismatch = 0;
    int repmatch = 0;
    int Ns = 0;
    int Qgapcount = 0;
    int Qgapbases = 0;
    int Tgapcount = 0;
    int Tgapbases = 0;
    char strand = '+';
    String Qname = null;
    int Qsize = 0;
    int Qstart = 0;
    int Qend = 0;
    String Tname = null;
    int Tsize = 0;
    int Tstart = 0;
    int Tend = 0;
    int blockcount = 0;
    ArrayList<Integer> blockSizes = new ArrayList<Integer>(8);
    ArrayList<Integer> qStarts = new ArrayList<Integer>(8);
    ArrayList<Integer> tStarts = new ArrayList<Integer>(8);
    // not strictly in the psl format, but can be included here
    byte[] sequence;
    byte[] qualities;

    /**
     * generic constructor that doesn't change any of the default values
     */
    public PslEntry() {
    }

    /**
     *
     * @param entry
     *
     * A long string, tab separated, output from blat
     *
     */
    public PslEntry(String entry) {

        // parse the long string into tab separated tokens
        String[] tokens = entry.split("\t");
        if (tokens.length < 21) {
            System.out.println("Invalid psl entry. Less than 21 columns");
            return;
        }

        // convert the tokens into values for the PslEntry
        match = Integer.parseInt(tokens[0]);
        mismatch = Integer.parseInt(tokens[1]);
        repmatch = Integer.parseInt(tokens[2]);
        Ns = Integer.parseInt(tokens[3]);
        Qgapcount = Integer.parseInt(tokens[4]);
        Qgapbases = Integer.parseInt(tokens[5]);
        Tgapcount = Integer.parseInt(tokens[6]);
        Tgapbases = Integer.parseInt(tokens[7]);
        strand = tokens[8].charAt(0);
        Qname = tokens[9];
        Qsize = Integer.parseInt(tokens[10]);
        Qstart = Integer.parseInt(tokens[11]);
        Qend = Integer.parseInt(tokens[12]);
        Tname = tokens[13];
        Tsize = Integer.parseInt(tokens[14]);
        Tstart = Integer.parseInt(tokens[15]);
        Tend = Integer.parseInt(tokens[16]);
        blockcount = Integer.parseInt(tokens[17]);

        String[] bstokens = tokens[18].split(",");
        String[] qstokens = tokens[19].split(",");
        String[] tstokens = tokens[20].split(",");
        for (int i = 0; i < bstokens.length; i++) {
            blockSizes.add(Integer.parseInt(bstokens[i]));
            qStarts.add(Integer.parseInt(qstokens[i]));
            tStarts.add(Integer.parseInt(tstokens[i]));
        }

        // trim the ending of the query name if it contains /1 or /2
        int qnamelen = Qname.length();
        if (Qname.substring(qnamelen - 2, qnamelen).equals("/1") || Qname.substring(qnamelen - 2, qnamelen).equals("/2")) {
            Qname = Qname.substring(0, qnamelen - 2);
        }

    }

    public String getEntryString() {

        StringBuilder sb = new StringBuilder(256);
        sb.append(match).append("\t");
        sb.append(mismatch).append("\t");
        sb.append(repmatch).append("\t");
        sb.append(Ns).append("\t");
        sb.append(Qgapcount).append("\t");
        sb.append(Qgapbases).append("\t");
        sb.append(Tgapcount).append("\t");
        sb.append(Tgapbases).append("\t");
        sb.append(strand).append("\t");
        sb.append(Qname).append("\t");
        sb.append(Qsize).append("\t");
        sb.append(Qstart).append("\t");
        sb.append(Qend).append("\t");
        sb.append(Tname).append("\t");
        sb.append(Tsize).append("\t");
        sb.append(Tstart).append("\t");
        sb.append(Tend).append("\t");
        sb.append(blockcount).append("\t");
        for (int i = 0; i < blockSizes.size(); i++) {
            sb.append(blockSizes.get(i)).append(",");
        }
        sb.append("\t");
        for (int i = 0; i < qStarts.size(); i++) {
            sb.append(qStarts.get(i)).append(",");
        }
        sb.append("\t");
        for (int i = 0; i < tStarts.size(); i++) {
            sb.append(tStarts.get(i)).append(",");
        }

        return sb.toString();
    }

    public String makeCigar() {
        StringBuilder cigar = new StringBuilder(8);

        // keep track of the current position along the query
        int querypos = this.qStarts.get(0);
        //int transcriptpos = entry.tStarts.get(0);
        int numblocks = this.blockSizes.size();

        // check if there is initial trimming       
        if (querypos > 0) {
            cigar.append(querypos).append("S");
        }

        // append the middle blocks
        for (int i = 1; i < numblocks; i++) {
            int nowqpos = this.qStarts.get(i);
            int nowtpos = this.tStarts.get(i);
            int prevtpos = this.tStarts.get(i - 1);
            int prevqpos = this.qStarts.get(i - 1);
            int prevblock = this.blockSizes.get(i - 1);

            int nowinsert = nowqpos - prevblock - prevqpos;
            int nowdelete = nowtpos - prevtpos - prevblock;

            // append the match cigar, which is the same always
            cigar.append(prevblock).append("M");

            // then append the insert and delete cigars
            if (nowinsert > 0) {
                // it's an insertion                
                cigar.append(nowinsert).append("I");
            }
            if (nowdelete > 0) {
                // there is a deletion                                                
                cigar.append(nowdelete).append("D");
            }
            querypos = nowqpos;
        }

        // append the last matching subsequence cigar
        cigar.append(this.blockSizes.get(numblocks - 1)).append("M");
        querypos += this.blockSizes.get(numblocks - 1);

        // check if there is final position clipping
        if (querypos < this.Qsize) {
            cigar.append(this.Qsize - querypos).append("S");
        }

        return cigar.toString();
    }
}

/**
 * comparator for Psl, used when sorting Psl entries, in decreasing score order.
 *
 *
 * @author tomasz
 */
class PslComparator implements Comparator {

    @Override
    public int compare(Object o1, Object o2) {
        int match1 = ((PslEntry) o1).match;
        int match2 = ((PslEntry) o2).match;
        if (match1 < match2) {
            return 1;
        } else if (match1 > match2) {
            return -1;
        } else {
            return 0;
        }
    }
}
