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

/**
 * Class that will hold locations of all subsequences. Class scans a long
 * sequence, considers all subsequences of given length. Computes an integer for
 * each subsequence as the sum of the bytes, then records location in a hashmap
 * using the sum of bytes as the key.
 *
 * The index stores the 0-based coordinates.
 *
 * @author tkonopka
 */
class ThesaurusEarIndex {

    // record the reference sequence string identifier and sequence
    private final String seqname;
    private final byte[] seq;
    // will also hold an index type object that will help with alignment
    //private final ArrayList<Integer>[] index;
    private final int[][] index;
    // arrays will store values for each base at each position
    // can be used later to build keys/hashcodes for 
    private final int[] Avals;
    private final int[] Tvals;
    private final int[] Cvals;
    private final int[] Gvals;
    private final int earlength;

    /**
     * get a hashcode for a subsequence [from, to)
     *
     * @param seq
     *
     * possibly long sequence
     *
     * @param from
     *
     * index of first base (included)
     * @param to
     *
     * index of last base (not included)
     *
     * @return
     *
     * a number that can be used as a hashcode
     *
     */
    public final int getCode(byte[] seq, int from, int to) {

        if (to - from > earlength) {
            return -1;
        }

        int ans = 0;
        for (int i = from; i < to; i++) {
            switch (seq[i]) {
                case 'A':
                    ans += Avals[i - from];
                    break;
                case 'T':
                    ans += Tvals[i - from];
                    break;
                case 'G':
                    ans += Gvals[i - from];
                    break;
                case 'C':
                    ans += Cvals[i - from];
                    break;
                default:
                    return -1;
            }
        }

        return ans;
    }

    private void buildATCGcodes(int earlength) {
        int Alarge = 0 << (2 * earlength - 2);
        int Tlarge = 1 << (2 * earlength - 2);
        int Clarge = 2 << (2 * earlength - 2);
        int Glarge = 3 << (2 * earlength - 2);

        for (int j = 0; j < earlength; j++) {
            Avals[j] = Alarge >> (j * 2);
            Tvals[j] = Tlarge >> (j * 2);
            Cvals[j] = Clarge >> (j * 2);
            Gvals[j] = Glarge >> (j * 2);
        }
    }

    /**
     * constructs an index for a sequence.
     *
     * @param seq
     * @param earlength
     *
     */
    public ThesaurusEarIndex(String seqname, byte[] seq, int earlength) {

        this.seqname = seqname;
        this.earlength = earlength;
        this.seq = seq;

        // build a data structure for the index
        int seqlen = seq.length;
        int indexsize = 4;
        for (int i = 1; i < earlength; i++) {
            indexsize = indexsize * 4;
        }

        ArrayList<Integer>[] tempindex = new ArrayList[indexsize];
        for (int i = 0; i < indexsize; i++) {
            tempindex[i] = new ArrayList<Integer>(seqlen / indexsize);
        }

        // make codes for all bases at all positions
        Avals = new int[earlength];
        Tvals = new int[earlength];
        Cvals = new int[earlength];
        Gvals = new int[earlength];
        buildATCGcodes(earlength);

        // now calculate the keys for all subsequences
        for (int i = 0; i < seqlen - earlength; i++) {
            int nowkey = getCode(seq, i, i + earlength);

            // add ear into the index (but not if ear contains N characters)
            if (nowkey > -1) {
                tempindex[nowkey].add(i);
            }
        }

        // now copy the information into primitive arrays.
        // They will be looked at multiple times, so perhaps this is worthwhile
        index = new int[indexsize][];
        for (int i = 0; i < indexsize; i++) {
            ArrayList<Integer> temp = tempindex[i];
            int tempsize = temp.size();
            index[i] = new int[tempsize];
            for (int k = 0; k < tempsize; k++) {
                index[i][k] = temp.get(k);
            }
            tempindex[i].clear();
        }

        // suggest to clean up now, because big arrays have been created and are no longer needed
        tempindex = null;
        System.gc();
    }

    public int getEarlength() {
        return earlength;
    }

    public String getSeqname() {
        return seqname;
    }

    /**
     * Prints out the status of the index
     */
    public void summarizeIndex() {
        for (int j = 0; j < index.length; j++) {
            int[] aa = index[j];
            System.out.print("Index for key " + j + ":\t");
            if (aa.length < 10) {
                for (int i = 0; i < aa.length; i++) {
                    System.out.print(aa[i] + "\t");
                }
                System.out.println();
            } else {
                for (int i = 0; i < 5; i++) {
                    System.out.print(aa[i] + "\t");
                }
                System.out.print("...\t");
                for (int i = aa.length - 5; i < aa.length; i++) {
                    System.out.print(aa[i] + "\t");
                }
                System.out.println();
            }
        }
    }

    /**
     * Prints out mapping between letter codes (e.g. AATC) and hashcodes
     *
     */
    public void printAllCodes() {
        byte[] alphabet = new byte[4];
        alphabet[0] = 'A';
        alphabet[1] = 'T';
        alphabet[2] = 'C';
        alphabet[3] = 'G';
        showCode("", alphabet);
    }

    /**
     * recursive function builds all possible substrings and outputs their codes
     *
     * @param s
     * @param alphabet
     */
    private void showCode(String s, byte[] alphabet) {
        // exit condition
        if (s.length() == earlength) {
            System.out.println(s + "\t" + getCode(s.getBytes(), 0, earlength));
            return;
        }
        for (int i = 0; i < alphabet.length; i++) {
            showCode(s + (char) alphabet[i], alphabet);
        }
    }

    /**
     * functions tries to match a short read sequence against the indexed
     * sequence. This looks at the first and last positions in the read and then
     * looks at the intermediate bases to check for agreement
     *
     * @param newseq
     *
     * a short read
     *
     * @return
     *
     * a list with integers following schema: alignment start position, mapping
     * quality, start position, mapping quality, etc.
     */
    public ArrayList<Integer> alignSequence(byte[] readseq, int maxmismatches) {

        int readlen = readseq.length;
        ArrayList<Integer> ans = new ArrayList<Integer>();

        if (readlen < 2 * earlength) {
            return ans;
        }

        // look at the start/end ear sequences
        int startcode = getCode(readseq, 0, earlength);
        int endcode = getCode(readseq, readlen - earlength, readlen);
        if (startcode < 0 || endcode < 0) {
            return ans;
        }
        int[] startpos = index[startcode];
        int[] endpos = index[endcode];
        int numstarts = startpos.length;
        int numends = endpos.length;
        if (numstarts < 1 || numends < 1) {
            return ans;
        }

        // walk through the positions, keep track of current indexes
        int starti = 0, endi = 0;
        int nowstart = -1;
        int nowend = -1;

        // consider all the possible start positions
        while (starti < numstarts) {

            nowstart = startpos[starti];
            //System.out.println("Looking at start " + nowstart);

            if (endi >= numends) {
                break;
            }

            // consider all the possible end positions
            int magicend = nowstart + readlen - earlength;
            while (nowend <= magicend && endi < numends) {
                nowend = endpos[endi];
                //System.out.println("Looking at end " + nowend);
                if (nowend == magicend) {
                    // this start and end combination matches readlength                    
                    int mapqual = getMappingQuality(readseq, nowstart, maxmismatches);
                    if (mapqual > 0) {
                        ans.add(nowstart);
                        ans.add(mapqual);
                    }
                }

                // look at the next index
                if (nowend <= magicend) {
                    endi++;
                }
            }

            // look at next start index
            starti++;
        }

        return ans;

    }

    /**
     * compares a read sequence to a fragment of the reference sequence
     *
     * @param readseq
     *
     * short sequence
     *
     * @param from
     *
     * starting index of reference sequence
     *
     * @param maxmismatches
     *
     * maximum number of mismatches to allow
     *
     * @return
     *
     * the number of bases of the read that match exactly against the reference
     *
     */
    public int getMappingQuality(byte[] readseq, int from, int maxmismatches) {
        // count the number of mismatches
        int nowmm = 0;
        int readlen = readseq.length;
        // loop over the readsequence bases and compare to the reference
        for (int i = 0; i < readlen; i++) {
            if (readseq[i] != seq[from + i]) {
                nowmm++;
                if (nowmm > maxmismatches) {
                    return -1;
                }
            }
        }
        return readlen - nowmm;
    }
}
