/*
 * Copyright 2013-2015 Tomasz Konopka.
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.genome.GenomeInfo;

/**
 * A "buffered reader" for thesaurus tables. Assumes thesaurus file is sorted,
 * allows users to query for entries overlapping loci.
 *
 * Similar to BamRegionsMap
 *
 * all coordinates in the thesaurus file and in this map are 1-based.
 *
 * @author tkonopka
 */
class ThesaurusRegionsMap {

    // information about the genome
    private final GenomeInfo ginfo;
    // the buffered reader providing access to the reader
    private final BufferedReader thesaurusreader;
    // interval representation of the thesaurus regions
    private final HashMap<Integer, ArrayList<ThesaurusEntry>> regions;
    // bookkeeping variables
    private String thesline = null;
    private ThesaurusEntry thesentry = null;
    private String curchr = null;
    private final int bufferlength;
    private int uptoposition = 0;

    /**
     * Creates and preps a reader of thesaurus intervals.
     *
     * This can called with a null argument. This should make the regions map
     * always to return no intervals for any variant.
     *
     * @param thesfile
     *
     * file with thesaurus entries
     *
     * @param ginfo
     *
     * Information about the genome (chromosome names, order, lengths)
     *
     * @param bufferlength
     *
     * when reading through the thesaurus entries, the reader will return
     * entries covering the requested position and entries that begin/end within
     * this distance for the requested position.
     *
     *
     * @throws IOException
     *
     */
    public ThesaurusRegionsMap(File thesfile, GenomeInfo ginfo, int bufferlength) throws IOException {
        this.ginfo = ginfo;
        this.bufferlength = bufferlength;
        if (thesfile != null) {
            thesaurusreader = BufferedReaderMaker.makeBufferedReader(thesfile);
            startLoad();
        } else {
            thesaurusreader = null;
        }
        regions = new HashMap<Integer, ArrayList<ThesaurusEntry>>(128);
    }

    private void startLoad() throws IOException {

        // read from the thesaurus buffer and skip the comment lines
        thesline = thesaurusreader.readLine();
        while (thesline != null && thesline.startsWith("#")) {
            thesline = thesaurusreader.readLine();
        }
        if (thesline == null) {
            return;
        }

        // at this stage, the thesline should be the header
        if (!thesline.startsWith("Align.chr")) {
            System.out.println("This is not a thesaurus file. Table must start with column Align.chr.");
            return;
        }
        thesline = thesaurusreader.readLine();
        thesentry = new ThesaurusEntry(thesline, ginfo);
        if (thesentry.alignChrIndex >= 0) {
            curchr = ginfo.getChrName(thesentry.alignChrIndex);
        }
    }

    /**
     * prints to system.out
     *
     * print a summary of the current status of the regions map.
     *
     */
    public void printCurrent() {
        System.out.println("Current uptoposition: " + uptoposition);
        System.out.println("Current regions map size: " + getNumEntries());

        System.out.println("Current status:");
        for (Map.Entry<Integer, ArrayList<ThesaurusEntry>> entry : regions.entrySet()) {
            ArrayList<ThesaurusEntry> map1 = entry.getValue();
            for (int i = 0; i < map1.size(); i++) {
                System.out.print(map1.get(i).toString());
            }
        }

    }

    private int getNumEntries() {
        int count = 0;
        for (Map.Entry<Integer, ArrayList<ThesaurusEntry>> entry : regions.entrySet()) {
            ArrayList<ThesaurusEntry> map1 = entry.getValue();
            count += map1.size();
        }
        return count;
    }

    /**
     * Provided the object was loadedUpto() a position, this function can return
     * an array of thesaurus entries that overlap with the given interval.
     *
     *
     * @param minpos
     *
     * starting position of interval
     *
     * @param maxpos
     *
     * maximum position of interval
     *
     *
     * @return
     *
     * array of entries overlapping with the interval, or null if no entries are
     * loaded.
     *
     */
    private ThesaurusEntry[] getEntries(int minpos, int maxpos) {

        // simple exit if no regions available
        if (isEmpty()) {
            return null;
        }

        // figure out how many regions are current stored in memory
        int count = getNumEntries();

        // make an array of ThesaurusEntries here
        ArrayList<ThesaurusEntry> temp = new ArrayList<ThesaurusEntry>(count);
        for (Map.Entry<Integer, ArrayList<ThesaurusEntry>> entry : regions.entrySet()) {
            ArrayList<ThesaurusEntry> map1 = entry.getValue();
            for (int i = 0; i < map1.size(); i++) {
                ThesaurusEntry hereentry = map1.get(i);
                if (ThesaurusEntry.overlap(minpos, maxpos, hereentry.alignStart, hereentry.alignEnd)) {
                    temp.add(hereentry);
                }
            }
        }

        if (temp.isEmpty()) {
            return null;
        }
        ThesaurusEntry[] ans = new ThesaurusEntry[temp.size()];
        for (int i = 0; i < temp.size(); i++) {
            ans[i] = temp.get(i);
        }
        return ans;
    }

    /**
     *
     * @return
     *
     * true if at current moment the regions map does not contain any items
     *
     */
    public boolean isEmpty() {
        return regions.isEmpty();
    }

    /**
     * prints to system.out
     *
     * prints the cached line of the thesaurus set
     *
     */
    public void printCached() {
        System.out.println("Cached:");
        if (thesline != null) {
            System.out.println(thesline);
        } else {
            System.out.println("thesline is null");
        }
        if (thesentry != null) {
            System.out.print(thesentry.toString());
        } else {
            System.out.println("thesentry is null");
        }
    }

    /**
     * Add an element into the nested regions hashmap.
     *
     * @param alignEnd
     * @param thesline
     * @param entry
     */
    private void addEntry(ThesaurusEntry entry) {
        ArrayList<ThesaurusEntry> temp = regions.get(entry.alignEnd);
        if (temp == null) {
            temp = new ArrayList<ThesaurusEntry>(4);
            regions.put(entry.alignEnd, temp);
        }
        temp.add(entry);
    }

    /**
     * reads lines from the thesaurus reader and returns when next line
     * describes a desired chromosome.
     *
     * Use this to ensure that the thesaurus is ready to process the chromosome.
     *
     * If current chromsome is "chr2" and user asks for "chr5", this function
     * will read and skip all lines describing to "chr3" and "chr4"
     *
     * @param chr
     * @throws IOException
     */
    private void skipToChr(String chr) throws IOException {

        if (thesline == null || thesentry == null) {
            return;
        }

        String chrtab = chr + "\t";
        int chrindex = ginfo.getChrIndex(chr);

        // check if need to get anything at all
        // if the cached line is chr5 and function is requested to skip to chr4, the function should not do anything
        if (thesentry.alignChrIndex > chrindex) {
            // pretend like the skipping has happened, but don't read anything from the file
            // or change the cached objects
            curchr = chr;
            uptoposition = 0;
            return;
        }

        // if reached here, 
        // read lines through the thesaurus until the desired chromosome appears in thesline
        // or until the chromosome in the thesline is further in the genome than the skipped-to chr
        while (!thesline.startsWith(chrtab) && ginfo.getChrIndex(thesline.split("\t", 2)[0]) < chrindex) {
            // if in here, the chromosome does not match, read next line of the thesaurus
            thesline = thesaurusreader.readLine();
            if (thesline == null) {
                thesentry = null;
                curchr = null;
                return;
            }
        }

        // if reached here, it means the thesaurus skipped to the correct chromosome
        thesentry = new ThesaurusEntry(thesline, ginfo);
        uptoposition = 0;
        curchr = thesentry.getAlignChr();
    }

    /**
     * Reads lines from the thesaurus. Those that overlap the position are added
     * to the map. Those that do not overlap with the position are just
     * skipped/ignored.
     *
     * This function uses a partial initialization trick to avoid setting up
     * full thesaurus entries unless they will actually be required to record.
     *
     * @param position
     * @throws IOException
     */
    private boolean loadUpToPosition(int positionMin, int positionMax) throws IOException {

        int curchrindex = ginfo.getChrIndex(curchr);

        // avoid loading anything unless the we are looking at the right chromosome
        // (this is necessary for when a chromosome is skipped in the thesaurus file)
        if (thesentry.alignChrIndex != curchrindex) {
            return true;
        }

        while (thesentry.isOK() && thesentry.alignChrIndex == curchrindex && thesentry.alignStart <= positionMax) {
            // maybe add the thesaurus entry, or maybe not
            // if an entry does not include the current position, the entry is not added
            if (thesentry.alignEnd >= positionMin) {
                thesentry.finishInit();
                addEntry(thesentry);
            }
            thesline = thesaurusreader.readLine();
            thesentry = new ThesaurusEntry(thesline, ginfo, false);
        }

        thesentry.finishInit();
        return true;
    }

    /**
     * the purpose of this class is to parse a thesaurus string in a partial
     * manner
     */
    class ThesaurusEntryLight {

        String thesline;
        String alignChr;
        int alignStart;

        public ThesaurusEntryLight(String thesline) {
            this.thesline = thesline;
        }
    }

    /**
     * This is a lookup/walk/search function with complex behavior.
     *
     * When the function is called, the ThesaurusRegionsMap object will assume
     * that the user is no longer interested in any regions with earlier genomic
     * coordinates than those specified.
     *
     * If the thesaurus information for the specified position is not loaded
     * yet, the function will load the data.
     *
     * Finally, the function outputs all the regions that contain the indicated
     * locus
     *
     *
     * @param chr
     *
     *
     * @param position
     * @param aroundbuffer
     *
     * a positive integer - interpret as read length
     *
     * @return
     *
     * a list of thesaurus entries that overlap the position chr:position with a
     * window maring of aroundbuffer i.e. overlap an interval
     * chr:(position-aroundbuffer) to chr:(position+aroundbuffer)
     *
     */
    public ThesaurusEntry[] lookup(String chr, int position, int aroundbuffer) throws IOException {

        if (thesaurusreader == null) {
            return null;
        }

        if (chr.equals(curchr)) {
            // clean up all thesaurus regions up to the locus   
            for (int i = uptoposition - bufferlength; i < position - bufferlength && !regions.isEmpty(); i++) {
                regions.remove(i);
            }
        } else {
            // change of chromosomes, so get rid of all regions on the current chromosome
            regions.clear();
            skipToChr(chr);
            // like in BamRegionsMap, clean up memory if the skip is successful
            if (chr.equals(curchr) && ginfo.getChrLength(chr) > 1e7) {
                System.gc();
            }
        }

        // it is possible thesaurus runs out of chromosomes while user still is trying to lookup things.
        // just abort here and report no entries
        if (curchr == null) {
            return null;
        }

        // make sure regions encompassing the position are loaded
        uptoposition = position;
        if (!loadUpToPosition(position - bufferlength, position + bufferlength)) {
            System.out.println("How did we get here?");
            return null;
        }

        return getEntries(position - aroundbuffer, position + aroundbuffer);
    }

    public void close() throws IOException {
        thesaurusreader.close();
    }
}
