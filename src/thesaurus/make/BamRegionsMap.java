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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import jsequtils.genome.GenomeInfo;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 * A "buffered reader" for bam files. Assumes bam file is sorted, allows users
 * to query for reads aligning onto loci.
 *
 * Similar to ThesaurusRegionsMap.
 *
 *
 * @author tkonopka
 */
class BamRegionsMap {

    private final GenomeInfo ginfo;
    // minimum mappqing quality to consider
    private final int minmapqual;
    // the buffered reader providing access to the reader
    private final SAMFileReader bamreader;
    private final SAMRecordIterator bamiterator;
    // interval representation of the thesaurus regions
    private final HashMap<Integer, ArrayList<SAMRecord>> regions;
    // bookkeeping variables
    private SAMRecord record = null;
    private String curchr = null;
    private int uptoposition = 0;

    public BamRegionsMap(File bamfile, GenomeInfo ginfo, String validate, int minmapqual) throws IOException {

        this.ginfo = ginfo;
        this.minmapqual = minmapqual;
        if (bamfile != null) {
            bamreader = new SAMFileReader(bamfile);
            if (validate.equals("STRICT")) {
                bamreader.setValidationStringency(SAMFileReader.ValidationStringency.STRICT);
            } else if (validate.equals("LENIENT")) {
                bamreader.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
            } else {
                bamreader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            }

            bamiterator = bamreader.iterator();
            if (bamiterator.hasNext()) {
                record = bamiterator.next();
                curchr = record.getReferenceName();
            }

        } else {
            bamreader = null;
            bamiterator = null;
        }

        regions = new HashMap<Integer, ArrayList<SAMRecord>>(128);
    }

    /**
     * make a copy of the current state of the regions object, into a primitive
     * array
     *
     * @return
     */
    public SAMRecord[] getRegions() {

        // simple exit if no regions available
        if (isEmpty()) {
            return null;
        }

        // figure out how many regions are current stored in memory
        int count = getRegionsCount();

        // make an array of ThesaurusEntries here
        SAMRecord[] ans = new SAMRecord[count];
        int index = 0;
        for (Map.Entry<Integer, ArrayList<SAMRecord>> entry : regions.entrySet()) {
            ArrayList<SAMRecord> map1 = entry.getValue();
            for (int i = 0; i < map1.size(); i++) {
                ans[index] = map1.get(i);
                index++;
            }
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

    private int getRegionsCount() {
        int count = 0;
        for (Map.Entry<Integer, ArrayList<SAMRecord>> entry : regions.entrySet()) {
            count += entry.getValue().size();
        }
        return count;
    }

    /**
     * reads lines from the bamfile reader and returns when next line describes
     * a desired chromosome.
     *
     * Use this to ensure that readiness to process the chromosome.
     *
     * If current chromsome is "chr2" and user asks for "chr5", this function
     * will read and skip all lines describing to "chr3" and "chr4"
     *
     * @param chr
     * @throws IOException
     */
    private void skipToChr(String chr) throws IOException {

        // *****************
        // should this function somehow keep track of all chromosomes read so far? 
        // *****************
        if (record == null) {
            return;
        }

        int chrindex = ginfo.getChrIndex(chr);

        // check if need to get anything at all
        // if the cached line is chr5 and function is requested to skip to chr4, the function should not do anything
        if (ginfo.getChrIndex(record.getReferenceName()) > chrindex) {
            // pretend like the skipping has happened, but don't read anything from the file
            // or change the cached objects
            curchr = chr;
            uptoposition = 0;
            return;
        }

        // if reached here, read lines through the bam file 
        // until the desired chromosome appears        
        while (!record.getReferenceName().equals(chr) && ginfo.getChrIndex(record.getReferenceName()) < chrindex) {
            // if in here, the chromosome does not match, read next record
            if (bamiterator.hasNext()) {
                record = bamiterator.next();
            } else {
                record = null;
                return;
            }
        }

        //System.out.println(record.getReferenceName() + ":" + record.getAlignmentStart() + " " + record.getReadName());
        uptoposition = 0;
        curchr = record.getReferenceName();
    }

    /**
     * Add an element into the nested regions hashmap.
     *
     * @param alignEnd
     * @param thesline
     * @param entry
     */
    private void addEntry(SAMRecord entry) {
        ArrayList<SAMRecord> temp = regions.get(entry.getAlignmentEnd());
        if (temp == null) {
            temp = new ArrayList<SAMRecord>(8);
            regions.put(entry.getAlignmentEnd(), temp);
        }        
        temp.add(entry);
    }

    /**
     * reads lines from the thesaurus
     *
     * @param position
     * @throws IOException
     */
    private boolean loadUpToPosition(int position) throws IOException {

        if (record == null) {
            return false;
        }

        // avoid loading anything unless the we are looking at the right chromosome
        // (this is necessary for when a chromosome is skipped in bam file)
        if (!record.getReferenceName().equals(curchr)) {
            return true;
        }

        uptoposition = position;
        //System.out.println("Bam loading up to : " + position);

        while (record != null && record.getReferenceName().equals(curchr) && record.getAlignmentStart() <= position) {
            // maybe add the thesaurus entry, or maybe not
            // if an entry does not include the current position, the entry is not added
            if (record.getAlignmentEnd() >= position
                    && record.getMappingQuality() >= minmapqual
                    && !record.getNotPrimaryAlignmentFlag()
                    && !record.getDuplicateReadFlag()
                    && !record.getReadUnmappedFlag()) {
                addEntry(record);
            }

            if (bamiterator.hasNext()) {
                record = bamiterator.next();
            } else {
                record = null;
            }
        }

        //System.out.println("Bam loading up to : " + position + " done");
        return true;
    }

    /**
     * This is a lookup/walk/search function with complex behavior.
     *
     * When the function is called, the BamRegionsMap object will assume that
     * the user is no longer interested in any regions with earlier genomic
     * coordinates than those specified.
     *
     * If the information for the specified position is not loaded yet, the
     * function will load the data.
     *
     * Finally, the function outputs all the SAMRecords that contain the
     * indicated locus.
     *
     *
     * @param chr
     * @param position
     * @return
     */
    public SAMRecord[] lookup(String chr, int position) throws IOException {

        if (bamreader == null) {
            return null;
        }

        if (chr.equals(curchr)) {
            // clean up all map regions up to the locus
            for (int i = uptoposition; i < position && !regions.isEmpty(); i++) {
                regions.remove(i);
            }
        } else {
            // change of chromosomes, so get rid of all regions on the current chromosome
            regions.clear();
            skipToChr(chr);
            // if the skip is successful, i.e. curchr is updated
            // clean up some memory after the old chromosome
            if (chr.equals(curchr)) {
                System.gc();
            }
        }

        // make sure regions encompassing the position are loaded
        if (!loadUpToPosition(position)) {
            if (uptoposition == 0) {
                System.out.println("Thesaurus does not contain information about chromosome: " + chr);
                return null;
            }
            return null;
        }

        return getRegions();
    }

    public void close() throws IOException {
        if (bamreader != null) {
            bamreader.close();
        }
    }
}
