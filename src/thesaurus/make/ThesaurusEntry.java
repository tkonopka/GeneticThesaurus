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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import jsequtils.genome.GenomeInfo;
import jsequtils.sequence.SequenceComplementer;
import jsequtils.sequence.SequenceMap;
import net.sf.samtools.SAMRecord;

/**
 * Container class to store one line of the thesaurus.
 *
 * Coordinates always follow a 1-based coordinate system.
 *
 * @author tkonopka
 */
class ThesaurusEntry {

    // GenomeInfo - allows storing chromosomes as integers
    private final GenomeInfo ginfo;
    // description of the region
    // String originChr, alignChr;
    int originChrIndex, alignChrIndex;
    int originStart, alignStart, originEnd, alignEnd;
    int penalty = -1;
    char alignStrand = '+';
    // description of the anchor point
    private ArrayList<ThesaurusAnchor> anchors;
    boolean ok = false;
    // used for comparing anchors, essentially sorts by alignment position (an integer)
    static ThesaurusAnchorComparator tac = new ThesaurusAnchorComparator();
    // parsed tokens is used during partial loading from strings
    private String[] parsedtokens = null;

    /**
     * Equivalent to calling constructor ThesaurusEntry(String, boolean) with a
     * full parse of the input.
     *
     * @param line
     */
    public ThesaurusEntry(String line, GenomeInfo ginfo) {
        this(line, ginfo, true);
    }

    /**
     * Construct an entry by parsing a line from a file. The line should be
     * formatted analogously to the output of toString(). The constructor will
     * set all fields of the ThesaurusEntry.
     *
     * @param line
     *
     * long string from thesaurus file
     *
     * @param fullparse
     *
     * if true, this constructor will full create a thesaurus entry.
     *
     * if false, loading will be faster but partial - run completeInit() to
     * finish parsing the anchors.
     *
     */
    public ThesaurusEntry(String line, GenomeInfo ginfo, boolean fullparse) {

        this.ginfo = ginfo;

        // parse the line
        String[] tokens = null;
        if (line != null) {
            tokens = line.split("\t");
        }

        // fill with dummy values if the input line is invalid
        if (line == null || tokens.length < 14) {
            originChrIndex = -1;
            alignChrIndex = -1;
            originStart = 0;
            originEnd = 0;
            alignStart = 0;
            alignEnd = 0;
            anchors = new ArrayList<ThesaurusAnchor>(2);
            return;
        }

        //alignChr = tokens[0];
        alignChrIndex = ginfo.getChrIndex(tokens[0]);
        alignStart = Integer.parseInt(tokens[1]);
        alignEnd = Integer.parseInt(tokens[2]);
        //originChr = tokens[3];
        originChrIndex = ginfo.getChrIndex(tokens[3]);
        originStart = Integer.parseInt(tokens[4]);
        originEnd = Integer.parseInt(tokens[5]);
        penalty = Integer.parseInt(tokens[6]);
        alignStrand = tokens[7].charAt(0);

        this.parsedtokens = tokens;
        if (fullparse) {
            finishInit();
        }

        ok = true;
    }

    public String getAlignChr() {
        return ginfo.getChrName(alignChrIndex);
    }

    public String getOriginChr() {
        return ginfo.getChrName(originChrIndex);
    }

    /**
     * call this function after an initialization with constructor
     * ThesaurusEntry(String, boolean=false)
     *
     */
    public final synchronized void finishInit() {
        if (parsedtokens == null) {
            return;
        }

        String alignPosition = parsedtokens[8];
        String alignRefBase = parsedtokens[9];
        String alignAltBase = parsedtokens[10];
        String originPosition = parsedtokens[11];
        String originRefBase = parsedtokens[12];
        String originAltBase = parsedtokens[13];

        anchors = toAnchorArray(alignPosition, originPosition, alignRefBase, originRefBase, alignAltBase, originAltBase);        
        if (penalty != anchors.size()) {
            String line = getAlignChr() + "\t" + this.alignStart + "\t" + this.alignEnd + "\t" + getOriginChr() + "\t" + this.originStart + "\t" + this.originEnd;
            System.out.println("Malformed ThesaurusEntry line (anchor and mismatch fields do not match):\n" + line);
        }

        // reset the line field so that the user cannot re-modifiy the anchor fields in this way.
        parsedtokens = null;
    }

    /**
     * This is a defensive copy constructor
     *
     * @param entry
     */
    public ThesaurusEntry(ThesaurusEntry entry) {
        ginfo = entry.ginfo;        
        originChrIndex = entry.originChrIndex;        
        alignChrIndex = entry.alignChrIndex;
        originStart = entry.originStart;
        originEnd = entry.originEnd;
        alignStart = entry.alignStart;
        alignEnd = entry.alignEnd;
        penalty = entry.penalty;
        alignStrand = entry.alignStrand;

        anchors = new ArrayList<ThesaurusAnchor>(entry.anchors.size());
        for (int i = 0; i < entry.anchors.size(); i++) {
            anchors.add(new ThesaurusAnchor(entry.anchors.get(i)));
        }

        ok = entry.ok;
    }

    /**
     * Construct an entry by parsing a SAM record for a read.
     *
     * The read name should be of the format, e.g. T100;chr1:10001-10100 i.e.
     * LABEL;CHROMOSOME;START-END
     *
     * Warning: this constructor will only set the following fields: originChr,
     * originStart, originEnd, alignChr, alignStart, alignEnd. Other fields must
     * be set via the separate set methods.
     *
     * @param record
     */
    public ThesaurusEntry(SAMRecord record, GenomeInfo ginfo) {
        this.ginfo = ginfo;

        // parse the read name to get the origin location        
        String[] readelements = record.getReadName().split(";|:|-");
        //originChr = readelements[1];
        originChrIndex = ginfo.getChrIndex(readelements[1]);
        originStart = Integer.parseInt(readelements[2]);
        originEnd = Integer.parseInt(readelements[3]);

        // get the alignment location
        //alignChr = record.getReferenceName();
        alignChrIndex = ginfo.getChrIndex(record.getReferenceName());
        alignStart = record.getAlignmentStart();
        alignEnd = record.getAlignmentEnd();

        if (record.getReadNegativeStrandFlag()) {
            alignStrand = '-';
        } else {
            alignStrand = '+';
        }

        anchors = new ArrayList<ThesaurusAnchor>(2);
        // the anchor array will not be properly formed, i.e. any mismatches will not
        // be recorded either in the anchors or in the penalty field.

        ok = true;
    }

    /**
     *
     * @return
     *
     * the status of this object as self-described
     *
     */
    public boolean isOK() {
        return ok;
    }

    /**
     * checks if the align and origin positions are the same.
     *
     * @return
     *
     * true if align and origin are the same.
     *
     */
    public boolean isTrivial() {
        //if (alignStart == originStart && alignEnd == originEnd && alignChr.equals(originChr) && alignStrand == '+') {
        //    return true;
        //}
        if (alignStart == originStart && alignEnd == originEnd && alignChrIndex == originChrIndex && alignStrand == '+') {
            return true;
        }
        return false;
    }

    /**
     *
     * @param Astart
     * @param Aend
     * @param Bstart
     * @param Bend
     * @return
     *
     * true if intervals [Astart, Aend] and [Bstart, Bend] have overlap.
     */
    public static boolean overlap(int Astart, int Aend, int Bstart, int Bend) {
        if (Astart <= Bstart) {
            if (Aend >= Bstart) {
                return true;
            } else {
                return false;
            }
        } else {
            if (Astart <= Bend) {
                return true;
            } else {
                return false;
            }

        }
    }

    /**
     *
     * @return
     *
     * true if the align and origin regions overlap
     */
    public boolean isAlignOriginOverlapping() {
        //if (alignChr.equals(originChr)) {
        //    return overlap(alignStart, alignEnd, originStart, originEnd);
        //} else {
        //    return false;
        //}

        if (alignChrIndex == originChrIndex) {
            return overlap(alignStart, alignEnd, originStart, originEnd);
        } else {
            return false;
        }
    }

    /**
     * converts information stored in the thesaurus into an array to easily
     * access individual anchors.
     *
     * @return
     */
    private static ArrayList<ThesaurusAnchor> toAnchorArray(String alignPosition, String originPosition,
            String alignRefBase, String originRefBase, String alignAltBase, String originAltBase) {

        if (alignPosition.equals("NA")) {
            return new ArrayList<ThesaurusAnchor>(1);
        }

        String[] alignPos = alignPosition.split(";");
        String[] alignRef = alignRefBase.split(";");
        String[] alignAlt = alignAltBase.split(";");
        String[] originPos = originPosition.split(";");
        String[] originRef = originRefBase.split(";");
        String[] originAlt = originAltBase.split(";");

        ArrayList<ThesaurusAnchor> array = new ArrayList<ThesaurusAnchor>(alignPos.length);

        for (int i = 0; i < alignPos.length; i++) {
            ThesaurusAnchor ta = new ThesaurusAnchor();
            ta.alignPosition = Integer.parseInt(alignPos[i]);
            ta.originPosition = Integer.parseInt(originPos[i]);
            ta.alignRef = alignRef[i].charAt(0);
            ta.alignAlt = alignAlt[i].charAt(0);
            ta.originRef = originRef[i].charAt(0);
            ta.originAlt = originAlt[i].charAt(0);
            array.add(ta);
        }
        return array;
    }

    public void addAnchor(ThesaurusAnchor anchor) {
        anchors.add(anchor);
        Collections.sort(anchors, tac);
        this.penalty = anchors.size();
    }

    public ThesaurusAnchor getAnchor(int index) {
        return anchors.get(index);
    }

    /**
     *
     * @return
     *
     * the number of anchors declared in this entry
     *
     */
    public int getNumAnchors() {
        return anchors.size();
    }

    /**
     * combines anchor points from "this" and "entry"
     *
     * @param entry
     * @return
     *
     * a list with unique anchor points
     *
     */
    private ArrayList<ThesaurusAnchor> getUniqueJointAnchors(ThesaurusEntry entryA, ThesaurusEntry entryB) {
        ArrayList<ThesaurusAnchor> anchorsA = entryA.anchors;
        ArrayList<ThesaurusAnchor> anchorsB = entryB.anchors;

        // copy the existing anchors into one array
        int aAs = anchorsA.size();
        int aBs = anchorsB.size();
        int aas = aAs + aBs;
        ThesaurusAnchor[] allanchors = new ThesaurusAnchor[aas];
        for (int i = 0; i < aAs; i++) {
            allanchors[i] = anchorsA.get(i);
        }
        for (int i = 0; i < aBs; i++) {
            allanchors[i + aAs] = anchorsB.get(i);
        }
        Arrays.sort(allanchors, tac);

        // make sure that all the anchors are unique, i.e. not equal
        ArrayList<ThesaurusAnchor> uniqueanchors = new ArrayList<ThesaurusAnchor>(aas);
        if (aas > 0) {
            uniqueanchors.add(allanchors[0]);
            for (int i = 1; i < aas; i++) {
                if (!allanchors[i].equals(allanchors[i - 1])) {
                    uniqueanchors.add(allanchors[i]);
                }
            }
        }


        return uniqueanchors;
    }

    /**
     * use the anchor array to create a string summarizing all the Anchor align
     * positions
     *
     * @param anchors
     * @return
     */
    private static String makeAnchorAlignPosition(ArrayList<ThesaurusAnchor> anchors) {
        if (anchors.isEmpty()) {
            return "NA";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < anchors.size(); i++) {
            sb.append(";").append(anchors.get(i).alignPosition);
        }
        return sb.toString().substring(1);
    }

    /**
     * use the anchor array to create a string summarizing all the anchor origin
     * positions
     *
     * @param anchors
     * @return
     */
    private static String makeAnchorOriginPosition(ArrayList<ThesaurusAnchor> anchors) {
        if (anchors.isEmpty()) {
            return "NA";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < anchors.size(); i++) {
            sb.append(";").append(anchors.get(i).originPosition);
        }
        return sb.toString().substring(1);
    }

    /**
     * use the anchor array to create a string summarizing all the anchor origin
     * reference base
     *
     * @param anchors
     * @return
     */
    private static String makeAnchorOriginRefBase(ArrayList<ThesaurusAnchor> anchors) {
        if (anchors.isEmpty()) {
            return "N";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < anchors.size(); i++) {
            sb.append(";").append(anchors.get(i).originRef);
        }
        return sb.toString().substring(1);
    }

    /**
     * use the anchor array to create a string summarizing all the anchor origin
     * alternative base
     *
     * @param anchors
     * @return
     */
    private static String makeAnchorOriginAltBase(ArrayList<ThesaurusAnchor> anchors) {
        if (anchors.isEmpty()) {
            return "N";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < anchors.size(); i++) {
            sb.append(";").append(anchors.get(i).originAlt);
        }
        return sb.toString().substring(1);
    }

    /**
     * use the anchor array to create a string summarizing all the anchor align
     * reference base
     *
     * @param anchors
     * @return
     */
    private static String makeAnchorAlignRefBase(ArrayList<ThesaurusAnchor> anchors) {
        if (anchors.isEmpty()) {
            return "N";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < anchors.size(); i++) {
            sb.append(";").append(anchors.get(i).alignRef);
        }
        return sb.toString().substring(1);
    }

    /**
     * use the anchor array to create a string summarizing all the anchor align
     * alternative base
     *
     * @param anchors
     * @return
     */
    private static String makeAnchorAlignAltBase(ArrayList<ThesaurusAnchor> anchors) {
        if (anchors.isEmpty()) {
            return "N";
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < anchors.size(); i++) {
            sb.append(";").append(anchors.get(i).alignAlt);
        }
        return sb.toString().substring(1);
    }

    /**
     * similar to mergeWith with number of mismatches, but this specifies a
     * constant error rate over the whole segment
     *
     * @param entry
     * @param errorrate
     * @return
     */
    public boolean mergeWith(ThesaurusEntry entry) {

        // Check that the chromosomes and strands match
        //if (alignStrand != entry.alignStrand || !alignChr.equals(entry.alignChr)
        //        || !originChr.equals(entry.originChr)) {
        //    return false;
        //}
        if (alignStrand != entry.alignStrand || alignChrIndex != entry.alignChrIndex
                || originChrIndex != entry.originChrIndex) {
            return false;
        }

        // Check that the alignments overlap        
        if (!overlap(originStart, originEnd, entry.originStart, entry.originEnd)) {
            return false;
        }
        if (!overlap(alignStart, alignEnd, entry.alignStart, entry.alignEnd)) {
            return false;
        }

        // make sure that offsets from start/end are correct so that joined segments are "in-phase"
        if (alignStrand == '+') {
            if (alignStart - entry.alignStart != originStart - entry.originStart) {
                return false;
            }
        } else {
            if (entry.alignEnd - alignEnd != originStart - entry.originStart) {
                return false;
            }
        }

        // require align and origin intervals to be the same after the merge
        // (This avoid bad results when two similar regions have long stretches
        // of polyX of slightly different lengths
        int newAlignStart = Math.min(alignStart, entry.alignStart);
        int newAlignEnd = Math.max(alignEnd, entry.alignEnd);
        int newOriginStart = Math.min(originStart, entry.originStart);
        int newOriginEnd = Math.max(originEnd, entry.originEnd);

        // check that length of region is the same        
        if (newOriginEnd - newOriginStart != newAlignEnd - newAlignStart) {
            return false;
        }

        // take care of a perverse case where (AAGGAA)GGAA is linked to AAGG(AAGGAA) on same chromosome
        // wherein after a merge the output would be linking a region to itself
        //if (newAlignStart == newOriginStart && newAlignEnd == newOriginEnd && alignChr.equals(originChr)) {
        //    return false;
        //}
        if (newAlignStart == newOriginStart && newAlignEnd == newOriginEnd && alignChrIndex == originChrIndex) {
            return false;
        }

        // if perfect overlap already, then do not do anything but just report as merge
        if (newAlignStart == alignStart && newAlignEnd == alignEnd) {
            return true;
        }

        // in this block avoid merging regions that overlap, but orientiation is not right
        // eg.   XXXXX------------XXXXX
        //           XXXXX----XXXXX
        //if (alignStrand == '+') {
        //    if (newAlignStart - newOriginStart != alignStart - originStart) {
        //        return false;
        //    }
        //} else {
        //    if (newAlignStart - newOriginStart == alignStart - originStart) {
        //        return false;
        //    }
        //}

        // if reached here, then accept the merge.
        // make a joint list of anchors here.
        anchors = getUniqueJointAnchors(this, entry);
        // change the coordinates in align and origin labels                        
        alignStart = newAlignStart;
        alignEnd = newAlignEnd;
        originStart = newOriginStart;
        originEnd = newOriginEnd;

        // perhaps change the anchor points
        this.penalty = anchors.size();

        return true;
    }

    /**
     * Helper function calculates how many anchors are stuck right at the edge
     * of a thesaurus entry.
     *
     * @param maxinarow
     *
     * dont' need to count beyond so many positions
     *
     * @param left
     *
     * set to true to look at positions near the alignment start and false to
     * look at positions near the alignment end
     *
     *
     * @return
     */
    private int getInARow(int maxinarow, boolean left) {
        int inarow = 0;
        int anchorssize = anchors.size();
        if (left) {
            for (int i = 0; i < maxinarow && i < anchorssize; i++) {
                if (anchors.get(i).alignPosition + i == alignStart + i) {
                    inarow++;
                } else {
                    break;
                }
            }
        } else {
            for (int i = 0; i < maxinarow && i < anchorssize; i++) {
                if (anchors.get(anchorssize - i - 1).alignPosition - i == this.alignEnd - i) {
                    inarow++;
                } else {
                    break;
                }
            }
        }
        return inarow;
    }

    /**
     * an entry shows similarity between two regions, but the true region of
     * similarity may be longer.
     *
     * This function looks at adjacent sequence and makes a thesaurus region
     * wider to the left.
     *
     * @param seqmap
     * @param errorrate
     *
     * maximum allowed error rate (anchors/segmentlength) in the entry
     *
     * @param maxinarow
     *
     * when extending, only allow so many anchors to appear together at the edge
     * of a segment
     *
     *
     */
    void extendLeft(SequenceMap seqmap, double errorrate, int maxinarow, int minalignstart) {
        // don't do anything if the sequancemap object is null
        if (seqmap == null) {
            return;
        }

        String alignChr = ginfo.getChrName(alignChrIndex);
        String originChr = ginfo.getChrName(originChrIndex);

        int originChrLength = seqmap.getChrLength(originChr);


        // count the number of anchors that are hugging the edge 
        // (prevents extension making too many anchors on the edge if the function is called multiple times)
        int inarow = getInARow(maxinarow, true);
        int anchorssize = anchors.size();

        // keep track of new anchors, their positions and number
        ArrayList<ThesaurusAnchor> newanchors = new ArrayList<ThesaurusAnchor>();
        int numnewanchors = 0;

        if (alignStrand == '+') {
            // origin and align reads are on same strand 

            // try to extend on the left               
            while (alignStart > minalignstart && originStart > 1 && inarow < maxinarow) {

                byte Aleft = seqmap.getSequenceBase1(alignChr, alignStart - 1);
                byte Oleft = seqmap.getSequenceBase1(originChr, originStart - 1);
                if (Aleft == Oleft) {
                    alignStart--;
                    originStart--;
                    inarow = 0;
                } else {
                    inarow++;
                    // there is a mismatch, either abort now or add mismatch to entry                    
                    if ((double) (anchorssize + numnewanchors + 1) >= errorrate * (1 + alignEnd - alignStart)) {
                        break;
                    } else {
                        // insert an anchor point and continue
                        alignStart--;
                        originStart--;
                        ThesaurusAnchor newanchor = new ThesaurusAnchor();
                        newanchor.alignPosition = alignStart;
                        newanchor.originPosition = originStart;
                        newanchor.alignRef = (char) Aleft;
                        newanchor.alignAlt = (char) Oleft;
                        newanchor.originRef = (char) Oleft;
                        newanchor.originAlt = (char) Aleft;
                        // add the new anchor to the beginning of the list                                                
                        newanchors.add(newanchor);
                        numnewanchors++;
                    }
                }
            }

        } else {
            // origin and align reads are on opposite strands
            // try to extend on the left             
            while (alignStart > minalignstart && originEnd < originChrLength && inarow < maxinarow) {
                byte Aleft = seqmap.getSequenceBase1(alignChr, alignStart - 1);
                byte Oright = seqmap.getSequenceBase1(originChr, originEnd + 1);
                if (Aleft == SequenceComplementer.complement(Oright)) {
                    alignStart--;
                    originEnd++;
                    inarow = 0;
                } else {
                    inarow++;
                    // there is a mismatch, either abort now or add mismatch to entry                    
                    if ((double) (anchorssize + numnewanchors + 1) >= errorrate * (1 + alignEnd - alignStart)) {
                        break;
                    } else {
                        // insert an anchor point and continue
                        alignStart--;
                        originEnd++;
                        ThesaurusAnchor newanchor = new ThesaurusAnchor();
                        newanchor.alignPosition = alignStart;
                        newanchor.originPosition = originEnd;
                        newanchor.alignRef = (char) Aleft;
                        newanchor.alignAlt = (char) SequenceComplementer.complement(Oright);
                        newanchor.originRef = (char) Oright;
                        newanchor.originAlt = (char) SequenceComplementer.complement(Aleft);
                        // add the new anchor to the beginning of the list                                                
                        newanchors.add(newanchor);
                        numnewanchors++;
                    }
                }
            }
        }

        if (!newanchors.isEmpty()) {
            Collections.reverse(newanchors);
            anchors.addAll(0, newanchors);
        }

        // finish by making sure the penalty field reflects any added anchors
        penalty = anchors.size();
    }

    /**
     * Similar to extendLeft, but this one shift boundary of alignment region to
     * the right
     *
     * @param seqmap
     * @param maxpenalty
     */
    void extendRight(SequenceMap seqmap, double errorrate, int maxinarow) {
        // don't do anything if the sequancemap object is null
        if (seqmap == null) {
            return;
        }

        String alignChr = ginfo.getChrName(alignChrIndex);
        String originChr = ginfo.getChrName(originChrIndex);

        int alignChrLength = seqmap.getChrLength(alignChr);
        int originChrLength = seqmap.getChrLength(originChr);

        // count the number of anchors that are hugging the edge 
        // (prevents extension making too many anchors on the edge if the function is called multiple times)
        int inarow = getInARow(maxinarow, false);
        int anchorssize = anchors.size();
        int numnewanchors = 0;

        if (alignStrand == '+') {
            // origin and align reads are on same strand 
            // try to extend on the right   
            while (alignEnd < alignChrLength && originEnd < originChrLength && inarow < maxinarow) {
                byte Aright = seqmap.getSequenceBase1(alignChr, alignEnd + 1);
                byte Oright = seqmap.getSequenceBase1(originChr, originEnd + 1);
                if (Aright == Oright) {
                    alignEnd++;
                    originEnd++;
                    inarow = 0;
                } else {
                    inarow++;
                    // there is a mismatch, either abort now or add mismatch to entry                    
                    if ((double) (anchorssize + numnewanchors + 1) >= errorrate * (1 + alignEnd - alignStart)) {
                        break;
                    } else {
                        // insert an anchor point and continue
                        alignEnd++;
                        originEnd++;
                        ThesaurusAnchor newanchor = new ThesaurusAnchor();
                        newanchor.alignPosition = alignEnd;
                        newanchor.originPosition = originEnd;
                        newanchor.alignRef = (char) Aright;
                        newanchor.alignAlt = (char) Oright;
                        newanchor.originRef = (char) Oright;
                        newanchor.originAlt = (char) Aright;
                        // add the new anchor at the end of the list                        
                        anchors.add(newanchor);
                        numnewanchors++;
                    }
                }
            }
        } else {
            // origin and align reads are on opposite strands            
            // try to extend on the right               
            while (alignEnd < alignChrLength && originStart > 1 && inarow < maxinarow) {
                byte Aright = seqmap.getSequenceBase1(alignChr, alignEnd + 1);
                byte Oleft = seqmap.getSequenceBase1(originChr, originStart - 1);
                if (Aright == SequenceComplementer.complement(Oleft)) {
                    alignEnd++;
                    originStart--;
                    inarow = 0;
                } else {
                    inarow++;
                    if ((double) (anchorssize + numnewanchors + 1) >= errorrate * (1 + alignEnd - alignStart)) {
                        break;
                    } else {
                        // insert an anchor point and continue
                        alignEnd++;
                        originStart--;
                        ThesaurusAnchor newanchor = new ThesaurusAnchor();
                        newanchor.alignPosition = alignEnd;
                        newanchor.originPosition = originStart;
                        newanchor.alignRef = (char) Aright;
                        newanchor.alignAlt = (char) SequenceComplementer.complement(Oleft);
                        newanchor.originRef = (char) Oleft;
                        newanchor.originAlt = (char) SequenceComplementer.complement(Aright);
                        // add the new anchor at the end of the list                        
                        anchors.add(newanchor);
                        numnewanchors++;
                    }
                }
            }
        }

        // finish by making sure the penalty field reflects any added anchors
        penalty = anchors.size();

    }

    /**
     *
     * @return
     *
     * a string that can be used as a table header. It explains the order of
     * values output by the toString() function.
     *
     */
    public static String getHeader() {
        return "Align.chr\tAlign.start\tAlign.end\t"
                + "Origin.chr\tOrigin.start\tOrigin.end\t"
                + "Mismatches\tAlign.strand\t"
                + "Align.position\tAlign.refBase\tAlign.altBase\t"
                + "Origin.position\tOrigin.refBase\tOrigin.altBase\n";
    }

    @Override
    public String toString() {
        return getAlignChr() + "\t" + alignStart + "\t" + alignEnd + "\t"
                + getOriginChr() + "\t" + originStart + "\t" + originEnd + "\t"
                + penalty + "\t" + alignStrand + "\t" + makeAnchorAlignPosition(anchors)
                + "\t" + makeAnchorAlignRefBase(anchors)
                + "\t" + makeAnchorAlignAltBase(anchors)
                + "\t" + makeAnchorOriginPosition(anchors)
                + "\t" + makeAnchorOriginRefBase(anchors)
                + "\t" + makeAnchorOriginAltBase(anchors) + "\n";
    }
}

/**
 * Comparator that will allow sorting of Thesaurus entries by their alignment
 * position. (Contrast with ThesaurusEntryOriginComparator)
 *
 * @author tkonopka
 */
class ThesaurusEntryAlignComparator implements Comparator {
    
    @Override
    public int compare(Object o1, Object o2) {
        ThesaurusEntry oi1 = (ThesaurusEntry) o1;
        ThesaurusEntry oi2 = (ThesaurusEntry) o2;

        int alignChr1 = oi1.alignChrIndex;
        int alignChr2 = oi2.alignChrIndex;


        if (alignChr1 < alignChr2) {
            return -1;
        } else if (alignChr1 > alignChr2) {
            return 1;
        }

        if (oi1.alignStart < oi2.alignStart) {
            return -1;
        } else if (oi1.alignStart > oi2.alignStart) {
            return 1;
        } else {
            // align positions are equal, so check for origin positions            

            int originChr1 = oi1.originChrIndex;
            int originChr2 = oi2.originChrIndex;
            
            if (originChr1 < originChr2) {
                return -1;
            } else if (originChr1 > originChr2) {
                return 1;
            }

            if (oi1.originStart < oi2.originStart) {
                return -1;
            } else if (oi1.originStart > oi2.originStart) {
                return 1;
            } else {
                return 0;
            }

        }

    }
}

/**
 * Comparator that will allow sorting of Thesaurus entries by their alignment
 * position and by their origin position. If all goes well, entries that are
 * overlapping should be close to each other.
 *
 * @author tkonopka
 */
class ThesaurusEntryMergingComparator implements Comparator {

    //private final HashMap<String, Integer> chrorder;

    //public ThesaurusEntryMergingComparator(GenomeInfo ginfo) {
    //    chrorder = new HashMap<String, Integer>(ginfo.getNumChromosomes() * 2);
    //    for (int i = 0; i < ginfo.getNumChromosomes(); i++) {
    //        chrorder.put(ginfo.getChrName(i), ginfo.getChrIndex(ginfo.getChrName(i)));
    //    }
    //}    
    
    @Override
    public int compare(Object o1, Object o2) {
        ThesaurusEntry oi1 = (ThesaurusEntry) o1;
        ThesaurusEntry oi2 = (ThesaurusEntry) o2;

        // always arrange same strand items first (will faciliate merging)
        if (oi1.alignStrand != oi2.alignStrand) {
            if (oi1.alignStrand == '+') {
                return -1;
            } else {
                return 1;
            }
        }

        // after the strand, use chromosomes, arrange them in order of the genome
        int alignChr1 = oi1.alignChrIndex;
        int alignChr2 = oi2.alignChrIndex;
        
        if (alignChr1 < alignChr2) {
            return -1;
        } else if (alignChr1 > alignChr2) {
            return 1;
        }

        // after chromosome, calculate distance between origin and alignment
        // segments that can be merged will have the same distance, i.e. will be able to merge
        int x1 = oi1.alignStart - oi1.originStart;
        int x2 = oi2.alignStart - oi2.originStart;

        if (x1 < x2) {
            return -1;
        } else if (x1 > x2) {
            return 1;
        } else {
            return 0;
        }


    }
}

/**
 * Container storing anchor positions for thesaurus entries.
 *
 * @author tkonopka
 */
class ThesaurusAnchor {

    int alignPosition, originPosition;
    char alignRef, alignAlt, originRef, originAlt;

    public ThesaurusAnchor() {
    }

    public ThesaurusAnchor(ThesaurusAnchor anchor) {
        alignPosition = anchor.alignPosition;
        originPosition = anchor.originPosition;
        alignRef = anchor.alignRef;
        alignAlt = anchor.alignAlt;
        originRef = anchor.originRef;
        originAlt = anchor.originAlt;
    }

    public boolean equals(ThesaurusAnchor ta) {
        if (alignPosition != ta.alignPosition || originPosition != ta.originPosition) {
            return false;
        }
        return (alignRef == ta.alignRef && alignAlt == ta.alignAlt
                && originRef == ta.originRef && originAlt == ta.originAlt);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("align: ").append(alignPosition).append(alignRef).append(alignAlt).append("\torigin: ").append(originPosition).append(originRef).append(originAlt);
        return sb.toString();
    }
}

class ThesaurusAnchorComparator implements Comparator {

    @Override
    public int compare(Object o1, Object o2) {
        ThesaurusAnchor ta1 = (ThesaurusAnchor) o1;
        ThesaurusAnchor ta2 = (ThesaurusAnchor) o2;

        if (ta1.alignPosition < ta2.alignPosition) {
            return -1;
        } else if (ta1.alignPosition > ta2.alignPosition) {
            return 1;
        } else {
            return 0;
        }
    }
}