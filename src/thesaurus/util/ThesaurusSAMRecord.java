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
package thesaurus.util;

import java.util.ArrayList;
import java.util.Arrays;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import static net.sf.samtools.CigarOperator.D;
import static net.sf.samtools.CigarOperator.I;
import static net.sf.samtools.CigarOperator.M;
import static net.sf.samtools.CigarOperator.N;
import static net.sf.samtools.CigarOperator.S;
import net.sf.samtools.SAMRecord;

/**
 * A wrapper for SAMRecord which provides some extra functionality. It allows
 * user to extract a base given a genomic position, which is not directly
 * possible from the SAMRecord. It also enables counting mismatches and
 * modifying the record (soft clipping ends)
 *
 * @author tkonopka
 */
public class ThesaurusSAMRecord {

    // the SAMrecord that this object is based on
    final private SAMRecord record;
    // shortcuts/processed versions of information stored in the record
    private CigarElement[] cigarelements;
    private int cigarsize;
    final private byte[] bases;
    private int alignStart;
    private int alignEnd;
    // additional information
    private int mismatches = -1;

    /**
     * basic constructor that copies a SAMrecord and computes minimal things
     *
     * @param record
     */
    public ThesaurusSAMRecord(SAMRecord record) {
        this.record = record;
        Cigar cigar = record.getCigar();
        cigarsize = cigar.numCigarElements();
        cigarelements = new CigarElement[cigarsize];
        for (int i = 0; i < cigarsize; i++) {
            cigarelements[i] = cigar.getCigarElement(i);
        }
        bases = record.getReadBases();
        alignStart = record.getAlignmentStart();
        alignEnd = record.getAlignmentEnd();
    }

    public ThesaurusSAMRecord(SAMRecord record, byte[] referencesequence, int cliplength) {

        // first create the basic record
        this(record);

        // the use reference sequence information to soft clip read ends
        int oldStart = alignStart;
        int oldEnd = alignEnd;
        //System.out.println("QWE1 " + record.getReadName() + " " + record.getAlignmentStart());
        softClipEnds(referencesequence, cliplength);
        //System.out.println("QWE2");
        int newStart = alignStart;
        int newEnd = alignEnd;
        //System.out.println(newStart + " " + newEnd);

        // if clipping has occurred, get a subsequence of the referencesequence that matches the current alignment positions        
        byte[] subrefseq = referencesequence;
        if (newStart != oldStart || oldEnd != newEnd) {
            // some clipping has occurred            
            int offset = newStart - oldStart;
            //System.out.println("offset is " + offset);
            subrefseq = Arrays.copyOfRange(referencesequence, offset, offset + (newEnd - newStart + 1));
            //System.out.println("after copy");
        }
        //System.out.println("QSE3");
        // calculate the number of mismatches
        this.mismatches = calcMismatches(subrefseq, newStart);
        //System.out.println("QSE4 mismatches " + mismatches);
        //System.out.println();
    }

    public int getNumMismatches() {
        return mismatches;
    }

    /**
     * Get a base in a read aligned onto a given position. e.g. if a read is
     * aligned at chr1:101-200, the user can query position=101 and obtain the
     * first base in the read, or position=104 and obtain the fourth base in the
     * read
     *
     * @param position
     * @return
     */
    public byte getBaseAtGenomicPosition(int position) {
        int temp = getPositionInRecord(position);
        if (temp < 0) {
            return 'N';
        } else {
            return bases[temp];
        }
    }

    /**
     * helper function that converts a genomic coordinate into an index relative
     * to the read sequence.
     *
     * @param position
     *
     * genomic position to look for (e.g. base 10000 on chr 3)
     *
     * @return
     *
     * the index in the bases array that corresponds to the requested genomic
     * position.
     *
     */
    public int getPositionInRecord(int position) {

        int chrpos = alignStart;
        int ans = 0;

        for (int i = 0; i < cigarsize; i++) {
            CigarElement ce = cigarelements[i];
            int celen = ce.getLength();

            switch (ce.getOperator()) {
                case M:
                    ans += celen;
                    chrpos += celen;
                    int temp = chrpos - position;
                    if (temp >= 0) {
                        return ans - temp;
                    }
                    break;
                case I:
                    ans += celen;
                    break;
                case D:
                    chrpos += celen;
                    break;
                case S:
                    ans += celen;
                    break;
                case N:
                    chrpos += celen;
                    break;
                default:
                    // Anything here? I think I have all the options covered
                    System.out.println("record: " + record.getReadName());
                    System.out.println("problematic cigar: " + record.getCigarString());
                    break;
            }
        }
        return -1;
    }

    /**
     * passes on the request to the SAMRecord
     *
     * @return
     */
    public byte[] getReadBases() {
        return record.getReadBases();
    }

    public int getAlignmentStart() {
        return this.alignStart;
    }

    public int getAlignmentEnd() {
        return this.alignEnd;
    }

    public String getChr() {
        return record.getReferenceName();
    }

    public String getReadName() {
        return record.getReadName();
    }

    /**
     *
     * @return
     *
     * an array of cigar elements as defined in the record.
     *
     * DO not change the values of these elements as the changes will not be
     * consistent any more with the record itself.
     *
     */
    public CigarElement[] getCigarElements() {
        return cigarelements;
    }

    /**
     * This block is almost-verbatim borrowed from the Bamformatics BamfoRecord
     * class. (Author: Tomasz Konopka, Licence: Apache 2.0)
     */
    private int[] getBasePositions() {
        int[] positions = new int[bases.length];
        int nowpos = alignStart;
        int nowindex = 0;
        for (int i = 0; i < cigarsize; i++) {
            CigarElement ce = cigarelements[i];
            int celen = ce.getLength();
            switch (ce.getOperator()) {
                case M:
                    for (int k = 0; k < celen; k++) {
                        positions[nowindex] = nowpos;
                        nowpos++;
                        nowindex++;
                    }
                    break;
                case D:
                    nowpos += celen;
                    break;
                case N:
                    nowpos += celen;
                    break;
                case S:
                    for (int j = 0; j < celen; j++) {
                        positions[nowindex] = -1;
                        nowindex++;
                    }
                    break;
                case I:
                    for (int j = 0; j < celen; j++) {
                        positions[nowindex] = -1;
                        nowindex++;
                    }
                    break;
                default:
                    break;
            }
        }
        return positions;
    }

    /**
     *
     * @return
     *
     * true if the cigar has operators such as IDNS close to each other, e.g.
     * 10M 3D 4I 10M
     *
     */
    private boolean isComplexCigar() {
        if (cigarsize == 1) {
            CigarOperator cone = cigarelements[0].getOperator();
            if (cone != CigarOperator.M) {
                return true;
            }
        } else if (cigarsize == 2) {
            CigarOperator cone = cigarelements[0].getOperator();
            CigarOperator ctwo = cigarelements[1].getOperator();
            if (cone != CigarOperator.M && ctwo != CigarOperator.M) {
                return true;
            }
        } else if (cigarsize > 2) {
            for (int i = 1; i < cigarsize - 1; i++) {
                CigarOperator co = cigarelements[i].getOperator();
                if (co == CigarOperator.D || co == CigarOperator.S || co == CigarOperator.I || co == CigarOperator.N) {
                    CigarOperator cprev = cigarelements[i - 1].getOperator();
                    CigarOperator cnext = cigarelements[i + 1].getOperator();
                    if (cprev != CigarOperator.M || cnext != CigarOperator.M) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Calculate the number of mismatches between the read sequence and the
     * reference sequence
     *
     * @param referencesequence
     *
     * section of sequence of the reference genome sequence
     *
     * @param refstart
     *
     * since referencesequence is a fragment, refstart denotes the position of
     * the first base in the referencesequence array.
     *
     *
     *
     */
    private int calcMismatches(byte[] referencesequence, int refstart) {
        if (isComplexCigar()) {
            return 0;
        }
        int count = this.getNumIndelOrSplice();

        // get position information about all bases        
        int[] positions = getBasePositions();
        int poslen = positions.length;

        for (int i = 0; i < poslen; i++) {
            int nowpos = positions[i];
            if (nowpos >= 0) {
                if (referencesequence[nowpos - refstart] != bases[i]) {
                    count++;
                }
            }
        }
        return count;
    }

    /**
     * Clip a record from both sides. This function changes a **copy** of the
     * original read alignStart, alignEnd, and cigar representations. The
     * original read stored within this object remains unchanged.
     *
     * @param referencesequence
     *
     * a byte array holding reference sequence. This should be equal in length
     * to the read length
     *
     * @param cliplength
     *
     * length to clip - requires these many consequent matches at read ends to
     * avoid clipping.
     *
     * @param cliplength
     */
    private void softClipEnds(byte[] referencesequence, int cliplength) {

        int readlength = bases.length;
        // avoid clipping reads that are too short
        if (readlength < cliplength * 2 || cliplength < 1) {
            return;
        }

        // avoid processing cigars that have adjacent non-M characters, e.g. 10M(3D4I)10M
        if (isComplexCigar()) {
            return;
        }

        // get position information about all bases        
        int[] positions = getBasePositions();

        // identify first good base on the left
        int leftindex = 0;
        boolean keepgoing = true;
        int readindex = 0;
        int referenceindex = 0;
        // need to keep count of the number of matches and mismatches        
        int rowmatches = 0; // count the number of matches in a row
        while (keepgoing && readindex < readlength) {
            if (positions[readindex] > 0) {
                if (bases[readindex] == referencesequence[referenceindex]) {
                    rowmatches++;
                } else {
                    rowmatches = 0;
                }
                // advance on the reference array by either one or by the offset from
                // previous position (offset may be >1 if read has splice/insertion
                if (readindex == 0) {
                    referenceindex++;
                } else {
                    if (positions[readindex - 1] > 0) {
                        referenceindex += positions[readindex] - positions[readindex - 1];
                    } else {
                        referenceindex++;
                    }
                }
            } else {
                // if there is an insertion/clip, so reset match count
                rowmatches = 0;
            }
            // exit condition is when 
            if (rowmatches >= cliplength) {
                keepgoing = false;
                leftindex = readindex - cliplength + 1;
            } else {
                readindex++;
            }
        }

        // identify first good base from the right edge  
        int rightindex = readlength;
        readindex = readlength - 1;
        referenceindex = referencesequence.length - 1;
        rowmatches = 0;
        keepgoing = true;
        while (keepgoing && readindex >= 0) {
            if (positions[readindex] > 0) {
                if (bases[readindex] == referencesequence[referenceindex]) {
                    rowmatches++;
                } else {
                    rowmatches = 0;
                }
                if (readindex == readlength - 1) {
                    referenceindex--;
                } else {
                    if (positions[readindex + 1] > 0) {
                        referenceindex -= (positions[readindex + 1] - positions[readindex]);
                    } else {
                        referenceindex--;
                    }
                }
            } else {
                rowmatches = 0;
            }
            if (rowmatches >= cliplength) {
                keepgoing = false;
                rightindex = readindex + cliplength - 1;
            } else {
                readindex--;
            }
        }

        // rightindex is equal to readlength when a read is messed up that rightindex loop
        // above finishes scanning the whole read but does not find a hit. 
        // Here handle this by clipping the whole read except for one base.
        if (rightindex == readlength) {
            rightindex = leftindex;
        }

        // don't do anything if no clipping will be required
        if (leftindex == 0 && rightindex == readlength - 1) {
            return;
        }
        // don't do anything if left and right indexes are weird
        // does this happen? I don't know... but safe to include this here
        if (leftindex > rightindex) {
            return;
        }

        // change the alignment start/end positions because of clipping
        if (positions[leftindex]>=alignStart) {
            alignStart = positions[leftindex];
        }
        if (positions[rightindex]>=alignStart) {
            alignEnd = positions[rightindex];
        }
                
        // change the position array to introduce clipping
        for (int i = 0; i < leftindex; i++) {
            positions[i] = -1;
        }
        for (int i = readlength - 1; i > rightindex; i--) {
            positions[i] = -1;
        }

        // make a new cigar using the positions array - add the new elements into a collections array
        ArrayList<CigarElement> tempcigar = new ArrayList<CigarElement>();
        int nowindex = 0;
        int previndex = 0;
        while (nowindex <= readlength) {
            if (nowindex == 0) {
                // skip this case
            } else if (nowindex == readlength) {
                // here must close the shop
                if (positions[nowindex - 1] > 0) {
                    tempcigar.add(new CigarElement(nowindex - previndex, CigarOperator.M));
                    previndex = nowindex; // but will no longer be needed
                } else {
                    tempcigar.add(new CigarElement(nowindex - previndex, CigarOperator.S));
                    previndex = nowindex; // but will no longer be needed
                }
            } else {
                if (positions[nowindex] > 0) {
                    if (positions[nowindex - 1] > 0) {
                        int temp = positions[nowindex] - positions[nowindex - 1];
                        if (temp > 1) {
                            // there is missing space between read-adjacent bases, so there must be a deletion
                            // perhaps record a match region, and then an insertion
                            if (nowindex - previndex > 0) {
                                tempcigar.add(new CigarElement(nowindex - previndex, CigarOperator.M));
                            }
                            tempcigar.add(new CigarElement(temp - 1, CigarOperator.D));
                            previndex = nowindex;
                        }
                    } else {
                        // current base is aligned, but previous base is handing (i.e. inserted)
                        tempcigar.add(new CigarElement(nowindex - previndex, CigarOperator.I));
                        previndex = nowindex;
                    }
                } else {
                    if (positions[nowindex - 1] > 0) {
                        // this is the first of clip or insertion
                        tempcigar.add(new CigarElement(nowindex - previndex, CigarOperator.M));
                        previndex = nowindex;
                    }
                }
            }
            nowindex++;
        }

        // convert the first and last cigar elements from insertions to softclips
        if (tempcigar.size() > 1) {
            CigarElement ce = tempcigar.get(0);
            if (ce.getOperator() == CigarOperator.I) {
                tempcigar.set(0, new CigarElement(ce.getLength(), CigarOperator.S));
            }
            ce = tempcigar.get(tempcigar.size() - 1);
            if (ce.getOperator() == CigarOperator.I) {
                tempcigar.set(tempcigar.size() - 1, new CigarElement(ce.getLength(), CigarOperator.S));
            }
        }

        // conver the arraylist into a simple array
        cigarelements = new CigarElement[tempcigar.size()];
        for (int i = 0; i < cigarelements.length; i++) {
            CigarElement tempce = tempcigar.get(i);
            cigarelements[i] = new CigarElement(tempce.getLength(), tempce.getOperator());
        }
        cigarsize = cigarelements.length;
        if (cigarsize == 0) {
            System.out.println("=== cigar zero size? " + this.getReadName() + " " + leftindex + " " + rightindex);
        }        
    }

    /**
     *
     * @return
     *
     * true if (possible modified) cigar holds an insertion or a deletion
     *
     */
    public boolean hasIndelOrSplice() {
        //return getNumIndelOrSplice()>0;
        int cigarlen = cigarelements.length;
        for (int i = 0; i < cigarlen; i++) {
            CigarElement ce = cigarelements[i];
            switch (ce.getOperator()) {
                case I:
                    return true;
                case D:
                    return true;
                case N:
                    return true;
                default:
            }
        }
        return false;
    }

    /**
     *
     * @return
     *
     * the number of insertion, deletion, or splice events in the cigar
     *
     */
    public int getNumIndelOrSplice() {
        int cigarlen = cigarelements.length;
        int count = 0;
        for (int i = 0; i < cigarlen; i++) {
            CigarElement ce = cigarelements[i];
            switch (ce.getOperator()) {
                case I:
                    count++;
                    break;
                case D:
                    count++;
                    break;
                case N:
                    count++;
                    break;
                default:
            }
        }
        return count;
    }

    public String getCigar() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < cigarsize; i++) {
            sb.append(cigarelements[i].getLength()).append(cigarelements[i].getOperator().toString());
        }
        return sb.toString();
    }

    /**
     *
     * @return
     *
     * mate information as reported in the SAM record
     */
    public String getMateChr() {
        return record.getMateReferenceName();
    }

    /**
     *
     * @return
     *
     * mate information recorded in the SAM record
     *
     */
    public int getMateStart() {
        return record.getMateAlignmentStart();
    }

    /**
     *
     * @return
     *
     * number of bases in the read
     *
     */
    public int getReadLength() {
        return bases.length;
    }

    public boolean hasMate() {
        return record.getReadPairedFlag();
    }

    public int getReferenceIndex() {
        return record.getReferenceIndex();
    }

    public int getMateReferenceIndex() {
        return record.getMateReferenceIndex();
    }
}
