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

import thesaurus.util.SNVPosition;
import thesaurus.util.ThesaurusSAMRecord;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.sequence.SequenceComplementer;
import jsequtils.variants.VCFEntry;
import jsequtils.variants.VCFEntrySet;

/**
 * This is an important class for variant filtering, i.e. determing whether a
 * variant can be linked with another genomic position via the thesaurus.
 *
 * Calculates the alternate loci based for a variant given reads that overlap
 * the variant and a set of thesaurus entries.
 *
 * @author tkonopka
 */
class ThesaurusSynonyms {

    // inputs to the class
    // the variant to be studied
    private final VCFEntry varentry;
    private final String[] varalt;
    private final GenomeInfo ginfo;
    private final GenomePositionComparator gcomp;
    // regions that overlap with the varaint
    private final ThesaurusEntry[] thesentries;
    // all variants called together with varentry
    private final VCFEntrySet calledvariants;

    /**
     * Define and find synonymous location for a variant. The constructor
     * collects the variant, thesaurus regions, and bam evidence. It also
     * computes the synonymous variants. Those variant can then be read out
     * using the public functions.
     *
     *
     * @param thesentries
     *
     * vector of thesaurus entries overlapping with the variant.
     *
     *
     * @param varentry
     *
     * the variant in question. The objective of the class is to characterize
     * synonymous locations for this variant.
     *
     * @param gcomp
     *
     * A comparator for variants. This is used to sort synonyms according to
     * their order in the genome.
     *
     * @param calledvariants
     *
     * an object holding all variants called together with this one VcfEntry
     *
     * (Used for computing forgiveness factors)
     *
     */
    public ThesaurusSynonyms(ThesaurusEntry[] thesentries,
            VCFEntry varentry, GenomeInfo ginfo,
            VCFEntrySet calledvariants) {
        this.thesentries = thesentries;
        this.varentry = varentry;
        this.varalt = varentry.getAlt().split(",");
        this.ginfo = ginfo;
        this.gcomp = new GenomePositionComparator();
        this.calledvariants = calledvariants;
    }

    /**
     * Call this function after object initiation to find synonymous loci. This
     * will compute the synonymous positions using the variant position and the
     * thesaurus information.
     *
     * @param bamrecords
     *
     * Set this to null if synonyms are to be detected based only on variant
     * position, not location of variant within reads.
     *
     * @param hitcount
     *
     * minumum number of reads that should support a thesaurus link.
     *
     * @param hitproportion
     *
     * minumum proportion of reads carrying a variant that should be consistent
     * with a thesaurus link.
     *
     * @param tolerance
     *
     * Function allows a certain number of mismatches
     *
     * @param maxtolerance
     *
     * Allowed mismatches are increased if align/origin regions contain
     * variants.
     *
     * @param many
     *
     * Function tries to get synonyms with given conditions. If there are more
     * synonyms than "many", the function tries again with stricter settings
     *
     * @return list of synonyms (loci)
     *
     */
    public ArrayList<SNVPosition> findVariants(ThesaurusSAMRecord[] tbamrecords,
            int hitcount, double hitproportion, int tolerance, int maxtolerance, 
            int many, BitSet chrbitset) {

        // the answer, i.e. the locations synonymous with varentry will be calculated    
        // and stored in the synonyms array
        ArrayList<SNVPosition> synonyms;

        if (tbamrecords == null) {            
            synonyms = computeAllSynonymousPositionsNoBam();
        } else {            
            // use precomputed information to get a list of synonymous loci
            synonyms = computeAllSynonymousPositions(tbamrecords, hitcount, hitproportion,
                    tolerance, maxtolerance, maxtolerance, chrbitset);            
            if (synonyms != null && synonyms.size() >= many && maxtolerance - 1 >= 0) {                
                synonyms = computeAllSynonymousPositions(tbamrecords, hitcount, hitproportion,
                        tolerance, Math.max(tolerance, maxtolerance - 1), maxtolerance, chrbitset);
            }            
        }
        return synonyms;
    }

    /**
     * Used from outside the class to get a string of all synonyms (for vtf
     * file)
     *
     * @param synonyms
     * @return
     */
    public static String chainSynonyms(ArrayList<SNVPosition> synonyms, GenomeInfo ginfo) {
        StringBuilder sb = new StringBuilder();
        int slen = synonyms.size();
        for (int i = 0; i < slen; i++) {
            SNVPosition nowentry = synonyms.get(i);
            sb.append("\t").append(nowentry.getChr(ginfo)).append(":").append(nowentry.getPosition());
        }
        return sb.toString();
    }

    /**
     *
     * @param varentry
     * @param thesentries
     * @return
     *
     * an array with two elements. First element will contain an integer with
     * the number of synonyms.
     *
     * Second element will contain all the synonyms listed in a string, tab
     * separated
     *
     *
     *
     */
    private ArrayList<SNVPosition> computeAllSynonymousPositionsNoBam() {

        ArrayList<SNVPosition> synonyms = new ArrayList<SNVPosition>();
        if (thesentries == null) {
            return synonyms;
        }
        int slen = thesentries.length;
        if (slen == 0) {
            return synonyms;
        }

        // get a first look at all the synonyms
        SNVPosition[] tempsynonyms = new SNVPosition[slen];
        for (int i = 0; i < slen; i++) {
            tempsynonyms[i] = getSynonymousLocus(thesentries[i]);
        }

        // sort the arrays, this will enable checking for duplicates
        Arrays.sort(tempsynonyms, gcomp);

        // look through the synonyms and record only the unique positions
        SNVPosition cursynonym = tempsynonyms[0];
        synonyms.add(cursynonym);
        for (int i = 1; i < slen; i++) {
            SNVPosition nowsynonym = tempsynonyms[i];
            if (gcomp.compare(cursynonym, nowsynonym) != 0) {
                synonyms.add(nowsynonym);
                cursynonym = nowsynonym;
            }
        }
        return synonyms;
    }

    private int sumBooleanArray(boolean[] aa) {
        int ans = 0;
        for (int i = 0; i < aa.length; i++) {
            if (aa[i]) {
                ans++;
            }
        }
        return ans;
    }

    /**
     * Use information on the one variant together with thesentries and
     * bamrecords to compute possible alternate locations for the variant
     *
     * At finish, the synonms arraylist is filled with synonyms for the base
     * variant.
     *
     * @param tbamrecords
     *
     * alignment records that cover the variant on interest
     *
     * @param hitcount
     *
     * minimum number of reads that should support a thesaurus link
     *
     * @param hitproportion
     *
     * require a certain fraction of reads carrying the variant to be consistent
     * with the thesaurus link
     *
     * @param tolerance
     *
     * allow a certain number of mismatches/errors (not counting anchors)
     *
     * @param maxtolerance
     *
     * when regions nearby variant or nearby synonym have called variants, the
     * tolerance is increased.
     *
     * @return
     *
     * a list of synonymous loci for the variant.
     *
     * The list can be empty if there are not detected links.
     *
     * The output can be null if the class cannot classify this variant because
     * the reads are too complex (too many mismatches, or containing indels,
     * etc.)
     *
     */
    private ArrayList<SNVPosition> computeAllSynonymousPositions(ThesaurusSAMRecord[] tbamrecords,
            int hitcount, double hitproportion, int tolerance, int maxtolerance, int maxerrors, BitSet chrbitset) {

        // easy cases first
        if (thesentries == null) {
            return new ArrayList<SNVPosition>(2);
        }
        int slen = thesentries.length;
        if (slen == 0) {
            return new ArrayList<SNVPosition>(2);
        }

        // find the bases at the variant locus
        // e.g. if all reads show a variant A->C, the array should hold all Cs.        
        byte[] basesAtPosition = computeBasesAtVariantPosition(tbamrecords);
        int numWithVariant = getNumReadsWithVariant(basesAtPosition);

        // require that a certain number of reads agree on a thesaurus link
        // set to 0 to have no threshold, or higher for more stringency        
        int minAmbiguous = Math.max(hitcount, (int) Math.floor(hitproportion * numWithVariant));

        // check that there are sufficiently many reads that will meet the requirement
        // of no indels, mismatches etc.
        boolean[] withVarNotDifficult = getNonDifficultReads(basesAtPosition, tbamrecords, maxerrors);
        int withVarNotDifficultCount = sumBooleanArray(withVarNotDifficult);
        if (withVarNotDifficultCount < minAmbiguous) {
            return null;
        }
        // check that mates are in thesaurus range
        boolean[] withMateInThesaurus = getMateInThesaurusReads(tbamrecords, chrbitset);

        // get adjustment of tolerances per read using number of variants at align position
        int[] alignTolerances = getAlignToleranceFactors(tbamrecords, tolerance, maxtolerance);

        // look at the thesaurus entries and the synonyms
        ArrayList<SNVPosition> tempsynonyms = new ArrayList<SNVPosition>(slen);
        for (int i = 0; i < slen; i++) {
            ThesaurusEntry nowthesentry = thesentries[i];
            // get the reads that support the variant
            boolean[] supporting = getSupportingReads(nowthesentry, basesAtPosition, withVarNotDifficult, withMateInThesaurus,
                    alignTolerances, tbamrecords, maxtolerance);
            int numInBam = sumBooleanArray(supporting);

            // require that at least a certain number of reads are consistent with a locus
            if (numInBam >= minAmbiguous) {
                if (isPositionAnAnchor(nowthesentry)) {
                    // the variant position matches an anchor. Other anchors could be variants.
                    ArrayList<SNVPosition> anchorsynonyms = getSynonymousLociByAnchor(nowthesentry, tbamrecords, supporting);
                    for (int k = 0; k < anchorsynonyms.size(); k++) {
                        tempsynonyms.add(anchorsynonyms.get(k));
                    }
                } else {
                    tempsynonyms.add(getSynonymousLocus(nowthesentry));
                }
            }
        }

        // perhaps none of the sites passed
        if (tempsynonyms.isEmpty()) {
            return new ArrayList<SNVPosition>(2);
        }

        // sort the arrays, this will enable checking for duplicates
        Collections.sort(tempsynonyms, gcomp);

        // look through the synonyms and record only the unique positions        
        slen = tempsynonyms.size();
        ArrayList<SNVPosition> synonyms = new ArrayList<SNVPosition>(slen);
        SNVPosition lastsynonym = tempsynonyms.get(0);
        if (gcomp.compare(lastsynonym, varentry) != 0) {
            synonyms.add(lastsynonym);
        }
        for (int i = 1; i < slen; i++) {
            SNVPosition nowsynonym = tempsynonyms.get(i);
            if (gcomp.compare(lastsynonym, nowsynonym) != 0 && gcomp.compare(nowsynonym, varentry) != 0) {
                synonyms.add(nowsynonym);
                lastsynonym = nowsynonym;
            }
        }
        return synonyms;
    }

    /**
     * This function detects a special case where a variant is an anchor point
     * in a thesaurus entry. Such cases are handled differently later.
     *
     * @param entry
     * @return
     */
    private boolean isPositionAnAnchor(ThesaurusEntry entry) {
        int numanchors = entry.getNumAnchors();
        for (int i = 0; i < numanchors; i++) {
            if (entry.getAnchor(i).alignPosition == varentry.getPosition()) {
                if (entry.getAnchor(i).alignAlt == varentry.getAlt().charAt(0)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * compute the number of reads that show the variant
     *
     * @param basesAtPosition
     * @return
     */
    private int getNumReadsWithVariant(byte[] basesAtPosition) {
        int ans = 0;
        int bsize = basesAtPosition.length;
        for (int i = 0; i < bsize; i++) {
            if (containsVariant(basesAtPosition[i], varalt)) {
                ans++;
            }
        }
        return ans;
    }

    /**
     *
     * get markers that indicate whether reads are not difficult to process.
     *
     *
     *
     * @param basesAtPosition
     * @param tbamrecords
     * @param tolerance
     * @param maxtolerance
     * @return
     */
    private boolean[] getNonDifficultReads(byte[] basesAtPosition, ThesaurusSAMRecord[] tbamrecords, int maxtolerance) {

        boolean[] supporting = new boolean[tbamrecords.length];

        int bsize = basesAtPosition.length;
        for (int i = 0; i < bsize; i++) {
            if (!tbamrecords[i].hasIndelOrSplice() && containsVariant(basesAtPosition[i], varalt)
                    && tbamrecords[i].getNumMismatches() <= (maxtolerance + 1)) {
                supporting[i] = true;
            }
        }

        return supporting;
    }

    private boolean[] getMateInThesaurusReads(ThesaurusSAMRecord[] tbamrecords, BitSet chrbitset) {

        boolean[] supporting = new boolean[tbamrecords.length];

        int bsize = tbamrecords.length;
        for (int i = 0; i < bsize; i++) {
            ThesaurusSAMRecord record = tbamrecords[i];
            int nowchr = record.getReferenceIndex();
            int matechr = record.getMateReferenceIndex();

            if (nowchr != matechr) {
                // when the mate is unmapped or on different chromosome, return true
                // that way the other function will use the read for synonym finding.
                supporting[i] = true;
            } else {
                // when the mate is mapped, make sure one of the mate overlaps with a set bit in the chrbitset
                int matestart = record.getMateStart() - 1;
                int mateend = record.getMateStart() + record.getReadLength() - 1;
                mateend = Math.min(mateend, chrbitset.size());
                supporting[i] = chrbitset.get(matestart) || chrbitset.get(mateend);
            }
        }

        return supporting;
    }

    /**
     * This function uses the information given in the arguments, plus it also
     * relies on parameters defined for the class: callvariants.
     *
     * @param thesentry
     *
     * @return
     *
     * the number of reads that overlap with the thesaurus range and that
     * display the variant described by this instantiation of the class
     *
     */
    private boolean[] getSupportingReads(ThesaurusEntry thesentry,
            byte[] basesAtPosition, boolean[] process, boolean[] mateok, int[] tolerances, ThesaurusSAMRecord[] tbamrecords,
            int maxtolerance) {

        boolean[] supporting = new boolean[tbamrecords.length];
        boolean thesentryoverlaps = thesentry.isAlignOriginOverlapping();

        // loop over the reads
        int bsize = basesAtPosition.length;
        for (int i = 0; i < bsize; i++) {
            // the record must fully be in the range of the thesaurus entry
            if (isRecordInThesaurusRange(tbamrecords[i], thesentry)) {

                // look at records that can be processed (no indels, mismatches, etc.)
                // and further where the records overlap with thesaurus range                
                if (process[i] && isMateInThesaurusRange(tbamrecords[i], thesentry) && (mateok[i] || thesentryoverlaps)) {
                    // if the entry has an anchor, make sure that the read holds the anchor
                    int herepenalties = getNumAnchorsInRecordRange(tbamrecords[i], thesentry);
                    if (herepenalties > 0) {
                        // from here need to check if the number of mismatches makes this a reasonable
                        // candidate for an alternate site. The tolerate factor will be adjusted
                        // by the number of called variants at the alignment and putative origin regions.
                        int nowtolerance = tolerances[i];
                        if (maxtolerance > nowtolerance) {
                            nowtolerance = computeOriginToleranceFactor(thesentry, tbamrecords[i], nowtolerance, maxtolerance);
                        }
                        if (herepenalties > nowtolerance) {
                            if (countAnchorsInRecord(tbamrecords[i], thesentry) >= (herepenalties - nowtolerance)) {
                                supporting[i] = true;
                            }
                        } else {
                            supporting[i] = true;
                        }
                    } else {
                        supporting[i] = true;
                    }
                }
            }
        }
        return supporting;
    }

    /**
     * adjusts the a tolerance factor by the number of variants called in the
     * record
     *
     * Actually, this is not simplified to give an array of constants. Consider
     * removing altogether?
     *
     * @param tbamrecords
     *
     * array of reads
     *
     * @param tolerance
     *
     * the basic tolerance level
     *
     * @param maxtolerance
     *
     * maximum tolerance values
     *
     *
     * @return
     *
     * an array of integers, one number per read in the array
     *
     * numbers will always be greater or equal to tolerance and less than or
     * equal to maxtolerance
     *
     */
    private int[] getAlignToleranceFactors(ThesaurusSAMRecord[] tbamrecords, int tolerance, int maxtolerance) {
        int tbamlen = tbamrecords.length;
        int[] ans = new int[tbamlen];

        for (int i = 0; i < tbamlen; i++) {
            ans[i] = tolerance;
        }

        return ans;
    }

    /**
     * Scans the thesaurus entry in the interval covered by the read. returns
     * the number of anchors that are in that range.
     *
     * This function does not look at the actual bases recorded within the read.
     *
     * @param record
     * @param thesentry
     * @return
     */
    private int getNumAnchorsInRecordRange(ThesaurusSAMRecord record, ThesaurusEntry thesentry) {

        // find out where the read is aligned
        int recstart = record.getAlignmentStart();
        int recend = record.getAlignmentEnd();

        // compare read coordinates with all the anchors, and count overlaps
        int numanchors = thesentry.getNumAnchors();
        int ans = 0;
        for (int i = 0; i < numanchors; i++) {
            ThesaurusAnchor anchor = thesentry.getAnchor(i);
            if (anchor.alignPosition >= recstart && anchor.alignPosition <= recend) {
                ans++;
            }
        }

        return ans;
    }

    /**
     * Count the number of anchors declared in the thesaurus entry that the
     * record actually contains
     *
     * @param record
     * @param thesentry
     * @return
     */
    private int countAnchorsInRecord(ThesaurusSAMRecord record, ThesaurusEntry thesentry) {

        int numanchors = thesentry.penalty;
        byte[] bases = record.getReadBases();
        int baseslen = bases.length;

        int count = 0;
        for (int i = 0; i < numanchors; i++) {
            ThesaurusAnchor anch = thesentry.getAnchor(i);
            int nowpos = anch.alignPosition;
            byte nowalt = (byte) anch.alignAlt;

            int posinrec = record.getPositionInRecord(nowpos);

            if (posinrec >= 0 && posinrec < baseslen) {
                if (bases[posinrec] == nowalt) {
                    count++;
                }
            }
        }

        return count;
    }

    /**
     * The Synonyms object has access to the calledvariants object, which is a
     * list of variants. This function looks through those variants and returns
     * the number of items present in the given range.
     *
     * @param thesentry
     * @param record
     * @param tolerance
     * @param maxtolerance
     * @return
     *
     * returns a value that is at least equal to the tolerance integer, plus the
     * numbers of variants called in either the origin or alternative. If this
     * number exceed maxtolerance, the function returns maxtolerance instead.
     *
     */
    private int computeOriginToleranceFactor(ThesaurusEntry thesentry, ThesaurusSAMRecord record, int tolerance, int maxtolerance) {

        int rs = record.getAlignmentStart();
        int re = record.getAlignmentEnd();
        int ts = thesentry.alignStart;

        // to look up called variants near the possible alternate locus, need to compute the 
        // the possible alternate alignment start and end positions
        int altStart = 0, altEnd = 0;
        if (thesentry.alignStrand == '+') {
            altStart = thesentry.originStart + (rs - ts);
            altEnd = altStart + (re - rs);
        } else {
            altEnd = thesentry.originEnd - (rs - ts);
            altStart = altEnd - (re - rs);
        }

        int numOnOriginInterval = calledvariants.getNumberInInterval(thesentry.getOriginChr(), altStart, altEnd);

        return Math.min(maxtolerance, tolerance + numOnOriginInterval);
    }

    /**
     * checks if the read is entirely in the range described by the thesaurus.
     *
     *
     *
     * @param record
     * @param thesentry
     * @return
     */
    private boolean isRecordInThesaurusRange(ThesaurusSAMRecord record, ThesaurusEntry thesentry) {
        int rs = record.getAlignmentStart();
        int re = record.getAlignmentEnd();
        int ts = thesentry.alignStart;
        int te = thesentry.alignEnd;

        return (rs >= ts && rs <= te && re >= ts && re <= te);
    }

    /**
     * checks if the mate of this read is in the range described by the
     * thesaurus
     *
     * since not all the mate alignment information is available in a SAMRecord
     * or ThesaurusSAMRecord this function is approximate (see code)
     *
     * @param record
     * @param thesentry
     * @return
     */
    private boolean isMateInThesaurusRange(ThesaurusSAMRecord record, ThesaurusEntry thesentry) {

        if (5 > 3) {
            return true;
        }

        // for unpaired data or mate-unmapped reads, return true
        // this will make the rest of the class determine synonyms for the read
        if (!record.hasMate()) {
            return true;
        }

        // match chromosomes
        int nowchr = record.getReferenceIndex();
        int matechr = record.getMateReferenceIndex();

        if (nowchr != matechr) {
            // when the mate is unmapped or on different chromosome, return true
            // that way the other function will use the read for synonym finding.
            return true;
        }

        int reclen = record.getReadLength();
        // get start position of the mate
        int matestart = record.getMateStart();
        // the mate end position is not easily available
        // just estimate it using the readlength
        int mateend = matestart + reclen;

        // for mates, maybe don't be too strict about boundaries, so arbitrarily extend the 
        // the thesaurus regions by half a readlength       

        // return true if mate overlaps the align position
        int ts = thesentry.alignStart - (reclen / 2);
        int te = thesentry.alignEnd + (reclen / 2);
        if (matestart >= ts && matestart <= te && mateend >= ts && mateend <= te) {
            return true;
        }

        // return true also if the mate overlaps the origin position
        int os = thesentry.originStart - (reclen / 2);
        int oe = thesentry.originEnd + (reclen / 2);        
        if (thesentry.alignChrIndex == thesentry.originChrIndex && (matestart >= os && matestart <= oe && mateend >= os && mateend <= oe)) {
            return true;
        }

        // the default thing to do is to reject the link
        return false;
    }

    /**
     * calculate basesAtPosition
     *
     * @param bamrecords
     *
     * @return
     *
     */
    private byte[] computeBasesAtVariantPosition(ThesaurusSAMRecord[] tbamrecords) {
        int blen = tbamrecords.length;
        byte[] ans = new byte[blen];
        for (int i = 0; i < blen; i++) {
            int temp = tbamrecords[i].getPositionInRecord(varentry.getPosition());
            if (temp < 0) {
                ans[i] = 'N';
            } else {
                ans[i] = tbamrecords[i].getReadBases()[temp];
            }
        }
        return ans;
    }

    /**
     * compare the base found in a record with a selection of detected variants.
     * If one of the variants matches the base, the return value is true.
     *
     * @param baseInRecord
     *
     * @param alt
     *
     * @return
     *
     */
    private boolean containsVariant(byte baseInRecord, String[] alt) {

        if (baseInRecord == 'N') {
            return false;
        }

        int altsize = alt.length;
        for (int i = 0; i < altsize; i++) {
            if (alt[i].length() == 1 && baseInRecord == (byte) alt[i].charAt(0)) {
                return true;
            }
        }

        return false;
    }

    /**
     *
     * @param varentry
     * @param thesentry
     * @return
     */
    private SNVPosition getSynonymousLocus(ThesaurusEntry thesentry) {

        int offset = 0;
        if (thesentry.alignStrand == '+') {
            offset = varentry.getPosition() - thesentry.alignStart;
        } else if (thesentry.alignStrand == '-') {
            offset = thesentry.alignEnd - varentry.getPosition();
        }

        char nowref, nowalt;

        // record if the synonym is expected to be of the same kind as the
        // the base varentry, or its complement, based on strand item in the thesaurus entry
        if (thesentry.alignStrand == '+') {
            nowalt = varentry.getAlt().charAt(0);
            nowref = varentry.getRef().charAt(0);
        } else {
            nowalt = SequenceComplementer.complement(varentry.getAlt().charAt(0));
            nowref = SequenceComplementer.complement(varentry.getRef().charAt(0));
        }
        SNVPosition entry = new SNVPosition(thesentry.getOriginChr(), thesentry.originStart + offset, nowref, nowalt, ginfo);

        return entry;
    }

    /**
     *
     * given an entry and records, produce a list of anchors that are covered by
     * the reads.
     *
     *
     * @param thesentry
     * @param records
     * @param supporting
     *
     * only considers a subset of reads for which the supporting flag is set to
     * true
     *
     * @return list of loci
     *
     *
     */
    private ArrayList<SNVPosition> getSynonymousLociByAnchor(ThesaurusEntry thesentry, ThesaurusSAMRecord[] records, boolean[] supporting) {

        ArrayList<SNVPosition> ans = new ArrayList<SNVPosition>();

        // get interval defined by those reads that support the variant
        int recstart = Integer.MAX_VALUE;
        int recend = Integer.MIN_VALUE;
        for (int i = 0; i < records.length; i++) {
            if (supporting[i]) {
                recstart = Math.min(recstart, records[i].getAlignmentStart());
                recend = Math.max(recend, records[i].getAlignmentEnd());
            }
        }

        // count the anchors in range
        int anchorsInRange = 0;
        int numanchors = thesentry.getNumAnchors();
        for (int i = 0; i < numanchors; i++) {
            ThesaurusAnchor anchor = thesentry.getAnchor(i);
            if (anchor.alignPosition >= recstart && anchor.alignPosition <= recend) {
                anchorsInRange++;
            }
        }

        // compare read coordinates with anchors
        for (int i = 0; i < numanchors; i++) {
            ThesaurusAnchor anchor = thesentry.getAnchor(i);
            if (anchor.alignPosition >= recstart && anchor.alignPosition <= recend) {
                String newchr = thesentry.getOriginChr();
                int newpos = anchor.originPosition;
                char newref, newalt;
                if (varentry.getPosition() == anchor.alignPosition && anchorsInRange == 1) {
                    // the variant is at this anchor and there are no other anchors                    
                    // Record the synonym location, but use the variant's declared ref and alt bases
                    if (thesentry.alignStrand == '+') {
                        newref = anchor.originRef;
                        newalt = varentry.getAlt().charAt(0);
                    } else {
                        newref = anchor.originRef;
                        newalt = SequenceComplementer.complement(varentry.getAlt().charAt(0));
                    }
                } else {
                    // generic case
                    if (thesentry.alignStrand == '+') {
                        newalt = anchor.originAlt;
                        newref = anchor.originRef;
                    } else {
                        newalt = anchor.originAlt;
                        newref = anchor.originRef;
                    }
                }
                SNVPosition newsynonym = new SNVPosition(newchr, newpos, newref, newalt, ginfo);
                if (newsynonym.getAlt() != newsynonym.getRef()) {
                    ans.add(newsynonym);
                }
            }
        }
        return ans;
    }
}
