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

import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePositionInterface;

/**
 *
 * Extends SNVPosition class in that it stores counts for ATCGN bases at a
 * locus.
 *
 * Internally, the counts are actually stored as doubles, not ints.
 *
 * @author tkonopka
 */
public class SNVPositionDetails extends SNVPosition implements GenomePositionInterface {

    long countA, countT, countC, countG, countN;

    public SNVPositionDetails(SNVPosition entry) {
        super(entry);
        countA = 0L;
        countT = 0L;
        countC = 0L;
        countG = 0L;
        countN = 0L;
    }

    @Override
    public String toString(GenomeInfo ginfo) {
        StringBuilder sb = new StringBuilder();
        sb.append(super.getChr(ginfo)).append("\t").append(super.getPosition())
                .append("\t").append(ref)
                .append("\t").append(alt)
                .append("\t").append(countA)
                .append("\t").append(countT)
                .append("\t").append(countC)
                .append("\t").append(countG)
                .append("\t").append(countN);
        return sb.toString();
    }

    public String getHeader() {
        return "chr\tposition\tA\tT\tC\tG\tN\n";
    }

    public void incrementA() {
        countA++;
    }

    public void incrementT() {
        countT++;
    }

    public void incrementC() {
        countC++;
    }

    public void incrementG() {
        countG++;
    }

    public void incrementN() {
        countN++;
    }

    /**
     * add counts collected at the given SNVPosition to the current object.
     *
     * The addition is only done when the two SNVPositionDetails objects are
     * consistent, i.e. describe the same variant or a complement variant.
     *
     *
     * @param dd
     */
    public void incrementCounts(SNVPositionDetails dd) {
        boolean ok = false;
        boolean complement = false;

        // check if the variants at the two positions are consistent or complements
        byte refbyte = (byte) this.ref;
        byte ddrefbyte = (byte) dd.ref;
        byte altbyte = (byte) this.alt;
        byte ddaltbyte = (byte) dd.alt;

        if (refbyte == ddrefbyte && altbyte == ddaltbyte) {
            ok = true;
        } else {
            if (refbyte == complement(ddrefbyte) && altbyte == complement(ddaltbyte)) {
                ok = true;
                complement = true;
            }
        }

        if (ok) {
            if (complement) {
                countA += dd.countT;
                countT += dd.countA;
                countC += dd.countG;
                countG += dd.countC;
            } else {
                countA += dd.countA;
                countT += dd.countT;
                countC += dd.countC;
                countG += dd.countG;
            }
            countN += dd.countN;
        }
    }

    public long getCountATCG() {
        return countA + countT + countC + countG;
    }

    public long getCountAT() {
        return countA + countT;
    }

    public long getCountCG() {
        return countC + countG;
    }

    public long getCountA() {
        return countA;
    }

    public long getCountT() {
        return countT;
    }

    public long getCountC() {
        return countC;
    }

    public long getCountG() {
        return countG;
    }

    public long getCountN() {
        return countN;
    }

    public long getCountRef() {
        byte refbyte = (byte) this.ref;
        return getCountFromByte(refbyte);
    }

    public long getCountAlt() {
        byte altbyte = (byte) this.alt;
        return getCountFromByte(altbyte);
    }

    private long getCountFromByte(byte b) {
        switch (b) {
            case 'A':
                return countA;
            case 'T':
                return countT;
            case 'C':
                return countC;
            case 'G':
                return countG;
            default:
                return countN;
        }
    }

    private byte complement(byte b) {
        switch (b) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            default:
                return 'N';
        }
    }
}
