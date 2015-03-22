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

import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePosition;
import jsequtils.variants.VCFEntry;

/**
 * A minimal class holding a genomic position (GenomicPosition) together with
 * variant information (ref/alt characters). This class works only for
 * single-nucleotide substitutions
 *
 * @author tkonopka
 */
public class SNVPosition extends GenomePosition {

    char ref, alt;

    /**
     *
     * @param entry
     */
    public SNVPosition(VCFEntry entry, GenomeInfo ginfo) {
        super(entry.getChr(), entry.getPosition(), ginfo);
        ref = entry.getRef().charAt(0);
        alt = entry.getAlt().charAt(0);
    }

    /**
     * defensive copy constructor
     *
     * @param entry
     */
    public SNVPosition(SNVPosition entry) {
        super(entry.getChrIndex(), entry.getPosition());
        this.ref = entry.getRef();
        this.alt = entry.getAlt();
    }

    /**
     * ab-initio constructor
     *
     * @param chr
     * @param position
     * @param ref
     * @param alt
     */
    public SNVPosition(String chr, int position, char ref, char alt, GenomeInfo ginfo) {
        super(chr, position, ginfo);
        this.ref = ref;
        this.alt = alt;
    }

    public char getRef() {
        return ref;
    }

    public char getAlt() {
        return alt;
    }

    /**
     * gives a text-based description of the SNV. Uses a representation of chromosomes
     * that is based on a numeric index, integer based position, and ref/alt pair. 
     * The chromosome and position is intentionally separated by ; rather than : 
     * as a reminder of the numeric index.
     * @return 
     */
    @Override
    public String toString() {
        return super.getChrIndex() + ";" + super.getPosition() + ":" + ref + ":" + alt;
    }

    @Override
    public String toString(GenomeInfo ginfo) {
        return super.getChr(ginfo) + ":" + super.getPosition() + ":" + ref + ":" + alt;
    }
}
