/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thesaurus.util;

import java.util.Comparator;

/**
 * a comparator for SNVPosition objects, i.e. genomic positions together with single-
 * nucleotide substitutions.
 * 
 * @author tkonopka
 */
public class SNVPositionComparator implements Comparator {

    /**
     * 
     * @param o1
     * @param o2
     * @return 
     * 
     * First set of comparisons are made for chromosome and position
     * 
     * negative number if o1 is before o2. 
     * (i.e. if chromosome o1 is before chromosome o2, or if position is smaller on same chromosome)
     *
     * positive num if o2 is before o1.
     * 
     * when the two position are exactly the same, the order of alt alleles is 
     * the natural char order
     * 
     * 
     */
    @Override
    public int compare(Object o1, Object o2) {
        SNVPosition e1 = (SNVPosition) o1;
        SNVPosition e2 = (SNVPosition) o2;

        int c1 = e1.getChrIndex();
        int c2 = e2.getChrIndex();
                
        // if the chromosomes are different, use the preset order
        if (c1 == -1 || c2 == -1) {
            if (c1 == -1) {
                return -1;
            } else {
                return 1;
            }
        } else if (c1 != c2) {
            return c1 - c2;
        }

        // if the chromosomes are the same compare the positions      
        int posdifference = e1.getPosition() - e2.getPosition();
        if (posdifference!=0) {
            return posdifference;
        }             
        return e1.getAlt()-e2.getAlt();
    }
    
}
