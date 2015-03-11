/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thesaurus.util;

import jsequtils.genome.GenomePositionInterface;

/**
 * Interface that promises to return chr, position, refbase, altbase
 * 
 * @author tkonopka
 */
public interface SNVPositionInterface extends GenomePositionInterface {
    public char getRef();
    public char getAlt();
}
