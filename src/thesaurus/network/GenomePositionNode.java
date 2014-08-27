/*
 * Copyright 2014 Tomasz Konopka.
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
package thesaurus.network;

import java.util.ArrayList;
import java.util.Collections;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePosition;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.genome.GenomePositionInterface;

/**
 * Meant to work with GenomePositionNetwork. Stores a node of a graph. The node
 * represents a genome position and the edges are Thesarus links to other
 * positions.
 *
 *
 * @author tkonopka
 */
public class GenomePositionNode implements GenomePositionInterface {

    private final GenomePosition gpos;
    private boolean called;
    // list of synonyms
    private final ArrayList<GenomePosition> synonyms = new ArrayList<GenomePosition>();
    // list of priority synonyms, a.k.a. synonyms among the called variants
    private final ArrayList<GenomePosition> prioritysynonyms = new ArrayList<GenomePosition>();

    public GenomePositionNode(GenomePosition gpos, boolean called) {
        this.gpos = gpos;
        this.called = called;
    }

    public boolean hasNeighbor(GenomePosition neighbor, GenomePositionComparator gcomp) {
        return Collections.binarySearch(synonyms, neighbor, gcomp) >= 0;
    }

    public boolean hasPriorityNeighbor(GenomePosition neighbor, GenomePositionComparator gcomp) {
        return Collections.binarySearch(prioritysynonyms, neighbor, gcomp) > 0;
    }

    public boolean isCalled() {
        return called;
    }

    public void setCalled(boolean called) {
        this.called = called;
    }

    /**
     * Adds links to neighbors
     *
     * @param neighbors
     * @param gcomp
     */
    public void addNeighbors(ArrayList<GenomePosition> neighbors, GenomePositionComparator gcomp) {
        int nsize = neighbors.size();
        ArrayList<GenomePosition> newneighbors = new ArrayList<GenomePosition>(nsize);
        for (int i = 0; i < nsize; i++) {
            if (!hasNeighbor(neighbors.get(i), gcomp)) {
                newneighbors.add(neighbors.get(i));
            }
        }
        this.synonyms.addAll(newneighbors);
        Collections.sort(synonyms, gcomp);
    }

    public int getNumNeighbors() {
        return synonyms.size();
    }

    public int getNumPriorityNeighbors() {
        return prioritysynonyms.size();
    }

    public GenomePosition getNeighbor(int index) {
        return synonyms.get(index);
    }

    public GenomePosition getPriorityNeighbor(int index) {
        return prioritysynonyms.get(index);
    }

    public void addNeighbor(GenomePosition neighbor, GenomePositionComparator gcomp) {
        if (!hasNeighbor(neighbor, gcomp)) {
            synonyms.add(neighbor);
            Collections.sort(synonyms, gcomp);
        }
    }

    public void addPriorityNeighbor(GenomePosition neighbor, GenomePositionComparator gcomp) {
        if (!hasPriorityNeighbor(neighbor, gcomp)) {
            prioritysynonyms.add(neighbor);
            Collections.sort(synonyms, gcomp);
        }
    }

    @Override
    public String getChr(GenomeInfo ginfo) {
        return gpos.getChr(ginfo);
    }

    @Override
    public int getPosition() {
        return gpos.getPosition();
    }

    public String getGenomePositionString(GenomeInfo ginfo) {
        return gpos.toString(ginfo);
    }

    @Override
    public String toString(GenomeInfo ginfo) {
        StringBuilder sb = new StringBuilder();
        int ssize = synonyms.size();
        int pssize = prioritysynonyms.size();
        sb.append(gpos.toString(ginfo)).append("\t").append(called).append("\t").append(ssize).append("\t").append(pssize);
        return sb.toString();
    }

    @Override
    public int getChrIndex() {
        return gpos.getChrIndex();
    }
}
