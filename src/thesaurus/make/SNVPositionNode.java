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
package thesaurus.make;

import thesaurus.util.SNVPosition;
import java.util.ArrayList;
import java.util.Collections;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.genome.GenomePositionInterface;

/**
 *
 * @author tkonopka
 */
class SNVPositionNode implements GenomePositionInterface {

    final SNVPosition snv;
    private boolean called;
    // list of synonyms
    private ArrayList<SNVPositionNode> links = null;

    public SNVPositionNode(SNVPosition snv, boolean called) {
        this.snv = snv;
        this.called = called;
    }

    public boolean hasNeighbor(SNVPositionNode neighbor, GenomePositionComparator gcomp) {
        return Collections.binarySearch(links, neighbor, gcomp) >= 0;
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
    void setNeighbors(ArrayList<SNVPositionNode> neighbors, GenomePositionComparator gcomp) {
        if (links != null) {
            System.out.println("[Warning] position: " + snv.getPosition()
                    + "; ignored attempt to set thesaurus links (complex variant?)");
            return;
        }
        links = new ArrayList<SNVPositionNode>(neighbors.size());
        links.addAll(neighbors);
        Collections.sort(links, gcomp);
    }

    public int getNumNeighbors() {
        if (links == null) {
            return 0;
        }
        return links.size();
    }

    public SNVPositionNode getNeighbor(int index) {
        return links.get(index);
    }

    @Override
    public String getChr(GenomeInfo ginfo) {
        return snv.getChr(ginfo);
    }

    @Override
    public int getChrIndex() {
        return snv.getChrIndex();
    }

    @Override
    public int getPosition() {
        return snv.getPosition();
    }

    public String getSNVPositionString() {
        return snv.toString();
    }

    public String getSNVPositionString(GenomeInfo ginfo) {
        return snv.toString(ginfo);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        int ssize = links.size();
        sb.append(snv.toString()).append("\t").append(called).append("\t").append(ssize);
        return sb.toString();
    }

    @Override
    public String toString(GenomeInfo ginfo) {
        StringBuilder sb = new StringBuilder();
        int ssize = links.size();
        sb.append(snv.toString(ginfo)).append("\t").append(called).append("\t").append(ssize);
        return sb.toString();
    }
}
