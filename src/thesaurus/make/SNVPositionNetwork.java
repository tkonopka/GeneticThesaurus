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
import thesaurus.util.VCFEntrySet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePositionComparator;

/**
 *
 * @author tkonopka
 */
public class SNVPositionNetwork {

    private final GenomePositionComparator gcomp;
    private final GenomeInfo ginfo;
    private final HashMap<String, SNVPositionNode> network;
    private boolean isok = false;

    public SNVPositionNetwork(VCFEntrySet calledvars, GenomeInfo ginfo) {
        // save the comparator
        this.ginfo = ginfo;
        this.gcomp = new GenomePositionComparator();
        int numcalled = calledvars.size();
        network = new HashMap<String, SNVPositionNode>(2*numcalled);
        for (int i = 0; i < numcalled; i++) {            
            SNVPositionNode newnode = new SNVPositionNode(new SNVPosition(calledvars.getVariant(i), ginfo), true);            
            String nodestring = newnode.getSNVPositionString();            
            network.put(nodestring, newnode);            
        }
    }

    public boolean isIsok() {
        return isok;
    }

    public void addLinks(SNVPosition from, ArrayList<SNVPosition> to) {
        String fromstring = from.toString();

        // get from object
        SNVPositionNode fromnode = network.get(fromstring);
        if (fromnode == null) {
            System.out.println("From node does not exist");
            return;
        }

        // check if all to objects exist - create them if necessary.
        int tolength = to.size();
        ArrayList<SNVPositionNode> tonodes = new ArrayList<SNVPositionNode>(to.size());
        for (int i = 0; i < tolength; i++) {
            String tostring = to.get(i).toString();
            SNVPositionNode tonode = network.get(tostring);
            if (tonode == null) {
                tonode = new SNVPositionNode(new SNVPosition(to.get(i)), false);
                network.put(tostring, tonode);                
            }
            tonodes.add(tonode);
        }
        
        // add link information 
        fromnode.setNeighbors(tonodes, gcomp);

    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, SNVPositionNode> entry : network.entrySet()) {
            sb.append(entry.getValue().toString()).append("\n");
        }
        return sb.toString();
    }

    public int size() {
        return network.size();
    }
    
    public ArrayList<SNVPosition> getAllNodes() {
        ArrayList<SNVPosition> ans = new ArrayList<SNVPosition>(network.size());
        for (Map.Entry<String, SNVPositionNode> entry: network.entrySet()) {
            SNVPositionNode nownode = entry.getValue();
            ans.add(nownode.snv);
        }
        Collections.sort(ans, gcomp);
        return ans;
    }
    
    public ArrayList<SNVPosition> getNeighborNodes(SNVPosition from) {        
        SNVPositionNode nownode = network.get(from.toString());
        int numnei = nownode.getNumNeighbors();        
        ArrayList<SNVPosition> ans = new ArrayList<SNVPosition>(numnei);
        for (int i=0; i<numnei; i++) {
            SNVPositionNode nownei = nownode.getNeighbor(i);
            ans.add(new SNVPosition(nownei.snv));
        }
        return ans;        
    }
}
