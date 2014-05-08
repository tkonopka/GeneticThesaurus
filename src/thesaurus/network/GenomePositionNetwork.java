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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomePosition;
import jsequtils.genome.GenomePositionComparator;
import jsequtils.variants.VcfEntry;

/**
 * A data structure that holds a graph of interconnected genomic positions. In contrast 
 * to a standard/simple graph, this one provides some functions are that specific to 
 * a Thesaurus graph (determine the number of clusters among called variants)
 * 
 * 
 * @author tkonopka
 */
public class GenomePositionNetwork {

    private GenomeInfo ginfo = null;
    private GenomePositionComparator gcomp = new GenomePositionComparator();
    private HashMap<String, GenomePositionNode> network = null;
    private boolean isok = false;

    public GenomePositionNetwork(File vcffile, File vtffile, File genome) {

        // load information about genome                 
        try {
            ginfo = new GenomeInfo(genome);                        
        } catch (Exception ex) {
            System.out.println("Failed to load genome info: " + ex.getMessage());
            return;
        }

        // load called positions from vcf file
        //System.out.println("Loading called variants");
        ArrayList<GenomePosition> vcfpos = new ArrayList<GenomePosition>(4096);
        try {
            BufferedReader br = BufferedReaderMaker.makeBufferedReader(vcffile);
            String s;
            while ((s = br.readLine()) != null) {
                if (!s.startsWith("#")) {
                    VcfEntry ve = new VcfEntry(s);
                    if (!ve.isIndel()) {
                        String[] stokens = s.split("\t", 4);
                        GenomePosition gpos = new GenomePosition(stokens[0], Integer.parseInt(stokens[1]), ginfo);
                        vcfpos.add(gpos);
                    }
                }
            }
            br.close();
        } catch (Exception ex) {
            System.out.println("Error loading from vcf file: " + ex.getMessage());
            return;
        }

        // create a network, first add the called nodes
        //System.out.println("Setting up initial network");
        network = new HashMap<String, GenomePositionNode>(vcfpos.size() * 2);
        for (int i = 0; i < vcfpos.size(); i++) {
            GenomePositionNode node = new GenomePositionNode(vcfpos.get(i), true);
            network.put(node.getGenomePositionString(ginfo), node);
        }

        // load the alternate position from vtf file
        //System.out.println("Adding neighbors");
        ArrayList<GenomePosition> vtfentries = new ArrayList<GenomePosition>(4096);
        try {
            BufferedReader br = BufferedReaderMaker.makeBufferedReader(vtffile);
            String s;
            while ((s = br.readLine()) != null) {
                if (!s.startsWith("#")) {
                    String[] stokens = s.split("\t|:");
                    GenomePosition gpos = new GenomePosition(stokens[0], Integer.parseInt(stokens[1]), ginfo);                    
                    ArrayList<GenomePosition> temp = new ArrayList<GenomePosition>(stokens.length / 2);
                    for (int i = 2; i < stokens.length; i += 2) {
                        GenomePosition gpos2 = new GenomePosition(stokens[i], Integer.parseInt(stokens[i + 1]), ginfo);                        
                        temp.add(gpos2);
                    }
                    network.get(gpos.toString(ginfo)).addNeighbors(temp, gcomp);
                }
            }
            br.close();
        } catch (Exception ex) {
            System.out.println("Error loading from vtf file: " + ex.getMessage());
            return;
        }

        //System.out.println("Connecting called variants");        
        int count = 0;
        for (Map.Entry<String, GenomePositionNode> entry : network.entrySet()) {
            GenomePositionNode node = entry.getValue();
            count += connectCalled(entry.getValue());
        }
        //System.out.println("done init " + count);

        isok = true;
    }

    public boolean isIsok() {
        return isok;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Map.Entry<String, GenomePositionNode> entry : network.entrySet()) {
            sb.append(entry.getValue().toString(ginfo)).append("\n");
        }
        return sb.toString();
    }

    /**
     * If node X has neighbors A, B, C and that C is called variant, make sure
     * that X records C as a priority neighbor and C records X as a priority
     * neighbor.
     *
     * @param node
     */
    private int connectCalled(GenomePositionNode node) {

        int numneighbors = node.getNumNeighbors();
        int count = 0;
        for (int i = 0; i < numneighbors; i++) {
            GenomePosition nowpos = node.getNeighbor(i);
            GenomePositionNode nownei = network.get(nowpos.toString(ginfo));
            if (nownei != null) {                
                node.addPriorityNeighbor(new GenomePosition(nownei.getChrIndex(), nownei.getPosition()), gcomp);
                nownei.addPriorityNeighbor(new GenomePosition(node.getChrIndex(), node.getPosition()), gcomp);
                count++;
            }
        }
        return count;
    }

    /**
     * counts distinct values in a hashmap by adding them into another hashmap
     * and getting the size of the resulting hashamp.
     *
     * @param clusterids
     * @return
     */
    private int countDistinctIds(HashMap<String, Integer> clusterids) {
        HashMap<Integer, Boolean> gotit = new HashMap<Integer, Boolean>(clusterids.size() * 2);
        for (Map.Entry<String, Integer> entry : clusterids.entrySet()) {
            gotit.put(entry.getValue(), true);
        }
        return gotit.size();
    }

    /**
     *
     * @param clusterids
     * @return
     *
     * the maximal cluster id assigned in the hashmap
     *
     */
    private boolean reduceClusterIds(HashMap<String, Integer> clusterids) {

        // count the number of distinct ids
        int beforecount = countDistinctIds(clusterids);

        for (Map.Entry<String, GenomePositionNode> entry : network.entrySet()) {
            GenomePositionNode node = entry.getValue();
            int numneighbors = node.getNumPriorityNeighbors();
            int nowcluster = clusterids.get(node.getGenomePositionString(ginfo));
            for (int i = 0; i < numneighbors; i++) {
                int neighborcluster = clusterids.get(node.getPriorityNeighbor(i).toString(ginfo));
                nowcluster = Math.min(nowcluster, neighborcluster);
            }
            // assign the minimal cluster number to all the neighbors
            clusterids.put(node.getGenomePositionString(ginfo), nowcluster);
            for (int i = 0; i < numneighbors; i++) {
                clusterids.put(node.getPriorityNeighbor(i).toString(ginfo), nowcluster);
            }
        }

        // find the maximal number assigned in the hashmap
        int aftercount = countDistinctIds(clusterids);

        return aftercount != beforecount;
    }

    /**
     * If two
     *
     * @return
     */
    @SuppressWarnings("empty-statement")
    public int countDisconnectedComponents() {


        // set up a table assigning cluster indexes and initialize each node to a separate cluster
        HashMap<String, Integer> clusters = new HashMap<String, Integer>(network.size() * 2);
        int count = 0;
        for (Map.Entry<String, GenomePositionNode> entry : network.entrySet()) {
            clusters.put(entry.getKey(), count);
            count++;
        }

        // loop reducing the cluster ids until they stabilize
        while (reduceClusterIds(clusters));

        return countDistinctIds(clusters);
    }

    @SuppressWarnings("empty-statement")
    public void outputClusterIds(OutputStream outstream) throws IOException {

        // set up a table assigning cluster indexes and initialize each node to a separate cluster
        HashMap<String, Integer> clusters = new HashMap<String, Integer>(network.size() * 2);
        int count = 0;
        for (Map.Entry<String, GenomePositionNode> entry : network.entrySet()) {
            clusters.put(entry.getKey(), count);
            count++;
        }

        // loop reducing the cluster ids until they stabilize
        while (reduceClusterIds(clusters));

        // transfer from hashmap into array and sort it to get correct order of entries
        GenomePositionNode[] gpnarray = new GenomePositionNode[network.size()];
        int i = 0;
        for (Map.Entry<String, GenomePositionNode> entry : network.entrySet()) {
            GenomePositionNode node = entry.getValue();
            gpnarray[i] = node;
            i++;
        }
        Arrays.sort(gpnarray, gcomp);

        // output a header and cluster ids
        outstream.write("chr\tposition\tcluster\n".getBytes());
        for (int j = 0; j < gpnarray.length; j++) {
            GenomePositionNode node = gpnarray[j];
            String temp = node.getGenomePositionString(ginfo);
            String s = node.getChr(ginfo) + "\t" + node.getPosition() + "\t" + clusters.get(temp) + "\n";
            outstream.write(s.getBytes());
        }

    }

    public int size() {
        return network.size();
    }
}
