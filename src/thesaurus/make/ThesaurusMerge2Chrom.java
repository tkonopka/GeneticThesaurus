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
package thesaurus.make;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Merges multiple bam files, splitting the output into one file per chromosome
 *
 * @author tkonopka
 *
 */
public class ThesaurusMerge2Chrom implements Runnable {

    private boolean isok = false;
    private String prefix = null;
    private ArrayList<File> infiles = new ArrayList<File>();
    private HashMap<String, SAMFileWriter> chrombams = new HashMap<String, SAMFileWriter>();

    private void printMerge2ChromHelp() {
        System.out.println("GeneticThesaurus merge2chrom: merge bam files");
        System.out.println();
        System.out.println("Usage: outprefix infile1 infile2 ...");
        System.out.println();
    }

    public ThesaurusMerge2Chrom(String[] args) {
        if (args == null || args.length == 0) {
            printMerge2ChromHelp();
            return;
        }

        // get the input files
        prefix = args[0];
        for (int i = 1; i < args.length; i++) {
            infiles.add(new File(args[i]));
        }

        isok = true;

        // check that all files are readable
        for (int i = 0; i < infiles.size(); i++) {
            File nowfile = infiles.get(i);
            if (!nowfile.canRead()) {
                System.out.println("cannot read file: " + nowfile.getAbsolutePath());
                isok = false;
            }
        }

    }

    @Override
    public void run() {
        if (!isok || infiles.isEmpty()) {
            return;
        }

        // first make lots of bamf files
        setupChromBamWriters(infiles.get(0));

        // process each input file in turn
        for (int i = 0; i < infiles.size(); i++) {
            splitIntoChroms(infiles.get(i));
        }

        // close all SAMwriter
        closeChromBamWriters();
    }

    /**
     * closes all the writers defined in chrombams
     *
     */
    private void closeChromBamWriters() {
        for (Map.Entry<String, SAMFileWriter> entry : chrombams.entrySet()) {
            SAMFileWriter nowwriter = entry.getValue();
            nowwriter.close();
        }
    }

    /**
     * Creates a collection of SAMFileWriters, one for each chromosome.
     *
     * The definition of the chromosomes is taken from the header of a template
     * bamfile.
     *
     * @param template
     */
    private void setupChromBamWriters(File template) {

        SAMFileReader inputSam = new SAMFileReader(template);
        SAMFileHeader inputHeader = inputSam.getFileHeader();
        ArrayList<SAMSequenceRecord> samchroms;
        samchroms = new ArrayList<SAMSequenceRecord>(inputSam.getFileHeader().getSequenceDictionary().getSequences());

        // create another list of just chromosome names
        ArrayList<String> chroms = new ArrayList<String>();
        for (int i = 0; i < samchroms.size(); i++) {
            chroms.add(samchroms.get(i).getSequenceName());
        }
        chroms.add("unmapped");

        // create a filewriter for each chromosome
        for (int i = 0; i < chroms.size(); i++) {
            String nowchrom = chroms.get(i);
            SAMFileHeader nowheader = inputHeader.clone();
            nowheader.addComment("Merged data using GeneticThesaurus merge2chrom: chromosome " + nowchrom);
            SAMFileWriter nowWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
                    nowheader, true, new File(prefix + "." + nowchrom + ".bam"));
            chrombams.put(nowchrom, nowWriter);
        }

    }

    /**
     * copy reads from one bam file and put them into files by chromosome.
     *
     * @param nowbam
     */
    private void splitIntoChroms(File nowbam) {
        SAMFileReader inputSam = new SAMFileReader(nowbam);

        // write each aligned record into its separate file, or into "unmapped"
        for (final SAMRecord record : inputSam) {
            String nowchrom = record.getReferenceName();
            SAMFileWriter nowwriter = chrombams.get(nowchrom);
            if (nowwriter == null) {
                nowwriter = chrombams.get("unmapped");
            }
            nowwriter.addAlignment(record);
        }

        inputSam.close();
    }
}
