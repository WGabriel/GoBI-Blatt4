import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

import net.sf.samtools.*;

public class Runner {
    public static void main(String[] args) {
        File gtf = null;
        File bam = null;
        File output = null;

        if (args.length != 6) {
            System.err.println("Required 6 parameters. Found: " + args.length);
            String consoleParameters = "";
            for (int i = 0; i < args.length; i++) {
                consoleParameters += " args[" + i + "]=" + args[i];
            }
            System.err.println("ConsoleParameters:" + consoleParameters);
        } else {
            // Console parameters are correct!
            for (int i = 0; i < args.length; i = i + 2) {
                switch (args[i]) {
                    case "-gtf":
                        gtf = new File(args[i + 1]);
                        continue;
                    case "-bam":
                        bam = new File(args[i + 1]);
                        continue;
                    case "-o":
                        output = new File(args[i + 1]);
                        continue;
                }
                // System.out.println("|User Input| "+i+args[i]);
            }
        }
        System.out.println("|UserInput|\ngtf:\t" + gtf.getAbsolutePath() + "\nbam:\t" + bam.getAbsolutePath()
                + "\noutput:\t" + output.getAbsolutePath());

        HashMap<String, Gene> genes = parseGene(gtf);

        SAMFileReader sam_reader = new SAMFileReader(bam, false);
        sam_reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> it = sam_reader.iterator();

        // Key = readName, Value = SAMRecord
        HashMap<String, SAMRecord> firstOfPairSAM = new HashMap<String, SAMRecord>();
        HashMap<String, SAMRecord> secondOfPairSAM = new HashMap<String, SAMRecord>();

        // int printToConsole = 5;
        // printToConsole = Integer.MAX_VALUE;
        int counterSAMRecordsInFile = 0;
        while (it.hasNext()) {
            SAMRecord sr = it.next();
            counterSAMRecordsInFile++;
            if (sr.getNotPrimaryAlignmentFlag()) {
                // Ignore, if read is secondary
                // System.err.println("Ignored: " + sr.getReadName() + " (secondary)");
            } else if (sr.getReadUnmappedFlag()) {
                // Ignore, if supplementary or unmapped
                // System.err.println("Ignored: " + sr.getReadName() + " (supplementary or
                // unmapped)");
            } else if (sr.getMateUnmappedFlag()) {
                // Ignore, if mate unmapped
                // System.err.println("Ignored: " + sr.getReadName() + " (mate unmapped)");
            } else if (sr.getReadNegativeStrandFlag() == sr.getMateNegativeStrandFlag()) {
                // Ignore, if mate is on same strand
                // System.err.println("Ignored: " + sr.getReadName() + " (mate is on same
                // strand)");
            } else if (!sr.getReferenceName().equals(sr.getMateReferenceName())) {
                // Ignore, if mate is not on same chromosome
                // System.err.println("Ignored: " + sr.getReadName() + " (mate is on not on same
                // chromosome)");
            } else {
                if (sr.getFirstOfPairFlag()) {
                    firstOfPairSAM.put(sr.getReadName(), sr);
                } else {
                    secondOfPairSAM.put(sr.getReadName(), sr);
                }
                // ---print SAMs to console---
                // if (printToConsole > 0) {
                // printToConsole--;
                // System.out.println("######################");
                // System.out.println("readName: " + sr.getReadName());
                // System.out.println("getNotPrimaryAlignmentFlag: " +
                // sr.getNotPrimaryAlignmentFlag());
                // System.out.println("getReferenceName: " + sr.getReferenceName());
                // System.out.println("getReadPairedFlag: " + sr.getReadPairedFlag());
                // System.out.println("getReadNegativeStrandFlag: " +
                // sr.getReadNegativeStrandFlag());
                // System.out.println("getFirstOfPairFlag: " + sr.getFirstOfPairFlag());
                // System.out.println("getMateUnmappedFlag: " + sr.getMateUnmappedFlag());
                // System.out.println("getMateReferenceName: " + sr.getMateReferenceName());
                // System.out.println("getMateAlignmentStart:" + sr.getMateAlignmentStart());
                // System.out.println("getReferenceIndex:" + sr.getReferenceIndex());
                // System.out.println("getReadString:" + sr.getReadString());
                // System.out.println("getSAMString:" + sr.getSAMString());
                // System.out.println("start: " + sr.getAlignmentStart() + " end: " +
                // sr.getAlignmentEnd());
                // for (AlignmentBlock ab : sr.getAlignmentBlocks()) {
                // int ref_s = ab.getReferenceStart();
                // int read_s = ab.getReadStart();
                // int ref_end = ref_s + ab.getLength();
                // int read_end = read_s + ab.getLength();
                // System.out.println(
                // "AlignmentBlock reference: " + ref_s + "-" + ref_end + " read: " + read_s +
                // "-" + read_end);
                // }
                // ---END: print SAMs to console---
                // }
            }
        }

        System.out.println("FirstOfPair size: " + firstOfPairSAM.size());
        System.out.println("SecondOfPair size: " + secondOfPairSAM.size());

        // Key = exon_id, Value = exon
        HashMap<String, Exon> skippedExons = new HashMap<String, Exon>();
        for (Gene g : genes.values()) {
            // parse transcripts to arrayList
            ArrayList<Transcript> transcripts = new ArrayList<Transcript>(g.transcripts.values());

            // Error checking
            if (transcripts.size() != 2) {
                System.err.println("No more/less than two transcripts allowed! (" + transcripts.size() + ")");
            }

            // converting Treeset to ArrayList for iteration
            ArrayList<Exon> exons_tr1 = new ArrayList<Exon>(transcripts.get(0).exons);
            ArrayList<Exon> exons_tr2 = new ArrayList<Exon>(transcripts.get(1).exons);

            // iterate through FIRST transcript
            for (int i = 0; i < exons_tr1.size() - 1; i++) {
                boolean leftBoundaryExists = false;
                boolean rightBoundaryExists = false;
                for (int j = 0; j < exons_tr2.size(); j++) {
                    // Check in SECOND transcript if there's a exon with left boundary
                    if (exons_tr2.get(j).start == exons_tr1.get(i).start
                            && exons_tr2.get(j).end == exons_tr1.get(i).end) {
                        leftBoundaryExists = true;
                    }
                    // Check in SECOND transcript if there's a exon with right boundary
                    if (exons_tr2.get(j).start == exons_tr1.get(i + 1).start
                            && exons_tr2.get(j).end == exons_tr1.get(i + 1).end) {
                        rightBoundaryExists = true;
                    }
                }
                // check if there exists an Exon in the middle
                if (leftBoundaryExists && rightBoundaryExists) {
                    // first & last exons cannot be skipped exons FRAGLICH
                    for (int j = 0; j < exons_tr2.size(); j++) {
                        if (exons_tr2.get(j).start > exons_tr1.get(i).end
                                && exons_tr2.get(j).end < exons_tr1.get(i + 1).start) {
                            skippedExons.put(exons_tr2.get(j).exon_id, new Exon(exons_tr2.get(j)));
                        }
                    }
                }
            }

            // iterate through SECOND transcript
            for (int i = 0; i < exons_tr2.size() - 1; i++) {
                boolean leftBoundaryExists = false;
                boolean rightBoundaryExists = false;
                for (int j = 0; j < exons_tr1.size(); j++) {
                    // Check in SECOND transcript if there's a exon with left boundary
                    if (exons_tr1.get(j).start == exons_tr2.get(i).start
                            && exons_tr1.get(j).end == exons_tr2.get(i).end) {
                        leftBoundaryExists = true;
                    }
                    // Check in SECOND transcript if there's a exon with right boundary
                    if (exons_tr1.get(j).start == exons_tr2.get(i + 1).start
                            && exons_tr1.get(j).end == exons_tr2.get(i + 1).end) {
                        rightBoundaryExists = true;
                    }
                }
                // check if there exists an Exon in the middle
                if (leftBoundaryExists && rightBoundaryExists) {
                    // first & last exons cannot be skipped exons //FRAGLICH
                    for (int j = 0; j < exons_tr1.size(); j++) {
                        if (exons_tr1.get(j).start > exons_tr2.get(i).end
                                && exons_tr1.get(j).end < exons_tr2.get(i + 1).start) {
                            skippedExons.put(exons_tr1.get(j).exon_id, new Exon(exons_tr1.get(j)));
                        }
                    }
                }
            }
        }
        System.out.println("Found a total of " + skippedExons.size() + " skippedExons.");

        ArrayList<ResultLine> result = new ArrayList<ResultLine>();
        // evtl LinkedHashSet
        for (Exon exon : skippedExons.values()) {
            int num_incl_reads = 0;
            for (SAMRecord sr : firstOfPairSAM.values()) {
                // Check, if exon and samRecord are on the same chromosome
                if (sr.getReferenceName().equals(exon.chr)) {
                    // System.out.println("SAME CHR: Exon " + exon.exon_id + " (" + exon.chr + ")
                    // SAMRecord "
                    // + sr.getReadName() + " (" + sr.getReferenceName() + ")");
                    boolean numInclRaised = false;
                    // num_incl_reads++, if any alignmentblock of samRecord is within exon
                    for (AlignmentBlock ab : sr.getAlignmentBlocks()) {
                        int ref_start = ab.getReferenceStart();
                        int ref_end = ref_start + ab.getLength();
                        if (ref_start <= exon.end && ref_end >= exon.start) {
                            num_incl_reads++;
                            numInclRaised = true;
                            break;
                        }
                    }
                    // check on other SAM-Pair, if num_incl_reads was not already raised
                    if (!numInclRaised) {
                        for (AlignmentBlock ab : secondOfPairSAM.get(sr.getReadName()).getAlignmentBlocks()) {
                            int ref_start = ab.getReferenceStart();
                            int ref_end = ref_start + ab.getLength();
                            if (ref_start <= exon.end && ref_end >= exon.start) {
                                num_incl_reads++;
                                break;
                            }
                        }
                    }

                } else {
                    // continue with next SAMRecord
                    continue;
                }
            }

            int num_excl_reads = 0;
            int num_total_reads = num_incl_reads + num_excl_reads;
            double psi = num_total_reads == 0 ? 0 : num_incl_reads / num_total_reads;

            ResultLine currLine = new ResultLine(exon.gene_id, exon.exon_id, exon.start, exon.end, num_incl_reads,
                    num_excl_reads, num_total_reads, psi);
            result.add(currLine);
        }

        System.out.println("Begin: write  bam.annot.");
        try {
            // Create or overwrite new file
            PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(output, false)));
            out.println("gene\texon\tnum_incl_reads\tnum_excl_reads\tnum_total_reads\tpsi");
            for (ResultLine rl : result) {
                // exon.gene_id + "\t" + exon.exon_id + "\t" + (exon.start + 1) + "-" +
                // (exon.end + 1)

                // out.println(rl.gene_id + "\t" + rl.exon_id + "\t" + rl.exon_start + "-" +
                // rl.exon_end + "\t"
                // + rl.num_incl_reads + "\t" + rl.num_excl_reads + "\t" + rl.num_total_reads +
                // "\t" + rl.psi);
                out.println(rl.gene_id + "\t" + rl.exon_start + "-" + rl.exon_end + "\t" + rl.num_incl_reads + "\t"
                        + rl.num_excl_reads + "\t" + rl.num_total_reads + "\t" + rl.psi);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished: write output file");

        sam_reader.close();
        System.out.println("End of main method.");
    }

    public static HashMap<String, Gene> parseGene(File gtfInput) {
        // Key is the always the respective ID
        HashMap<String, Gene> allGenes = new HashMap<String, Gene>();
        // First String = gene_id of CDS
        int exonCounter = 0;
        int transcriptCounter = 0;
        int geneCounter = 0;
        // read gtfInput-file
        try {
            BufferedReader br = new BufferedReader(new FileReader(gtfInput));
            String line = null;
            int linecounter = 0; // for error printing only
            System.out.println("Begin: Parsing gtf file.");
            while ((line = br.readLine()) != null) {
                linecounter++;
                // ignore comments (beginning with "#")
                if (!line.startsWith("#")) {
                    // For every line in input
                    String[] tabSeparated = line.split("\\t");

                    String seqname = tabSeparated[0];
                    // String source = tabSeparated[1];
                    String feature = tabSeparated[2];
                    String start = tabSeparated[3];
                    String end = tabSeparated[4];
                    // String score = tabSeparated[5];
                    String strand = tabSeparated[6];
                    // String frame = tabSeparated[7];
                    String attribute = tabSeparated[8];

                    // -------For lines, which are exons:-------
                    if (feature.equals("exon")) {
                        exonCounter++;
                        // parameters needed to construct new Exon()
                        String exon_id = "";
                        String exon_number = "";
                        String transcript_id = "";
                        String transcript_name = "";
                        String gene_id = "";
                        String gene_name = "";

                        // gather parameters from String "attribute"
                        String[] attributeSeparated = attribute.split(";");
                        // search in attributeSeparated for parameters
                        for (int i = 0; i < attributeSeparated.length; i++) {
                            if (attributeSeparated[i].contains("exon_id")) {
                                // get only value between quotation marks
                                exon_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1,
                                        attributeSeparated[i].lastIndexOf("\""));
                                continue;
                            } else if (attributeSeparated[i].contains("transcript_id")) {
                                transcript_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1,
                                        attributeSeparated[i].lastIndexOf("\""));
                                continue;
                            } else if (attributeSeparated[i].contains("transcript_name")) {
                                transcript_name = attributeSeparated[i].substring(
                                        attributeSeparated[i].indexOf("\"") + 1,
                                        attributeSeparated[i].lastIndexOf("\""));
                                continue;
                            } else if (attributeSeparated[i].contains("exon_number")) {
                                exon_number = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1,
                                        attributeSeparated[i].lastIndexOf("\""));
                                continue;
                            } else if (attributeSeparated[i].contains("gene_id")) {
                                gene_id = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1,
                                        attributeSeparated[i].lastIndexOf("\""));
                                continue;
                            } else if (attributeSeparated[i].contains("gene_name")) {
                                gene_name = attributeSeparated[i].substring(attributeSeparated[i].indexOf("\"") + 1,
                                        attributeSeparated[i].lastIndexOf("\""));
                                continue;
                            }
                        }
                        // Construct exon and add to exonMap
                        Exon e = new Exon(seqname, Integer.parseInt(start), Integer.parseInt(end) + 1, strand, exon_id,
                                Integer.parseInt(exon_number), transcript_id, transcript_name, gene_id, gene_name);
                        // Check, if all exon values are actually filled
                        // Note: e.gene_name can be empty, is not checked
                        if (e.start == 0 || e.end == 0 || e.chr.isEmpty() || e.strand.isEmpty() || e.exon_id.isEmpty()
                                || e.exon_number == 0 || e.transcript_id.isEmpty() || e.gene_id.isEmpty()) {
                            System.err.println("Exon in line " + linecounter + " has an empty value!");
                        }
                        // exonCounter++;

                        // Create data structure allGenes, contains transcripts,
                        // contains exons
                        if (!allGenes.containsKey(e.gene_id)) {
                            // add a new gene to allGenes, which first requires
                            // a set of transcripts, which requires a
                            // set of exons
                            TreeSet<Exon> exons = new TreeSet<Exon>(new Exon());
                            exons.add(new Exon(e));

                            HashMap<String, Transcript> transcripts = new HashMap<String, Transcript>();
                            transcripts.put(e.transcript_id, new Transcript(e.transcript_id, exons));

                            allGenes.put(e.gene_id, new Gene(e.gene_id, transcripts));
                            geneCounter++;
                            transcriptCounter++;
                        } else {
                            // check, if existing gene already contains the
                            // found transcript
                            Gene currentGene = allGenes.get(e.gene_id);
                            if (!currentGene.transcripts.containsKey(e.transcript_id)) {
                                // add a new transcript to transcripts, which
                                // first requires a set of exons
                                TreeSet<Exon> exons = new TreeSet<Exon>(new Exon());
                                exons.add(new Exon(e));

                                currentGene.transcripts.put(new String(e.transcript_id),
                                        new Transcript(e.transcript_id, exons));
                                // System.out.println("New transcript " +
                                // t.transcript_id + " created and exon " +
                                // e.exon_id
                                // + " added.");

                                transcriptCounter++;
                            } else {
                                // If transcript already exists, just add new
                                // exon to this transcript
                                currentGene.transcripts.get(e.transcript_id).exons.add(new Exon(e));
                                // System.out.println("Exon " + e.exon_id + "
                                // added to existing transcript "
                                // +
                                // allTranscripts.get(e.transcript_id).transcript_id);
                            }
                        }
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        // Print Transcript Information to console
        System.out.println("Finished: Parsing file.");
        System.out.println("exonCounter: " + exonCounter);
        System.out.println("transcriptCounter:" + transcriptCounter);
        System.out.println("geneCounter:" + geneCounter);
        return allGenes;
    }

}