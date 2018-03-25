
public class ResultLine implements Comparable<ResultLine> {
    String gene_id;
    String exon_id;
    int exon_start;
    int exon_end;
    int num_incl_reads;
    int num_excl_reads;
    int num_total_reads;
    double psi;

    public ResultLine(String gene_id, String exon_id, int exon_start, int exon_end, int num_incl_reads,
                      int num_excl_reads, int num_total_reads, double psi) {
        this.gene_id = gene_id;
        this.exon_id = exon_id;
        this.exon_start = exon_start;
        this.exon_end = exon_end;
        this.num_incl_reads = num_incl_reads;
        this.num_excl_reads = num_excl_reads;
        this.num_total_reads = num_total_reads;
        this.psi = psi;

    }

    @Override
    public int compareTo(ResultLine rl) {
        // Sorts GeneRegions ascending to their start-points
        if (this.gene_id.compareTo(rl.gene_id) > 0) {
            return 1;
        } else if (this.gene_id.compareTo(rl.gene_id) < 0) {
            return -1;
        } else {
            return 0;
        }

    }

}
