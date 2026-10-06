process VSEARCH_MAP_OTUS {
    label 'process_high'
    publishDir "${params.outdir}/vsearch_clustering", mode: 'copy'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    path closed_centroids
    path denovo_centroids
    path all_samples
    val cutoff

    output:
    path "final_rep_seqs.fasta", emit: rep_seqs
    path "final_otu_table.shared", emit: shared

    script:
    def id_float = 1.0 - (cutoff as Float)
    """
    cat \({closed_centroids}\){denovo_centroids} > final_rep_seqs.fasta

    # Mapeia amostras originais nas referências finais para gerar a contagem
    vsearch --usearch_global ${all_samples} \\
            --db final_rep_seqs.fasta \\
            --id ${id_float} \\
            --otutabout vsearch_otu_table.txt \\
            --threads ${task.cpus}

    # TRADUTOR MOTHUR -> .shared
    awk -v label="${cutoff}" '
    BEGIN { FS="\\t"; OFS="\\t" }
    NR==1 {
        num_samples = NF - 1
        for (i=2; i<=NF; i++) samples[i-1] = \$i
    }
    NR>1 {
        num_otus++
        otus[num_otus] = \$1
        for (i=2; i<=NF; i++) counts[i-1, num_otus] = \$i
    }
    END {
        printf "label\\tGroup\\tnumOtus"
        for (o=1; o<=num_otus; o++) printf "\\t%s", otus[o]
        printf "\\n"
        
        for (s=1; s<=num_samples; s++) {
            printf "%s\\t%s\\t%d", label, samples[s], num_otus
            for (o=1; o<=num_otus; o++) printf "\\t%s", counts[s, o]
            printf "\\n"
        }
    }' vsearch_otu_table.txt > final_otu_table.shared
    """
}