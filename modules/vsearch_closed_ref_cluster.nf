process VSEARCH_CLOSED_REF_CLUSTER {
    label 'process_high'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    path global_derep
    path reference_db
    val cutoff

    output:
    path "closed_centroids.fasta", emit: centroids
    path "unmapped.fasta", emit: unmapped

    script:
    def id_float = 1.0 - (cutoff as Float)
    """
    vsearch --usearch_global ${global_derep} \\
            --db ${reference_db} \\
            --id ${id_float} \\
            --centroids closed_centroids.fasta \\
            --notmatched unmapped.fasta \\
            --sizein --sizeout \\
            --fasta_width 0 \\
            --threads ${task.cpus}
    """
}