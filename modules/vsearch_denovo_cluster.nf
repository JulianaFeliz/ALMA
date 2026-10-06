process VSEARCH_DENOVO_CLUSTER {
    label 'process_high'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    path unmapped
    val cutoff

    output:
    path "denovo_centroids.fasta", emit: centroids

    script:
    def id_float = 1.0 - (cutoff as Float)
    """
    if [ -s ${unmapped} ]; then
        vsearch --cluster_size ${unmapped} \\
                --id ${id_float} \\
                --centroids denovo_centroids.fasta \\
                --sizein --sizeout \\
                --fasta_width 0 \\
                --threads ${task.cpus}
    else
        touch denovo_centroids.fasta
    fi
    """
}