process VSEARCH_GLOBAL_DEREP {
    label 'process_medium'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    path clean_fastas

    output:
    // Solta o arquivo junto com as amostras originais (precisaremos para contar as OTUs depois)
    path "all_samples.fasta", emit: all_samples
    // Solta as sequências desreplicadas para o Open Reference
    path "global_derep.fasta", emit: global_derep

    script:
    """
    cat *.glomero.fasta > all_samples.fasta

    vsearch --derep_fulllength all_samples.fasta \\
            --sizein --sizeout \\
            --fasta_width 0 \\
            --fastaout global_derep.fasta \\
            --threads ${task.cpus}
    """
}