process VSEARCH_DEREP {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.derep.fasta"), emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    vsearch --derep_fulllength ${fasta} \\
            --sizeout \\
            --minuniquesize 2 \\
            --fastaout ${meta.id}.derep.fasta \\
            --threads ${task.cpus}
    """
}