process VSEARCH_FILTER {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    tuple val(meta), path(fastq)
    val min_len
    val max_len
    val max_ns // Antigo max_ambig

    output:
    tuple val(meta), path("${meta.id}.filtered.fasta"), emit: fasta 
    path "versions.yml", emit: versions

    script:
    """
    vsearch --fastq_filter ${fastq} \\
            --fastq_maxee 1.0 \\
            --fastq_minlen ${min_len} \\
            --fastq_maxlen ${max_len} \\
            --fastq_maxns ${max_ns} \\
            --fastaout ${meta.id}.filtered.fasta \\
            --threads ${task.cpus}
    """
}