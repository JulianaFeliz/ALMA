process VSEARCH_MERGE {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("${meta.id}.merged.fastq"), emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    vsearch --fastq_mergepairs ${reads_1} \\
            --reverse ${reads_2} \\
            --fastqout ${meta.id}.merged.fastq \\
            --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}