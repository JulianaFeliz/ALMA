process MOTHUR_COUNT_SEQS {
    tag "$meta.id"
    label 'process_low'
    memory '32 GB'
    cpus '1'
    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'
    
    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(names)

    output:
    tuple val(meta), path("${meta.id}.count_table"), emit: count
    path "versions.yml", emit: versions

    script:
    """
    # O Nextflow entrega o arquivo do passo anterior na variavel \${names}. 
    # Vamos apenas checar o final do nome desse arquivo.

    if [[ "${names}" == *.count_table ]]; then
        # Se ja for a tabela de contagem, copia pro formato de saida exigido
        cp "${names}" "${meta.id}.count_table"
    else
        # Se for formato antigo (.names), pede pro Mothur converter
        mothur "#count.seqs(name=${names})"
        mv *.count_table "${meta.id}.count_table"
    fi

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
