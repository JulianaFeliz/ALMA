process MOTHUR_COUNT_SEQS {
    tag "$meta.id"
    label 'process_low'

    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'

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
