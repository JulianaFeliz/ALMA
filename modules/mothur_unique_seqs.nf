process MOTHUR_UNIQUE_SEQS {
    tag "$meta.id"
    label 'process_low'

    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    // Trocamos o nome para "files" porque pode chegar 1 arquivo ou 2 arquivos!
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.unique.fasta"), emit: fasta
    // Emitimos a nova count_table (ou names)
    tuple val(meta), path("${meta.id}.count_table"), emit: names
    path "versions.yml", emit: versions

    script:
    """
    # 1. Pega os nomes exatos dos arquivos que o Nextflow jogou na pasta
    MY_FASTA=\$(ls *.fasta | head -n 1)
    
    # 2. Verifica se veio uma tabela de contagem junto no cano
    if ls *.count_table 1> /dev/null 2>&1; then
        MY_COUNT=\$(ls *.count_table | head -n 1)
        # Roda o comando separando certinho Fasta e Count!
        mothur "#unique.seqs(fasta=\$MY_FASTA, count=\$MY_COUNT)"
    else
        # Se não tiver tabela (ex: no UNIQUE_SEQS_2), roda o comando antigo
        mothur "#unique.seqs(fasta=\$MY_FASTA)"
    fi

    # 3. Renomear o Fasta gerado (o Mothur sempre adiciona .unique no meio)
    mv *.unique.fasta "${meta.id}.unique.fasta"

    # 4. Renomear a tabela gerada
    if ls *.unique.count_table 1> /dev/null 2>&1; then
        mv *.unique.count_table "${meta.id}.count_table"
    elif ls *.names 1> /dev/null 2>&1; then
        mv *.names "${meta.id}.count_table"
    fi

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
