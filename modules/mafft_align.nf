process MAFFT_ALIGN {
    tag "$meta.id"
        
    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1'

    input:
    tuple val(meta), path(query_seqs)
    path reference_seqs

    output:
    tuple val(meta), path("${meta.id}.aligned.fasta"), emit: alignment
    path "versions.yml", emit: versions
    
    script:
    """
    export TMPDIR=\$PWD
    export TEMP=\$PWD
    export TMP=\$PWD

    # 1. Limpeza e Renomeação (Versão Segura para Nextflow)
    # Remove retornos de carro, linhas vazias e caracteres especiais
    # Em seguida, adiciona o ID da amostra antes do nome original
    cat ${query_seqs} | tr -d '\\r' | sed '/^\$/d' | sed 's/|/_/g' | sed 's/[[:blank:]]/_/g' | tr -s '_' | sed "s/>/>${meta.id}_/g" > clean_query.fasta

    # 2. Limpar a Referencia
    cat ${reference_seqs} | tr -d '\\r' | sed '/^\$/d' > clean_ref.fasta

    # 3. ALINHAMENTO
    mafft --thread ${task.cpus} \\
          --anysymbol \\
          --memsavetree \\
          --add clean_query.fasta \\
          --reorder \\
          clean_ref.fasta > ${meta.id}.aligned.fasta

    # 4. Versao
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//')
    END_VERSIONS
    """

}
