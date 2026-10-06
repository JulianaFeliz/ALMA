process VSEARCH_CHIMERA {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input: 
    tuple val(meta), path(fasta)

    output:
    // Emitimos SÓ O PATH para facilitar a fusão das amostras no próximo passo
    path "${meta.id}.clean.fasta", emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    # 1. Remove quimeras (de novo)
    vsearch --uchime_denovo ${fasta} \\
            --nonchimeras nonchimeric.fasta \\
            --sizein \\
            --fasta_width 0 \\
            --threads ${task.cpus}

    # 2. O truque de Mestre: Adiciona o nome da amostra (meta.id) dentro do FASTA
    # Exemplo: '>1' vira '>Amostra1_1' 
    sed "s/^>/>\({meta.id}_/" nonchimeric.fasta >\){meta.id}.clean.fasta
    """
}