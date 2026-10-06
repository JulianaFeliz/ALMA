process MAFFT_ALIGN {
    label 'process_high'
    tag "Global_OTUs" 
        
    // Link do Galaxy Project mantido
    container 'https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1'

    input:
    path query_seqs     // Recebe VSEARCH_MAP_OTUS.out.rep_seqs
    path reference_seqs // params.reference_mafrax

    output:
    path "aligned_otus.fasta", emit: alignment
    path "versions.yml", emit: versions
    
    script:
    """
    export TMPDIR=\$PWD
    export TEMP=\$PWD
    export TMP=\$PWD

    # 1. Limpeza Global
    cat ${query_seqs} | tr -d '\\r' | sed '/^\$/d' | sed 's/|/_/g' | sed 's/[[:blank:]]/_/g' | tr -s '_' > clean_query.fasta

    # 2. Limpar a Referência
    cat ${reference_seqs} | tr -d '\\r' | sed '/^\$/d' > clean_ref.fasta

    # 3. ALINHAMENTO
    # --keeplength é vital: Força o MAFFT a não alterar o comprimento do alinhamento de referência,
    # garantindo que o RAxML consiga sobrepor os resultados perfeitamente na sua árvore!
    mafft --thread ${task.cpus} \\
          --anysymbol \\
          --memsavetree \\
          --keeplength \\
          --add clean_query.fasta \\
          --reorder \\
          clean_ref.fasta > aligned_otus.fasta

    # 4. Versão
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//')
    END_VERSIONS
    """
}