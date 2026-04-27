process MOTHUR_UNIQUE_SEQS {
    tag "$meta.id"
    label 'process_low'
    
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'
    
    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path("${meta.id}.unique.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.{names,count_table}"), emit: names
    path "versions.yml", emit: versions
    
    script:
    """
    # Run Mothur unique.seqs
    mothur "#unique.seqs(fasta=${fasta})"
    
    
    # 2. Renomeia o FASTA gerado (seja ele good.unique da etapa 1 ou filter.unique da etapa 2)
    for f in *.unique.fasta; do
        if [ -f "\$f" ]; then
            mv "\$f" "${meta.id}.unique.fasta"
        fi
    done

    # 3. Renomeia a tabela gerada (cobre mothur novo e velho)
    if ls *.count_table 2>/dev/null; then
        mv *.count_table "${meta.id}.count_table" || mv *.count_table "${meta.id}.names"
    elif ls *.names 2>/dev/null; then
        mv *.names "${meta.id}.names"
    fi

    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
