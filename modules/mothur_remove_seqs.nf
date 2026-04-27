process MOTHUR_REMOVE_SEQS {
    tag "${meta.id}"
    label 'process_low'
    
    conda 'bioconda::mothur=1.48.0'
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'
    
    input:
    tuple val(meta), path(fasta), path(accnos), path(count)
    
    output:
    tuple val(meta), path("*.pick.fasta"), emit: fasta
    tuple val(meta), path("*.pick.count_table"), emit: count
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    mothur "#remove.seqs(fasta=${fasta}, accnos=${accnos}, count=${count})"
    
    # Rename output files
    mv *.pick.fasta ${prefix}.pick.fasta || true
    mv *.pick.count_table ${prefix}.pick.count_table || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v | head -n 1 | cut -d'=' -f2 | cut -d' ' -f1)
    END_VERSIONS
    """
}
