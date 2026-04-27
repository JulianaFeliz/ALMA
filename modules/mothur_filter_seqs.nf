process MOTHUR_FILTER_SEQS {
    tag "${meta.id}"
    label 'process_low'
    
    conda 'bioconda::mothur=1.48.0'
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'
    
    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path("*.filter.fasta"), emit: filtered
    tuple val(meta), path("*.filter"), emit: filter
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    mothur "#filter.seqs(fasta=${fasta}, vertical=T, processors=${task.cpus})"
    
    # Rename output files
    mv *.filter.fasta ${prefix}.filter.fasta || true
    mv *.filter ${prefix}.filter || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v | head -n 1 | cut -d'=' -f2 | cut -d' ' -f1)
    END_VERSIONS
    """
}
