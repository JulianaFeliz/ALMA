process MOTHUR_ALIGN_SEQS {
    tag "${meta.id}"
    label 'process_medium'
    
    conda 'bioconda::mothur=1.48.0'
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'
    
    input:
    tuple val(meta), path(fasta)
    path reference
    
    output:
    tuple val(meta), path("*.align"), emit: aligned
    tuple val(meta), path("*.align.report"), emit: report
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    mothur "#align.seqs(fasta=${fasta}, reference=${reference}, processors=${task.cpus})"
    
    # Rename output files
    mv *align_report ${prefix}.align.report 2>/dev/null || true
    mv *.align.report ${prefix}.align.report 2>/dev/null || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v | head -n 1 | cut -d'=' -f2 | cut -d' ' -f1)
    END_VERSIONS
    """
}
