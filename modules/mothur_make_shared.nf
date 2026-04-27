process MOTHUR_MAKE_SHARED {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/mothur/shared", mode: 'copy'
    
    conda 'bioconda::mothur=1.48.0'
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'
    
    input:
    tuple val(meta), path(list), path(count)
    
    output:
    tuple val(meta), path("*.shared"), emit: shared
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    mothur "#make.shared(list=${list}, count=${count})"
    
    # Rename output files
    mv *.shared ${prefix}.shared || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v | head -n 1 | cut -d'=' -f2 | cut -d' ' -f1)
    END_VERSIONS
    """
}
