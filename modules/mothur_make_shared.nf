process MOTHUR_MAKE_SHARED {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/mothur/shared", mode: 'copy'
    
    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'
    
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
