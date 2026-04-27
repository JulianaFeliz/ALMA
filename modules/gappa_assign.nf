process GAPPA_ASSIGN {
    tag "${meta.id}"
    label 'process_low'

    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/gappa:0.8.5--h077b44d_3'

    publishDir "${params.outdir}/gappa_results", mode: 'copy'    

    input:
    tuple val(meta), path(jplace)
    path taxonomy_file
    val lwr_threshold

    output:
    tuple val(meta), path("*.taxonomy.tsv"), emit: taxonomy
    tuple val(meta), path("*.profile.tsv"), emit: profile, optional: true
    tuple val(meta), path("*.per_query.tsv"), emit: per_query, optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    """
    gappa examine assign \\
        --jplace-path ${jplace} \\
        --taxon-file ${taxonomy_file} \\
        --per-query-results \\
        --krona \\
        --sativa \\
        --best-hit \\
        --distribution-ratio ${lwr_threshold} \\
        --out-dir ./

    # Rename output files with consistent prefix
    if [ -f "profile.tsv" ]; then
        mv profile.tsv ${prefix}.taxonomy.tsv
    elif [ -f "${prefix}.profile.tsv" ]; then
        cp ${prefix}.profile.tsv ${prefix}.taxonomy.tsv
    fi

    if [ -f "per_query.tsv" ]; then
        mv per_query.tsv ${prefix}.per_query.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gappa: \$(gappa --version 2>&1 | head -n 1 | sed 's/gappa version: //g')
    END_VERSIONS
    """
}
