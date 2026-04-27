process MOTHUR_CLASSIFY_SEQS {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}/mothur/taxonomy", mode: 'copy'

    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'
    
    input:
    tuple val(meta), path(fasta), path(count)
    path reference
    path taxonomy
    val cutoff

    output:
    tuple val(meta), path("${meta.id}.taxonomy"), emit: taxonomy
    tuple val(meta), path("${meta.id}.tax.summary"), emit: summary
    path "versions.yml", emit: versions

    script:
    def cutoff_value = cutoff ?: params.cutoff ?: 80
    """
    # Run Mothur classify.seqs
    mothur "#classify.seqs(fasta=${fasta}, count=${count}, reference=${reference}, taxonomy=${taxonomy}, cutoff=${cutoff_value}, processors=${task.cpus})"

    # Rename outputs to expected format
    # Mothur creates files like: filename.taxonomy and filename.tax.summary
    for file in *.taxonomy; do
        [ -f "\$file" ] && mv "\$file" ${meta.id}.taxonomy
    done

    for file in *.tax.summary; do
        [ -f "\$file" ] && mv "\$file" ${meta.id}.tax.summary
    done

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
