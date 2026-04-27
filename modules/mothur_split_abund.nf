process MOTHUR_SPLIT_ABUND {
    tag "$meta.id"
    label 'process_low'

    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    // Takes the synchronized package straight from PRE_CLUSTER
    tuple val(meta), path(fasta), path(count)

    output:
    // Emits the cleaned, abundant sequences
    tuple val(meta), path("${meta.id}.abund.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.abund.count_table"), emit: count
    path "versions.yml", emit: versions

    script:
    """
    # Franco's exact standard: cutoff=1 to remove singletons
    mothur "#split.abund(fasta=${fasta}, count=${count}, cutoff=1, accnos=true)"

    # Bulletproof renaming: grab the '.abund' files Mothur generates
    MY_FASTA=\$(ls *.abund.fasta 2>/dev/null | head -n 1)
    MY_COUNT=\$(ls *.abund.count_table 2>/dev/null | head -n 1)

    if [ ! -z "\$MY_FASTA" ]; then
        mv "\$MY_FASTA" "${meta.id}.abund.fasta"
    fi

    if [ ! -z "\$MY_COUNT" ]; then
        mv "\$MY_COUNT" "${meta.id}.abund.count_table"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
