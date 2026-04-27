process MOTHUR_GET_LINEAGE {
    tag "$meta.id"
    label 'process_low'

    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(count)
    tuple val(meta), path(taxonomy)

    output:
    tuple val(meta), path("${meta.id}.unique.pick.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.pick.count_table"), emit: count
    tuple val(meta), path("${meta.id}.pick.taxonomy"), emit: taxonomy
    path "versions.yml", emit: versions

    script:
    """
    # Run Mothur get.lineage to extract Glomeromycota sequences
    mothur "#get.lineage(fasta=${fasta}, count=${count}, \\
            taxonomy=${taxonomy}, taxon=Glomeromycota)"

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
