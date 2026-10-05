process MOTHUR_PRE_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    tuple val(meta), path(fasta), path(count)

    output:
    tuple val(meta), path("${meta.id}.precluster.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.precluster.count_table"), emit: count
    tuple val(meta), path("${meta.id}.precluster.map"), emit: map, optional: true
    path "versions.yml", emit: versions

    script:
    //  'precluster_diffs = 3' in your nextflow.config, it uses 3.
    def diffs = params.precluster_diffs ?: 3 
    
    """
    mothur "#pre.cluster(fasta=${fasta}, count=${count}, diffs=${diffs}, processors=${task.cpus})"

    # Bulletproof renaming - only rename if the file actually exists!
    if ls *.precluster.fasta 1> /dev/null 2>&1; then
        mv *.precluster.fasta "${meta.id}.precluster.fasta"
    fi

    if ls *.precluster.count_table 1> /dev/null 2>&1; then
        mv *.precluster.count_table "${meta.id}.precluster.count_table"
    fi

    if ls *.precluster.*.map 1> /dev/null 2>&1; then
        mv *.precluster.*.map "${meta.id}.precluster.map"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
