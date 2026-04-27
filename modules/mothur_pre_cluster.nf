process MOTHUR_PRE_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'
    
    conda 'bioconda::mothur=1.48.0'
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'
    
    input:
    tuple val(meta), path(fasta), path (count)
    
    output:
    tuple val(meta), path("*.precluster.fasta"), emit: fasta
    tuple val(meta), path("*.precluster.count_table"), emit: count
    tuple val(meta), path("*.precluster.map"), emit: map, optional: true
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def diffs = 2  // Allow 2 differences for pre-clustering
    """
    mothur "#pre.cluster(fasta=${fasta}, count=${count}, diffs=${diffs}, processors=${task.cpus})"
    
    # Rename output files
    mv *.precluster.fasta ${prefix}.precluster.fasta || true
    mv *.precluster.count_table ${prefix}.precluster.count_table || true
    mv *.precluster.*.map ${prefix}.precluster.map || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v | head -n 1 | cut -d'=' -f2 | cut -d' ' -f1)
    END_VERSIONS
    """
}
