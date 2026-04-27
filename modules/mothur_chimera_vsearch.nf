process MOTHUR_CHIMERA_VSEARCH {
    tag "${meta.id}"
    label 'process_medium'

    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(count)

    output:
    tuple val(meta), path("${meta.id}.denovo.vsearch.chimeras"), emit: chimeras
    tuple val(meta), path("${meta.id}.denovo.vsearch.accnos"), emit: accnos
    path "versions.yml", emit: versions

    script:
    """
    mothur "#chimera.vsearch(fasta=${fasta}, count=${count}, dereplicate=T, processors=${task.cpus})"

    # Rename output files - only if needed
    if [ -f *.denovo.vsearch.chimeras ] && [ ! -f ${meta.id}.denovo.vsearch.chimeras ]; then
        mv *.denovo.vsearch.chimeras ${meta.id}.denovo.vsearch.chimeras
    fi
    
    if [ -f *.denovo.vsearch.accnos ] && [ ! -f ${meta.id}.denovo.vsearch.accnos ]; then
        mv *.denovo.vsearch.accnos ${meta.id}.denovo.vsearch.accnos
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v 2>&1 | head -n 1 | sed 's/.*v\\.//; s/ .*//')
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//')
    END_VERSIONS
    """
}
