process MOTHUR_CHIMERA_VSEARCH {
    tag "$meta.id"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    // 💅 DIVA FIX: One single pair of hands to catch the VIP package!
    tuple val(meta), path(fasta), path(count)

    output:
    tuple val(meta), path("${meta.id}.denovo.vsearch.chimeras"), emit: chimeras, optional: true
    tuple val(meta), path("${meta.id}.denovo.vsearch.accnos"), emit: accnos
    path "versions.yml", emit: versions

    script:
    """
    mothur "#chimera.vsearch(fasta=${fasta}, count=${count}, dereplicate=T, processors=${task.cpus})"

    # Rename chimeras output if it exists
    if ls *.denovo.vsearch.chimeras 1> /dev/null 2>&1; then
        mv *.denovo.vsearch.chimeras "${meta.id}.denovo.vsearch.chimeras"
    fi

    # DIVA ANTI-CRASH MECHANISM:
    # If no chimeras were found, Mothur doesn't create an accnos file.
    # We create an empty one so Nextflow doesn't panic and crash.
    if ls *.denovo.vsearch.accnos 1> /dev/null 2>&1; then
        mv *.denovo.vsearch.accnos "${meta.id}.denovo.vsearch.accnos"
    else
        touch "${meta.id}.denovo.vsearch.accnos"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
