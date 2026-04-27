process MOTHUR_FILTER_SEQS {
    tag "$meta.id"
    label 'process_low'

    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("${meta.id}.filtered.*"), emit: filtered
    path "versions.yml", emit: versions

    script:
    """
    # 1. Safely find the alignment file
    MY_ALN=""
    for f in *.align *.fasta; do
        if [ -f "\$f" ]; then
            MY_ALN="\$f"
            break
        fi
    done

    # 2. Run Mothur
    mothur "#filter.seqs(fasta=\$MY_ALN, vertical=T, trump=., processors=${task.cpus})"

    # 3. Safely rename the filtered output
    for f in *.filter.fasta *.filter.align; do
        if [ -f "\$f" ]; then
            mv "\$f" "${meta.id}.filtered.fasta"
            break
        fi
    done

    # 4. Safely copy the count table to travel with it
    for f in *.count_table; do
        if [ -f "\$f" ]; then
            cp "\$f" "${meta.id}.filtered.count_table"
            break
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
