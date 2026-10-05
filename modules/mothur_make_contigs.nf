process MOTHUR_MAKE_CONTIGS {
    tag "$meta.id"
    label 'process_medium'

    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    // AGORA SIM! Emitindo Fasta e Groups na mesma tupla
    tuple val(meta), path("${meta.id}.contigs.fasta"), path("${meta.id}.contigs.count_table"), emit: contigs
    tuple val(meta), path("${meta.id}.contigs.report"), emit: report
    path "versions.yml", emit: versions

    script:
    """
    # Create file list for Mothur
    echo "${meta.id}\t${fastq_1}\t${fastq_2}" > ${meta.id}.files

    # Run Mothur make.contigs
    mothur "#make.contigs(file=${meta.id}.files, processors=${task.cpus})"

    # 1. Renomear o FASTA
    if [ -f "${meta.id}.trim.contigs.fasta" ]; then
        mv "${meta.id}.trim.contigs.fasta" "${meta.id}.contigs.fasta"
    else
        echo "ERROR: ${meta.id}.trim.contigs.fasta not found!"
        exit 1
    fi

    # 2. Renomear o GROUPS (O Mothur gera ele com base no nome do .files)
    if [ -f "${meta.id}.groups" ]; then
        mv "${meta.id}.groups" "${meta.id}.contigs.groups"
    elif [ -f "${meta.id}.contigs.count_table" ]; then
        echo "Count table file already has correct name."
    else
        echo "ERROR: Count table file not found!"
        ls -la *.count_table* || echo "No count files found"
        exit 1
    fi

    # 3. Renomear o REPORT
    if [ -f "${meta.id}.contigs_report" ]; then
        mv "${meta.id}.contigs_report" "${meta.id}.contigs.report"
    elif [ -f "${meta.id}.contigs.report" ]; then
        echo "Report file already has correct name."
    else
        echo "ERROR: Report file not found!"
        ls -la *report* || echo "No report files found"
        exit 1
    fi

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
