process MOTHUR_MAKE_CONTIGS {
    tag "$meta.id"
    label 'process_medium'
    
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'
    
    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    
    output:
    tuple val(meta), path("${meta.id}.contigs.fasta"), emit: contigs
    tuple val(meta), path("${meta.id}.contigs.report"), emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    # Create file list for Mothur
    echo "${meta.id}\t${fastq_1}\t${fastq_2}" > ${meta.id}.files
    
    # Run Mothur make.contigs
    mothur "#make.contigs(file=${meta.id}.files, processors=${task.cpus})"
    
    # Rename output files to match expected names

    if [ -f "${meta.id}.trim.contigs.fasta" ]; then
        mv "${meta.id}.trim.contigs.fasta" "${meta.id}.contigs.fasta"
    else
        echo "ERROR: ${meta.id}.trim.contigs.fasta not found!"
        exit 1
    fi    


    # 2. Renomear o REPORT (de _contigs_report ou .contigs_report para .contigs.report)
    if [ -f "${meta.id}.contigs_report" ]; then
        mv "${meta.id}.contigs_report" "${meta.id}.contigs.report"
    elif [ -f "${meta.id}.contigs.report" ]; then
        echo "Report file already has correct name."
    else
        echo "ERROR: Report file not found!"
        ls -la *report* || echo "No report files found"
        exit 1
    fi

    # 3. Verificação de segurança para o Nextflow não travar
    if [ ! -f "${meta.id}.contigs.fasta" ]; then
        echo "ERROR: Fasta file not found!"
        ls -la
        exit 1
    fi
    
    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
   """
}
