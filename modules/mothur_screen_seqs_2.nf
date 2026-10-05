process MOTHUR_SCREEN_SEQS_2 {
    tag "$meta.id"
    label 'process_low'
    time '12h'
    cpus 6
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    tuple val(meta), path(aligned_fasta), path(count_table)
    val start_pos
    val end_pos

    output:
    tuple val(meta), path("${meta.id}.good.*"), emit: seqs
    tuple val(meta), path("${meta.id}.bad.accnos"), emit: bad_seqs, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Roda o Mothur exigindo que as sequências comecem e terminem nas posições certas
    mothur "#screen.seqs(fasta=${aligned_fasta}, count=${count_table}, start=${start_pos}, end=${end_pos}, processors=${task.cpus})"

    # Renomear saídas
    if [ -f *.good.align ] && [ ! -f ${meta.id}.good.align ]; then
        mv *.good.align ${meta.id}.good.align
    fi
    
    # O Fasta normal caso a extensão seja fasta em vez de align
    if [ -f *.good.fasta ] && [ ! -f ${meta.id}.good.fasta ]; then
        mv *.good.fasta ${meta.id}.good.fasta
    fi

    if [ -f *.good.count_table ] && [ ! -f ${meta.id}.good.count_table ]; then
        mv *.good.count_table ${meta.id}.good.count_table
    fi

    if [ -f *.bad.accnos ] && [ ! -f ${meta.id}.bad.accnos ]; then
        mv *.bad.accnos ${meta.id}.bad.accnos
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
