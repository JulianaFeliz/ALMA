process MOTHUR_SCREEN_SEQS {
    tag "$meta.id"
    label 'process_low'
    time '12h'
    cpus 6
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'

    input:
    
    tuple val(meta), path(fasta), path(count_file) 
    val min_length
    val max_length
    val max_ambig
    val max_homop

    output:
    
    tuple val(meta), path("${meta.id}.good.*"), emit: seqs 
    tuple val(meta), path("${meta.id}.bad.accnos"), emit: bad_seqs, optional: true
    path "versions.yml", emit: versions

    script:
    """
    mothur "#screen.seqs(fasta=${fasta}, count=${count_file}, maxambig=${max_ambig}, minlength=${min_length}, maxlength=${max_length}, maxhomop=${max_homop}, processors=${task.cpus})"

    # Renomear saídas
    if [ -f *.good.fasta ] && [ ! -f ${meta.id}.good.fasta ]; then
        mv *.good.fasta ${meta.id}.good.fasta
    fi

    
    # Se for count_table
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
