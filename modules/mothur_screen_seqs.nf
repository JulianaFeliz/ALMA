process MOTHUR_SCREEN_SEQS {
    tag "$meta.id"
    label 'process_low'

    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'

    input:
    tuple val(meta), path(fasta)
    val min_length
    val max_length
    val max_ambig
    val max_homop

    output:
    tuple val(meta), path("${meta.id}.good.*"), emit: seqs
    tuple val(meta), path("${meta.id}.bad.accnos"), emit: bad_seqs, optional: true
    path "versions.yml", emit: versions

    script:
    def length_args = (min_length != null && max_length != null) ?
        "minlength=${min_length}, maxlength=${max_length}" : ""
    def ambig_args = (max_ambig != null) ? "maxambig=${max_ambig}" : ""
    def homop_args = (max_homop != null) ? "maxhomop=${max_homop}" : ""

    def all_args = [length_args, ambig_args, homop_args]
        .findAll { arg -> arg != "" }
        .join(", ")

    """
    # Run Mothur screen.seqs
    mothur "#screen.seqs(fasta=${fasta}, ${all_args}, processors=${task.cpus})"

    # Rename outputs ONLY if the name doesn't match expected pattern
    # Handle both .fasta and .align extensions
    if [ -f *.good.fasta ] && [ ! -f ${meta.id}.good.fasta ]; then
        mv *.good.fasta ${meta.id}.good.fasta
    fi
    
    if [ -f *.good.align ] && [ ! -f ${meta.id}.good.align ]; then
        mv *.good.align ${meta.id}.good.align
    fi

    if [ -f *.bad.accnos ] && [ ! -f ${meta.id}.bad.accnos ]; then
        mv *.bad.accnos ${meta.id}.bad.accnos
    fi

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur --version 2>&1 | sed 's/^.*v\\.//; s/\\..*\$//')
    END_VERSIONS
    """
}
