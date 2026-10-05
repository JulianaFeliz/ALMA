process MOTHUR_GET_OTUREP {
    tag "${meta.id}"
    label 'process_medium'
    
    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/mothur:1.48.0--hb64bf22_1'
    
    input:
    tuple val(meta), path(list), path(fasta), path(count)
    
    output:
    tuple val(meta), path("*.rep.fasta"), emit: rep_seqs
    tuple val(meta), path("*.rep.names"), emit: rep_names, optional: true
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    """
    mothur "#get.oturep(list=${list}, fasta=${fasta}, count=${count}, method=abundance)"
    
    # Rename output files
    mv *.rep.fasta ${prefix}.rep.fasta || true
    mv *.rep.names ${prefix}.rep.names || true
    
    # Mark the OTU from sample
    sed -i "s/^>.*/&|${meta.id}/" *.fasta
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v | head -n 1 | cut -d'=' -f2 | cut -d' ' -f1)
    END_VERSIONS
    """
}
