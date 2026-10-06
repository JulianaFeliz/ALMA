process GAPPA_ASSIGN {
    tag "Global_Taxonomy"
    label 'process_low'

    // Link do Galaxy Project
    container 'https://depot.galaxyproject.org/singularity/gappa:0.8.5--h077b44d_3'

    publishDir "${params.outdir}/gappa_results", mode: 'copy'    

    input:
    path jplace          // Recebe global_placement.jplace direto do RAXML
    path taxonomy_file   // params.taxonomy_ALMA
    val lwr_threshold    // params.gappa_lwr_threshold

    output:
    path "global.taxonomy.tsv", emit: taxonomy
    path "global.per_query.tsv", emit: per_query
    path "global.krona.txt", optional: true
    path "versions.yml", emit: versions

    script:
    """
    # Adicionada a flag --lwr-cutoff que estava faltando no seu comando original!
    gappa examine assign \\
        --jplace-path ${jplace} \\
        --taxon-file ${taxonomy_file} \\
        --per-query-results \\
        --krona \\
        --sativa \\
        --best-hit \\
        --lwr-cutoff ${lwr_threshold} \\
        --consensus-thresh 0.5 \\
        --resolve-missing-paths \\
        --out-dir ./

    # Renomear os arquivos para o padrão global de forma elegante
    if [ -f "profile.tsv" ]; then
        mv profile.tsv global.taxonomy.tsv
    fi

    if [ -f "per_query.tsv" ]; then
        mv per_query.tsv global.per_query.tsv
    fi
    
    if [ -f "krona.txt" ]; then
        mv krona.txt global.krona.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gappa: \$(gappa --version 2>&1 | head -n 1 | sed 's/gappa version: //g')
    END_VERSIONS
    """
}