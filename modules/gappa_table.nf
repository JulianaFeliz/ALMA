process GAPPA_TABLE {
    label 'process_low'

    publishDir "${params.outdir}/gappa_results", mode: 'copy'
    container 'https://depot.galaxyproject.org/singularity/python:3.9--1'
    containerOptions '--no-home -e'

    input:
    path per_query_tsvs
    path shared_files

    output:
    path "ALMA_OTU_abundance.tsv"

   
    script:
    """
    python3 \$(command -v gappa_to_table.py)
    echo "Forçando a criacao da tabela nova!"
    """
   
}
