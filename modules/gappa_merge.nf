process GAPPA_MERGE {
    label 'process_low'

    // Joga a árvore final direto na sua pasta de resultados!
    publishDir "${params.outdir}/gappa_results", mode: 'copy'

    // Container do Gappa
    container 'https://depot.galaxyproject.org/singularity/gappa:0.8.5--h077b44d_3'

    input:
    path jplace_files

    output:
    path "ALMA_with_samples.jplace.gz", emit: merged_jplace

    script:
    """
    # 1. Cria uma pastinha limpa só para receber a saída do Gappa
    mkdir -p saida_gappa

    # 2. Roda a fusão e joga o resultado lá dentro
    gappa edit merge \\
        --jplace-path ./*.jplace \\
        --out-dir saida_gappa/ \\
        --file-prefix ALMA_TMP
    
    # 3. Pega o arquivo de lá (seja qual for o nome que o Gappa deu), renomeia pro nome bonito e traz pra cá
    mv saida_gappa/ALMA_TMP* ./ALMA_with_samples.jplace
    
    # 4. Zipa a árvore com sucesso garantido!
    gzip ALMA_with_samples.jplace
    """
}
