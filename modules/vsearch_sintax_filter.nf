process VSEARCH_SINTAX_FILTER {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/vsearch:2.22.1--h1b792b2_0'

    input:
    tuple val(meta), path(fasta)
    path sintax_db // O seu SINTAX.udb

    output:
    // Emitimos o FASTA filtrado, pronto para o agrupamento!
    path "${meta.id}.glomero.fasta", emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    # 1. Roda o SINTAX contra o Eukaryome
    # O cutoff 0.8 garante 80% de confiança para cravar o nome
    vsearch --sintax ${fasta} \\
            -db ${sintax_db} \\
            -tabbedout sintax_results.txt \\
            -sintax_cutoff 0.8 \\
            --threads ${task.cpus}

    # 2. Extrai APENAS os IDs que contêm a palavra "Glomeromycota"
    # O cut -f1 pega a primeira coluna (que é o nome da sequência)
    grep "Glomeromycota" sintax_results.txt | cut -f1 > glomero_ids.txt

    # 3. TRUQUE DE MESTRE: Usa AWK para filtrar o seu FASTA
    # Ele lê os IDs aprovados e cria um FASTA novo só com os Glomeros!
    awk 'NR==FNR{ids[">"\$1]; next} /^>/{f=(\$1 in ids)} f' glomero_ids.txt \({fasta} >\){meta.id}.glomero.fasta
    """
}