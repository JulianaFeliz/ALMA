// --- MODULES/RAXML_EPA.NF ---

process RAXML_EVALUATE {
    tag "Global_Tree"
    cpus 8
    memory '16 GB'

    publishDir "${params.outdir}/raxml_results", mode: 'copy'
    
    // Bolha 1: Exclusiva do RAxML-NG
    container 'https://depot.galaxyproject.org/singularity/raxml-ng:1.2.2--h6747034_1'

    input:
    path msa          // O alinhamento único gerado pelo MAFFT
    path tree         // A árvore do ALMAdb
    path reference_aln // O alinhamento de referência do ALMAdb

    output:
    path "ref_eval.raxml.bestModel", emit: best_model
    path "ref_only.fasta", emit: ref_only
    path "query_only.fasta", emit: query_only

    script:
    """
    # 1. Separar Queries e Referências de forma universal e segura
    # A) Extrai apenas os nomes dos cabeçalhos do seu alinhamento de referência
    grep "^>" ${reference_aln} | cut -d ' ' -f 1 > ref_headers.txt

    # B) Usa AWK para varrer o MSA do MAFFT. O que NÃO estiver na lista de referências, é OTU (Query).
    # Este script funciona perfeitamente mesmo que as sequências FASTA tenham múltiplas linhas.
    awk 'NR==FNR{refs[\$1]; next} /^>/{keep= !(\$1 in refs)} keep' ref_headers.txt ${msa} > query_only.fasta

    # C) A referência (ALMAdb) já está perfeitamente alinhada e com as dimensões corretas da árvore.
    cp ${reference_aln} ref_only.fasta

    # 2. Avaliar parâmetros do modelo na árvore de referência
    # Como as referências são as mesmas, ele avalia super rápido.
    raxml-ng --evaluate \\
             --msa ref_only.fasta \\
             --tree ${tree} \\
             --model GTR+G \\
             --prefix ref_eval \\
             --threads ${task.cpus}
    """
}

process EPA_NG_PLACEMENT {
    tag "Global_Placement"
    cpus 8
    memory '32 GB' // Aumentamos a RAM, pois EPA-ng exige memória para o projeto global

    publishDir "${params.outdir}/raxml_results", mode: 'copy'

    // Bolha 2: Exclusiva do EPA-ng
    container 'https://depot.galaxyproject.org/singularity/epa-ng:0.3.8--h9a82719_1'

    input:
    path bestModel
    path ref_only
    path query_only
    path tree

    output:
    path "global_placement.jplace", emit: jplace
    path "versions.yml", emit: versions

    script:
    """
    # 3. EPA-ng Placement (Lendo os dados globais que saíram da Bolha 1)
    epa-ng --tree ${tree} \\
           --ref-msa ${ref_only} \\
           --query ${query_only} \\
           --model ${bestModel} \\
           --outdir . \\
           --redo \\
           --threads ${task.cpus}
    
    # 4. Finalizar renomeando para o padrão global
    mv epa_result.jplace global_placement.jplace

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epa-ng: \$(epa-ng --version 2>&1 | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}

// 5. 
workflow RAXML_EPA {
    take:
        msa
        tree
        ref_aln

    main:
        // Avalia o modelo baseado na referência
        RAXML_EVALUATE(msa, tree, ref_aln)
        
        // Faz o posicionamento das OTUs globais na árvore
        EPA_NG_PLACEMENT(
            RAXML_EVALUATE.out.best_model,
            RAXML_EVALUATE.out.ref_only,
            RAXML_EVALUATE.out.query_only,
            tree
        )

    emit:
        jplace = EPA_NG_PLACEMENT.out.jplace
        versions = EPA_NG_PLACEMENT.out.versions
}