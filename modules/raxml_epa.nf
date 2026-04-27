// --- MODULES/RAXML_EPA.NF ---

process RAXML_EVALUATE {
    tag "$meta.id"
    cpus 8
    memory '16 GB'

    publishDir "${params.outdir}/raxml_results", mode: 'copy'
    
    // Bolha 1: Exclusiva do RAxML-NG
    container 'https://depot.galaxyproject.org/singularity/raxml-ng:1.2.2--h6747034_1'

    input:
    tuple val(meta), path(msa)
    path tree
    path reference_aln

    output:
    // Nós já passamos os arquivos separados adiante para ganhar tempo!
    tuple val(meta), path("ref_eval.raxml.bestModel"), path("ref_only.fasta"), path("query_only.fasta"), emit: prep_data

    script:
    """
    # 1. Separar Queries e Referências (Agora buscando pelo ID da Amostra em vez de Q_)
    awk '/^>${meta.id}_/{p=1} /^>/ && !/^>${meta.id}_/{p=0} p{print}' ${msa} > query_only.fasta
    awk 'BEGIN{p=1} /^>${meta.id}_/{p=0} /^>/ && !/^>${meta.id}_/{p=1} p{print}' ${msa} > ref_only.fasta

    # 2. Avaliar parâmetros do modelo na árvore de referência
    raxml-ng --evaluate \\
             --msa ref_only.fasta \\
             --tree ${tree} \\
             --model GTR+G \\
             --prefix ref_eval \\
             --threads ${task.cpus}
    """
}

process EPA_NG_PLACEMENT {
    tag "$meta.id"
    cpus 8
    memory '16 GB'

    // Bolha 2: Exclusiva do EPA-ng
    container 'https://depot.galaxyproject.org/singularity/epa-ng:0.3.8--h9a82719_1'

    input:
    tuple val(meta), path(bestModel), path(ref_only), path(query_only)
    path tree

    output:
    tuple val(meta), path("${meta.id}.jplace"), emit: jplace
    path "versions.yml", emit: versions

    script:
    """
    # 3. EPA-ng Placement (Lendo os dados que saíram da Bolha 1)
    epa-ng --tree ${tree} \\
           --ref-msa ${ref_only} \\
           --query ${query_only} \\
           --model ${bestModel} \\
           --outdir . \\
           --redo \\
           --threads ${task.cpus}
    
    # 4. Finalizar
    mv epa_result.jplace ${meta.id}.jplace

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epa-ng: \$(epa-ng --version 2>&1 | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}

// 5. A Mágica: Este bloco engana o main.nf e faz ele achar que é um processo só!
workflow RAXML_EPA {
    take:
        msa_ch
        tree
        ref_aln

    main:
        RAXML_EVALUATE(msa_ch, tree, ref_aln)
        EPA_NG_PLACEMENT(RAXML_EVALUATE.out.prep_data, tree)

    emit:
        jplace = EPA_NG_PLACEMENT.out.jplace
        versions = EPA_NG_PLACEMENT.out.versions
}
