process RAXML_EPA {
    tag "$meta.id"
    cpus 8
    memory '16 GB'

    container 'https://depot.galaxyproject.org/singularity/raxml-ng:1.2.2--h6747034_1'


    input:
    tuple val(meta), path(msa)
    path tree
    path reference_aln

    output:
    tuple val(meta), path("${meta.id}.jplace"), emit: jplace
    path "versions.yml", emit: versions

    
    script:
    """
    # 1. Separate Queries and remove gaps (preserving headers)
    awk '/^>Q_/{p=1} /^>/ && !/^>Q_/{p=0} p{print}' ${msa} > query_only.fasta
    awk 'BEGIN{p=1} /^>Q_/{p=0} /^>/ && !/^>Q_/{p=1} p{print}' ${msa} > ref_only.fasta

    echo "Check - Reference sequences: \$(grep -c '^>' ref_only.fasta)"
    echo "Check - Query sequences: \$(grep -c '^>' query_only.fasta)"

    # 2. Avaliar parâmetros do modelo na árvore de referência
    
    raxml-ng --evaluate \\
             --msa ref_only.fasta \\
             --tree ${tree} \\
             --model GTR+G+I \\
             --prefix ref_eval \\
             --threads ${task.cpus}

    # 4. EPA-ng Placement using the optimized model
    epa-ng --tree ${tree} \\
           --ref-msa ref_only.fasta \\
           --query query_only.fasta \\
           --model ref_eval.raxml.bestModel \\
           --outdir . \\
           --redo \\
           --threads ${task.cpus}
    # 5. Finaliz
    mv epa_result.jplace ${meta.id}.jplace

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epa-ng: \$(epa-ng --version 2>&1 | head -n 1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
    
}
