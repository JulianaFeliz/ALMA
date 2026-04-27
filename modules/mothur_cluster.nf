process MOTHUR_CLUSTER {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/mothur/cluster", mode: 'copy'

    conda 'bioconda::mothur=1.48.0'
    container '/home/jrbfelix/ALMA/glomeromycota-pipeline/mothur_v1.48.sif'

    input:
    tuple val(meta), path(fasta), path(count_table)
    val cutoff
    val method

    output:
    tuple val(meta), path("*.list"), emit: list
    tuple val(meta), path("*.rabund"), emit: rabund, optional: true
    tuple val(meta), path("*.sabund"), emit: sabund, optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    """
    mothur "#cluster(fasta=${fasta}, count=${count_table}, method=${method}, cutoff=${cutoff}, processors=${task.cpus})" || true

    # Renomeia os arquivos de saida. Se o Mothur nao gerar rabund/sabund, o 'touch' cria um arquivo vazio pro Nextflow parar de encher o saco.
    mv *.list ${prefix}.list || true
    mv *.rabund ${prefix}.rabund 2>/dev/null || touch ${prefix}.rabund
    mv *.sabund ${prefix}.sabund 2>/dev/null || touch ${prefix}.sabund

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mothur: \$(mothur -v | head -n 1 | cut -d'=' -f2 | cut -d' ' -f1)
    END_VERSIONS
    """
}
