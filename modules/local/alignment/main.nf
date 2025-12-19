process ALIGNMENT_local {

    scratch true
    publishDir "${params.runfolderDir}/../results/14_alignment/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(r1), path(r2), val(genome_index)

    output:
    tuple val(sample), path("${sample}.Aligned.*"), emit: aligned

    script:
    """
    STAR --genomeDir ${genome_index} \
        --readFilesIn ${r1} ${r2}\
        --readFilesCommand zcat \
        --runThreadN 8 \
        --alignIntronMax 5 \
        --outFilterMismatchNoverReadLmax 0.15 \
        --outFilterScoreMin 27 \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax 1 \
        --scoreDelOpen -10 --scoreInsOpen -10 \
        --outSAMtype BAM SortedByCoordinate\
        --outFileNamePrefix ${sample}.

    """

}