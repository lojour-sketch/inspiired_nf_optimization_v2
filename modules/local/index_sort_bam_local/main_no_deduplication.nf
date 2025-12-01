process INDEX_SORT_BAM_local {

    publishDir "${params.runfolderDir}/../results/15_index_sort_bam_singleend/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(aligned)

    output:
    tuple val(sample), path("${sample}_nameSorted.bam"), emit: sorted

    script:
    """
    #we need to index the bam file
    samtools index -b ${sample}.Aligned.sortedByCoord.out.bam

    # and sort the BAM by read name
    samtools sort -n -o ${sample}_nameSorted.bam ${sample}.Aligned.sortedByCoord.out.bam
    """

}