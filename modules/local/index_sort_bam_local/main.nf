process INDEX_SORT_BAM_local {

    publishDir "${params.outdir}/14_index_sort_bam/${params.projectName}", mode: 'copy', overwrite: true

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