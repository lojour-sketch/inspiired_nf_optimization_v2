process INDEX_SORT_BAM_local {

    publishDir "${params.runfolderDir}/../results/15_index_sort_bam_singleend/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(aligned1_r1), path(aligned_r2)

    output:
    tuple val(sample), path("${sample}_R1_nameSorted.bam"), path("${sample}_R2_nameSorted.bam"), emit: sorted

    script:
    """
    #we need to index the bam file
    samtools index -b ${sample}.R1.Aligned.sortedByCoord.out.bam
    samtools index -b ${sample}.R2.Aligned.sortedByCoord.out.bam

    # and sort the BAM by read name
    samtools sort -n -o ${sample}_R1_nameSorted.bam ${sample}.R1.Aligned.sortedByCoord.out.bam
    samtools sort -n -o ${sample}_R2_nameSorted.bam ${sample}.R2.Aligned.sortedByCoord.out.bam
    """

}