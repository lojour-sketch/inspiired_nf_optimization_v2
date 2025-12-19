process ALIGNMENT_SINGLE_local {

    scratch true
    publishDir "${params.runfolderDir}/../results/14_alignment/${params.projectName}", mode: 'symlink', overwrite: true
    memory '40GB'

    input:
    tuple val(sample), path(r1), path(r2), val(genome_index)

    output:
    tuple val(sample), path("${sample}.R1.Aligned.sortedByCoord.out.bam"), path("${sample}.R2.Aligned.sortedByCoord.out.bam"), emit: aligned // we should have two alignment files, one for each read

    script:

    """
    #we will perform two STAR runs, one for the R1 and one for the R2

    #unzip r1 and r2 if they are compressed, and assign to r1/r2 variable if not
    if [[ ${r1} == *.gz ]]; then
        gunzip -c ${r1} > ${sample}.R1.fastq
        r1=${sample}.R1.fastq
    else
        r1=${r1}
    fi

    if [[ ${r2} == *.gz ]]; then
        gunzip -c ${r2} > ${sample}.R2.fastq
        r2=${sample}.R2.fastq
    else
        r2=${r2}
    fi

    STAR --genomeDir ${genome_index} \
        --readFilesIn \$r1 \
        --runThreadN 1 \
        --alignIntronMax 5 \
        --outFilterMismatchNoverReadLmax 0.15 \
        --alignEndsType Local \
        --outFilterMultimapNmax 1 \
        --scoreDelOpen -10 --scoreInsOpen -10 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}.R1.

    STAR --genomeDir ${genome_index} \
        --readFilesIn \$r2 \
        --runThreadN 1 \
        --alignIntronMax 5 \
        --outFilterMismatchNoverReadLmax 0.15 \
        --alignEndsType Local \
        --outFilterMultimapNmax 1 \
        --scoreDelOpen -10 --scoreInsOpen -10 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}.R2.

    """

}