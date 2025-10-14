process ALIGNMENT_local {

    scratch true
    publishDir '/home/lrenteria/inspiired_nf/results/14_alignment', mode: 'symlink', overwrite: true
    memory '40GB'

    input:
    tuple val(sample), path(r1), path(r2), path(keys), val(genome_index)

    output:
    tuple val(sample), path("${sample}.Aligned.sam.R1.*"), path("${sample}.Aligned.sam.R2.*"), path(keys), emit: aligned

    script:
    // we want to align the reads independently, so in single-end mode, and later join them with the keys
    """
    # we will run STAR twice
    STAR --genomeDir ${genome_index} \
        --readFilesIn ${r1} \
        --readFilesCommand zcat \
        --runThreadN 8 \
        --alignIntronMax 5 \
        --outFilterMismatchNoverReadLmax 0.15 \
        --outFilterScoreMin 27 \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax 1 \
        --scoreDelOpen -10 --scoreInsOpen -10 \
        --outSAMtype SAM \
        --outFileNamePrefix ${sample}.Aligned.sam.R1.

    STAR --genomeDir ${genome_index} \
        --readFilesIn ${r2} \
        --readFilesCommand zcat \
        --runThreadN 8 \
        --alignIntronMax 5 \
        --outFilterMismatchNoverReadLmax 0.15 \
        --outFilterScoreMin 27 \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax 1 \
        --scoreDelOpen -10 --scoreInsOpen -10 \
        --outSAMtype SAM \
        --outFileNamePrefix ${sample}.Aligned.sam.R2.
    """

}