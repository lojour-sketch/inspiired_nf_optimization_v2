process SHORTREMOVE_local{

    publishDir '../../../results/9_short_seq_removal', mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path('*.short_removed_R1.fastq.gz'), path('*.short_removed_R2.fastq.gz')

    script:
    """
    seqkit seq -m 20 ${read1} -o ${sample}.short_removed_R1.fastq.gz
    seqkit seq -m 20 ${read2} -o ${sample}.short_removed_R2.fastq.gz
    """

}