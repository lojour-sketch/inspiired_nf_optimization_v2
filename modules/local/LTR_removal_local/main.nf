process LTRremove {

    publishDir '/home/lrenteria/inspiired_nf/results/7_LTR_presence/7.2_LTRremovalR2', mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2), path(log)

    output:
    tuple val(sample), path("*.ltr_removed_R1.fastq.gz"), path("*.ltr_removed_R2.fastq.gz")

    script:
    def LTRbit = "GAAAATCTCTAGCA"
    """
    #echo "Sample: ${sample}" > ${sample}.ltr_filtered.log
    zcat ${read2} | awk 'NR % 4 == 2 {gsub(/${LTRbit}/,"")} {print}' | gzip > ${sample}.ltr_removed_R2.fastq.gz

    zcat ${read1} | gzip > ${sample}.ltr_removed_R1.fastq.gz
    """

}