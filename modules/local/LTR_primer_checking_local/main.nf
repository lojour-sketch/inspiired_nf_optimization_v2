process LTRchecking {

    publishDir '/home/lrenteria/inspiired_nf/results/7_LTR_presence/7.1_LTRpresencecheck', mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.id), path("${meta.id}.ltr_filtered_R1.fastq.gz"), path("${meta.id}.ltr_filtered_R2.fastq.gz"), path("${meta.id}.ltr_checked.log")


    script:
    def LTRbit = "GAAAATCTCTAGCA"
    def read1 = reads[0]
    def read2 = reads[1]

    """
    echo "Sample: ${meta.id}" > ${meta.id}.ltr_checked.log
    echo -n " Total R2 reads: " >> ${meta.id}.ltr_checked.log
    zcat ${read2} | awk 'NR % 4 == 2' | wc -l >> ${meta.id}.ltr_checked.log

    seqkit grep -s -i -p "${LTRbit}" "${read2}" -o "${meta.id}.ltr_filtered_R2.fastq.gz"
    seqkit seq -n "${meta.id}.ltr_filtered_R2.fastq.gz" | sed 's/[ \\t].*//; s/\\/\\([12]\\)\$//' > "${meta.id}.keep_ids.txt"

    echo -n "Reads with LTRbit: " >> ${meta.id}.ltr_checked.log
    wc -l < ${meta.id}.keep_ids.txt >> ${meta.id}.ltr_checked.log

    seqkit grep -f ${meta.id}.keep_ids.txt ${read1} -o ${meta.id}.ltr_filtered_R1.fastq.gz
    """


}