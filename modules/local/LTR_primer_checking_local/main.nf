process LTRchecking_seqkit_local {

    publishDir "${params.outdir}/8_LTR_presence/${params.projectName}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads), val(primer), val(ltrbit), val(largeLTRFrag), val(project), val(mingDNA)

    output:
    tuple val(meta),
          path("${meta}.ltr_filtered_R1.fastq.gz"),
          path("${meta}.ltr_filtered_R2.fastq.gz"),
          val(primer), val(ltrbit), val(largeLTRFrag), val(project), val(mingDNA), emit: reads

    script:
    def prefixLen = primer.size() + ltrbit.size()
    def startindex = prefixLen + 1
    def primer_ltr = primer + ltrbit
    def out_r1 = "${meta}.ltr_filtered_R1.fastq.gz"
    def out_r2 = "${meta}.ltr_filtered_R2.fastq.gz"
    def range = "${startindex}:-1"
    """
    # keep only the r2 reads that start with primer+LTR
    seqkit grep -s -r -p "^${primer_ltr}" "${reads[1]}" -o temp_R2_kept_${meta}.fastq.gz

    # Extract the IDs of kept reads
    seqkit seq -n temp_R2_kept_${meta}.fastq.gz > keep_ids_${meta}.txt

    # Extract paired reads from R1, after changing header names to match R1 headers
    sed 's/2:N:0/1:N:0/' keep_ids_${meta}.txt > ids_R1_${meta}.txt
    seqkit grep -n -f ids_R1_${meta}.txt ${reads[0]} -o "${out_r1}"

    # We remove the primer+LTR sequence from the beginning of the reads
    seqkit subseq -r "${range}" temp_R2_kept_${meta}.fastq.gz -o "${out_r2}"

    #remove temporal files
    rm temp_R2_kept_${meta}.fastq.gz keep_ids_${meta}.txt ids_R1_${meta}.txt

    """
}
