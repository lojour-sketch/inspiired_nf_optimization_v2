process LTRchecking_seqkit {

    memory '60 GB'  // much less than 50 GB, seqkit is very light
    cpus 2

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRFrag), val(mingDNA), val(idx)

    output:
    tuple val(meta),
          path("${meta}.ltr_filtered_R1.chunk${idx}.fastq.gz"),
          path("${meta}.ltr_filtered_R2.chunk${idx}.fastq.gz"),
          val(primer), val(ltrbit), val(largeLTRFrag), val(mingDNA), val(idx), emit: reads

    publishDir '../../../results/8_LTR_presence', mode: 'copy', overwrite: true

    script:
    def prefixLen = primer.size() + ltrbit.size()
    def startindex = prefixLen + 1
    def primer_ltr = primer + ltrbit
    def out_r1 = "${meta}.ltr_filtered_R1.chunk${idx}.fastq.gz"
    def out_r2 = "${meta}.ltr_filtered_R2.chunk${idx}.fastq.gz"
    def range = "${startindex}:-1"
    """
    # keep only the r2 reads that start with primer+LTR
    seqkit grep -s -r -p "^${primer_ltr}" "${read2}" -o temp_R2_kept.fastq.gz

    # Extract the IDs of kept reads
    seqkit seq -n temp_R2_kept.fastq.gz > keep_ids.txt

    # Extract paired reads from R1, after changing header names to match R1 headers
    sed 's/2:N:0/1:N:0/' keep_ids.txt > ids_R1.txt
    seqkit grep -n -f ids_R1.txt ${read1} | gzip > "${out_r1}"

    # We remove the primer+LTR sequence from the beginning of the reads
    seqkit subseq -r "${range}" temp_R2_kept.fastq.gz | gzip > "${out_r2}"

    # cleanup
    rm temp_R2_kept.fastq.gz keep_ids.txt ids_R1.txt

    """
}
