process LTRchecking_seqkit_local {

    cpus 6

    publishDir "${params.runfolderDir}/../results/8_LTR_presence/${params.projectName}", mode: 'symlink', overwrite: true

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

    # we don't want to remove the paired r1 reads because we will do single end alignment
    # so we will just change the name of the input file to the output name
    cp "${reads[0]}" "${out_r1}"

    # We remove the primer+LTR sequence from the beginning of the reads
    seqkit subseq -r "${range}" temp_R2_kept_${meta}.fastq.gz -o "${out_r2}"

    #remove temporal files
    rm temp_R2_kept_${meta}.fastq.gz keep_ids_${meta}.txt

    # MANUAL COPY FOR DEBUGGING
    # cp "${meta}.ltr_filtered_R1.fastq.gz" /home/lrenteria/inspiired_nf/results/8_LTR_presence/
    # cp "${meta}.ltr_filtered_R2.fastq.gz" /home/lrenteria/inspiired_nf/results/8_LTR_presence/
    """
}
