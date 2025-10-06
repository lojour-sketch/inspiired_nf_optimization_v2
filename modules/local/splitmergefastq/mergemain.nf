process mergeChunks {

    input:
    tuple val(meta), path(filtered_r1_chunks), path(filtered_r2_chunks), val(outdir)


    publishDir { outdir }, mode: 'copy', overwrite: true

    output:
    tuple val(meta), path("${meta}.ltr_filtered_R1.fastq.gz"), path("${meta}.ltr_filtered_R2.fastq.gz"), emit: merged_reads

    script:
    """
    # Merge filtered R1 fastq chunks
    zcat ${filtered_r1_chunks} | gzip > ${meta}.ltr_filtered_R1.fastq.gz

    # Merge filtered R2 fastq chunks
    zcat ${filtered_r2_chunks} | gzip > ${meta}.ltr_filtered_R2.fastq.gz
    """
}