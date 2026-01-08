process FINDVECTOR_local {

    publishDir "${params.runfolderDir}/../results/10_findvector/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(project), val(mingDNA)
    path(vector_fasta)

    output:
    tuple val(sample), path("${sample}.R1_*.fastq.gz"), path("${sample}.R2_*.fastq.gz"), emit: ch_vector_removed

    script:
    if (!vector_fasta.toString().endsWith('null'))
        """
        find_vector.R "${sample}" "${read1}" "${read2}" "${primer}${ltrbit}" "${vector_fasta}"

        """

    else
        """
        # If when: is false, just copy the input reads to the output
        cp ${read1} ${sample}.R1_skipped_vector.fastq.gz
        cp ${read2} ${sample}.R2_skipped_vector.fastq.gz

        """

}