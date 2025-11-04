process FINDVECTOR_local {

    memory '20 GB'
    publishDir "${params.runfolderDir}/../results/10_findvector", mode: 'symlink', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(project), val(mingDNA)
    path(vector_fasta)

    output:
    tuple val(meta), path("${meta}.R1_vector_removed.fastq.gz"), path("${meta}.R2_vector_removed.fastq.gz"), emit: ch_vector_removed

    script:
    """
    find_vector_2.R "${meta}" "${read1}" "${read2}" "${primer}${ltrbit}" "${vector_fasta}"

    """

}