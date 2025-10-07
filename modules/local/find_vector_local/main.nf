process FINDVECTOR_local {

    conda "${params.runfolderDir}/../modules/local/find_vector/environment.yml"
    publishDir "${params.runfolderDir}/../results/10_findvector", mode: 'symlink', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)
    tuple val(meta), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)
    path(vector_fasta)

    output:
    tuple val(meta), path("${meta}.R1_vector_removed.fastq.gz"), path("${meta}.R2_vector_removed.fastq.gz"), emit: ch_vector_removed

    script:
    """
    find_vector.R "${meta}" "${read1}" "${read2}" "${primer}${ltrbit}" "${vector_fasta}"

    """

}