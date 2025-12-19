process DEREPLICATE_local {

    publishDir "${params.runfolderDir}/../results/12_dereplicated_reads/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}.R1_unique_by_sequence.fastq.gz"), path("${sample}.R2_unique_by_sequence.fastq.gz"), path("${sample}_keys.RData"), emit: ch_dereplicated

    script:
    """
    dereplicate.R "${sample}" "${read1}" "${read2}"

    #compressing the output files
    gzip ${sample}.R1_unique_by_sequence.fastq
    gzip ${sample}.R2_unique_by_sequence.fastq
    """

}