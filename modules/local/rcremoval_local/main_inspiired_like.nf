process RCremoval_inspiired_local {

    cpus 12
    cache 'deep'
    memory '100GB'

    publishDir "${params.runfolderDir}/../results/9_reverse_complement_removal_inspiired", mode: 'symlink', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(project), val(mingDNA), val(unique_linker), val(common_linker)

    output:
    tuple val(meta), path("${meta}.rc_removed_R1.fastq.gz"), path("${meta}.rc_removed_R2.fastq.gz"), emit: reads
    tuple val(meta), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)

    script:
    """
    rc_removal_inspiired.R "${meta}" "${read1}" "${read2}" "${largeLTRfrag}" "${common_linker}"
    """
}