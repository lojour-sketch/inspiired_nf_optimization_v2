process CREATE_KEYS_local {

    publishDir "${params.runfolderDir}/../results/15_keys/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(R1_vector_removed), path(R2_vector_removed), path(vqNames), path(r1_zipped), path(r2_zipped)

    output:
    tuple val(sample), path("${sample}_R1_numeric.fastq.gz"), path("${sample}_R2_numeric.fastq.gz"), path("${sample}_keys.RData"), emit: keys

    script:
    """
    #unzip r1 and r2 if they are compressed
    if [[ ${r1_zipped} == *.gz ]]; then
        gunzip -c ${r1_zipped} > ${sample}.R1.fastq
        r1=${sample}.R1.fastq
    fi

    if [[ ${r2_zipped} == *.gz ]]; then
        gunzip -c ${r2_zipped} > ${sample}.R2.fastq
        r2=${sample}.R2.fastq
    fi

    create_keys.R \$r1 \$r2 ${sample}

    """
}