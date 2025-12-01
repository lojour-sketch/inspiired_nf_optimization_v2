process DEMUXING_FASTQ_local {

    publishDir "${params.runfolderDir}/../results/1_demuxed_undetermined/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple path(FASTQfolderDir), path(samplesheet)

    output:
    tuple path("*.R1.fq.gz"), path("*.R2.fq.gz"), emit: fastq

    script:
    """
    fqtk demux \\
        --inputs ${params.FASTQfolderDir}/Undetermined_S0_L001_R1_001.fastq.gz ${params.FASTQfolderDir}/Undetermined_S0_L001_I1_001.fastq.gz ${params.FASTQfolderDir}/Undetermined_S0_L001_R2_001.fastq.gz \\
        --read-structures 173T 12B 138T \\
        --sample-metadata ${params.demuxSampleSheet} \\
        --output .
    """
}
