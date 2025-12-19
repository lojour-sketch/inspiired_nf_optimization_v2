process DEMUXING_FASTQ_inspiired_like {

    publishDir "${params.runfolderDir}/../results/1_demuxed_undetermined_inspiired_like/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(samples), val(indices), path(FASTQfolderDir), path(samplesheet)

    output:
    tuple path("*_R1.fastq.gz"), path("*_R2.fastq.gz"), emit: fastq

    script:
    def samples_str = samples.join(',')
    def indices_str = indices.join(',')

    """
    demux_inspiired_like.R ${params.FASTQfolderDir} ${params.SampleSheet} ${indices_str} ${samples_str}
    """
}