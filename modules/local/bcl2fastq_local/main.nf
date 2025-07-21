process BCL2FASTQ_local {
    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("${run_dir}/demuxed/alloCART/*_R*_001.fastq.gz")        , emit: fastq
    tuple val(meta), path("${run_dir}/demuxed/alloCART/*_I*_001.fastq.gz")       , optional:true, emit: fastq_idx
    tuple val(meta), path("${run_dir}/demuxed/Undetermined_S0_R*_001.fastq.gz")  , optional:true, emit: undetermined
    tuple val(meta), path("${run_dir}/demuxed/Undetermined_S0_I*_001.fastq.gz")  , optional:true, emit: undetermined_idx
    tuple val(meta), path("${run_dir}/demuxed/Reports")                             , emit: reports
    tuple val(meta), path("${run_dir}/demuxed/Stats")                               , emit: stats
    tuple val(meta), path("InterOp/*.bin")                       , emit: interop
    path("versions.yml")                                         , emit: versions

    script: 
    """

    # changed some parameters that differ from Patxi's script
    bcl2fastq \
        -R ${run_dir}/RUN329 \
        -o ${run_dir}/demuxed \
        -r 25 -p 25 -w 25 \
        --use-bases-mask I20Y159,I12,Y143 \
        --create-fastq-for-index-reads \
        --barcode-mismatches 2,2 \
        --no-lane-splitting   

    cp -r ${run_dir}/RUN329/InterOp .

    #logging bcl2fastq version used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    bcl2fastq: \$(bcl2fastq -V 2>&1 | grep -m 1 bcl2fastq | sed 's/^.*bcl2fastq v//')
    END_VERSIONS

    """


}