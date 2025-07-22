process BCL2FASTQ_local {

    publishDir '/home/lrenteria/inspiired_nf/results/1_demuxed'

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("results/**_S[1-9]*_R?_00?.fastq.gz")          , emit: fastq
    tuple val(meta), path("results/**_S[1-9]*_I?_00?.fastq.gz")          , optional:true, emit: fastq_idx
    tuple val(meta), path("results/**Undetermined_S0*_R?_00?.fastq.gz")  , optional:true, emit: undetermined
    tuple val(meta), path("results/**Undetermined_S0*_I?_00?.fastq.gz")  , optional:true, emit: undetermined_idx
    tuple val(meta), path("results/Reports")                             , emit: reports
    tuple val(meta), path("results/Stats")                               , emit: stats
    tuple val(meta), path("InterOp/*.bin")                       , emit: interop
    path("versions.yml")                                         , emit: versions

    script: 
    """
    # changed some parameters that differ from Patxi's script
    bcl2fastq \\
        --runfolder-dir ${run_dir} \\
        --output-dir results \\
        --no-lane-splitting \\
        --barcode-mismatches 2,2 \\
        --create-fastq-for-index-reads \\
        -r 25 \\
        -p 25 \\
        -w 25 \\
        --barcode-mismatches 2 \\
        #         --processing-threads ${task.cpus}
        #         --sample-sheet ${samplesheet}
        #         --use-bases-mask I20Y159,I12,Y143 \\

    cp -r ${run_dir}/InterOp .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    bcl2fastq: \$(bcl2fastq -V 2>&1 | grep -m 1 bcl2fastq | sed 's/^.*bcl2fastq v//')
    END_VERSIONS

    """


}