process BCL2FASTQ_local {

    //we want to save bcl2fastq output in our results/1_demuxed folder 
    publishDir "${params.runfolderDir}/../results/1_demuxed/${params.projectName}", pattern: 'results/**/*', mode: 'symlink', overwrite: true
    publishDir "${params.runfolderDir}/../results/1_demuxed/${params.projectName}", pattern: 'InterOp/*.bin', mode: 'symlink', overwrite: true

    input:
    tuple val(sample), val(primer), val(ltrbit), val(largeLTRFrag), val(project), val(mingDNA), val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("results/${project}/*/*_R*_001.fastq.gz")        , emit: fastq
    tuple val(meta), path("results/${project}/*/*_I*_001.fastq.gz")       , optional:true, emit: fastq_idx
    tuple val(meta), path("results/Undetermined_S0_R*_001.fastq.gz")  , optional:true, emit: undetermined
    tuple val(meta), path("results/Undetermined_S0_I*_001.fastq.gz")  , optional:true, emit: undetermined_idx
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
        --use-bases-mask I20Y159,I12,Y143 \\
        #         --processing-threads ${task.cpus}
        #         --sample-sheet ${samplesheet}
                 
    
    cp -r ${run_dir}/InterOp .

    #logging bcl2fastq version used
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    bcl2fastq: \$(bcl2fastq -V 2>&1 | grep -m 1 bcl2fastq | sed 's/^.*bcl2fastq v//')
    END_VERSIONS

    """


}