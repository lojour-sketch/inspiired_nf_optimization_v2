process BCL2FASTQ {
    //tag {"$meta.lane" ? "$meta.id"+"."+"$meta.lane" : "$meta.id" }
    label 'process_high'

    //container "nf-core/bcl2fastq:2.20.0.422"
    container "./bcl2fastq.sif"
    
    
    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("output/**_S[1-9]*_R?_00?.fastq.gz")          , emit: fastq
    tuple val(meta), path("output/**_S[1-9]*_I?_00?.fastq.gz")          , optional:true, emit: fastq_idx
    tuple val(meta), path("output/**Undetermined_S0*_R?_00?.fastq.gz")  , optional:true, emit: undetermined
    tuple val(meta), path("output/**Undetermined_S0*_I?_00?.fastq.gz")  , optional:true, emit: undetermined_idx
    tuple val(meta), path("output/Reports")                             , emit: reports
    tuple val(meta), path("output/Stats")                               , emit: stats
    tuple val(meta), path("InterOp/*.bin")                       , emit: interop
    path("versions.yml")                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BCL2FASTQ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def input_tar = run_dir.toString().endsWith(".tar.gz") ? true : false
    def input_dir = input_tar ? run_dir.toString() - '.tar.gz' : run_dir

    """
    echo "input directory: ${input_dir}"
    echo "run directory: ${run_dir}"
    echo "PWD below!"
    pwd
    echo "ls -l 220401_M03894_0295_000000000-KB998/RunInfo.xml:"
    ls -l 220401_M03894_0295_000000000-KB998/RunInfo.xml
    echo "Absolute path of input_dir:"
    readlink -f ${input_dir}
    pwd .


    if [ ! -d ${input_dir} ]; then
        mkdir -p ${input_dir}
    fi

    if ${input_tar}; then
        ## Ensures --strip-components only applied when top level of tar contents is a directory
        ## If just files or multiple directories, place all in $input_dir

        if [[ \$(tar -taf ${run_dir} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
            tar \\
                -C $input_dir --strip-components 1 \\
                -xavf \\
                $args2 \\
                $run_dir \\
                $args3
        else
            tar \\
                -C $input_dir \\
                -xavf \\
                $args2 \\
                $run_dir \\
                $args3
        fi
    fi

    cd ${run_dir} # added this line to be able to access the data

    # changed some parameters to match Patxi's script
    bcl2fastq \\
        --runfolder-dir ${input_dir} \\
        --output-dir results \\
        --use-bases-mask I20Y159,I12,Y143 \\
        --no-lane-splitting \\
        --barcode-mismatches 2,2 \\
        --create-fastq-for-index-reads \\
        -r 25 \\
        -p 25 \\
        -w 25 \\
        --sample-sheet ${samplesheet}
       # --processing-threads ${task.cpus}

    cp -r ${input_dir}/InterOp .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcl2fastq: \$(bcl2fastq -V 2>&1 | grep -m 1 bcl2fastq | sed 's/^.*bcl2fastq v//')
    END_VERSIONS
    """
}
