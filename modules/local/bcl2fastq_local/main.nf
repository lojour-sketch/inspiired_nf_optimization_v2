process BCL2FASTQ_local {

    //we want to save bcl2fastq output in our results/1_demuxed folder 
    publishDir "${params.outdir}/1_demuxed/${params.projectName}", pattern: 'results/**/*', mode: 'copy', overwrite: true
    publishDir "${params.outdir}/1_demuxed/${params.projectName}", pattern: 'InterOp/*.bin', mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(primer), val(ltrbit), val(largeLTRFrag), val(project), val(mingDNA), val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("results/*/*/*_R*_001.fastq.gz")        , emit: fastq
    tuple val(meta), path("results/*/*/*_I*_001.fastq.gz")       , optional:true, emit: fastq_idx
    tuple val(meta), path("results/Undetermined_S0_R*_001.fastq.gz")  , optional:true, emit: undetermined
    tuple val(meta), path("results/Undetermined_S0_I*_001.fastq.gz")  , optional:true, emit: undetermined_idx
    tuple val(meta), path("results/Reports")                             , emit: reports
    tuple val(meta), path("results/Stats")                               , emit: stats
    tuple val(meta), path("InterOp/*.bin")                       , emit: interop

    shell: 
    '''

    bcl2fastq \\
        --runfolder-dir !{run_dir} \\
        --output-dir results \\
        --no-lane-splitting \\
        --barcode-mismatches 2,2 \\
        --create-fastq-for-index-reads \\
        -r 25 \\
        -p 25 \\
        -w 25 \\
        --use-bases-mask !{params.bcl_bases_mask ?: 'I20Y159,I12,Y143'} \\
        --sample-sheet !{samplesheet} \\

                 
    
    cp -r !{run_dir}/InterOp .


    # Route every declared project, including bcl2fastq's flat project layout.
    for project_dir in results/*; do
        [ -d "$project_dir" ] || continue
        for file in "$project_dir"/*.fastq.gz; do
            [ -f "$file" ] || continue
            sample_id=$(basename "$file" | sed 's/_S[0-9]*.*//')
            mkdir -p "$project_dir/$sample_id"
            mv "$file" "$project_dir/$sample_id/"
        done
    done

    '''


}
