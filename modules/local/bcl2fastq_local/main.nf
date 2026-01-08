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
        --use-bases-mask I20Y159,I12,Y143 \\
        --sample-sheet !{samplesheet} \\

                 
    
    cp -r !{run_dir}/InterOp .


    if [ -d "results/!{project}" ]; then
        cd "results/!{project}"
        
        # Check if there are fastq.gz files directly in this directory
        if ls *.fastq.gz 1> /dev/null 2>&1; then
            echo "Files found directly in project dir - reorganizing into subdirectories..."
            
            for file in *.fastq.gz; do
                # Extract sample ID (everything before _S[0-9])
                sample_id=$(echo "$file" | sed 's/_S[0-9]*.*//')
                
                echo "Moving $file to $sample_id/"
                mkdir -p "$sample_id"
                mv "$file" "$sample_id/"
            done
            
            echo "Reorganization complete"
        
            # Show final structure
            echo "Final structure:"
            ls -la
        else
            echo "Files already in subdirectories - no reorganization needed"
        fi
    else
        echo "Warning: results/!{project} directory not found"
    fi

    '''


}