process DEMUXING_FASTQ_local {

    debug true

    publishDir "${params.runfolderDir}/../results/1_demuxed_undetermined/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple path(FASTQfolderDir), path(samplesheet), val(readStructure)

    output:
    tuple path("*.R1.fq.gz"), path("*.R2.fq.gz"), emit: fastq

    script:
    // Get the absolute path that Nextflow staged
    def fastqDir = FASTQfolderDir.toString()
    """
    
    # Use find to get absolute paths (follows symlinks with -L)
    R1_FILES=\$(find -L ${FASTQfolderDir} -maxdepth 1 -name "Undetermined_*_R1*.fastq.gz" -type f)
    I1_FILES=\$(find -L ${FASTQfolderDir} -maxdepth 1 -name "Undetermined_*_I1*.fastq.gz" -type f)
    R2_FILES=\$(find -L ${FASTQfolderDir} -maxdepth 1 -name "Undetermined_*_R2*.fastq.gz" -type f)


    fqtk demux \\
        --inputs \$R1_FILES \$I1_FILES \$R2_FILES \\
        --read-structures ${readStructure} \\
        --sample-metadata ${samplesheet} \\
        --output .
    """
}
