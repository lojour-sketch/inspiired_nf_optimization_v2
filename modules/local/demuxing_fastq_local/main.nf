process DEMUXING_FASTQ_local {

    debug true

    publishDir "${params.outdir}/1_demuxed_undetermined/${params.projectName}", mode: 'copy', overwrite: true

    input:
    tuple path(FASTQfolderDir), path(samplesheet), val(readStructure)

    output:
    tuple path("*.R1.fq.gz"), path("*.R2.fq.gz"), emit: fastq

    script:
    // Get the absolute path that Nextflow staged
    def fastqDir = FASTQfolderDir.toString()
    """
    
    set -euo pipefail
    mapfile -t R1_FILES < <(find -L "${FASTQfolderDir}" -maxdepth 1 -name "Undetermined_*_R1*.fastq.gz" -type f | sort)
    mapfile -t I1_FILES < <(find -L "${FASTQfolderDir}" -maxdepth 1 -name "Undetermined_*_I1*.fastq.gz" -type f | sort)
    mapfile -t R2_FILES < <(find -L "${FASTQfolderDir}" -maxdepth 1 -name "Undetermined_*_R2*.fastq.gz" -type f | sort)
    if [[ \${#R1_FILES[@]} -ne 1 || \${#I1_FILES[@]} -ne 1 || \${#R2_FILES[@]} -ne 1 ]]; then
        echo "Provide one synchronized R1/I1/R2 triplet; concatenate corresponding lanes in the same order before running." >&2
        exit 1
    fi

    fqtk demux \\
        --inputs "\${R1_FILES[0]}" "\${I1_FILES[0]}" "\${R2_FILES[0]}" \\
        --read-structures ${readStructure} \\
        --sample-metadata ${samplesheet} \\
        --output .
    """
}
