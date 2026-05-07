process RCremoval_inspiired_local {

    cache 'lenient'

    publishDir "${params.outdir}/9_reverse_complement_removal/${params.projectName}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(project), val(mingDNA), val(unique_linker), val(common_linker)

    output:
    tuple val(meta), path("${meta}.rc_removed_R1.fastq.gz"), path("${meta}.rc_removed_R2.fastq.gz"), emit: reads

    script:
    """
    # Check if files are empty first
    r1_lines=\$(zcat ${read1} | wc -l)
    r2_lines=\$(zcat ${read2} | wc -l)
    
    echo "R1 has \$r1_lines lines"
    echo "R2 has \$r2_lines lines"
    
    if [ \$r1_lines -eq 0 ] || [ \$r2_lines -eq 0 ]; then
        echo "Empty input files detected, creating empty outputs"
        touch ${meta}.rc_removed_R1.fastq
        touch ${meta}.rc_removed_R2.fastq
        gzip ${meta}.rc_removed_R1.fastq
        gzip ${meta}.rc_removed_R2.fastq
    else
        # Decompress files to avoid R/ShortRead issues with compressed files
        gunzip -c ${read1} > temp_r1.fastq
        gunzip -c ${read2} > temp_r2.fastq
        
        # Run R script with uncompressed files
        rc_removal_inspiired.R "${meta}" "temp_r1.fastq" "temp_r2.fastq" "${largeLTRfrag}" "${common_linker}"
        
        # Cleanup temporary files
        rm -f temp_r1.fastq temp_r2.fastq
    fi
    """
}