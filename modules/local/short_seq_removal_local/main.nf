process SHORTREMOVE_local{

    publishDir "${params.runfolderDir}/../results/11_short_remove", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)

    output:
    tuple val(sample), path("${sample}.short_removed_R1.paired.fastq.gz"), path("${sample}.short_removed_R2.paired.fastq.gz"), emit: reads

    script:
    """
    seqkit seq -m ${mingDNA} ${read1} -o ${sample}.short_removed_R1.fastq
    seqkit seq -m ${mingDNA} ${read2} -o ${sample}.short_removed_R2.fastq

    # Check the filtered file sizes and content
    echo "Filtered file sizes:"
    ls -lh ${sample}.short_removed_R1.fastq ${sample}.short_removed_R2.fastq

    echo "Number of reads in filtered files:"
    seqkit stats ${sample}.short_removed_R1.fastq ${sample}.short_removed_R2.fastq

    # Check if files are empty
    if [[ ! -s "${sample}.short_removed_R1.fastq" ]] || [[ ! -s "${sample}.short_removed_R2.fastq" ]]; then
        echo "ERROR: One or both filtered files are empty"
        echo "This suggests the mingDNA value (${mingDNA}) might be too high"
        exit 1
    fi
    
    # Try to pair with verbose output
    echo "Attempting to pair files..."

    # as we have removed different reads from each file, we synchronize the files again
    seqkit pair -1 ${sample}.short_removed_R1.fastq -2 ${sample}.short_removed_R2.fastq

    # Check if pairing was successful
    if [[ -f "${sample}.short_removed_R1.paired.fastq" ]]; then
        echo "Pairing successful, compressing files..."
        rm -f ${sample}.short_removed_R*.paired.fastq.gz
        gzip ${sample}.short_removed_R1.paired.fastq
        gzip ${sample}.short_removed_R2.paired.fastq
    else
        echo "ERROR: Pairing failed - no paired files created"
        exit 1
    fi
    """

}