process DEDUPLICATE_local {
    publishDir "${params.runfolderDir}/../results/15_deduplicated/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(files)

    output:
    tuple val(sample), path("${sample}_deduped_nameSorted.bam"),  emit: deduped

    script:
    """

    #we need to index the bam file
    samtools index -b ${files[0]} 
    

    #default umi-separator: _
    #default extract-umi-method: read_id, take the UMI from the read name
    #default unmapped-reads: discard
    umi_tools dedup \
                -I ${files[0]} \
                -S ${sample}_deduped.bam \
                --output-stats=${sample}_dedup \
                --method=directional \
                --chimeric-pairs=discard \
                --unpaired-reads=discard \
                --paired 

    
    # Sort the deduplicated BAM by read name
    samtools sort -n -o ${sample}_deduped_nameSorted.bam ${sample}_deduped.bam

    """
}