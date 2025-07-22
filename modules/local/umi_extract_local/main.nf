process UMI_EXTRACT_LOCAL {
   
    publishDir '/home/lrenteria/inspiired_nf/results/2_extractedumi'
   
    input:
    tuple val(sample_id), val(linker1), val(linker2), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.umi_R1.fastq.gz"), path("${sample_id}.umi_R2.fastq.gz"), path("${sample_id}.umi_extract.log")

    script:
    def (r1, r2) = reads
    //umi is 12nts long, but in the first sample is 13 nt long. it is followed by the common linker. 
    def bc_pattern = "(?P<umi_1>[ATCGN]{12,13})(?P<discard_2>${linker2})"
    """
    umi_tools extract \
        --extract-method=regex \
        --bc-pattern="${bc_pattern}" \
        -I ${r1} \
        --read2-in ${r2} \
        -S ${sample_id}.umi_R1.fastq.gz \
        --read2-out ${sample_id}.umi_R2.fastq.gz \
        --log=${sample_id}.umi_extract.log
    """
}