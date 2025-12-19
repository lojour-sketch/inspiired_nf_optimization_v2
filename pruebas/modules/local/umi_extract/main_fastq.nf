process UMIEXTRACT_local {

    publishDir "${params.runfolderDir}/../results/2_extractedumi/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample_id), val(linker1), val(linker2), path(reads), val(was_modified)

    output:
    tuple val(sample_id), path("${sample_id}.umi_R1.fastq.gz"), path("${sample_id}.umi_R2.fastq.gz"), path("${sample_id}.umi_extract.log")


    script:
    def (r1, r2) = reads
    //in order for this script to be available for every size of umi, we will add one nt to the umis that have a sample linker longer than the mean length
    def umi_length = was_modified ? 13 : 12
    def bc_pattern = "(?P<cell_1>${linker1})(?P<umi_1>[ATCGN]{${umi_length}})(?P<cell_2>${linker2})"

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