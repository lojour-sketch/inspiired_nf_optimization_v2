process RCremoval_cutadapt {

    publishDir '/home/lrenteria/inspiired_nf/results/9_reverse_complement_removal_cutadapt', mode: 'symlink', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)
    tuple val(sample), val(unique_linker), val(common_linker)

    output:
    tuple val(meta), path("*.rc_removed_R1.fastq.gz"), path("*.rc_removed_R2.fastq.gz"), emit: reads
    tuple val(meta), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)

    script:
    def err_ltr = 3 / largeLTRfrag.size()
    def err_primer = 3 / primer.size()

    //we want to ensure Groovy takes it as integer
    def min_overlap_ltr = (int) Math.floor(largeLTRfrag.size() / 2) 
    def min_overlap_primer = (int) Math.floor(primer.size() / 2)

    //we will create a reverse complement function in groovy to create rc variables
    def rc(String dna) {
        def complement = [
            'A':'T', 'T':'A',
            'C':'G', 'G':'C'
        ]
        return dna.reverse().collect{
            base -> complement[base] ?: base 
        }.join()
    }
    """
    #!/bin/bash

    R1_out="${meta}.rc_removed_R1.fastq.gz"
    R2_out="${meta}.rc_removed_R2.fastq.gz"


    cutadapt \
        -j 8 \
        -g "${rc(largeLTRfrag)}" -e ${err_ltr} -O ${min_overlap_ltr} \
        -G "${rc(common_linker)}" -e ${err_primer} -O ${min_overlap_primer} \
        --minimum-length ${mingDNA} \
        -o "\$R1_out" -p "\$R2_out" \
        "${read1}" "${read2}"
    """
}
