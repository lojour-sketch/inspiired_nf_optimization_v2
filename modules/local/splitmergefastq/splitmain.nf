process splitFastq {

    memory '32 GB'

    input:
    tuple val(meta), path(reads), val(primer), val(ltrbit), val(largeLTRFrag), val(mingDNA)

    output:
    tuple val(meta), path("chunks/*_1_N.part_*.fq.gz"), path("chunks/*_2_N.part_*.fq.gz"), val(primer), val(ltrbit), val(largeLTRFrag), val(mingDNA)

    script: 
    // Make chunks directory
    """
    mkdir -p chunks
    seqkit split -p 4 -O chunks --quiet ${reads[0]}
    seqkit split -p 4 -O chunks --quiet ${reads[1]}
    """
}