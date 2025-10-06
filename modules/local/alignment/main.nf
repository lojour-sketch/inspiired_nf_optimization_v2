process ALIGNMENT_local {

    publishDir '/home/lrenteria/inspiired_nf/results/10_alignment', mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(r1), path(r2), path(keys)
    path genome_index

    output:
    tuple val(sample), path("*.Aligned.sam"), path(keys) emit: aligned

    script:
    """
    STAR --genomeDir ${genome_index} \
        --readFilesIn ${r1} ${r2} \
        --readFilesCommand cat \
        --runThreadN 8 \
        --alignIntronMax 5 \
        --outSAMtype SAM \
        --outFileNamePrefix ${sample}.Aligned.sam.
    """

}