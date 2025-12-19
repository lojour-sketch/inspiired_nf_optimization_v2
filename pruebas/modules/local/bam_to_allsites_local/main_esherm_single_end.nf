process BAM_TO_HITS_local {

    publishDir "${params.runfolderDir}/../results/16_allsites/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(aligned_r1), path(aligned_r2), val(minPctIdent), val(maxAlignStart), val(maxFragLength)

    output:
    tuple val(sample), path("${sample}_R1_hits.tsv"), path("${sample}_R2_hits.tsv"), emit: hits

    script:
    """

    bam_to_hits_file.py ${aligned_r1} ${aligned_r2} ${maxAlignStart} ${minPctIdent} ${maxFragLength} ${sample}

    """
}