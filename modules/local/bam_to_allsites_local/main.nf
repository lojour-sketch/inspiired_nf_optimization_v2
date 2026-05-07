process BAM_TO_ALLSITES_local {

    publishDir "${params.outdir}/15_allsites/${params.projectName}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(sorted), val(minPctIdent), val(maxAlignStart), val(maxFragLength)

    output:
    tuple val(sample), path("${sample}_allSites.tsv"), path("${sample}_individualInsertions_notPaired.tsv"), emit: allsites

    script:
    """

    bam_to_allsites_start_by_strand.py "${sorted}" "${sample}_allSites.tsv" "${maxAlignStart}" "${minPctIdent}" "${maxFragLength}" "${sample}"

    """
}