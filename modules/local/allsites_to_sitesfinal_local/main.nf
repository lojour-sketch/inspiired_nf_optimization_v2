process ALLSITES_TO_SITESFINAL_edited_grouping_local {

    publishDir "${params.outdir}/16_sitesfinal/${params.projectName}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(allsites), path(indivFile)

    output:
    tuple val(sample), path("${sample}_allsites.rds"), path("${sample}_sitesfinal.rds"), emit: sitesfinal

    script:
    """
    allsites_to_sitesfinal_edited_grouping.R "${sample}" "${allsites}"
    """


}