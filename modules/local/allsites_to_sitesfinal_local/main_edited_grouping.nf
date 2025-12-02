process ALLSITES_TO_SITESFINAL_edited_grouping_local {

    publishDir "${params.runfolderDir}/../results/17_sitesfinal_inspiired/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(allsites), path(tmpFile)

    output:
    tuple val(sample), path("${sample}_allsites_nostandard.rds"), path("${sample}_sitesfinal.rds"), emit: sitesfinal

    script:
    """
    allsites_to_sitesfinal_edited_grouping.R "${sample}" "${allsites}"
    """


}