process ALLSITES_TO_SITESFINAL_inspiired_local {

    publishDir "${params.runfolderDir}/../results/17_sitesfinal_inspiired/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(allsites), path(tmpFile)

    output:
    tuple val(sample), path("${sample}_allsites_nostandard.rds"), path("${sample}_sitesfinal.rds"), emit: sitesfinal

    script:
    """
    allsites_to_sitesfinal_Libe_inspiired.R "${sample}" "${allsites}"
    """


}