process SITESFINAL_TO_POINTS_local {

    publishDir "${params.runfolderDir}/../results/19_sitesfinal_to_points/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(allsites), path(sitesfinal), val(refGenome), path(refGenomeFile), val(refKnowngeneFile)

    output:
    tuple val(sample), path("fig_points_${sample}.pdf"), path("annotated_points_${sample}.xlsx"), emit: points

    script:
    """

    sitesfinal_to_points.R "${sample}" "${sample}_sitesfinal.rds" "${refGenomeFile}" "${refKnowngeneFile}"

    """

}