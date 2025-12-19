process ANNOTATE_SITES_local {

    publishDir "${params.runfolderDir}/../results/18_annotated/${params.projectName}", mode: 'symlink', overwrite: true
    

    input:
    tuple val(sample), path(allsites), path(sitesfinaltsv), val(refGenome), path(refGenomeFile), val(refKnowngeneFile)

    output:
    tuple val(sample), path("fig_annot_${sample}.pdf"), path("annotated_${sample}.xlsx"), emit: annotated

    script:
    """
    annotate_sites.R "${sample}" "${sample}_allsites.rds" "${sample}_sitesfinal.rds" "${refGenomeFile}" "${refKnowngeneFile}"
    """

}