process ANNOTATE_SITES_local {

    publishDir "${params.runfolderDir}/../results/18_annotated_esherm", mode: 'symlink', overwrite: true
    
    memory '40GB'

    input:
    tuple val(sample), path(allsites), path(allsites_nostandard), path(sitesfinaltsv), val(refGenome), path(refGenomeFile), val(refKnowngeneFile)

    output:
    tuple val(sample), path("fig_annot_${sample}.pdf"), path("annotated_${sample}.xlsx"), emit: annotated

    script:
    """

    annotate_sites.R "${sample}" "${sample}_sitesfinal.rds" "${sample}_allsites.rds" "${refGenomeFile}" "${refKnowngeneFile}"

    """

}