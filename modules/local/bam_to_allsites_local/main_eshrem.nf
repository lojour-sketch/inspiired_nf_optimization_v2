process BAM_TO_ALLSITES_local {

    publishDir "${params.runfolderDir}/../results/16_allsites/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(deduped), val(minPctIdent), val(maxAlignStart), val(maxFragLength)

    output:
    tuple val(sample), path("${sample}_allSites.tsv"), path("${sample}_tmpFile.tsv"),emit: allsites

    script:
    """
    #log the input parameters
    echo "minPctIdent: ${minPctIdent}"

    bam_to_allsites.py ${deduped} ${sample}_allSites.tsv ${maxAlignStart} ${minPctIdent} ${maxFragLength} ${sample}

    """
}