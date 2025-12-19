process INSPIIRED_HITS_TO_SITESFINAL_local {

    publishDir "${params.runfolderDir}/../results/17_sitesfinal_singleend_inspiired/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(hits1), path(hits2), path(numeric1), path(numeric2), path(keys)

    output:
    tuple val(sample), path("*.RData"), emit: sitesfinal

    script:
    """
    
    inspiired_hits_to_sitesfinal.R ${hits1} ${hits2} ${sample}_sitesfinal.RData ${sample} ${keys}


    """


}