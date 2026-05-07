process NORMALIZE_index_length_local {

    publishDir "${params.outdir}/00_normalized_index_length/${params.projectName}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(samplesheet), path(runfolder)

    output:
    tuple val(sample), path("*_normalized.csv"), path(runfolder), emit: normalized
    path("modified_samples.txt"), emit: modified_samples, optional: true

    script:
    """
    normalize_index_length.py --samplesheet ${samplesheet}
    """

}