process NORMALIZE_index_length_local {

    publishDir "${params.outdir}/00_normalized_index_length/${params.projectName}", pattern: '*_normalized.csv', mode: 'copy', overwrite: true
    publishDir "${params.outdir}/00_normalized_index_length/${params.projectName}", pattern: 'modified_samples.txt', mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(samplesheet), path(runfolder)

    output:
    tuple val(sample), path("*_normalized.csv"), path(runfolder), emit: normalized
    path("modified_samples.txt"), emit: modified_samples, optional: true

    script:
    def sourceHash = java.security.MessageDigest.getInstance('SHA-256').digest(file("${moduleDir}/../../../bin/normalize_index_length.py").bytes).encodeHex()
    def parserHash = java.security.MessageDigest.getInstance('SHA-256').digest(file("${moduleDir}/../../../bin/validate_inputs.py").bytes).encodeHex()
    """
    # source_sha256=${sourceHash}; parser_sha256=${parserHash}
    normalize_index_length.py --samplesheet ${samplesheet}
    """

}
