process CREATE_demux_samplesheet_local {

    publishDir "${params.outdir}/00_create_demux_samplesheet/${params.projectName}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(normalized_samplesheet), path(rundir)

    output:
    path("DemuxSampleSheet.tsv"), emit: demux_sheet

    script:
    def reverseHash = java.security.MessageDigest.getInstance('SHA-256').digest(file("${moduleDir}/../../../bin/create_demux_samplesheet_rev_comp_index2.py").bytes).encodeHex()
    def sourceHash = java.security.MessageDigest.getInstance('SHA-256').digest(file("${moduleDir}/../../../bin/create_demux_samplesheet.py").bytes).encodeHex()
    def parserHash = java.security.MessageDigest.getInstance('SHA-256').digest(file("${moduleDir}/../../../bin/validate_inputs.py").bytes).encodeHex()
    """
    # source_sha256=${sourceHash}; parser_sha256=${parserHash}
    # reverse_source_sha256=${reverseHash}
    if [ "${params.instrument}" == "MiSeq" ]
    then
        create_demux_samplesheet.py --samplesheet ${normalized_samplesheet}
    elif [ "${params.instrument}" == "NextSeq2000" ] || [ "${params.instrument}" == "NextSeq500" ]
    then
        create_demux_samplesheet_rev_comp_index2.py --samplesheet ${normalized_samplesheet}
    else
        echo "ERROR: Unsupported instrument: ${params.instrument}" >&2
        exit 1
    fi
    """
}
