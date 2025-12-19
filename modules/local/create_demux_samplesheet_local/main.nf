process CREATE_demux_samplesheet_local {

    publishDir "${params.runfolderDir}/../results/00_create_demux_samplesheet/${params.projectName}", mode: 'symlink', overwrite: true

    input:
    tuple val(sample), path(normalized_samplesheet), path(rundir)

    output:
    path("DemuxSampleSheet.tsv"), emit: demux_sheet

    script:
    """
    if [${params.instrument} == 'MiSeq']
    then
        create_demux_samplesheet.py --samplesheet ${normalized_samplesheet}
    elif  [${params.instrument} == 'NextSeq2000'] || [${params.instrument} == 'NextSeq500']
    then
        create_demux_samplesheet_rev_comp_index2.py --samplesheet ${normalized_samplesheet}
    fi

    """

}