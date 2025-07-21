include { BCL2FASTQ_local } from '/home/lrenteria/inspiired_nf/modules/local/bcl2fastq_local/main.nf'

workflow BCL2FASTQ_wfl {
    take:
    ch_bcl_input


    main:
    // Run bcl2fastq process
    BCL2FASTQ_local(ch_bcl_input)

    //save bcl2fastq output channels

    BCL2FASTQ_local.out.fastq           .collect().view { "Demultiplexed FASTQ files: ${it}" }.set { ch_demux_fastq }
    BCL2FASTQ_local.out.fastq_idx       .collect().view { "Demultiplexed FASTQ index files: ${it}" }.set { ch_fastq_idx }
    BCL2FASTQ_local.out.undetermined    .collect().view { "Undetermined FASTQ files: ${it}" }.set { ch_undetermined }
    BCL2FASTQ_local.out.undetermined_idx.collect().view { "Undetermined FASTQ index files: ${it}" }.set { ch_undetermined_idx }
    BCL2FASTQ_local.out.reports         .collect().view { "BCL2FASTQ reports: ${it}" }.set { ch_bcl_reports }
    BCL2FASTQ_local.out.stats           .collect().view { "BCL2FASTQ stats: ${it}" }.set { ch_bcl_stats }
    BCL2FASTQ_local.out.interop         .collect().view { "BCL2FASTQ interop files: ${it}" }.set { ch_bcl_interop }
    BCL2FASTQ_local.out.versions        .collect().view { "BCL2FASTQ versions: ${it}" }.set { ch_bcl_versions }


    emit:
    ch_demux_fastq
    ch_fastq_idx
    ch_undetermined
    ch_undetermined_idx
    ch_bcl_reports
    ch_bcl_stats
    ch_bcl_interop
    ch_bcl_versions
}