include { RCremoval_local } from '../../../modules/local/rcremoval_local/main'

workflow RCREMOVAL_wfl {
    take:
    ch_ltr_removed      //input in type: [sample, read1, read2]

    main:
    RCremoval_local(ch_ltr_removed)
    RCremoval_local.out.view{ "RC removed channel: ${it}"}.set{ ch_rc_removed }

    emit:
    ch_rc_removed       //
}