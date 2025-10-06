include { SHORTREMOVE_local } from '../../../modules/local/short_seq_removal_local/main'

workflow SHORTREMOVE_wfl{
    take:
    ch_rc_removed
    ch_primer_ltr

    main:

    SHORTREMOVE_local(ch_rc_removed, ch_primer_ltr)
    SHORTREMOVE_local.out.reads.set{ ch_short_removed }


    emit:
    ch_short_removed

}