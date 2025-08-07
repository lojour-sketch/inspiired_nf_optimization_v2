include { SHORTREMOVE_local } from '../../../modules/local/short_seq_removal_local/main'

workflow SHORTREMOVE_wfl{
    take:
    ch_rc_removed

    main:
    SHORTREMOVE_local(ch_rc_removed)
    SHORTREMOVE_local.out.view{ "Short sequences removed: ${it}" }.set{ ch_short_removed }

    emit:
    ch_short_removed

}