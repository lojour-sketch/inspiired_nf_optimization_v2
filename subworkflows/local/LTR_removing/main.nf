
include { LTRremove } from '../../../modules/local/LTR_removal_local/main'

workflow LTRREMOVING_wfl{
    take:
    ch_ltr_checked

    main:
    LTRremove(ch_ltr_checked)
    LTRremove.out.view{ "Reads with LTR removed ${it}" }.set{ ch_ltr_removed }

    emit:
    ch_ltr_removed
}