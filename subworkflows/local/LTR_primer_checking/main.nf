include { LTRchecking } from '../../../modules/local/LTR_primer_checking_local/main'

workflow LTRCHECKING_wfl {
    take:
    reads

    main:
    LTRchecking(reads)
    LTRchecking.out.view { "LTR checking output: ${it}" }.set{ ch_ltr_checked }



    emit:
    ch_ltr_checked      //the format is: [sample_id, LTRfilteredR1, LTRfilteredR2, LTRlog]

}