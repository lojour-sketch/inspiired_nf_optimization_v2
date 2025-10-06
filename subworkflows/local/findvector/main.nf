include { FINDVECTOR_local } from '../../../modules/local/find_vector/main'

workflow FINDVECTOR_wfl {

    take:
    ch_rc_removed
    ch_primer_ltr
    vectorfasta

    main:
    FINDVECTOR_local(ch_rc_removed, ch_primer_ltr, vectorfasta)
    FINDVECTOR_local.out.ch_vector_removed.set{ ch_vector_removed }

    emit:
    ch_vector_removed

}