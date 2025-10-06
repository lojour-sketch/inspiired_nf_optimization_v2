include { DEREPLICATE_local } from '../../../modules/local/dereplicate_local/main.nf'


workflow DEREPLICATE_wfl {
    take:
    ch_genomic_reads

    main:

    DEREPLICATE_local(ch_genomic_reads)
    DEREPLICATE_local.out.ch_dereplicated.set{ ch_dereplicated }

    emit:
    ch_dereplicated
}