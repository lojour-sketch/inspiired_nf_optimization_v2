include { RCremoval_inspiired } from '../../../modules/local/rcremoval_local/main_inspiired_like'
//include { RCremoval_cutadapt } from '../../../modules/local/rcremoval_local/main_cutadapt'

workflow RCREMOVAL_wfl {
    take:
    ch_ltr_removed      //input in type: [sample, read1, read2, primer, ltrbit, largeLTRFrag, mingDNA, idx]
    ch_linkers          //input type: [sample, unique_linker, common_linker]

    main:
    RCremoval_inspiired(ch_ltr_removed, ch_linkers) 

    // Debug the outputs
    RCremoval_inspiired.out.reads.map{ it -> tuple(it[0], it[1], it[2]) }

    RCremoval_inspiired.out.reads.set{ ch_rc_removed } // [sample_id, rc_removed_R1, rc_removed_R2]

    emit:
    ch_rc_removed      
}