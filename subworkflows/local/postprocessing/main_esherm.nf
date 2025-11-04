include { DEDUPLICATE_local } from '../../../modules/local/deduplicate_local/main'
include { BAM_TO_ALLSITES_local } from '../../../modules/local/bam_to_allsites_local/main_eshrem'
include { ALLSITES_TO_SITESFINAL_esherm_local } from '../../../modules/local/allsites_to_sitesfinal_local/main_esherm'
include { ANNOTATE_SITES_local } from '../../../modules/local/annotate_sites_local/main_esherm'

workflow POSTPROCESSING_wfl {
    take:
    ch_aligned
    ch_processing_params
    ch_refGenome

    main:

    DEDUPLICATE_local(ch_aligned)
    ch_deduped = DEDUPLICATE_local.out.deduped

    //we need to join the channels to be able to give it to the process as input
    ch_baminput = ch_deduped.merge(ch_processing_params)

    BAM_TO_ALLSITES_local(ch_baminput)
    ch_allsites = BAM_TO_ALLSITES_local.out.allsites

    ALLSITES_TO_SITESFINAL_esherm_local(ch_allsites)
    ch_sitesfinal = ALLSITES_TO_SITESFINAL_esherm_local.out.sitesfinal

    //we also need to join the channels to be able to give it to the process as input
    ch_annotinput = ch_sitesfinal.merge(ch_refGenome)

    ANNOTATE_SITES_local(ch_annotinput)
    ch_annotated = ANNOTATE_SITES_local.out.annotated

    emit:
    ch_annotated
}