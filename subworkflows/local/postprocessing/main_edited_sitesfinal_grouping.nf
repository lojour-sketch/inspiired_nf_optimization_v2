include { INDEX_SORT_BAM_local } from '../../../modules/local/index_sort_bam_local/main_no_deduplication'
include { BAM_TO_ALLSITES_local } from '../../../modules/local/bam_to_allsites_local/main_eshrem'
include { ALLSITES_TO_SITESFINAL_edited_grouping_local } from '../../../modules/local/allsites_to_sitesfinal_local/main_edited_grouping'
include { ANNOTATE_SITES_local } from '../../../modules/local/annotate_sites_local/main_inspiired'

workflow POSTPROCESSING_twice_wfl {
    take:
    ch_aligned
    ch_processing_params
    ch_refGenome

    main:

    INDEX_SORT_BAM_local(ch_aligned)
    ch_sorted = INDEX_SORT_BAM_local.out.sorted

    //we need to join the channels to be able to give it to the process as input
    ch_baminput = ch_sorted.merge(ch_processing_params)

    BAM_TO_ALLSITES_local(ch_baminput)
    ch_allsites = BAM_TO_ALLSITES_local.out.allsites

    ALLSITES_TO_SITESFINAL_edited_grouping_local(ch_allsites)
    ch_sitesfinal = ALLSITES_TO_SITESFINAL_edited_grouping_local.out.sitesfinal

    //we also need to join the channels to be able to give it to the process as input
    ch_annotinput = ch_sitesfinal.merge(ch_refGenome)

    ANNOTATE_SITES_local(ch_annotinput)
    ch_annotated = ANNOTATE_SITES_local.out.annotated

    emit:
    ch_annotated
}